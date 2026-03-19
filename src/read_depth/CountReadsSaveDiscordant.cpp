// compiling on RAP (to get libcurl + dependencies):
// sudo apt --yes install libdeflate-dev # (to get libdeflate headers)
// tar -zxf /mnt/project/lohdata/resources/htslib/htslib-1.18-libdeflate.tar.gz -C .
// g++ -O3 -Wall CountReadsSaveDiscordant.cpp -o countReadsSaveDiscordant -Ihtslib-1.18-libdeflate/ htslib-1.18-libdeflate/libhts.a -lcurl -ldeflate -lz -lpthread

// compiling on O2 (without curl):
// g++ -O3 -Wall CountReadsSaveDiscordant.cpp -o countReadsSaveDiscordant -I/n/data1/bwh/medicine/loh/ploh/external_software/htslib/htslib-1.18/ -L/n/data1/bwh/medicine/loh/ploh/external_software/htslib/htslib-1.18/ -L/n/data1/bwh/medicine/loh/ploh/external_software/libdeflate-1.18/ -Wl,-Bstatic -lhts -ldeflate -Wl,-Bdynamic -lz -lpthread

#include <iostream>
#include <vector>
#include <numeric>
#include <map>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <cassert>
#include <sys/time.h>
#include <htslib/sam.h>
#include <htslib/bgzf.h>

#include "RefMasks.hpp"

using namespace std;

struct QualFlag {
  unsigned char qual;
  unsigned short flag;
  QualFlag(unsigned char _qual=0, unsigned short _flag=0) : qual(_qual), flag(_flag) {}
};

typedef unsigned long long uint64;

inline uint64 makeKey(uint64 posLeft, uint64 posRight) {
  return (posLeft<<32)|posRight;
}
inline void incrementBinCount(uint16_t &binCount) {
  if (binCount < 65535) binCount++;
}

void writeDiscordant(FILE *foutDiscordant, const char *IDstr,
		     const map <uint64, QualFlag> &regionToQualFlag,
		     const vector <QualFlag> &posToQualFlag, const char *chrStr) {
  map < uint64, pair <QualFlag, QualFlag> > regionToQualFlagPair;
  for (map <uint64, QualFlag>::const_iterator it = regionToQualFlag.begin();
       it != regionToQualFlag.end(); it++) {
    int pos = (it->first)>>32, mpos = (it->first)&0xFFFFFFFF;
    QualFlag qf = it->second;
    QualFlag qfMate;
    map <uint64, QualFlag>::const_iterator itOpp = regionToQualFlag.find(makeKey(mpos, pos));
    if (itOpp != regionToQualFlag.end())
      qfMate = itOpp->second;
    else if (qf.flag & BAM_FSUPPLEMENTARY)
      qfMate = posToQualFlag[mpos]; // if SA and mate not in discordant map, look up by pos

    if (qfMate.qual) { // store region, flipping read/mate if necessary
      if (pos < mpos)
	regionToQualFlagPair[makeKey(pos, mpos)] = make_pair(qf, qfMate);
      else
	regionToQualFlagPair[makeKey(mpos, pos)] = make_pair(qfMate, qf);
    }
  }
  for (map < uint64, pair <QualFlag, QualFlag> >::iterator it = regionToQualFlagPair.begin();
       it != regionToQualFlagPair.end(); it++) {
    uint64 key = it->first;
    pair <QualFlag, QualFlag> qfs = it->second;
    fprintf(foutDiscordant, "%s\t%s\t%Ld\t%Ld\t%d\t%d\t%d\n",
	    IDstr, chrStr, (key>>32)+1, (key&0xFFFFFFFF)+1, qfs.first.flag, qfs.second.flag,
	    min(qfs.first.qual, qfs.second.qual));
  }
}

// assumes outputs have been zero-initialized
void processCram(uint32_t readsWithGCpct[MAX_GC_PCT+1], uint64_t &totReadsNoDups,
		 uint64_t &totReadsAllowDups, vector < vector <uint16_t> > &readsPerBin,
		 int binSize, const char *cramFile, int shift,
		 const vector < vector <uint16_t> > &GCmasks, FILE *foutDiscordant,
		 const char *IDstr) {

  // open cram file
  samFile *fin = hts_open(cramFile, "r"); assert(fin!=NULL);
  // only extract required fields
  hts_set_opt(fin, CRAM_OPT_REQUIRED_FIELDS,
	      SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_RNEXT | SAM_PNEXT);

  // find which tid's in cram header correspond to chr1-22,X,Y
  bam_hdr_t *hdr = sam_hdr_read(fin);
  int n_targets = hdr->n_targets;
  char *tid_to_chr = (char *) calloc(n_targets, 1); // map target IDs to 1-22,23=X,24=Y or else 0
  for (int tid = 0; tid < n_targets; tid++) {
    const char *tname = hdr->target_name[tid];
    if (strlen(tname) > 3 && tname[0]=='c' && tname[1]=='h' && tname[2]=='r') {
      if (strlen(tname) == 4) {
	if (isdigit(tname[3]))
	  tid_to_chr[tid] = tname[3]-'0';
	else if (tname[3]=='X')
	  tid_to_chr[tid] = 23;
	else if (tname[3]=='Y')
	  tid_to_chr[tid] = 24;
      }
      else if (strlen(tname) == 5 && isdigit(tname[3]) && isdigit(tname[4]))
	tid_to_chr[tid] = 10*(tname[3]-'0') + tname[4]-'0';
    }
  }

  bam1_t *aln = bam_init1();

  // declare discordant read-related variables and lookup tables
  const int minDiscordantMAPQ = 30, minPosDiff = 5000;
  map <uint64, QualFlag> regionToQualFlag; // discordant region (pos, mpos) -> (qual, flag)
  vector <QualFlag> posToQualFlag; // all qualifying (discordant or not) pos -> (qual, flag)
  int prev_tid = -1;

  // iterate through cram file
  int binIndex = 0, binStart = 0, binEnd = binSize;
  while (true) {
    int ret = sam_read1(fin, hdr, aln);
    if (ret == -1) break; // end of stream
    if (ret < -1) {
      fprintf(stderr, "ERROR: sam_read1 returned %d; truncated file?\n", ret);
      exit(1);
    }

    int chr0 = tid_to_chr[aln->core.tid]-1;

    if (aln->core.tid != prev_tid) { // starting new chromosome
      if (prev_tid != -1) // write collated discordant reads for previous chromosome
	writeDiscordant(foutDiscordant, IDstr, regionToQualFlag, posToQualFlag,
			hdr->target_name[prev_tid]);
      if (chr0 == -1) // not autosome or X/Y => done
	break;
      else { // reset discordant read lookup tables
	regionToQualFlag.clear();
	posToQualFlag.clear();
	posToQualFlag.resize(hdr->target_len[aln->core.tid]);
	prev_tid = aln->core.tid;
      }
    }

    if (!(aln->core.tid==aln->core.mtid)) // require read and mate to align to same chromosome
      continue;

    int flag = aln->core.flag;
    int qual = aln->core.qual;
    int pos = aln->core.pos; // leftmost position of alignment (zero-based coord)

    // update lookup tables for discordant read analysis
    if (qual >= minDiscordantMAPQ // high MAPQ, read+mate mapped, not dup/secondary/QCfail
	&& !(flag & (BAM_FUNMAP|BAM_FMUNMAP|BAM_FSECONDARY|BAM_FDUP|BAM_FQCFAIL))) {
      if (qual > posToQualFlag[pos].qual) // update best MAPQ observed at POS
	posToQualFlag[pos] = QualFlag(qual, flag);
      if (abs(aln->core.mpos - pos) >= minPosDiff) // discordant
	regionToQualFlag[makeKey(pos, aln->core.mpos)] = QualFlag(qual, flag);
    }

    // move on to read-counting

    if (flag&(BAM_FMREVERSE|BAM_FUNMAP|BAM_FMUNMAP|BAM_FSECONDARY|BAM_FQCFAIL|BAM_FSUPPLEMENTARY))
      continue; // require mate forward, read and mate mapped, not secondary/QCfail/supplementary
    if (!(flag & BAM_FREVERSE)) // require read reverse
      continue;
    
    uint16_t gcMasks = GCmasks[chr0][pos>>shift];
    if (gcMasks & ALL_MASKS_EXCEPT_GC) // filter if pos has N in window or is in lc|sv|cnvmask
      continue;

    if (!(binStart <= pos && pos < binEnd)) { // update binIndex, binEnd if necessary
      binIndex = pos / binSize;
      binStart = binIndex * binSize;
      binEnd = binStart + binSize;
    }

    if (!(flag & BAM_FDUP)) // if not duplicate read, update bin count
      incrementBinCount(readsPerBin[chr0][binIndex]); // update bin count if not overflowing

    if (noMasks(gcMasks)) { // if no mask bits set, update GC profile counts
      int gc = getGCpct(gcMasks);
      readsWithGCpct[gc] += !(flag & BAM_FDUP); // update GC count disallowing duplicate reads
      totReadsAllowDups++;
    }
  }

  for (int gc = 0; gc <= MAX_GC_PCT; gc++)
    totReadsNoDups += readsWithGCpct[gc];

  bam_destroy1(aln);
  bam_hdr_destroy(hdr);
  free(tid_to_chr);
  sam_close(fin);
}

void writePerBinCounts(BGZF *foutCounts, const char *IDstr, uint16_t masks,
		       const vector <float> &GCprofile, int binSize,
		       const vector < vector <uint16_t> > &readsPerBin) {
  int IDlen = strlen(IDstr);
  assert(bgzf_write(foutCounts, &IDlen, sizeof(IDlen)) > 0);
  assert(bgzf_write(foutCounts, IDstr, IDlen) > 0);
  assert(bgzf_write(foutCounts, &masks, sizeof(masks)) > 0);
  assert(bgzf_write(foutCounts, &GCprofile[0], sizeof(GCprofile[0])*GCprofile.size()) > 0);
  assert(bgzf_write(foutCounts, &binSize, sizeof(binSize)) > 0);
  for (int i = 0; i < NUM_CHROMS; i++) {
    int numBins = readsPerBin[i].size();
    assert(bgzf_write(foutCounts, &numBins, sizeof(numBins)) > 0);
    assert(bgzf_write(foutCounts, &readsPerBin[i][0], numBins * sizeof(readsPerBin[i][0])) > 0);
  }
}

int main(int argc, char *argv[]){

  if (argc != 6) {
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "- arg1 = ID\n");
    fprintf(stderr, "- arg2 = cram file or URL\n");
    fprintf(stderr, "- arg3 = binary file containing precomputed GC content and masks\n");
    fprintf(stderr, "- arg4 = bin size (e.g., 1000 for 1kb bins)\n");
    fprintf(stderr, "- arg5 = output prefix\n");
    exit(1);
  }

  struct timeval tv;
  gettimeofday(&tv, NULL); double tPrev = tv.tv_sec + tv.tv_usec*1e-6;

  const char *IDstr = argv[1];
  const char *cramFile = argv[2];
  const char *GCmasksFile = argv[3];
  int binSize; sscanf(argv[4], "%d", &binSize);
  const char *outPrefix = argv[5];
  const int READ_LEN = 150; // could set as input param, but only used for coverage estimation

  cout << "Counting reads for ID " << IDstr << " in " << binSize << "bp bins" << endl;

  cout << "Opening output files: " << outPrefix << ".{counts.bin.bgz,GCprofile.txt}" << endl;
  char buf[1000];
  sprintf(buf, "%s.counts.bin.bgz", outPrefix);
  BGZF *foutCounts = bgzf_open(buf, "w"); assert(foutCounts!=NULL);
  sprintf(buf, "%s.GCprofile.txt", outPrefix);
  FILE *foutGCprofile = fopen(buf, "w"); assert(foutGCprofile!=NULL);
  sprintf(buf, "%s.discordant.txt", outPrefix);
  FILE *foutDiscordant = fopen(buf, "w"); assert(foutDiscordant!=NULL);

  // load precomputed GC content and masks
  cout << "Loading precomputed GC content and masks: " << GCmasksFile << endl;
  int shift; // bit shift
  vector <int> chrLens(NUM_CHROMS); // chromosome lengths
  vector < vector <uint16_t> > GCmasks(NUM_CHROMS); // GCmasks vectors
  readGCmasks(GCmasksFile, shift, chrLens, GCmasks);

  // tabulate bps with each GC percent
  vector <uint32_t> bpsWithGCpct(MAX_GC_PCT+1);
  for (int i = 0; i < NUM_CHROMS; i++) {
    int numSegs = (chrLens[i]+(1<<shift)-1)>>shift;
    for (int j = 0; j < numSegs; j++)
      if (noMasks(GCmasks[i][j])) // segment is good
	bpsWithGCpct[getGCpct(GCmasks[i][j])] += 1<<shift;
  }

  uint32_t bpsInGCprofile = accumulate(bpsWithGCpct.begin(), bpsWithGCpct.end(), 0);
  cout << "Number of bps included in GC profile: " << bpsInGCprofile << endl;

  // allocate arrays for read counts
  uint32_t readsWithGCpct[MAX_GC_PCT+1]; // not allowing duplicate reads
  memset(readsWithGCpct, 0, sizeof(readsWithGCpct));
  uint64_t totReadsNoDups = 0, totReadsAllowDups = 0;
  vector < vector <uint16_t> > readsPerBin(NUM_CHROMS);
  for (int i = 0; i < NUM_CHROMS; i++)
    readsPerBin[i].resize((chrLens[i]+binSize-1) / binSize);

  // iterate through cram file
  processCram(readsWithGCpct, totReadsNoDups, totReadsAllowDups, readsPerBin, binSize, cramFile,
	      shift, GCmasks, foutDiscordant, IDstr);

  // write GC profile data
  fprintf(foutGCprofile, "%s", IDstr);
  fprintf(foutGCprofile, "\t%.3f", totReadsNoDups * 2.0 * READ_LEN / bpsInGCprofile); // coverage
  fprintf(foutGCprofile, "\t%.3f", totReadsAllowDups * 2.0 * READ_LEN / bpsInGCprofile);
  vector <float> GCprofile(MAX_GC_PCT+1);
  for (int gc = 0; gc <= MAX_GC_PCT; gc++) {
    GCprofile[gc] = readsWithGCpct[gc] / (double) (bpsWithGCpct[gc]+1e-9);
    fprintf(foutGCprofile, "\t%.4f", GCprofile[gc]);
  }
  fprintf(foutGCprofile, "\n");

  // write per-bin non-duplicate read counts
  writePerBinCounts(foutCounts, IDstr, ALL_MASKS_EXCEPT_GC, GCprofile, binSize, readsPerBin);
  
  gettimeofday(&tv, NULL); double tCur = tv.tv_sec + tv.tv_usec*1e-6;
  cout << "Processed " << IDstr << ": " << tCur-tPrev << " sec" << endl;

  bgzf_close(foutCounts);
  fclose(foutGCprofile);
  fclose(foutDiscordant);

  return 0;
}
