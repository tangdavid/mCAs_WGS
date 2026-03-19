// compiling on RAP (to get libcurl + dependencies):
// sudo apt --yes install libdeflate-dev # (to get libdeflate headers)
// tar -zxf /mnt/project/lohdata/resources/htslib/htslib-1.18-libdeflate.tar.gz -C .
// g++ -O3 -Wall ComputeDepthProfiles.cpp -o computeDepthProfiles -Ihtslib-1.18-libdeflate/ htslib-1.18-libdeflate/libhts.a -lcurl -ldeflate -lz -lpthread

// compiling on O2 (without curl):
// g++ -O3 -Wall ComputeDepthProfiles.cpp -o computeDepthProfiles -I/n/data1/bwh/medicine/loh/ploh/external_software/htslib/htslib-1.18/ -L/n/data1/bwh/medicine/loh/ploh/external_software/htslib/htslib-1.18/ -L/n/data1/bwh/medicine/loh/ploh/external_software/libdeflate-1.18/ -Wl,-Bstatic -lhts -ldeflate -Wl,-Bdynamic -lz -lpthread

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <set>
#include <algorithm>
#include <numeric>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cassert>
#include <htslib/bgzf.h>

#include "RefMasks.hpp"

using namespace std;

struct GCpctCoeffs {
  int c[MAX_GC_PCT+1];
};

struct State {
  double logP[3];
  char prev[3];
};

inline double sq(double x) { return x*x; }


vector < vector <GCpctCoeffs> > computeCoeffs(int shift, const vector <int> &chrLens,
					      const vector < vector <uint16_t> > &GCmasks,
					      int binSize, int masks, const char *maskBinsFile) {

  set <int> isMasked[24];
  if (maskBinsFile != NULL) {
    int ctrMasked = 0;
    cout << "Reading list of bins to mask: " << maskBinsFile << flush;
    ifstream finMaskBins(maskBinsFile); assert(finMaskBins);
    string chrStr; int binStart, binEnd;
    while (finMaskBins >> chrStr >> binStart >> binEnd) {
      assert(binEnd-binStart == binSize);
      int chr0 = parseChrStrToChr0(chrStr.c_str());
      isMasked[chr0].insert(binStart);
      ctrMasked++;
    }
    finMaskBins.close();
    cout << " [" << ctrMasked << " bins masked]" << endl;
  }

  vector < vector <GCpctCoeffs> > coeffs(NUM_CHROMS);
  for (int chr0 = 0; chr0 < NUM_CHROMS; chr0++) {
    int numBins = (chrLens[chr0]+binSize-1) / binSize;
    coeffs[chr0].resize(numBins);
    for (int b = 0; b < numBins; b++) {
      int binStart = b*binSize;
      if (isMasked[chr0].count(binStart)) continue;
      int binEnd = min(binStart+binSize, chrLens[chr0]);
      for (int posShift = binStart>>shift; (posShift<<shift) < binEnd; posShift++) {
	uint16_t gcMasks = GCmasks[chr0][posShift];
	if (!(gcMasks & masks)) // no mask bits set
	  coeffs[chr0][b].c[getGCpct(gcMasks)] +=
	    min(binEnd, (posShift+1)<<shift) - max(binStart, posShift<<shift);
      }
    }
  }
  return coeffs;
}

vector < vector <int> > computeUnmaskedBpPerBin(const vector < vector <GCpctCoeffs> > &coeffs) {
  vector < vector <int> > unmaskedBpPerBin(NUM_CHROMS);
  for (int chr0 = 0; chr0 < NUM_CHROMS; chr0++) {
    int numBins = coeffs[chr0].size();
    unmaskedBpPerBin[chr0].resize(numBins);
    for (int b = 0; b < numBins; b++) {
      for (int gc = 0; gc <= MAX_GC_PCT; gc++)
	unmaskedBpPerBin[chr0][b] += coeffs[chr0][b].c[gc];
    }
  }
  return unmaskedBpPerBin;
}

vector < vector <int> > computeChunksWriteHeader(ofstream &fout,
						 const vector < vector <int> > &unmaskedBpPerBin,
						 int binSize, double unmaskedBpPerChunk) {

  vector < vector <int> > chunkBounds(NUM_CHROMS);
  fout << "ID";
  for (int chr0 = 0; chr0 < NUM_CHROMS; chr0++) {
    int numBins = unmaskedBpPerBin[chr0].size();
    double unmaskedBpTot =
      accumulate(unmaskedBpPerBin[chr0].begin(), unmaskedBpPerBin[chr0].end(), 0);
    int numChunks = (int) (unmaskedBpTot / unmaskedBpPerChunk);
    chunkBounds[chr0].resize(numChunks+1);
    chunkBounds[chr0][0] = 0;
    int bCur = 0, unmaskedBpCum = 0;
    for (int c = 0; c < numChunks; c++) {
      while (unmaskedBpCum < c * unmaskedBpTot / numChunks)
	unmaskedBpCum += unmaskedBpPerBin[chr0][bCur++];
      chunkBounds[chr0][c] = bCur;
    }
    chunkBounds[chr0][numChunks] = numBins;
    
    for (int c = 0; c < numChunks; c++) {
      fout << "\tchr";
      if (chr0<22) fout << chr0+1; else fout << (char) ('X'+chr0-22);
      fout << "_" << chunkBounds[chr0][c]*binSize << "_" << chunkBounds[chr0][c+1]*binSize;
    }
  }
  fout << endl;

  return chunkBounds;
}

int main(int argc, char *argv[]){

  if (argc != 5 && argc != 6) {
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "- arg1 = counts.bin.bgz file\n");
    fprintf(stderr, "- arg2 = binary file containing precomputed GC content and masks\n");
    fprintf(stderr, "- arg3 = unmasked bp per chunk (e.g., 1e5 or 1e6 for 100kb or 1Mb)\n");
    fprintf(stderr, "- arg4 = output prefix -> .{CNVs,profiles}.txt\n");
    fprintf(stderr, "- (optional) arg5 = bed file of bins to mask\n");
    exit(1);
  }

  const char *countsBinFile = argv[1];
  const char *GCmasksFile = argv[2];
  double unmaskedBpPerChunk; sscanf(argv[3], "%lf", &unmaskedBpPerChunk);
  const string outPrefix(argv[4]);
  const char *maskBinsFile = argc==5 ? NULL : argv[5];

  // load precomputed GC content and masks
  cout << "Loading precomputed GC content and masks: " << GCmasksFile << endl;
  int shift; // bit shift
  vector <int> chrLens(NUM_CHROMS); // chromosome lengths
  vector < vector <uint16_t> > GCmasks(NUM_CHROMS); // GCmasks vectors
  readGCmasks(GCmasksFile, shift, chrLens, GCmasks);
  
  cout << "Opening output files: " << outPrefix << ".{CNVs,profiles}.txt" << endl;
  FILE *foutCNVs = fopen((outPrefix + ".CNVs.txt").c_str(), "w"); assert(foutCNVs != NULL);
  ofstream foutProfiles((outPrefix + ".profiles.txt").c_str()); assert(foutProfiles);
  foutProfiles << fixed << setprecision(3);

  // for each chr, for each bin, compute number of positions with each GC percent
  vector < vector <GCpctCoeffs> > coeffs;
  vector < vector <int> > unmaskedBpPerBin;
  vector < vector <int> > chunkBounds; // bin [start,end) indices for chunks of ~unmaskedBpPerChunk
  uint16_t masksPrev = (uint16_t) -1;

  // load read counts per bin
  cout << "Processing read counts: " << countsBinFile << endl;
  BGZF *finCounts = bgzf_open(countsBinFile, "r"); assert(finCounts!=NULL);
  int IDlen;
  vector <float> GCprofile(MAX_GC_PCT+1);
  vector <uint16_t> readsPerBin; // buffer to store unzipped read counts
  vector <float> expPerBin;
  vector <float> ratioPerBin;
  vector <State> dpTable;
  vector <char> cnInBin;
  while (bgzf_read(finCounts, &IDlen, sizeof(IDlen)) == sizeof(IDlen)) {
    char IDstr[IDlen+1]; assert(bgzf_read(finCounts, IDstr, IDlen) > 0); IDstr[IDlen] = '\0';
    uint16_t masks; assert(bgzf_read(finCounts, &masks, sizeof(masks)) > 0); // masks used
    assert(bgzf_read(finCounts, &GCprofile[0], sizeof(GCprofile[0])*GCprofile.size()) > 0);
    int binSize; assert(bgzf_read(finCounts, &binSize, sizeof(binSize)) > 0);
    if (masks != masksPrev) { // compute coeffs and define chunks once for first ID
      assert(masksPrev == (uint16_t) -1);
      masksPrev = masks;
      coeffs = computeCoeffs(shift, chrLens, GCmasks, binSize, masks, maskBinsFile);
      unmaskedBpPerBin = computeUnmaskedBpPerBin(coeffs);
      chunkBounds = computeChunksWriteHeader(foutProfiles, unmaskedBpPerBin, binSize,
					     unmaskedBpPerChunk);
    }
    foutProfiles << IDstr;
    for (int chr0 = 0; chr0 < NUM_CHROMS; chr0++) {
      int numBins; assert(bgzf_read(finCounts, &numBins, sizeof(numBins)) > 0);
      readsPerBin.resize(numBins);
      assert(bgzf_read(finCounts, &readsPerBin[0], numBins * sizeof(readsPerBin[0])) > 0);
      expPerBin.resize(numBins);
      memset(&expPerBin[0], 0, numBins*sizeof(expPerBin[0]));
      ratioPerBin.resize(numBins);
      dpTable.resize(numBins+2);
      dpTable[0].logP[0] = dpTable[0].logP[2] = -1e9;
      cnInBin.resize(numBins);
      memset(&cnInBin[0], 0, numBins*sizeof(cnInBin[0]));

      int binsForMedian = 0;
      const int PAR1_LEN = 2800000; // rough length of PAR1 to mask
      for (int b = 0; b < numBins; b++) {
	for (int gc = 0; gc <= MAX_GC_PCT; gc++)
	  expPerBin[b] += coeffs[chr0][b].c[gc] * GCprofile[gc];
	if (expPerBin[b] > 0 && (chr0+1 <= 22 || b*binSize > PAR1_LEN)) // skip PAR1
	  ratioPerBin[binsForMedian++] = readsPerBin[b] / expPerBin[b];
      }
      nth_element(ratioPerBin.begin(), ratioPerBin.begin() + binsForMedian/2,
		  ratioPerBin.begin() + binsForMedian);
      double medianRatio = ratioPerBin[binsForMedian/2];
      int ploidy;
      if (chr0+1 <= 22)
	ploidy = 2;
      else if (chr0+1 == 23) { // chrX
	if (medianRatio < 0.75) // probably male (ok if expanded mLOX; just won't call many CNVs)
	  ploidy = 1;
	else
	  ploidy = 2;
      }
      else { // chrY
	if (medianRatio < 0.1) // probably female (ok if expanded mLOY; just won't call CNVs)
	  ploidy = 0;
	else
	  ploidy = 1;
      }
      fprintf(foutCNVs, "%s\tchr", IDstr);
      if (chr0<22) fprintf(foutCNVs, "%d", chr0+1); else fprintf(foutCNVs, "%c", 'X'+chr0-22);
      fprintf(foutCNVs, "\tMEDIAN\t%.3f\tPLOIDY\t%d\n", medianRatio, ploidy);

      const double logPjump = log(0.001); // TODO: run using a few and take union?
      for (int b = 0; b <= numBins; b++) {
	for (int CN = 1; CN <= 3; CN++) { // 1=DEL, 2=noCNV, 3=DUP
	  double logPemit = 0;
	  if (b < numBins && expPerBin[b] > 0) {
	    double expGivenCN;
	    if (ploidy == 2) expGivenCN = expPerBin[b] * CN / 2.0 * medianRatio;
	    else if (ploidy == 1 && b*binSize > PAR1_LEN)
	      expGivenCN = expPerBin[b] * max((CN-1.0), 0.1) * medianRatio;
	    else expGivenCN = 1; // don't call CNVs (by setting expGivenCN to a constant)
	    logPemit = -0.5 * (log(expGivenCN) + sq(readsPerBin[b]-expGivenCN) / expGivenCN);
	  }
	  dpTable[b+1].logP[CN-1] = -1e9;
	  for (int prevCN = 1; prevCN <= 3; prevCN++) {
	    double logPtot = dpTable[b].logP[prevCN-1] + abs(CN-prevCN)*logPjump + logPemit;
	    if (dpTable[b+1].logP[CN-1] < logPtot) {
	      dpTable[b+1].logP[CN-1] = logPtot;
	      dpTable[b+1].prev[CN-1] = prevCN-1;
	    }
	  }
	}
      }

      { // backtrack
	int bTable = numBins+1, CN0 = 2-1;
	CN0 = dpTable[bTable--].prev[CN0]; // backtrack from CN=2 state at right end of HMM
	while (bTable > 0) {
	  cnInBin[bTable-1] = CN0+1;
	  CN0 = dpTable[bTable--].prev[CN0];
	}
      }

      // record and output DELs and DUPs
      for (int bStart = 0; bStart < numBins; bStart++)
	if (cnInBin[bStart] != 2) {
	  int bEnd = bStart+1;
	  while (bEnd < numBins && cnInBin[bEnd]==cnInBin[bEnd-1])
	    bEnd++;

	  double obsSum = 0, expSum = 0;
	  for (int b = bStart; b < bEnd; b++) {
	    obsSum += readsPerBin[b];
	    expSum += expPerBin[b];
	  }
	  
	  fprintf(foutCNVs, "%s\tchr", IDstr);
	  if (chr0<22) fprintf(foutCNVs, "%d", chr0+1); else fprintf(foutCNVs, "%c", 'X'+chr0-22);
	  fprintf(foutCNVs, "\t%d\t%d\t", bStart*binSize, bEnd*binSize);
	  if (cnInBin[bStart] == 1) fprintf(foutCNVs, "DEL"); else fprintf(foutCNVs, "DUP");
	  fprintf(foutCNVs, "\t%.3f\t%.0f\t%.1f\n", obsSum/expSum, obsSum, expSum);

	  bStart = bEnd-1;
	}

      // output read depth profile = per-chunk obs/exp (after masking HMM DUPs + DELs)
      for (int c = 0; c+1 < (int) chunkBounds[chr0].size(); c++) {
	double obsSum = 0, expSum = 0; int unmaskedBp = 0;
	for (int b = chunkBounds[chr0][c]; b < chunkBounds[chr0][c+1]; b++)
	  if (cnInBin[b] == 2) {
	    obsSum += readsPerBin[b];
	    expSum += expPerBin[b];
	    unmaskedBp += unmaskedBpPerBin[chr0][b];
	  }
	foutProfiles << "\t" << (unmaskedBp < 0.5*unmaskedBpPerChunk ? NAN : obsSum/expSum);
      }
    }
    foutProfiles << endl;
  }
  bgzf_close(finCounts);

  fclose(foutCNVs);
  foutProfiles.close();
  
  return 0;
}
