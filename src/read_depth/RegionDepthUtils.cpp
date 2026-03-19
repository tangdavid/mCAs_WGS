// TODO: create RegionDepthUtils.hpp; compile + link in Makefile

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <utility>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <htslib/kstring.h>
#include <htslib/bgzf.h>

using namespace std;

inline double crop(double x, double delta) {
  return x > delta ? delta : (x < -delta ? -delta : x);
}

struct ProfileRegion {
  int chr0, start, end;
  ProfileRegion(int _chr0, int _start, int _end) : chr0(_chr0), start(_start), end(_end) {}
};

struct RegionsPerChr {
  vector < pair <int, int> > regions[NUM_CHROMS];
};

double medianNonMissing(const vector <double> &vec, int &ctrNonMissing) {
  vector <double> vecNonMissing;
  for (int i = 0; i < (int) vec.size(); i++)
    if (!isnan(vec[i]))
      vecNonMissing.push_back(vec[i]);
  ctrNonMissing = vecNonMissing.size();
  if (ctrNonMissing == 0) return 0;
  sort(vecNonMissing.begin(), vecNonMissing.end());
  return 0.5 * (vecNonMissing[(ctrNonMissing-1)/2] + vecNonMissing[ctrNonMissing/2]);
}

// read CNVs and determine sex
map <string, RegionsPerChr> readCNVsexData(map <string, bool> &IDisMale, const char *CNVsFile) {

  map <string, RegionsPerChr> CNVs;
  IDisMale.clear();

  BGZF *finCNVs = bgzf_open(CNVsFile, "r"); assert(finCNVs!=NULL);
  kstring_t kstr = KS_INITIALIZE;
  int ctrCNV = 0, ctrMale = 0;
  while (true) {
    int ret = bgzf_getline(finCNVs, '\n', &kstr);
    if (ret == -1) break;
    if (ret < -1) {
      cerr << "ERROR: bgzf_getline returned " << ret << " while reading CNVs" << endl;
      exit(1);
    }
    const char *kstr_c_str = ks_c_str(&kstr);
    ks_tokaux_t aux;
    char *tok;
    int token_ctr = 0;
    string IDstr; int chr0 = 0, CNVstart = 0, CNVend;
    bool isMedianChrXratio = false, isMedianChrYratio = false; double medianChrXratio = 0;
    for (tok = kstrtok(kstr_c_str, "\t", &aux); tok; tok = kstrtok(0, 0, &aux), token_ctr++) {
      *(char *) aux.p = 0;
      if (token_ctr == 0)
	IDstr = tok;
      else if (token_ctr == 1)
	chr0 = parseChrStrToChr0(tok);
      else if (token_ctr == 2) {
	if (strcmp(tok, "MEDIAN")==0) {
	  if (chr0 == 23-1)
	    isMedianChrXratio = true;
	  else if (chr0 == 24-1)
	    isMedianChrYratio = true;
	  else
	    break;
	}
	else
	  assert(sscanf(tok, "%d", &CNVstart) == 1);
      }
      else if (token_ctr == 3) {
	if (isMedianChrXratio)
	  assert(sscanf(tok, "%lf", &medianChrXratio) == 1);
	else if (isMedianChrYratio) {
	  double medianChrYratio;
	  assert(sscanf(tok, "%lf", &medianChrYratio) == 1);
	  bool isMale = medianChrYratio > 0.05 && medianChrXratio < 0.75;;
	  IDisMale[IDstr] = isMale;
	  ctrMale += isMale;
	}
	else {
	  assert(sscanf(tok, "%d", &CNVend) == 1);
	  CNVs[IDstr].regions[chr0].push_back(make_pair(CNVstart, CNVend));
	  ctrCNV++;
	}
	break;
      }
    }
  }
  ks_free(&kstr);
  bgzf_close(finCNVs);
  cerr << "Read " << ctrCNV << " CNV calls" << endl;
  cerr << "Male proportion: " << ctrMale << " / " << CNVs.size() << endl;

  return CNVs;
}

// read profile regions and baselines; populate chrBinToProfileInd if requested (binSize!=0)
// - chrBinToProfileInd is assumed to be initialized to the right size (per chr0) if not NULL
vector <ProfileRegion> readProfileRegionsBaselines(vector <int> &profileChrStartInds,
						   vector <double> baselines[2], // female, male
						   vector <int> chrBinToProfileInd[],
						   int binSize, // 0 to skip chrBinToProfileInd
						   const char *PCbaselinesPrefix) {
  vector <ProfileRegion> profileRegions;
  profileChrStartInds = vector <int> (NUM_CHROMS+1);

  int profileLen = 0;
  for (int sex = 0; sex < 2; sex++) {
    char baselinesFile[strlen(PCbaselinesPrefix)+50];
    sprintf(baselinesFile, "%s.isMale%d.baselines.txt", PCbaselinesPrefix, sex);
    cerr << "Reading per-region depth baselines: " << baselinesFile << endl;
    ifstream finBaselines(baselinesFile); assert(finBaselines);
    string header; getline(finBaselines, header);
    if (sex == 0) { // parse header = list of chunk boundaries: chr#_##_##
      istringstream iss(header);
      string chr_start_end;
      while (iss >> chr_start_end) {
	size_t underPos = chr_start_end.find('_');
	assert(underPos != string::npos);
	int chr0 = parseChrStrToChr0(chr_start_end.substr(0, underPos).c_str());
	int start, end;
	assert(sscanf(chr_start_end.substr(underPos+1).c_str(), "%d_%d", &start, &end) == 2);
	if (binSize) {
	  assert(start % binSize == 0); assert(end % binSize == 0);
	  for (int b = start/binSize; b < end/binSize; b++)
	    chrBinToProfileInd[chr0][b] = profileRegions.size();
	}
	profileRegions.push_back(ProfileRegion(chr0, start, end));
	profileChrStartInds[chr0+1]++;
      }
      profileLen = profileRegions.size();
      cerr << "Read " << profileLen << " genomic chunks in profile file" << endl;
      for (int chr0 = 0; chr0 < NUM_CHROMS; chr0++)
	profileChrStartInds[chr0+1] += profileChrStartInds[chr0];
    }
    // read baselines
    baselines[sex].resize(profileLen);
    for (int i = 0; i < profileLen; i++)
      assert(finBaselines >> baselines[sex][i]);
    double x; assert(!(finBaselines >> x)); // check eof reached
    finBaselines.close();
    cerr << "Read baseline depths for " << baselines[sex].size() << " genomic chunks" << endl;
  }
  assert(baselines[0].size() == baselines[1].size());

  return profileRegions;
}

// read genomic wave PCs
void readGenomicWavePCs(vector < vector <double> > PCs[2], const char *PCbaselinesPrefix,
			int numPCs, int profileLen) {

  for (int sex = 0; sex < 2; sex++) {
    char PCsFile[strlen(PCbaselinesPrefix)+50];
    sprintf(PCsFile, "%s.isMale%d.PCs.txt", PCbaselinesPrefix, sex);
    cerr << "Reading " << numPCs << " PCs from " << PCsFile << endl;
    ifstream finPCs(PCsFile); assert(finPCs);
    string line;
    getline(finPCs, line); // ignore header
    for (int p = 0; p < numPCs; p++) {
      getline(finPCs, line);
      vector <double> pc(profileLen);
      istringstream iss(line);
      for (int i = 0; i < profileLen; i++)
	assert(iss >> pc[i]);
      iss.peek();
      assert(iss.eof());
      PCs[sex].push_back(pc);
    }
    finPCs.close();
  }
}

// read line of depth profile file; return IDstr (or "" if ID is not requested => don't parse)
bool readProfileLine(string &IDstr, vector <double> &depthProfile, BGZF *finDepthProfile,
		     const map <string, RegionsPerChr> *IDregionsPtr, int isHeader,
		     kstring_t *kstr_ptr) {

  int profileLen = depthProfile.size();
  assert(profileLen); // assume depthProfile has been initialized to the correct size
  
  // read line of profile file
  int ret = bgzf_getline(finDepthProfile, '\n', kstr_ptr);
  if (ret == -1) return false;
  if (ret < -1) {
    cerr << "ERROR: bgzf_getline returned " << ret << " while reading profiles" << endl;
    exit(1);
  }
  const char *kstr_c_str = ks_c_str(kstr_ptr);

  if (isHeader) { // if header line, check first token is "ID" and skip rest of line
    assert(strncmp(kstr_c_str, "ID\t", 3) == 0);
    return true;
  }
      
  // parse ID and depth profile data
  // ID = first token
  ks_tokaux_t aux;
  char *tok = kstrtok(kstr_c_str, "\t", &aux);
  *(char *) aux.p = 0;
  IDstr = tok;
  if (IDregionsPtr != NULL && IDregionsPtr->find(IDstr) == IDregionsPtr->end()) {
    IDstr = "";
    return true; // no region requested for sample => ignore rest of line
  }

  // depth profile data = remaining tokens
  int profileInd = 0;
  while ((tok = kstrtok(0, 0, &aux))) {
    *(char *) aux.p = 0;
    assert(profileInd < profileLen);
    if (!sscanf(tok, "%lf", &depthProfile[profileInd]))
      depthProfile[profileInd] = NAN; // set to missing if token doesn't parse
    profileInd++;
  }
  assert(profileInd == profileLen);

  return true;
}

// compute PC coeffs and combine to estimate genomic wave (after applying crop!)
vector <double> computeGenomicWave(const vector <double> &depthProfile,
				   const vector < vector <double> > &PCs,
				   const vector <double> &baselines, double readDepthRelCrop) {

  int profileLen = depthProfile.size();
  vector <double> wave(profileLen);

  for (int p = 0; p < (int) PCs.size(); p++) {
    const vector <double> &pc = PCs[p];
    double coeff = 0;
    for (int i = 0; i < profileLen; i++)
      if (!isnan(depthProfile[i]))
	coeff += pc[i] * (crop(depthProfile[i]/baselines[i] - 1, readDepthRelCrop));
    for (int i = 0; i < profileLen; i++)
      wave[i] += coeff * pc[i];
  }

  return wave;
}

// compute PC-adjusted depth per chromosome (median across profile regions)
vector <double> computeMedianDepthRatioPerChrom(const vector <double> &depthProfile,
						const vector <double> &baselines,
						const vector <double> &wave,
						const vector <int> &profileChrStartInds, bool adjPC) {

  vector <double> medianRatioPerChrom(NUM_CHROMS);

  vector <double> ratioPerProfileRegion;
  for (int chr0 = 0; chr0 < NUM_CHROMS; chr0++) {
    ratioPerProfileRegion.resize(profileChrStartInds[chr0+1] - profileChrStartInds[chr0]);
    for (int i = profileChrStartInds[chr0]; i < profileChrStartInds[chr0+1]; i++)
      ratioPerProfileRegion[i-profileChrStartInds[chr0]] =
	depthProfile[i] / (baselines[i] * (1+wave[i]*adjPC));
    int ctrNonMissing;
    medianRatioPerChrom[chr0] = medianNonMissing(ratioPerProfileRegion, ctrNonMissing);
  }
  return medianRatioPerChrom;
}

// compute final calibration (median of median-per-autosome)
double computeMedianAutosomeObsExp(const vector <double> &depthProfile,
				   const vector <double> &baselines, const vector <double> &wave,
				   const vector <int> &profileChrStartInds, bool adjPC) {

  vector <double> medianRatioPerChrom =
    computeMedianDepthRatioPerChrom(depthProfile, baselines, wave, profileChrStartInds, adjPC);

  medianRatioPerChrom.resize(NUM_CHROMS-2); // restrict to autosomes
  int ctrNonMissing;
  return medianNonMissing(medianRatioPerChrom, ctrNonMissing);
}
