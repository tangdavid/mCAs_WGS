#ifndef _REFMASKS_HPP
#define _REFMASKS_HPP

#include <vector>
#include <string>
#include <cstdio>
#include <cstring>
#include <cassert>

typedef unsigned short uint16_t;
typedef unsigned int uint32_t;



#define NUM_CHROMS 24 // chr1-22,X,Y
#define GC_WINDOW 400

// GC_with_masks 16-bit int data type:
// - first 9 bits = GC content of GC_WINDOW bases ending at current 0-based pos
// - next 3 bits = indicators of mask status within 0-based segment [pos, pos+(1<<shift))
// - next 1 bit = indicator of mask status anywhere in GC_WINDOW bases ending at current pos
#define GC_CONTENT_BITMASK 0x1ff
enum mask_fields {
  GCMASK_AT_SEG = 0x200,
  LCMASK_AT_SEG = 0x400,
  SVMASK_AT_SEG = 0x800,
  CNVMASK_AT_SEG = 0x1000,
  REFMASK_IN_WINDOW = 0x2000
};
#define getGCpct(gcMasks) ((((gcMasks)&GC_CONTENT_BITMASK)+2)>>2)
#define MAX_GC_PCT (getGCpct(GC_WINDOW))

#define ALL_MASKS_EXCEPT_GC (LCMASK_AT_SEG|SVMASK_AT_SEG|CNVMASK_AT_SEG|REFMASK_IN_WINDOW)
#define ALL_MASKS (GCMASK_AT_SEG|LCMASK_AT_SEG|SVMASK_AT_SEG|CNVMASK_AT_SEG|REFMASK_IN_WINDOW)
#define noMasks(gcMasks) (!((gcMasks)&ALL_MASKS))

// TODO: move to cpp file
void readGCmasks(const char *GCmasksFile, int &shift, std::vector <int> &chrLens,
		 std::vector < std::vector <unsigned short> > &GCmasks) {

  FILE *finGCmasks = fopen(GCmasksFile, "rb"); assert(finGCmasks!=NULL);
  assert(fread(&shift, sizeof(shift), 1, finGCmasks) == 1);
  assert(fread(&chrLens[0], sizeof(chrLens[0]), NUM_CHROMS, finGCmasks) == NUM_CHROMS);
  for (int i = 0; i < NUM_CHROMS; i++) {
    int numSegs = (chrLens[i]+(1<<shift)-1)>>shift;
    GCmasks[i].resize(numSegs);
    assert((int) fread(&GCmasks[i][0], sizeof(GCmasks[i][0]), numSegs, finGCmasks) == numSegs);
  }
  assert(fgetc(finGCmasks)==EOF);
  fclose(finGCmasks);
}

// TODO: move to cpp file
int parseChrStrToChr0(const char *chrStr) {
  int chr0 = 0;
  if (strcmp(chrStr, "chrX")==0)
    chr0 = 23-1;
  else if (strcmp(chrStr, "chrY")==0)
    chr0 = 24-1;
  else {
    assert(strncmp(chrStr, "chr", 3)==0);
    sscanf(chrStr, "chr%d", &chr0); chr0--;
    assert(chr0 >= 0 && chr0 < 22);
  }
  return chr0;
}

// TODO: move to cpp file
std::string chr0toString(int chr0) {
  assert(chr0 >= 0 && chr0 < 24);  
  if (chr0+1 <= 22) {
    char buf[6];
    sprintf(buf, "chr%d", chr0+1);
    return buf;
  }
  else if (chr0+1 == 23)
    return "chrX";
  else
    return "chrY";
}

#endif
