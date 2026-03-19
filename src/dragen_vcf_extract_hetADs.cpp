// compiling on RAP (to get libcurl + dependencies):
// sudo apt --yes install libdeflate-dev # (to get libdeflate headers)
// tar -zxf /mnt/project/lohdata/resources/htslib/htslib-1.18-libdeflate.tar.gz -C .
// g++ -O3 -Wall dragen_vcf_extract_hetADs.cpp -o ../bin/dragen_vcf_extract_hetADs -Ihtslib-1.18-libdeflate/ htslib-1.18-libdeflate/libhts.a -lcurl -ldeflate -lz -lpthread

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <sys/time.h>
#include <htslib/bgzf.h>

#define BUF_LEN 0x20000 // buffer size

int main(int argc, char *argv[]) {
  assert(argc==2);

  struct timeval t0, t1;
  gettimeofday(&t0, NULL);

  BGZF *fp = bgzf_open(argv[1], "r");
  if (fp == NULL) {
    fprintf(stderr, "ERROR: failed to open file: %s\n", argv[1]);
    exit(1);
  }

  char buf[BUF_LEN]; int buf_used = 0;
  long long bytes_tot = 0; int ctr_variants = 0, ctr_hets = 0;
  while (true) {
    int bytes_read = bgzf_read(fp, buf+buf_used, BUF_LEN-buf_used); // fill up buffer
    assert(bytes_read >= 0); // bgzf_read() returns -1 on error
    bool last_buf = false;
    if (bytes_read < BUF_LEN-buf_used) {
      //assert(bgzf_peek(fp) == -1); // check that EOF was reached if buffer didn't fill up -- unnecessary if bytes_read>=0 already checked?
      last_buf = true;
    }

    buf_used += bytes_read;
    bytes_tot += bytes_read;
    
    char *newline_ptr = buf-1, *newline_ptr_next;
    while ((newline_ptr_next = (char *) memchr(newline_ptr+1, '\n', buf_used-(newline_ptr+1-buf)))
	   != NULL) {

      // process line: [newline_ptr+1, newline_ptr_next)
      if (*(newline_ptr+1) != '#') { // not header line => check if GT is het
	ctr_variants++;
	char *tab_ptr_GT = (char *) memrchr(newline_ptr+1, '\t', newline_ptr_next-(newline_ptr+1));
	assert(tab_ptr_GT != NULL); // check that tab character was found
	if ((*(tab_ptr_GT+1)-'0') + (*(tab_ptr_GT+3)-'0') == 1) { // 0/1, 0|1, 1|0
	  ctr_hets++;
	  char *tab_ptr_FMT = (char *) memrchr(newline_ptr+1, '\t', tab_ptr_GT-(newline_ptr+1));
	  assert(tab_ptr_FMT != NULL); // check that tab character was found

	  *tab_ptr_GT = '\0'; // null-terminate end of FMT string
	  *newline_ptr_next = '\0'; // null-terminate end of GT string
	  char *token_FMT, *token_data;
	  char *remaining_FMT = tab_ptr_FMT+1, *remaining_data = tab_ptr_GT+1;
	  int AD0 = 0, AD1 = 0, GQ = 0;
	  while ((token_FMT = strtok_r(remaining_FMT, ":", &remaining_FMT)) &&
		 (token_data = strtok_r(remaining_data, ":", &remaining_data))) {
	    if (strcmp(token_FMT, "AD") == 0)
	      assert(sscanf(token_data, "%d,%d", &AD0, &AD1) == 2);
	    if (strcmp(token_FMT, "GQ") == 0)
	      assert(sscanf(token_data, "%d", &GQ) == 1);
	  }

	  if (GQ >= 20) {
	    char *tab_ptr_POS = (char *) memchr(newline_ptr+1, '\t', tab_ptr_FMT-(newline_ptr+1));
	    assert(tab_ptr_POS != NULL); // check that tab character was found
	    char *tab_ptr_ID = (char *) memchr(tab_ptr_POS+1, '\t', tab_ptr_FMT-(tab_ptr_POS+1));
	    assert(tab_ptr_ID != NULL); // check that tab character was found
	    char *tab_ptr_REF = (char *) memchr(tab_ptr_ID+1, '\t', tab_ptr_FMT-(tab_ptr_ID+1));
	    assert(tab_ptr_REF != NULL); // check that tab character was found
	    char *tab_ptr_ALT = (char *) memchr(tab_ptr_REF+1, '\t', tab_ptr_FMT-(tab_ptr_REF+1));
	    assert(tab_ptr_ALT != NULL); // check that tab character was found
	    char *tab_ptr_QUAL = (char *) memchr(tab_ptr_ALT+1, '\t', tab_ptr_FMT-(tab_ptr_ALT+1));
	    assert(tab_ptr_QUAL != NULL); // check that tab character was found
	    
	    fwrite(newline_ptr+1, 1, tab_ptr_ID - newline_ptr, stdout);
	    fwrite(tab_ptr_REF+1, 1, tab_ptr_QUAL - tab_ptr_REF, stdout);
	    printf("%d\t%d\n", AD0, AD1);
	  }
	}
      }
      newline_ptr = newline_ptr_next;
    }
    
    if (last_buf) {
      assert(newline_ptr == buf+buf_used-1); // check that file ended with newline
      break;
    }
    
    if (newline_ptr == buf-1) { // exit if newline was not found within buffer
      fprintf(stderr, "ERROR: line longer than %d\n", BUF_LEN);
      exit(1);
    }

    // move end of buffer to beginning
    memmove(buf, newline_ptr+1, buf_used - (newline_ptr+1-buf));
    buf_used = buf_used - (newline_ptr+1-buf);
  }

  bgzf_close(fp);

  gettimeofday(&t1, NULL);

  fprintf(stderr, "Decompressed %Ld bytes (%d variants, %d 0/1-hets) in %.3f sec\n",
	  bytes_tot, ctr_variants, ctr_hets, t1.tv_sec-t0.tv_sec + 1e-6*(t1.tv_usec-t0.tv_usec));
  
  return 0;
}
