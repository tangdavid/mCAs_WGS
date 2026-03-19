// compiling on RAP (to get libcurl + dependencies):
// sudo apt --yes install libdeflate-dev # (to get libdeflate headers)
// tar -zxf /mnt/project/lohdata/resources/htslib/htslib-1.18-libdeflate.tar.gz -C .
// g++ -O3 -Wall vcf_extract_common_AC_AN.cpp -o vcf_extract_common_AC_AN -Ihtslib-1.18-libdeflate/ htslib-1.18-libdeflate/libhts.a -lcurl -ldeflate -lz -lpthread

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <sys/time.h>
#include <htslib/bgzf.h>

#define BUF_LEN 0x10000 // buffer size

int main(int argc, char *argv[]) {
  assert(argc==2);

  const double min_MAF = 0.01;

  struct timeval t0, t1;
  gettimeofday(&t0, NULL);

  BGZF *fp = bgzf_open(argv[1], "r");
  assert(fp != NULL);
  /*
  const int threads = 2;
  bgzf_mt(fp, threads, 256); // adding a thread doesn't seem to help
  */
  char buf[BUF_LEN+1]; buf[BUF_LEN] = '\0';
  int bytes_read, buf_used = 0;
  long long bytes_tot = 0; int ctr_header_lines = 0, ctr_variants = 0, ctr_common = 0;
  while (true) {

    // fill up buffer
    bytes_read = bgzf_read(fp, buf+buf_used, BUF_LEN-buf_used);
    assert(bytes_read >= 0); // bgzf_read() returns -1 on error
    buf_used += bytes_read;
    bytes_tot += bytes_read;

    if (buf_used == 0) {
      //assert(bgzf_peek(fp) == -1); // should be at end of file -- unnecessary to check if bytes_read>=0 already checked?
      break;
    }

    // decide whether line should be printed
    bool print = false;
    if (buf[0] == '#') { // if header line, set print flag
      ctr_header_lines++;
      print = true;
    }
    else { // not header line => check AC and AN
      ctr_variants++;
      if (ctr_variants % 10000 == 0) { // print progress '.'
	fprintf(stderr, "variants processed: %d\n", ctr_variants);
	fflush(stderr);
      }
      char *tab_ptr = buf-1;
      for (int t = 0; t < 7; t++) { // seek to 7th tab = position before 8th field (INFO)
	tab_ptr = (char *) memchr(tab_ptr+1, '\t', (buf+buf_used)-(tab_ptr+1));
	assert(tab_ptr != NULL);
      }
      char *AC_ptr = strstr(tab_ptr+1, ";AC="); // search for INFO/AC
      if (AC_ptr != NULL)
	AC_ptr += 4;
      else {
	assert(strncmp(tab_ptr+1, "AC=", 3)==0);
	AC_ptr = tab_ptr+4;
      }
      int AC; int chars_read;
      assert(sscanf(AC_ptr, "%d%n", &AC, &chars_read));
      char *AN_ptr = strstr(tab_ptr+1, ";AN="); // search for INFO/AN
      assert(AN_ptr != NULL);
      AN_ptr += 4;
      int AN;
      assert(sscanf(AN_ptr, "%d%n", &AN, &chars_read));
      print = AN*min_MAF <= AC && AC <= AN*(1-min_MAF);
      if (print)
	ctr_common++;
    }

    // consume characters (and print if desired) until next newline
    while (true) {
      char *newline_ptr = (char *) memchr(buf, '\n', buf_used);
      if (newline_ptr == NULL) { // newline not yet found
	if (print)
	  fwrite(buf, 1, buf_used, stdout);
	// refill buffer
	bytes_read = bgzf_read(fp, buf, BUF_LEN);
	assert(bytes_read != 0);
	buf_used = bytes_read;
	bytes_tot += bytes_read;
      }
      else { // newline found
	if (print)
	  fwrite(buf, 1, newline_ptr+1-buf, stdout);
	// move end of buffer to beginning
	memmove(buf, newline_ptr+1, buf_used - (newline_ptr+1-buf));
	buf_used = buf_used - (newline_ptr+1-buf);
	break;
      }
    }
  }

  bgzf_close(fp);

  gettimeofday(&t1, NULL);

  fprintf(stderr, "\n");
  fprintf(stderr, "Bytes decompressed: %Ld\n", bytes_tot);
  fprintf(stderr, "Header lines: %d\n", ctr_header_lines);
  fprintf(stderr, "Variants: %d\n", ctr_variants);
  fprintf(stderr, "Common variants extracted: %d\n", ctr_common);
  fprintf(stderr, "Elapsed time: %.3f sec\n", t1.tv_sec-t0.tv_sec + 1e-6*(t1.tv_usec-t0.tv_usec));

  return 0;
}
