#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "common.h"
#include "datatype.h"

static void
usage ()
{
}

int
main (int argc, char **argv)
{
  char *mseedfile = NULL;
  data_t *data    = NULL;
  double sampleRate;
  uint64_t totalSamples;
  //double complex *filterResult;
  //double complex *freqResponse;
  float complex *filterResult;
  float complex *freqResponse;
  int nfft = 2000;
  float lowcut, highcut; /* low and high cutoff frequencies */
  int order;
  int passes;
  char *outputFile;
  const char *outputScript = "filter_result.m";
  int rv;
}
