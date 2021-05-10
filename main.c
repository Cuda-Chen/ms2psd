#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "autocorrelation.h"
#include "bandpass_filter.h"
#include "common.h"
#include "datatype.h"
#include "output2Octave.h"
#include "parse_miniSEED.h"
#include "spgram.h"

static void
usage ()
{
  printf ("Usage: ./ms2psd [c1] [c2] [order] [passes] [input] [output]");
  printf ("\n\nInput parameters:\n");
  printf ("c1: low cut frequency (Hz)\n");
  printf ("c2: high cut frequency (Hz)\n");
  printf ("order: filter order\n");
  printf ("passes: set '1' for forward filtering or set '2' for forward-backward filtering\n");
  printf ("input: a miniSEED seismic record\n");
  printf ("output: a text file containing filtered result (real and imaginary part)\n");
  printf ("\nOutput files and their format: \n");
  printf ("1. A MATLAB script to plot the input and filter result.\n");
  //printf ("2. A text file containing filtered result (real and imaginary part).\n");
  //printf ("3. A miniSEED file which contains filtered result (real part only).\n");
}

int
main (int argc, char **argv)
{
  char *mseedfile = NULL;
  data_t *data    = NULL;
  double sampleRate;
  uint64_t totalSamples;
  data_t *autocorr;
  float complex *filterResult;
  float complex *freqResponse;
  int nfft = 2000;
  float lowcut, highcut; /* low and high cutoff frequencies */
  int order;
  int passes;
  char *outputFile;
  const char *outputScript = "psd_result.m";
  float *psd;
  int rv;

  /* Simple argement parsing */
  if (argc != 7)
  {
    usage ();
    return 1;
  }
  lowcut     = atof (argv[1]);
  highcut    = atof (argv[2]);
  order      = atoi (argv[3]);
  passes     = atoi (argv[4]);
  mseedfile  = argv[5];
  outputFile = argv[6];

  /* Get data from input miniSEED file */
  rv = parse_miniSEED (mseedfile, &data, &sampleRate, &totalSamples);
  if (rv != 0)
  {
    return rv;
  }
  if (data == NULL)
  {
    printf ("Input data read unsuccessfully\n");
    return -1;
  }

  /* Apply auto-correlation */
  autocorr = (data_t *)malloc (sizeof (data_t) * totalSamples);
  autocorrelation_float (data, totalSamples, autocorr);
  if (autocorr == NULL)
  {
    printf ("Auto-correlation failed\n");
    return -1;
  }

  /* Filter the data with band pass filter */
  filterResult = (float complex *)malloc (sizeof (float complex) * totalSamples);
  freqResponse = (float complex *)malloc (sizeof (float complex) * nfft);
  bandpass_filter_float (autocorr, sampleRate, totalSamples, nfft,
                         lowcut, highcut, order, passes,
                         filterResult, freqResponse);
  if (filterResult == NULL || freqResponse == NULL)
  {
    printf ("Something wrong when applying band pass filter\n");
    return -1;
  }

  /* Calculate spectral periodogram */
  psd = (float *)malloc (sizeof (float) * nfft);
  spgram (filterResult, totalSamples, nfft, psd);
  if (psd == NULL)
  {
    printf ("spectral periodogram failed\n");
    return -1;
  }

  /* Save the filtered result to MATLAB script */
  output2Octave (outputScript, nfft, psd);

  /* Free alllocated objects */
  free (data);
  free (autocorr);
  free (filterResult);
  free (freqResponse);
  free (psd);

  return 0;
}
