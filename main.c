#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
//#include <assert.h>

#include "fftw3.h"

#include "autocorr.h"
#include "bandpass_filter.h"
#include "common.h"
#include "cosine_taper.h"
#include "datatype.h"
#include "detrend.h"
#include "fft.h"
#include "output2Octave.h"
#include "parse_miniSEED.h"
#include "spgram.h"
#include "parse_sacpz.h"
#include "freq_response.h"

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
#if 0
  data_t *autocorr;
  float complex *filterResult;
  float complex *freqResponse;
  int nfft = 2000;
#endif
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
  /*autocorr = (data_t *)malloc (sizeof (data_t) * totalSamples);
  autocorrelation_float (data, totalSamples, autocorr);
  if (autocorr == NULL)
  {
    printf ("Auto-correlation failed\n");
    return -1;
  }*/
#if 0
  autocorr = (data_t *)malloc (sizeof (data_t) * totalSamples * 2 - 1);
  autocorr_float (data, totalSamples, autocorr);
  if (autocorr == NULL)
  {
    printf ("Auto-correlation failed\n");
    return -1;
  }

  /* Filter the data with band pass filter */
  filterResult = (float complex *)malloc (sizeof (float complex) * totalSamples * 2 - 1);
  freqResponse = (float complex *)malloc (sizeof (float complex) * nfft);
  bandpass_filter_float (autocorr, sampleRate, totalSamples * 2 - 1, nfft,
                         lowcut, highcut, order, passes,
                         filterResult, freqResponse);
  if (filterResult == NULL || freqResponse == NULL)
  {
    printf ("Something wrong when applying band pass filter\n");
    return -1;
  }

  /* Calculate spectral periodogram */
  psd = (float *)malloc (sizeof (float) * nfft);
  spgram (filterResult, totalSamples * 2 - 1, nfft, psd);
  if (psd == NULL)
  {
    printf ("spectral periodogram failed\n");
    return -1;
  }

  /* Save the filtered result to MATLAB script */
  output2Octave (outputScript, nfft, psd);
#endif

  /* Method #2 */
  FILE *dataout = fopen ("dataout.txt", "w");
  for (int i = 0; i < totalSamples; i++)
  {
    fprintf (dataout, "%f%c", data[i], " \n"[i != totalSamples - 1]);
  }
  fclose (dataout);
  /* Detrend */
  data_t *detrended;
  detrend (data, (int)totalSamples, &detrended);
  FILE *detrendout = fopen ("detrendout.txt", "w");
  for (int i = 0; i < totalSamples; i++)
  {
    fprintf (detrendout, "%f%c", detrended[i], " \n"[i != totalSamples - 1]);
  }
  fclose (detrendout);
  /* First taper the signal with 5%-cosine-window */
  float *taperedSignal = (float *)malloc (sizeof (float) * totalSamples);
  cosineTaper (detrended, (int)totalSamples, 0.05, taperedSignal);
  /* Then execute FFT */
  double *tapered = (double *)malloc (sizeof (double) * totalSamples);
  for (int i = 0; i < totalSamples; i++)
  {
    tapered[i] = (double)taperedSignal[i];
  }
  double complex *fftResult;
  fft(tapered, totalSamples, &fftResult);
  /* Apply instrument response removal */
  /* Read SACPZ file */
  const char *sacpzfile = "foo";
  double complex *poles, *zeros;
  double npoles, nzeros;
  double constant;
  parse_sacpz(sacpzfile, 
              &poles, &npoles,
              &zeros, &nzeros,
              &constant);
  /* Get frequency response */
  double *freq;
  double complex *freqResponse;
  int flag = 0;
  get_freq_response(poles, npoles,
                    zeros, nzeros,
                    constant, sampleRate, totalSamples,
                    &freq, &freqResponse, flag);
  remove_response(fftResult, freqResponse, totalSamples);
  /* Apply PSD */

#if 0
  FILE *fftResult = fopen ("fft_result.txt", "w");
  if (fftResult == NULL)
  {
    printf ("Error opening!");
    return -1;
  }
  fftToFileHalf (tapered, totalSamples, sampleRate, fftResult);
  fclose (fftResult);

  double complex *fftout;
  fft (tapered, totalSamples, &fftout);
  double *amp, *phase;
  spectra (fftout, totalSamples, &amp, &phase);
  FILE *fftOut = fopen ("fft_out.txt", "w");
  for (int i = 0; i < totalSamples; i++)
  {
    fprintf (fftOut, "%lf %lf%c", amp[i], phase[i], " \n"[i != totalSamples - 1]);
  }
  fclose (fftOut);
#endif

  /* Free allocated objects */
  free (data);
#if 0
  free (autocorr);
  free (filterResult);
  free (freqResponse);
  free (psd);
#endif
  free (detrended);
  free (taperedSignal);
  free( tapered);
  free(fftResult);
  free(poles);
  free(zeros);
  free(freq);
  free(freqResponse);
#if 0
  free (fftout);
  free (amp);
  free (phase);
#endif

  return 0;
}
