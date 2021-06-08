#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "fftw3.h"

#include "fft.h"

/*
 * Create an array of evenly spaced numbers.
 * Input: malloced array `array`, min, max, n
 * Output: array `array` filled with values
 * Reference: https://gist.github.com/mortenpi/f20a93c8ed3ee7785e65
 */
static void
range (double *array, double min, double max, size_t n)
{
  double delta = (max - min) / (double)(n - 1);
  size_t i;
  for (i = 0; i < n; i++)
  {
    array[i] = min + i * delta;
  }
}

void
fftToFileHalf (double *data, uint64_t dataSamples, double sampleRate, FILE *fptr)
{
  fftw_plan fft;
  uint64_t i;
  double min = 0.0, max = sampleRate;
  uint64_t mid = (dataSamples % 2 == 0) ? (dataSamples / 2) : (dataSamples / 2 + 1);
  /* allocate memory */
  fftw_complex *in  = (fftw_complex *)fftw_malloc (sizeof (fftw_complex) * dataSamples);
  fftw_complex *out = (fftw_complex *)fftw_malloc (sizeof (fftw_complex) * dataSamples);
  double *freq      = (double *)malloc (sizeof (double) * dataSamples);

  /* prepare input data */
  for (i = 0; i < dataSamples; i++)
  {
    in[i][0] = data[i];
    in[i][1] = 0;
  }

  /* prepare frequency index and data */
  range (freq, min, max, dataSamples);

  /* Fourier transform and save result to `out` */
  fft = fftw_plan_dft_1d (dataSamples, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute (fft);
  fftw_destroy_plan (fft);
  for (i = 0; i < mid; i++)
  {
    fprintf (fptr, "%lf %lf %lf\n", freq[i], out[i][0], out[i][1]);
  }

  /* free allocated memory */
  fftw_free (in);
  fftw_free (out);
  free (freq);
}

void
fftToFile (double *data, uint64_t dataSamples, double sampleRate, FILE *fptr)
{
  fftw_plan fft;
  uint64_t i;
  double min = 0, max = sampleRate;
  /* allocate memory */
  fftw_complex *in  = (fftw_complex *)fftw_malloc (sizeof (fftw_complex) * dataSamples);
  fftw_complex *out = (fftw_complex *)fftw_malloc (sizeof (fftw_complex) * dataSamples);
  double *freq      = (double *)malloc (sizeof (double) * dataSamples);

  /* prepare input data */
  for (i = 0; i < dataSamples; i++)
  {
    in[i][0] = data[i];
    in[i][1] = 0;
  }

  /* prepare frequency index and data */
  range (freq, min, max, dataSamples);

  /* Fourier transform and save result to `out` */
  fft = fftw_plan_dft_1d (dataSamples, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute (fft);
  fftw_destroy_plan (fft);

  for (i = 0; i < dataSamples; i++)
  {
    fprintf (fptr, "%lf %lf %lf\n", freq[i], out[i][0], out[i][1]);
  }

  /* free allocated memory */
  fftw_free (in);
  fftw_free (out);
  free (freq);
}

/* Usage: Unit test function of FFT utility.
 * Return: None
 */
void
testFFT (double *data, fftw_complex *in, fftw_complex *out,
         fftw_complex *ref, uint64_t totalSamples)
{
  fftw_plan fft, ifft;
  uint64_t i;
  /* allocate memory */
  in  = (fftw_complex *)fftw_malloc (sizeof (fftw_complex) * totalSamples);
  out = (fftw_complex *)fftw_malloc (sizeof (fftw_complex) * totalSamples);
  ref = (fftw_complex *)fftw_malloc (sizeof (fftw_complex) * totalSamples);

  /* prepare input data */
  for (i = 0; i < totalSamples; i++)
  {
    in[i][0] = data[i];
    in[i][1] = 0;
  }

  /* Fourier transform and save result to `out` */
  fft = fftw_plan_dft_1d (totalSamples, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute (fft);
  fftw_destroy_plan (fft);

  /* inverse Fourier transform and save result to `ref` */
  printf ("\ninverse transform:\n");
  ifft = fftw_plan_dft_1d (totalSamples, out, ref, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute (ifft);
  /* normalize */
  for (i = 0; i < totalSamples; i++)
  {
    ref[i][0] *= 1. / totalSamples;
    ref[i][1] *= 1. / totalSamples;
  }
  for (i = 0; i < totalSamples; i++)
  {
    printf ("recover: %" PRId64 " %+9.5f %+9.5f I v.s. %+9.5f %+9.5f I\n",
            i, in[i][0], in[i][1], ref[i][0], ref[i][1]);
  }
  fftw_destroy_plan (ifft);

  /* cleanup plan */
  fftw_cleanup ();
}
