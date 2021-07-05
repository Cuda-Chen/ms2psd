#include <complex.h>
#include <stdio.h>
#include <stdlib.h>

#define PI 3.14159265358979323846

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

int
get_freq_response (double complex *poles, int npoles,
                   double complex *zeros, int nzeros,
                   double constant, double sampling_rate, int data_samples,
                   double **freq, double **freq_response,
                   int flag)
{
  double min = 0.0f;
  double max = sampling_rate;
  int i;

  /* Prepare frequency index and data */
  *freq = (double *)malloc (sizeof (double) * data_samples);
  if (*freq == NULL)
  {
    fprintf (stderr, "Allocate frequency index failed\n");
    return -1;
  }
  range (*freq, min, max, data_samples);

  /* Allocate frequency response array */
  *freq_response = (double *)malloc (sizeof (double) * data_samples);
  if (*freq_response == NULL)
  {
    fprintf (stderr, "Allocate frequency response failed\n");
    return -1;
  }

  /* Calculate frequency response */
  for (i = 0; i < data_samples; i++)
  {
    double complex numerator   = 1.0f + 0.0f * I;
    double complex denominator = 1.0f + 0.0f * I;
    double omega               = 2 * PI * (*freq)[i];
    double complex i_omega     = 0.0f + omega * I;
    int j;

    /* Calculate zeros */
    for (j = 0; j < nzeros; j++)
      numerator = numerator * (i_omega - zeros[j]);
    /* Calculate poles */
    for (j = 0; j < npoles; j++)
      denominator = denominator * (i_omega - poles[j]);

    double complex amp = numerator / denominator;
    /* phase */
    /* amplitude */
    (*freq_response)[i] = cabs (amp);
  }

  return 0;
}
