#include <complex.h>
#include <stdio.h>
#include <stdlib.h>

#include "range.h"

#define PI 3.14159265358979323846

int
get_freq_response (double complex *poles, int npoles,
                   double complex *zeros, int nzeros,
                   double constant, double sampling_rate, int data_samples,
                   double **freq, double complex **freq_response,
                   int flag)
{
  double k = constant;
  int i;
  /* Check user need displacement, velocity, or acceleration response */
  if (flag == 1)
    nzeros -= 1;
  else if (flag == 2)
    nzeros -= 2;

  /* Prepare frequency index and data */
  *freq = (double *)malloc (sizeof (double) * data_samples);
  if (*freq == NULL)
  {
    fprintf (stderr, "Allocate frequency index failed\n");
    return -1;
  }
  //range (*freq, sampling_rate, data_samples);
  double delta         = 1. / sampling_rate;
  double totalDuration = delta * data_samples;
  for (i = 0; i < data_samples / 2 + 1; i++)
  {
    (*freq)[i] = i / totalDuration;
  }
  for (i = data_samples / 2 + 1; i < data_samples; i++)
  {
    (*freq)[i] = -(data_samples - i) / totalDuration;
  }

  /* Allocate frequency response array */
  *freq_response = (double complex *)malloc (sizeof (double complex) * data_samples);
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

    (*freq_response)[i] = k * numerator / denominator;
    /* phase */
    /* amplitude */
  }

  return 0;
}

int
remove_response (double complex *data, double complex *freq_response, int data_samples)
{
  int i;
  data[0] = 0.0f + 0.0f * I;
  for (i = 1; i < data_samples; i++)
  {
    data[i] /= freq_response[i];
  }

  return 0;
}
