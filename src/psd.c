#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void
calculatePSD (float complex *in, int num_samples,
              double sampling_rate, double *psd)
{
  double delta = 1. / sampling_rate;

  int i;
  for (i = 0; i < num_samples; i++)
  {
    psd[i] = 2 * delta / num_samples * (pow (creal (in[i]), 2) + pow (cimag (in[i]), 2));
  }
}
