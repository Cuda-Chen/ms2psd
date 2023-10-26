#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void
calculatePSD (float complex *in, int num_samples,
              double sampling_rate, float *psd)
{
  float delta = 1. / sampling_rate;

  int i;
  for (i = 0; i < num_samples; i++)
  {
    psd[i] = 2 * delta / num_samples * (pow (crealf (in[i]), 2) + pow (cimagf (in[i]), 2));
  }
}
