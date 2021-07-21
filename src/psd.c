#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int
calculatePSD (double complex *in, int num_samples,
              double sampling_rate, double **psd)
{
  double delta = 1. / sampling_rate;
  *psd         = (double *)malloc (sizeof (double) * num_samples);
  if (psd == NULL)
  {
    fprintf (stderr, "cannot allocate PSD space\n");
    return -1;
  }

  int i;
  for (i = 0; i < num_samples; i++)
  {
    (*psd)[i] = 2 * delta / num_samples * (pow (creal (in[i]), 2) + pow (cimag (in[i]), 2));
  }

  return 0;
}
