#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "liquid.h"

void
spgram (float complex *data, uint64_t totalSamples, int nfft, float *psd)
{
  // options
  int windowType          = LIQUID_WINDOW_RCOSTAPER;
  unsigned int windowSize = 64;
  unsigned int delay      = (uint64_t)totalSamples;

  // create our spectral periodogram
  spgramcf q = spgramcf_create (nfft,
                                windowType,
                                windowSize,
                                delay);
#ifdef DEBUG
  spgramcf_print (q);
#endif

  unsigned int i;

  for (i = 0; i < totalSamples; i++)
  {
    spgramcf_push (q, data[i]);
  }

  spgramcf_get_psd (q, psd);

  spgramcf_destroy (q);
}

int
spectra (double complex *input, uint64_t totalSamples, double **amp, double **phase)
{
  // allocate amplitude and phase array
  *amp = (double *)malloc (sizeof (double) * totalSamples);
  if (amp == NULL)
  {
    fprintf (stderr, "Allocate amplitude array failed.\n");
    return -1;
  }
  *phase = (double *)malloc (sizeof (double) * totalSamples);
  if (phase == NULL)
  {
    fprintf (stderr, "Allocate phase array failed.\n");
    return -1;
  }

  uint64_t i;
  for (i = 0; i < totalSamples; i++)
  {
    *(*(amp) + i)   = sqrt (pow (creal (input[i]), 2) + pow (cimag (input[i]), 2));
    *(*(phase) + i) = atan2 (cimag (input[i]), creal (input[i]));
  }

  return 0;
}
