#include <complex.h>
#include <stdint.h>

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
                                LIQUID_WINDOW_RCOSTAPER,
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
