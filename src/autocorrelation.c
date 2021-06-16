#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "liquid.h"

#include "autocorrelation.h"
#include "datatype.h"

void
autocorrelation_float (data_t *data, uint64_t totalSamples, data_t *autoCorrelationResult)
{
  // options
  unsigned int sequenceLen = (unsigned int)totalSamples; // may overflow in the future
  unsigned int windowSize  = sequenceLen;
  unsigned int delay       = 8;
  int normalizeByEnergy    = 0; // normalize output by E{x^2}?

  // create autocorr object
  autocorr_rrrf q = autocorr_rrrf_create (windowSize, delay);
#ifdef DEBUG
  autocorr_rrrf_print (q);
#endif

  // Need padding to signal first?

  unsigned int i = 0;
  for (i = 0; i < totalSamples; i++)
  {
    autocorr_rrrf_push (q, data[i]);
    autocorr_rrrf_execute (q, &autoCorrelationResult[i]);

    if (normalizeByEnergy)
    {
      autoCorrelationResult[i] /= autocorr_rrrf_get_energy (q);
    }
  }

  FILE *temp = fopen ("autocorr_result.txt", "w");
  for (i = 0; i < totalSamples; i++)
  {
    fprintf (temp, "%d %f\n", i, autoCorrelationResult[i]);
  }
  fclose (temp);

  autocorr_rrrf_destroy (q);
}
