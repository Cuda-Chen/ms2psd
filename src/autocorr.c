#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "liquid.h"

#include "autocorrelation.h"
#include "datatype.h"

void
autocorr_float (data_t *data, uint64_t totalSamples, data_t *rxx)
{
  int i, k;
  int N = (int)totalSamples;
  for (i = 0; i <= 2 * N - 2; i++)
  {
    int j  = i - N + 1;
    rxx[i] = 0;
    for (k = 0; k <= N - 1; k++)
    {
      if ((k - j < 0) || (k - j >= N))
      {
        rxx[i] += 0.0;
      }
      else
      {
        rxx[i] += data[k] * data[k - j];
      }
    }
  }

  FILE *temp = fopen ("autocorr_result.txt", "w");
  for (i = 0; i < 2 * N - 1; i++)
  {
    fprintf (temp, "%d %f\n", i, rxx[i]);
  }
  fclose (temp);
}
