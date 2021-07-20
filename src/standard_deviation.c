#include <inttypes.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>

#include "standard_deviation.h"

double
calculateSD (double *data, uint64_t dataSize)
{
  double sum = 0.0, mean, SD = 0.0;
  uint64_t i;
  for (i = 0; i < dataSize; i++)
  {
    sum += data[i];
  }
  mean = sum / (double)dataSize;
  for (i = 0; i < dataSize; i++)
  {
    SD += pow (data[i] - mean, 2);
  }
  printf ("mean: %lf / %" PRId64 " = %lf\n", sum, dataSize, mean);
  return sqrt (SD / dataSize);
}

void
getMeanAndSD (float *data, uint64_t dataSize, float *_mean, float *_SD)
{
  float sum = 0.0, mean, SD = 0.0;
  uint64_t i;
  for (i = 0; i < dataSize; i++)
  {
    sum += data[i];
  }
  mean = sum / (float)dataSize;
  for (i = 0; i < dataSize; i++)
  {
    SD += pow (data[i] - mean, 2);
  }
#ifdef DEBUG
  printf ("mean: %lf / %" PRId64 " = %lf\n", sum, dataSize, mean);
#endif
  *_mean = round (mean * 100) / 100;
  *_SD   = round (sqrt (SD / dataSize) * 100) / 100;
}

void
testCalculateSD ()
{
  double test[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  printf ("Standard deviation of test array: %lf\n", calculateSD (test, 10));
}
