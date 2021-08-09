#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "range.h"

int
main ()
{
  double *left, *right;
  int freqLen;
  double sampleRate = 100.0f;
  int windowLength  = 90000;

  int rv = setLeftAndRightFreq (&left, &right, &freqLen, sampleRate, windowLength);
  assert (rv == 0);

  for (int i = 0; i < freqLen; i++)
  {
    printf ("%lf %lf %lf %lf \n", left[i], right[i], 1. / left[i], 1. / right[i]);
  }

  free (left);
  free (right);

  return 0;
}
