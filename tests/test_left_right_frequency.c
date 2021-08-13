#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "range.h"

int
main ()
{
  double *left, *right;
  int freqLen;
  double sampleRate = 100.0f;
  int windowLength  = 900;

  int rv = setLeftAndRightFreq (&left, &right, &freqLen, sampleRate, windowLength);
  assert (rv == 0);

  for (int i = 0; i < freqLen; i++)
  {
    double ts = 1. / left[i];
    double tl = 1. / right[i];
    double tc = sqrt (ts * tl);
    printf ("%lf %lf %lf %lf %lf\n", left[i], right[i], ts, tl, tc);
  }

  free (left);
  free (right);

  return 0;
}
