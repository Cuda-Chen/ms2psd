#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "range.h"

void
range (double *array, double sampleRate, int totalSamples)
{
  double delta         = 1. / sampleRate;
  double totalDuration = delta * totalSamples;
  int i;

  for (i = 0; i < totalSamples; i++)
    array[i] = i / totalDuration;
}

int
setLeftAndRightFreq (double **leftFreqs, double **rightFreqs, int *freqLen, double sampleRate, int windowLength)
{
  double nyquiestFreq = sampleRate / 2.0f;
  double octaveStep   = pow (2, 0.125);
  int i;
  int count = 0;

  /* short period (high frequency) */
  double fh = nyquiestFreq;
  double ts = 1. / fh;

  /* long period (low frequency) */
  double fl;
  double tl = 1.5f * ts;

  double tr = (double)windowLength / 5.0f; /* longest resolvable preiod */

  /* Count the length of left/right frequency length */
  while ((ts <= tr) && (tl <= tr))
  {
    count += 1;
    ts *= octaveStep;
    tl *= octaveStep;
  }
  *freqLen = count;

  *leftFreqs  = (double *)malloc (sizeof (double) * count);
  *rightFreqs = (double *)malloc (sizeof (double) * count);
  if (((*leftFreqs) == NULL) || ((*rightFreqs) == NULL))
  {
    fprintf (stderr, "Memory allocation failed of leftFreqs or rightFreqs\n");
    return -1;
  }

  /* Set the left and right frequency of each octave */
  ts = 1. / fh;
  tl = 1.5f * ts;
  for (i = 0; i < count; i++)
  {
    fh               = 1. / ts;
    fl               = 1. / tl;
    (*leftFreqs)[i]  = fh;
    (*rightFreqs)[i] = fl;
    ts *= octaveStep;
    tl *= octaveStep;
  }

  return 0;
}
