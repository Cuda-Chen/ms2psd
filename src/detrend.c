#include <stdio.h>
#include <stdlib.h>

#include "datatype.h"

int
detrend (data_t *input, int npts, data_t **output)
{
  /* allocate output array memory */
  *output = (data_t *)malloc (sizeof (data_t) * npts);
  if (*output == NULL)
  {
    fprintf (stderr, "Cannot allocate detrend output memory.");
    return -1;
  }

  int i;
  data_t sx  = 0.0,
         sy  = 0.0,
         sxy = 0.0,
         sx2 = 0.0;
  for (i = 0; i < npts; i++)
  {
    sx += i;
    sy += input[i];
    sxy += i * input[i];
    sx2 += i * i;
  }

  data_t slope    = (sx * sy - npts * sxy) / (sx * sy - npts * sx2);
  data_t constant = (sy - slope * sx) / (float)npts;

  for (i = 0; i < npts; i++)
  {
    *(*(output) + i) = input[i] - (constant + slope * i);
  }

  return 0;
}
