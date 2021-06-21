#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "datatype.h"
#include "detrend.h"

int
main ()
{
  int npoints   = 1000;
  data_t *input = (data_t *)malloc (sizeof (data_t) * npoints);
  data_t *output;
  srand (42);

  int i;
  for (i = 0; i < npoints; i++)
    input[i] = rand ();

  int result = detrend (input, npoints, &output);
  assert (result == 0);

  free (input);
  free (output);
}
