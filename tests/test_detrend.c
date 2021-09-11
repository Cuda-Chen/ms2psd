#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "datatype.h"
#include "detrend.h"

const data_t EPSILON = 1e-6;

int
main ()
{
  /*int npoints   = 1000;
  data_t *input = (data_t *)malloc (sizeof (data_t) * npoints);
  data_t *output = (data_t *)malloc(sizeof(data_t) * npoints);
  srand (42);

  int i;
  for (i = 0; i < npoints; i++)
    input[i] = rand ();*/

  const int npoints = 5;
  data_t input[npoints] = {1.1, 1.2, 1.3, 1.4, 1.5};
  data_t output[npoints];

  detrend (input, npoints, output);
  data_t mean = 0.0;
  for(int i = 0; i < npoints; i++)
    mean += output[i];
  mean /= (data_t)npoints;
  assert(mean < EPSILON);

  /*free (input);
  free (output);*/

  return 0;
}
