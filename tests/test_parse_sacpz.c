#include <complex.h>
#include <stdio.h>
#include <stdlib.h>

#include "parse_sacpz.h"

int
main ()
{
  const char *sacpzfile = "SAC_PZs_TW_NACB_BHZ__2007.254.07.25.20.0000_99999.9999.24.60.60.99999";
  double complex *poles;
  double complex *zeros;
  double constant;
  int nzeros, npoles;
  int result = parse_sacpz (sacpzfile, &poles, &npoles, &zeros, &nzeros, &constant);

  int i;
  printf ("zeros\n");
  for (i = 0; i < nzeros; i++)
    printf ("%lf %lf\n", creal (zeros[i]), cimag (zeros[i]));
  printf ("poles\n");
  for (i = 0; i < npoles; i++)
    printf ("%lf %lf\n", creal (poles[i]), cimag (poles[i]));
  printf ("constant\n");
  printf ("%lf\n", constant);

  free (poles);
  free (zeros);

  return 0;
}
