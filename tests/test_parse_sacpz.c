#include <complex.h>
#include <stdio.h>

#include "parse_sacpz.h"

int
main ()
{
  const char *sacpzfile = "SAC_PZs_TW_NACB_BHZ__2007.254.07.25.20.0000_99999.9999.24.60.60.99999";
  double complex *poles;
  double complex *zeros;
  double constant;
  int result = parse_sacpz (sacpzfile, &poles, &zeros, &constant);

  return 0;
}
