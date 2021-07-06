#include <complex.h>
#include <stdio.h>
#include <stdlib.h>

#include "parse_sacpz.h"

#define OUTPUT_FILENAME "response.py"

int
main ()
{
  const char *sacpzfile = "SAC_PZs_TW_NACB_BHZ__2007.254.07.25.20.0000_99999.9999.24.60.60.99999";
  double complex *poles;
  double complex *zeros;
  double constant;
  int nzeros, npoles;
  int result = parse_sacpz (sacpzfile, &poles, &npoles, &zeros, &nzeros, &constant);

  /* Output poles, zeros, and constant. */
  int i;
  printf ("zeros\n");
  for (i = 0; i < nzeros; i++)
    printf ("%lf %lf\n", creal (zeros[i]), cimag (zeros[i]));
  printf ("poles\n");
  for (i = 0; i < npoles; i++)
    printf ("%lf %lf\n", creal (poles[i]), cimag (poles[i]));
  printf ("constant\n");
  printf ("%lf\n", constant);

  /* Plot frequency and phase response of this instrument. */
  FILE *fp = fopen (OUTPUT_FILENAME, "w");
  fprintf (fp, "from scipy import signal\n");
  fprintf (fp, "import matplotlib.pyplot as plt\n");
  fprintf (fp, "\n");

  fprintf (fp, "z = [");
  for (i = 0; i < nzeros; i++)
  {
    fprintf (fp, "%lf+%lfj", creal (zeros[i]), cimag (zeros[i]));
    if (i != nzeros - 1)
      fprintf (fp, ",");
  }
  fprintf (fp, "]\n");
  fprintf (fp, "p = [");
  for (i = 0; i < npoles; i++)
  {
    fprintf (fp, "%lf+%lfj", creal (poles[i]), cimag (poles[i]));
    if (i != npoles - 1)
      fprintf (fp, ",");
  }
  fprintf (fp, "]\n");
  fprintf (fp, "k = %lf\n", constant);

  fclose (fp);

  free (poles);
  free (zeros);

  return 0;
}
