#include <complex.h>
#include <stdio.h>
#include <stdlib.h>

#define LEN 0x200

int
parse_sacpz (const char *sacpzfile, double complex **poles,
             double complex **zeros, double *constant)
{
  FILE *fp = fopen (sacpzfile, "r");
  if (!fp)
  {
    printf ("cannot open %s\n", sacpzfile);
    return -1;
  }

  char buffer[LEN];

  while (fgets (buffer, LEN, fp))
  {
    printf ("%s\n", buffer);
  }

  fclose (fp);

  return 0;
}
