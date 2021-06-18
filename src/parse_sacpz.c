#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LEN 0x200

/* Parse SAC PZ file. 
 * Referenced from here: https://github.com/savage13/sacpz2resp/blob/master/sacpz2resp.py#L76
 * */
int
parse_sacpz (const char *sacpzfile, double complex **poles,
             int *npoles,
             double complex **zeros,
             int *nzeros, double *constant)
{
  FILE *fp = fopen (sacpzfile, "r");
  if (!fp)
  {
    printf ("cannot open %s\n", sacpzfile);
    return -1;
  }

  char buffer[LEN];
  int read_zeros = 0, read_poles = 0;
  int i;
  int zeros_cnt = 0, poles_cnt = 0;

  while (fgets (buffer, LEN, fp))
  {
    char *pch = strtok (buffer, " ");
    if (strcmp (pch, "CONSTANT") == 0)
    {
      pch       = strtok (NULL, " ");
      *constant = atof (pch);
    }
    else if (strcmp (pch, "ZEROS") == 0)
    {
      pch     = strtok (NULL, " ");
      *nzeros = atoi (pch);
      *zeros  = (double complex *)malloc (sizeof (double complex) * (*nzeros));
      if ((*zeros) == NULL)
      {
        printf ("cannot allocate memery\n");
        return -1;
      }
      for (i = 0; i < *nzeros; i++)
        *((*zeros) + i) = 0.0 + 0.0 * I;
      read_zeros = 1;
      read_poles = 0;
    }
    else if (strcmp (pch, "POLES") == 0)
    {
      pch     = strtok (NULL, " ");
      *npoles = atoi (pch);
      *poles  = (double complex *)malloc (sizeof (double complex) * (*npoles));
      if ((*poles) == NULL)
      {
        printf ("cannot allocate memery\n");
        return -1;
      }
      for (i = 0; i < *npoles; i++)
        *((*poles) + i) = 0.0 + 0.0 * I;
      read_zeros = 0;
      read_poles = 1;
    }
    else
    {
      double real, imag;
      real = atof (pch);

      pch  = strtok (NULL, " ");
      imag = atof (pch);

      if (read_poles == 1)
      {
        *((*poles) + poles_cnt) = real + imag * I;
        poles_cnt++;
      }
      else if (read_zeros == 1)
      {
        *((*zeros) + zeros_cnt) = real + imag * I;
        zeros_cnt++;
      }
    }
  }

  fclose (fp);

  return 0;
}
