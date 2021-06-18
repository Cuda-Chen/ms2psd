#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
  int read_zeros = 0, read_poles = 0;
  int nzeros, npoles;

  while (fgets (buffer, LEN, fp))
  {
    char *pch = strtok(buffer, " ");
    if(strcmp(pch, "CONSTANT") == 0) {
        printf("CONSTANT ");
        printf("%s", pch);
        pch = strtok(NULL, " ");
        printf(" %s", pch);
    } else if(strcmp(pch, "ZEROS") == 0) {
        printf("ZEROS");
        printf("%s", pch);
        pch = strtok(NULL, " ");
        printf(" %s", pch);
    } else if(strcmp(pch, "POLES") == 0) {
        printf("POLES");
        printf("%s", pch);
        pch = strtok(NULL, " ");
        printf(" %s", pch);
    } else {
        printf("%s", pch);
        pch = strtok(NULL, " ");
        printf(" %s", pch);
    }
    printf("\n");
  }

  fclose (fp);

  return 0;
}
