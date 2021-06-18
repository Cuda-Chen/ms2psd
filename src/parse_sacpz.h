#ifndef PARSE_SACPZ_H
#define PARSE_SACPZ_H

#include <complex.h>

int parse_sacpz (const char *sacpzfile, double complex **poles,
                 int *npoles,
                 double complex **zeros,
                 int *nzeros, double *constant);

#endif
