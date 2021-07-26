#ifndef PARSE_MINISEED_H
#define PARSE_MINISEED_H

#include <stdint.h>

#include "libmseed.h"

#include "datatype.h"

int parse_miniSEED (const char *mseedfile, MS3Selections *selection, data_t **data, double *sampleRate, uint64_t *totalSamples);

#endif
