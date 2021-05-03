#ifndef PARSE_MINISEED_H
#define PARSE_MINISEED_H

#include <stdint.h>

#include "datatype.h"

int parse_miniSEED (const char *mseedfile, data_t **data, double *sampleRate, uint64_t *totalSamples);

#endif
