#ifndef PARSE_MINISEED_H
#define PARSE_MINISEED_H

#include <stdint.h>

#include "libmseed.h"

#include "datatype.h"

int parse_miniSEED_from_file (const char *mseedfile, MS3Selections *selection, data_t **data, double *sampleRate, uint64_t *totalSamples);
int parse_miniSEED_from_stream (char *buffer, uint64_t bufferlength, MS3Selections *selection, data_t **data, double *sampleRate, uint64_t *totalSamples);

#endif
