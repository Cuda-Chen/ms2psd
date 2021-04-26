#ifndef BANSPASS_FILTER_H
#define BANSPASS_FILTER_H

#include <complex.h>
#include <stdint.h>

#include "datatype.h"
#if 0
void bandpass_filter (double *data, double sampleRate, uint64_t totalSamples, int nfft,
                      double complex *filterResult, double complex *freqResponse);
#endif
void bandpass_filter_float (data_t *data, double sampleRate, uint64_t totalSamples, int nfft,
                            float lowcutFreq, float highcutFreq, int _order, int passes,
                            float complex *filterResult, float complex *freqResponse);

#endif
