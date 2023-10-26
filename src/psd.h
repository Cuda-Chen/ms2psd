#ifndef PSD_H
#define PSD_H

#include <complex.h>

void calculatePSD (float complex *in, int num_samples,
                   double sampling_rate, float *psd);

#endif
