#ifndef PSD_H
#define PSD_H

#include <complex.h>

int calculatePSD (double complex *in, int num_samples,
                  double sampling_rate, double **psd);

#endif
