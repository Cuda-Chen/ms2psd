#ifndef SPGRAM_H
#define SPGRAM_H

#include <complex.h>
#include <stdint.h>

void spgram (float complex *data, uint64_t totalSamples, int nfft, float *psd);

int spectra (double complex *input, uint64_t totalSamples, double **amp, double **phase);

#endif
