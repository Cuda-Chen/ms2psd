#ifndef FFT_H
#define FFT_H

#include "fftw3.h"

void fft (float *data, uint64_t dataSamples, double complex *output);
void ifft (double complex *data, uint64_t dataSamples, double *output);
void fftToFileHalf (double *data, uint64_t dataSamples, double sampleRate, FILE *fptr);
void fftToFile (double *data, uint64_t dataSamples, double sampleRate, FILE *fptr);
void testFFT (double *data, fftw_complex *in, fftw_complex *out,
              fftw_complex *ref, uint64_t totalSamples);

#endif
