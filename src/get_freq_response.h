#ifndef GET_FREQ_RESPONSE_H
#define GET_FREQ_RESPONSE_H

#include <complex.h>

int get_freq_response (double complex *poles, int npoles,
                       double complex *zeros, int nzeros,
                       double constant, double sampling_rate, int data_samples,
                       double **freq, double **freq_response, int flag);

#endif
