#ifndef GET_FREQ_RESPONSE_H
#define GET_FREQ_RESPONSE_H

#include <complex.h>

int get_freq_response (double complex *poles, int npoles,
                       double complex *zeros, int nzeros,
                       double constant, double sampling_rate, int data_samples,
                       double **freq, double complex **freq_response, int flag);
int remove_response (double complex *data, double complex *freq_response, int data_samples);

#endif
