#ifndef COSINE_TAPER_H
#define COSINE_TAPER_H

void cosineTaper (float *data, int n, float alpha, float *tapered);
void sacCosineTaper (double *freqs, int n, float f1, float f2, float f3, float f4, double sampling_rate, float *tapered);

#endif
