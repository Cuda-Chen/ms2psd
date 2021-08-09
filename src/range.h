#ifndef RANGE_H
#define RANGE_H

void range (double *array, double sampleRate, int totalSamples);
int setLeftAndRightFreq (double **leftFreqs, double **rightFreqs, int *freqLen, double sampleRate, int windowLength);

#endif
