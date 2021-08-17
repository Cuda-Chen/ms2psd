#ifndef PROCESS_TRACE_H
#define PROCESS_TRACE_H

int
getTraceProperties (const char *mseedfile, nstime_t *starttime, nstime_t *endtime, double *samplingRate, int *totalSegmentsOfHour);

int processTrace (const char *mseedfile,
                  float f1,
                  float f2,
                  float f3,
                  float f4,
                  int totype,
                  const char *sacpzfile,
                  const char *outputFile);

#endif
