#ifndef AUTOCORRELATION_H
#define AUTOCORRELATION_H

#include <stdint.h>

#include "datatype.h"

void autocorrelation_float (data_t *data, uint64_t totalSamples,
                            data_t **autoCorrelationResult);

#endif
