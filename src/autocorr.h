#ifndef AUTOCORR_H
#define AUTOCORR_H

#include <stdint.h>

#include "datatype.h"

void autocorr_float (data_t *data, uint64_t totalSamples,
                     data_t *autoCorrelationResult);

#endif
