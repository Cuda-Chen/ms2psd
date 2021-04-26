#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "liquid.h"

#include "bandpass_filter.h"
#include "datatype.h"
#if 0
void
bandpass_filter (double *data, double sampleRate, uint64_t totalSamples, int nfft,
                 double complex *filterResult, double complex *freqResponse)
{
  // options
  liquid_iirdes_filtertype ftype = LIQUID_IIRDES_BUTTER;
  liquid_iirdes_bandtype btype   = LIQUID_IIRDES_BANDPASS;
  liquid_iirdes_format format    = LIQUID_IIRDES_TF;
  unsigned int order             = 2;            // filter order
  float fc                       = 0.01f;        // cutoff frequency
  float f0                       = 0.1f;         // center frequency
  float Ap                       = 1.0f;         // pass-band ripple
  float As                       = 40.0f;        // stop-band attenuation
  uint64_t n                     = totalSamples; // number of samples

  // design filter from prototype
  iirfilt_crcf q = iirfilt_crcf_create_prototype (
      ftype, btype, format, order, fc, f0, Ap, As);
  iirfilt_crcf_print (q);

  uint64_t i;
  for (i = 0; i < n; i++)
  {
    // run filter
    iirfilt_crcf_execute (q, data[i], &filterResult[i]);
  }

  // compute frequency response
  int counter;
  for (counter = 0; counter < nfft; counter++)
  {
    double freq = 0.5f * (double)counter / (double)nfft;
    iirfilt_crcf_freqresponse (q, freq, &freqResponse[counter]);
  }

  // destroy filter object
  iirfilt_crcf_destroy (q);
}
#endif
void
bandpass_filter_float (data_t *data, double sampleRate, uint64_t totalSamples, int nfft,
                       float lowcutFreq, float highcutFreq, int _order, int passes,
                       float complex *filterResult, float complex *freqResponse)
{
  // options
  liquid_iirdes_filtertype ftype = LIQUID_IIRDES_BUTTER;
  liquid_iirdes_bandtype btype   = LIQUID_IIRDES_BANDPASS;
  liquid_iirdes_format format    = LIQUID_IIRDES_SOS;
  unsigned int order             = _order; // filter order
  int passesFlag                 = (passes == 2) ? 1 : 0;
  float fc                       = lowcutFreq / sampleRate;                      // cutoff frequency
  float f0                       = sqrt (lowcutFreq * highcutFreq) / sampleRate; // center frequency
  //float f0 = (lowcutFreq + highcutFreq) / 2 / sampleRate;
  float Ap                  = 1.0f;         // pass-band ripple
  float As                  = 40.0f;        // stop-band attenuation
  uint64_t n                = totalSamples; // number of samples
  float complex *filterTemp = (float complex *)malloc (sizeof (float complex) * n);

  printf ("lowcut: %f highcut: %f : center: %f\n", fc, highcutFreq / sampleRate, f0);

  // design filter from prototype
  iirfilt_crcf q = iirfilt_crcf_create_prototype (
      ftype, btype, format, order, fc, f0, Ap, As);
  iirfilt_crcf_print (q);

  //uint64_t i; // setting this result in segfault when backward filtering
  int i;
  for (i = 0; i < n; i++) // forward filtering
  {
    iirfilt_crcf_execute (q, data[i], &filterTemp[i]);
  }
  if (passesFlag)
  {
    iirfilt_crcf_reset (q);

    for (i = n - 1; i >= 0; i--) // backward filtering
    {
      iirfilt_crcf_execute (q, filterTemp[i], &filterResult[i]);
    }
  }
  else
  {
    for (i = 0; i < n; i++)
    {
      filterResult[i] = filterTemp[i];
    }
  }

  // compute frequency response
  int counter;
  for (counter = 0; counter < nfft; counter++)
  {
    double freq = 0.5f * (double)counter / (double)nfft;
    iirfilt_crcf_freqresponse (q, freq, &freqResponse[counter]);
  }

  // destroy filter object
  iirfilt_crcf_destroy (q);
}
