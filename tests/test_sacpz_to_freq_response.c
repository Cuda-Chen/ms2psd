#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "get_freq_response.h"
#include "parse_sacpz.h"

int
main ()
{
  const char *sacpzfile = "SAC_PZs_TW_NACB_BHZ__2007.254.07.25.20.0000_99999.9999.24.60.60.99999";
  double complex *poles;
  double complex *zeros;
  double constant;
  int nzeros, npoles;
  int result;

  result = parse_sacpz (sacpzfile, &poles, &npoles, &zeros, &nzeros, &constant);
  assert (result == 0);

  double *freq;
  double *freq_response;
  double sampling_rate = 100.0f;
  int data_samples     = 1024;
  result               = get_freq_response (poles, npoles, zeros, nzeros,
                              constant, sampling_rate, data_samples,
                              &freq, &freq_response, 0);
  assert (result == 0);
  int i;
  for (i = 0; i < data_samples; i++)
  {
    printf ("%lf %e\n", freq[i], freq_response[i]);
  }

  /* Free allocated objects */
  free (poles);
  free (zeros);
  free (freq);
  free (freq_response);

  return 0;
}
