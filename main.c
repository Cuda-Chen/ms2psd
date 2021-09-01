#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
//#include <assert.h>

#include "fftw3.h"

#include "autocorr.h"
#include "bandpass_filter.h"
#include "common.h"
#include "cosine_taper.h"
#include "datatype.h"
#include "detrend.h"
#include "fft.h"
#include "freq_response.h"
#include "output2Octave.h"
#include "parse_miniSEED.h"
#include "parse_sacpz.h"
#include "psd.h"
#include "standard_deviation.h"

#include "process_trace.h"

static void
usage ()
{
  printf ("Usage: ./ms2psd [f1] [f2] [f3] [f4] [totype] [input] [resp]");
  printf ("\n\nInput parameters:\n");
  printf ("f1, f2, f3, f4: four-corner frequencies (Hz)\n");
  printf ("totype: specify the following numbers for output waveform format:\n");
  printf ("        0: displacement\n");
  printf ("        1: velocity\n");
  printf ("        2: acceleration\n");
  printf ("input: input waveform. Should be miniSEED format\n");
  printf ("resp: response file in SACPZ format\n");
}

int
main (int argc, char **argv)
{
  char *mseedfile = NULL;

  float f1, f2, f3, f4; /* four-corner frequencies */
  int totype;           /* output waveform model, e.g., displacement, velocity, or acceleration */
  char *sacpzfile;
  int rv;

  /* Simple argement parsing */
  if (argc != 8)
  {
    usage ();
    return 1;
  }
  f1        = atof (argv[1]);
  f2        = atof (argv[2]);
  f3        = atoi (argv[3]);
  f4        = atoi (argv[4]);
  totype    = atoi (argv[5]);
  mseedfile = argv[6];
  sacpzfile = argv[7];

  rv = processTrace (mseedfile,
                     f1,
                     f2,
                     f3,
                     f4,
                     totype,
                     sacpzfile);

  return 0;
}
