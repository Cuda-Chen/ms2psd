#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "libmseed.h"

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

nstime_t NSECS = 1000000000;

/* 1-hour long segment properties */
int lengthOfOneHour  = 3600;
int overlapOfOneHour = 50;

/* 900 seconds, or 15-minute long segment properties */
int lengthOfSegment  = 900;
int overlapOfSegment = 75;

static int
getTraceProperties (const char *mseedfile, nstime_t *starttime, nstime_t *endtime, double *samplingRate)
{
  uint32_t flags = 0;
  int8_t verbose = 0;
  flags |= MSF_VALIDATECRC;
  //flags |= MSF_RECORDLIST;
  int rv;

  MS3TraceList *mstl = NULL;
  rv                 = ms3_readtracelist (&mstl, mseedfile, NULL, 0, flags, verbose);
  if (rv != 0)
    return rv;
  *starttime    = mstl->traces->earliest;
  *endtime      = mstl->last->latest;
  *samplingRate = mstl->traces->first->samprate;

#ifdef DEBUG
  char starttimestr[30];
  char endtimestr[30];
  ms_nstime2timestr (mstl->traces->earliest, starttimestr, ISOMONTHDAY, NANO);
  ms_nstime2timestr (mstl->last->latest, endtimestr, ISOMONTHDAY, NANO);
  printf ("%s - %s, %lf Hz\n", starttimestr, endtimestr, *samplingRate);
#endif

  /* Clean up */
  if (mstl)
    mstl3_free (&mstl, 0);

  return 0;
}

int
processTrace (const char *mseedfile,
              float f1,
              float f2,
              float f3,
              float f4,
              int totype,
              const char *sacpzfile,
              const char *outputFile)
{
  data_t *data = NULL;
  double sampleRate;
  uint64_t totalSamples;

  int rv;
  uint8_t pubversion = 0;

  /* Get the start time and end time of this trace */
  nstime_t starttimeOfTrace;
  nstime_t endtimeOfTrace;
  rv = getTraceProperties (mseedfile, &starttimeOfTrace, &endtimeOfTrace, &sampleRate);
  if (rv != 0)
  {
    fprintf (stderr, "Cannot open input miniSEED for getting start time and end time\n");
    return -1;
  }

  int nextTimeStamp;
  nstime_t nextTimeStampNS;
  nstime_t starttime;
  nstime_t endtime;

  /* Split trace to 1-hour long segment with 50% overlapping
   * for reducing processing time */

  /* Split 1-hour long segment to 15-minute long segment
   * with 75% overlapping for reducing data variance */
  nextTimeStamp                = lengthOfSegment - (lengthOfSegment * overlapOfSegment / 100);
  nextTimeStampNS              = nextTimeStamp * NSECS;
  starttime                    = starttimeOfTrace;
  endtime                      = starttime + ((nstime_t)lengthOfSegment * NSECS);
  char *sidpattern             = "*";
  MS3Selections *fooselections = NULL;

  int segments         = 0;
  nstime_t endtimeTick = endtime;
  while (endtimeTick <= endtimeOfTrace)
  {
    rv = ms3_addselect (&fooselections, sidpattern, starttime, endtime, pubversion);
    starttime += nextTimeStampNS;
    endtime += nextTimeStampNS;
    endtimeTick += nextTimeStampNS;
    segments++;
  }
  int psdBinWindowSize = sampleRate * lengthOfSegment;
  double *psdBin       = (double *)malloc (sizeof (double) * segments * psdBinWindowSize);
  for (int i = 0; i < segments * psdBinWindowSize; i++)
  {
    psdBin[i] = 0.0f;
  }

#ifdef DEBUG
  ms3_printselections (fooselections);
#endif

  if (fooselections)
    ms3_freeselections (fooselections);

  // reset
  starttime = starttimeOfTrace;
  endtime   = starttime + ((nstime_t)lengthOfSegment * NSECS);
  /* Get data from input miniSEED file */
  for (int a = 0; a < segments; a++)
  {
    MS3Selections *selection = NULL;
    rv                       = ms3_addselect (&selection, sidpattern, starttime, endtime, pubversion);
    rv                       = parse_miniSEED (mseedfile, selection, &data, &sampleRate, &totalSamples);
    if (rv != 0)
    {
      return rv;
    }
    if (data == NULL)
    {
      printf ("Input data read unsuccessfully\n");
      return -1;
    }

    /* Demean */
    float mean, std;
    getMeanAndSD (data, totalSamples, &mean, &std);
    for (int i = 0; i < totalSamples; i++)
    {
      data[i] -= mean;
    }
    /* Detrend */
    data_t *detrended;
    detrend (data, (int)totalSamples, &detrended);
    /* First taper the signal with 5%-cosine-window */
    float *taperedSignal = (float *)malloc (sizeof (float) * totalSamples);
    cosineTaper (detrended, (int)totalSamples, 0.05, taperedSignal);
    double *tapered = (double *)malloc (sizeof (double) * totalSamples);
    for (int i = 0; i < totalSamples; i++)
      tapered[i] = (double)taperedSignal[i];
    /* Then execute FFT */
    double complex *fftResult;
    fft (tapered, totalSamples, &fftResult);

    double delta         = 1. / sampleRate;
    double totalDuration = delta * totalSamples;

    /* instrument response removal */
    /* Read SACPZ file */
    //const char *sacpzfile = "./tests/SAC_PZs_TW_NACB_BHZ__2007.254.07.25.20.0000_99999.9999.24.60.60.99999";
    double complex *poles, *zeros;
    int npoles, nzeros;
    double constant;
    parse_sacpz (sacpzfile,
                 &poles, &npoles,
                 &zeros, &nzeros,
                 &constant);
    /* Get frequency response */
    double *freq;
    double complex *freqResponse;
    get_freq_response (poles, npoles,
                       zeros, nzeros,
                       constant, sampleRate, totalSamples,
                       &freq, &freqResponse, totype);

    /* Apply freqyency response removal */
    remove_response (fftResult, freqResponse, totalSamples);

    /* band-pass filtering to prevent overamplification */
    double *taper_window;
    sacCosineTaper (freq, totalSamples, f1, f2, f3, f4, sampleRate, &taper_window);
    for (int i = 0; i < totalSamples; i++)
    {
      fftResult[i] *= taper_window[i];
    }

    /* Apply iFFT */
    double *output;
    ifft (fftResult, totalSamples, &output);
    /* Output signal with instrument response removed */
    FILE *out = fopen (outputFile, "w");
    if (out == NULL)
    {
      fprintf (stderr, "Cannot output signal output file.\n");
      return -1;
    }
    for (int i = 0; i < totalSamples; i++)
    {
      fprintf (out, "%e\n", output[i]);
    }
    fclose (out);

    /* Get power spetral density (PSD) of this segment */
    double *psd;

    /* Though McMarana 2004 mentions you should divide delta-t for each frequency response,
   * we do not do such operation as this produces resaonable result
   * Yet I still reverse this block for future use
   */
    /*for(int i = 0; i < totalSamples; i++) {
      fftResult[i] /= (1. / sampleRate);
  }*/

    rv = calculatePSD (fftResult, totalSamples, sampleRate, &psd);
    if (rv != 0)
    {
      fprintf (stderr, "Something wrong in calculate PSD procedure,\n");
      return -1;
    }

    for (int i = 0; i < psdBinWindowSize; i++)
    {
      int psdBinIndex         = a * psdBinWindowSize + i;
      *(psdBin + psdBinIndex) = *(psd + i);
    }
#if 0
  /* Set unit to decibel (dB) */
  for (int i = 0; i < totalSamples; i++)
  {
    psd[i] = 10 * log10 (psd[i]);
  }
  FILE *psd_out = fopen ("psd_out.txt", "w");
  for (int i = 0; i < totalSamples; i++)
  {
    fprintf (psd_out, "%e %e\n", freq[i], psd[i]);
  }
  fclose (psd_out);
#endif

    /* Free allocated objects */
    free (data);
    free (detrended);
    free (taperedSignal);
    free (tapered);
    free (fftResult);
    free (poles);
    free (zeros);
    free (freq);
    free (freqResponse);
    free (taper_window);
    free (psd);
    free (output);

    starttime += nextTimeStampNS;
    endtime += nextTimeStampNS;
    if (selection)
      ms3_freeselections (selection);
  }

  /* PSD properties (min, max, mean, median) summary */
  double *psdMean = (double *)malloc (sizeof (double) * psdBinWindowSize);
  for (int i = 0; i < psdBinWindowSize; i++)
    psdMean[i] = 0.0f;
  /* Mean calculation */
  for (int i = 0; i < psdBinWindowSize; i++)
  {
    for (int j = 0; j < segments; j++)
    {
      int psdBinIndex = j * psdBinWindowSize + i;
      *(psdMean + i) += *(psdBin + psdBinIndex);
    }
  }
  for (int i = 0; i < psdBinWindowSize; i++)
    psdMean[i] /= (double)segments;

  /* Set unit to decibel (dB) */
  for (int i = 0; i < psdBinWindowSize; i++)
  {
    psdMean[i] = 10 * log10 (psdMean[i]);
  }

  /* Output PSD */
  double *estimatedFreqs = (double *)malloc (sizeof (double) * psdBinWindowSize);
  range (estimatedFreqs, sampleRate, psdBinWindowSize);
  FILE *psd_out = fopen ("psd_out.txt", "w");
  for (int i = 0; i < psdBinWindowSize; i++)
  {
    fprintf (psd_out, "%e %e\n", estimatedFreqs[i], psdMean[i]);
  }
  fclose (psd_out);
  // all segment PSD
#if 0
  for(int i = 0; i < segments; i++) {
      char psd_file[50];
      sprintf(psd_file, "psd_trace_%d.txt", i);
      FILE *out = fopen(psd_file, "w");
      for(int j = 0; j < psdBinWindowSize; j++) {
          int index = i * psdBinWindowSize + j;
          fprintf(out, "%e %e\n", estimatedFreqs[j], *(psdBin + index));
      }
      fclose(out);
  }
#endif

  free (psdBin);
  free (psdMean);
  free (estimatedFreqs);

  return 0;
}
