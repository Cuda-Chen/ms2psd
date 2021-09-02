#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

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
#include "range.h"
#include "standard_deviation.h"

nstime_t NSECS = 1000000000;

/* 1-hour long segment properties, 50% overlap*/
int lengthOfOneHour  = 3600;
int overlapOfOneHour = 50;

/* 900 seconds, or 15-minute long segment properties, 75% overlap */
int lengthOfSegment  = 900;
int overlapOfSegment = 75;

/* PDF properties */
int mindB = -200;
int maxdB = -50;

static int
compare (const void *a, const void *b)
{
  const double *da = (const double *)a;
  const double *db = (const double *)b;
  return (*da > *db) - (*da < *db);
}

static double
decibel (const double a)
{
  return 10 * log10 (a);
}

static int
binLocation (const double v, double start)
{
  return (int)(fabs (round (v)) - fabs (start));
}

int
getTraceProperties (const char *mseedfile, nstime_t *starttime, nstime_t *endtime, double *samplingRate, int *totalSegmentsOfOneHour)
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
  *totalSegmentsOfOneHour            = 0;
  int nextTimeStampOfSegments        = lengthOfOneHour - (lengthOfOneHour * overlapOfOneHour / 100);
  nstime_t nextTimeStampOfSegmentsNS = nextTimeStampOfSegments * NSECS;
  nstime_t start                     = mstl->traces->earliest;
  nstime_t end                       = start + ((nstime_t)lengthOfOneHour * NSECS);
  while (end <= mstl->last->latest)
  {
    start += nextTimeStampOfSegmentsNS;
    end += nextTimeStampOfSegmentsNS;
    (*totalSegmentsOfOneHour)++;
  }

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
              const char *sacpzfile)
{
  //data_t *data = NULL;
  double sampleRate;
  uint64_t totalSamples;

  int rv;
  char *sidpattern   = "*";
  uint8_t pubversion = 0;

  /* Get the start time and end time of this trace */
  nstime_t starttimeOfTrace;
  nstime_t endtimeOfTrace;
  int totalSegmentsOfOneHour = 1;
  rv                         = getTraceProperties (mseedfile, &starttimeOfTrace, &endtimeOfTrace, &sampleRate, &totalSegmentsOfOneHour);
  if (rv != 0)
  {
    fprintf (stderr, "Cannot open input miniSEED for getting start time and end time\n");
    return -1;
  }
  int nextTimeStampOfHours        = lengthOfOneHour - (lengthOfOneHour * overlapOfOneHour / 100);
  nstime_t nextTimeStampOfHoursNS = nextTimeStampOfHours * NSECS;
  nstime_t starttimeOfThisHour    = starttimeOfTrace;
  nstime_t endtimeOfThisHour      = starttimeOfThisHour + ((nstime_t)lengthOfOneHour * NSECS);
  //totalSegmentsOfOneHour          = 1;

  /* Set the left and right frequency limit of each octave used in dimension reduction part */
  double *leftFreqs, *rightFreqs;
  int freqLen; /* i.e., length of center periods */
  double smoothingWidthFactor = 1.25;
  rv                          = setLeftAndRightFreq (&leftFreqs, &rightFreqs, &freqLen, sampleRate, lengthOfSegment, smoothingWidthFactor);
  if (rv != 0)
  {
    return -1;
  }
  /* Set center frequency (geometric mean of low and high frequency) */
  double *centerPeriods = (double *)malloc (sizeof (double) * freqLen);
  for (int i = 0; i < freqLen; i++)
  {
    double ts        = 1. / leftFreqs[i];
    double tl        = 1. / rightFreqs[i];
    centerPeriods[i] = sqrt (ts * tl);
  }
  /* PDF statistics */
  int psdBinWindowSize                 = sampleRate * lengthOfSegment;
  double *psdBinReducedMeanAggerated   = (double *)malloc (sizeof (double) * totalSegmentsOfOneHour * freqLen);
  double *psdBinReducedMinAggerated    = (double *)malloc (sizeof (double) * totalSegmentsOfOneHour * freqLen);
  double *psdBinReducedMaxAggerated    = (double *)malloc (sizeof (double) * totalSegmentsOfOneHour * freqLen);
  double *psdBinReducedMedianAggerated = (double *)malloc (sizeof (double) * totalSegmentsOfOneHour * freqLen);
  double *estimatedFreqs               = (double *)malloc (sizeof (double) * psdBinWindowSize);
  range (estimatedFreqs, sampleRate, psdBinWindowSize);
  /* PSD properties (min, max, mean, median) summary statistics for 1-hour long segment */
  double *psdMin    = (double *)malloc (sizeof (double) * psdBinWindowSize);
  double *psdMax    = (double *)malloc (sizeof (double) * psdBinWindowSize);
  double *psdMean   = (double *)malloc (sizeof (double) * psdBinWindowSize);
  double *psdMedian = (double *)malloc (sizeof (double) * psdBinWindowSize);

  /* Read miniSEED file to buffer */
  FILE *fp;
  char *inputmseedBuffer = NULL;
  uint64_t inputmseedBufferLength;
  struct stat sb = {0};
  if (!(fp = fopen (mseedfile, "rb")))
  {
    fprintf (stderr, "Cannot open input miniSEED file\n");
    return -1;
  }
  if (fstat (fileno (fp), &sb))
  {
    fprintf (stderr, "Cannot stat input miniSEED file\n");
    return -1;
  }
  if (!(inputmseedBuffer = (char *)malloc (sb.st_size)))
  {
    fprintf (stderr, "Cannot allocate input miniSEED buffer\n");
    return -1;
  }
  if (1 != (fread (inputmseedBuffer, sb.st_size, 1, fp)))
  {
    fprintf (stderr, "Cannot read miniSEED file into buffer\n");
    return -1;
  }
  fclose (fp);
  inputmseedBufferLength = sb.st_size;

  /* Read SACPZ file */
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
  int samples = lengthOfSegment * sampleRate;
  get_freq_response (poles, npoles,
                     zeros, nzeros,
                     constant, sampleRate, samples,
                     &freq, &freqResponse, totype);

  /* Split trace to 1-hour long segment with 50% overlapping
   * for reducing processing time */
  for (int traceIdx = 0; traceIdx < totalSegmentsOfOneHour; traceIdx++)
  {
    /* Split 1-hour long segment to 15-minute long segment
   * with 75% overlapping for reducing data variance */
    int nextTimeStampOfSegments        = lengthOfSegment - (lengthOfSegment * overlapOfSegment / 100);
    nstime_t nextTimeStampOfSegmentsNS = nextTimeStampOfSegments * NSECS;
    nstime_t starttimeOfThisSegment    = starttimeOfThisHour;
    nstime_t endtimeOfThisSegment      = starttimeOfThisSegment + ((nstime_t)lengthOfSegment * NSECS);
    nstime_t endtimeTemp               = starttimeOfThisHour + ((nstime_t)lengthOfOneHour * NSECS);
#ifdef DEBUG
    MS3Selections *fooselections = NULL;
#endif

    /* Count number of available 15-minutes segments in this 1-hour long segment */
    int segments = 0;
    while (endtimeOfThisSegment <= endtimeTemp)
    {
#ifdef DEBUG
      rv = ms3_addselect (&fooselections, sidpattern, starttimeOfThisSegment, endtimeOfThisSegment, pubversion);
#endif
      starttimeOfThisSegment += nextTimeStampOfSegmentsNS;
      endtimeOfThisSegment += nextTimeStampOfSegmentsNS;
      segments++;
    }
    double *psdBin = (double *)malloc (sizeof (double) * segments * psdBinWindowSize);
    for (int i = 0; i < segments * psdBinWindowSize; i++)
    {
      psdBin[i] = 0.0f;
    }

#ifdef DEBUG
    ms3_printselections (fooselections);
    if (fooselections)
      ms3_freeselections (fooselections);
#endif

    /* Reset start time and end time to the first 15-minute long segment */
    starttimeOfThisSegment = starttimeOfThisHour;
    endtimeOfThisSegment   = starttimeOfThisSegment + ((nstime_t)lengthOfSegment * NSECS);
    /* Get data from input miniSEED file */
    for (int a = 0; a < segments; a++)
    {
      MS3Selections *selection = NULL;
      data_t *data_temp;
      rv = ms3_addselect (&selection, sidpattern, starttimeOfThisSegment, endtimeOfThisSegment, pubversion);
      if (rv != 0)
      {
        printf ("selection failed\n");
        return rv;
      }

      rv = parse_miniSEED_from_stream (inputmseedBuffer, inputmseedBufferLength, selection, &data_temp, &sampleRate, &totalSamples);
      if (rv != 0)
      {
        return rv;
      }
      if (data_temp == NULL)
      {
        printf ("Input data read unsuccessfully\n");
        return -1;
      }

      /* Adjust input data length to 15-minute long equivalent */
      int desiredSamples = (int)(lengthOfSegment * sampleRate);
      data_t *data       = (data_t *)malloc (sizeof (data_t) * desiredSamples);
      if (data == NULL)
      {
        printf ("Segment data array allocation failed\n");
        return -1;
      }
      /* Padding if input data less than 15-minute */
      if (totalSamples < desiredSamples)
      {
        for (int i = 0; i < totalSamples; i++)
        {
          data[i] = data_temp[i];
        }
        for (int i = totalSamples; i < desiredSamples; i++)
        {
          data[i] = 0.0f;
        }
      }
      else
      {
        for (int i = 0; i < desiredSamples; i++)
        {
          data[i] = data_temp[i];
        }
      }
      totalSamples = desiredSamples;

      /* Demean */
      float mean, std;
      getMeanAndSD (data, totalSamples, &mean, &std);
      for (int i = 0; i < totalSamples; i++)
      {
        data[i] -= mean;
      }
      /* Detrend */
      data_t *detrended = (data_t *)malloc (sizeof (data_t) * totalSamples);
      if (detrended == NULL)
      {
        fprintf (stderr, "Cannot allocate detrend output memory.\n");
        return -1;
      }
      detrend (data, (int)totalSamples, detrended);
      /* First taper the signal with 5%-cosine-window */
      float *taperedSignal = (float *)malloc (sizeof (float) * totalSamples);
      cosineTaper (detrended, (int)totalSamples, 0.05, taperedSignal);
      double *tapered = (double *)malloc (sizeof (double) * totalSamples);
      for (int i = 0; i < totalSamples; i++)
        tapered[i] = (double)taperedSignal[i];
      /* Then execute FFT */
      double complex *fftResult = (double complex *)malloc (sizeof (double complex) * totalSamples);
      fft (tapered, totalSamples, fftResult);

      /* instrument response removal */
      /* Apply freqyency response removal */
      remove_response (fftResult, freqResponse, totalSamples);

      /* band-pass filtering to prevent overamplification */
      double *taper_window = (double *)malloc (sizeof (double) * totalSamples);
      if (taper_window == NULL)
      {
        fprintf (stderr, "taper window allocation failed\n");
        return -1;
      }
      for (int i = 0; i < totalSamples; i++)
        taper_window[i] = 0.0f;
      sacCosineTaper (freq, totalSamples, f1, f2, f3, f4, sampleRate, taper_window);
      for (int i = 0; i < totalSamples; i++)
      {
        fftResult[i] *= taper_window[i];
      }

      /* Get power spetral density (PSD) of this segment */
      double *psd = (double *)malloc (sizeof (double) * totalSamples);
      if (psd == NULL)
      {
        fprintf (stderr, "cannot allocate PSD space\n");
        return -1;
      }

      /* Though McMarana 2004 mentions you should divide delta-t for each frequency response,
   * we do not do such operation as this produces resaonable result
   * Yet I still reverse this block for future use
   */
      /*for(int i = 0; i < totalSamples; i++) {
      fftResult[i] /= (1. / sampleRate);
  }*/

      calculatePSD (fftResult, totalSamples, sampleRate, psd);

      for (int i = 0; i < psdBinWindowSize; i++)
      {
        int psdBinIndex         = a * psdBinWindowSize + i;
        *(psdBin + psdBinIndex) = *(psd + i);
      }

      /* Free allocated objects */
      free (data_temp);
      free (data);
      free (detrended);
      free (taperedSignal);
      free (tapered);
      free (fftResult);

      free (taper_window);
      free (psd);

      starttimeOfThisSegment += nextTimeStampOfSegmentsNS;
      endtimeOfThisSegment += nextTimeStampOfSegmentsNS;
      if (selection)
        ms3_freeselections (selection);
    }

    /* PSD summary of this 1-hour long segment */
    /* Temporary array for PSD summary sorting */
    double *psdArr = (double *)malloc (sizeof (double) * segments);
    for (int i = 0; i < psdBinWindowSize; i++)
      psdMean[i] = 0.0f;
    for (int i = 0; i < psdBinWindowSize; i++)
    {
      for (int j = 0; j < segments; j++)
      {
        int psdBinIndex = j * psdBinWindowSize + i;
        *(psdMean + i) += *(psdBin + psdBinIndex);
        *(psdArr + j) = *(psdBin + psdBinIndex);
      }

      qsort (psdArr, segments, sizeof (psdArr[0]), compare);
      psdMin[i]    = psdArr[0];
      psdMax[i]    = psdArr[segments - 1];
      psdMedian[i] = (segments % 2 == 0) ? ((psdArr[segments / 2 - 1] + psdArr[segments / 2]) / 2.0) : (psdArr[segments / 2]);
    }
    /* Mean calculation */
    for (int i = 0; i < psdBinWindowSize; i++)
      psdMean[i] /= (double)segments;

    /* Set unit to decibel (dB) */
    for (int i = 0; i < psdBinWindowSize; i++)
    {
      psdMin[i]    = decibel (psdMin[i]);
      psdMax[i]    = decibel (psdMax[i]);
      psdMean[i]   = decibel (psdMean[i]);
      psdMedian[i] = decibel (psdMedian[i]);
    }

    /* Dimension reduction technique escribed in McMarana 2004 */
    /* Dimension reduction */
    double *psdBinReduced = (double *)malloc (sizeof (double) * segments * freqLen);
    for (int i = 0; i < segments * freqLen; i++)
    {
      psdBinReduced[i] = 0.0f;
    }
    for (int i = 0; i < segments; i++)
    {
      for (int j = 0; j < freqLen; j++)
      {
        int count              = 0;
        int psdBinReducedIndex = i * freqLen + j;
        for (int k = 0; k < psdBinWindowSize; k++)
        {
          if ((leftFreqs[j] >= estimatedFreqs[k]) && (estimatedFreqs[k] >= rightFreqs[j]))
          {
            int psdBinIndex = i * psdBinWindowSize + k;
            psdBinReduced[psdBinReducedIndex] += psdBin[psdBinIndex];
            count++;
          }
        }

        psdBinReduced[psdBinReducedIndex] /= (double)count;
      }
    }

    free (psdBin);

    free (psdArr);

    /* Statistics with dimension reduction */
    double *psdBinReducedMean   = (double *)malloc (sizeof (double) * freqLen);
    double *psdBinReducedMin    = (double *)malloc (sizeof (double) * freqLen);
    double *psdBinReducedMax    = (double *)malloc (sizeof (double) * freqLen);
    double *psdBinReducedMedian = (double *)malloc (sizeof (double) * freqLen);
    double *psdBinReducedArr    = (double *)malloc (sizeof (double) * segments);
    for (int i = 0; i < freqLen; i++)
    {
      psdBinReducedMean[i] = 0.0f;
    }
    for (int i = 0; i < freqLen; i++)
    {
      for (int j = 0; j < segments; j++)
      {
        int psdBinReducedIndex = j * freqLen + i;
        psdBinReducedMean[i] += psdBinReduced[psdBinReducedIndex];
        psdBinReducedArr[j] = psdBinReduced[psdBinReducedIndex];
      }
      qsort (psdBinReducedArr, segments, sizeof (psdBinReducedArr[0]), compare);
      psdBinReducedMin[i]    = psdBinReducedArr[0];
      psdBinReducedMax[i]    = psdBinReducedArr[segments - 1];
      psdBinReducedMedian[i] = (segments % 2 == 0) ? ((psdBinReducedArr[segments / 2 - 1] + psdBinReducedArr[segments / 2]) / 2.0) : (psdBinReducedArr[segments / 2]);
      psdBinReducedMean[i] /= (double)segments;

      psdBinReducedMean[i]   = decibel (psdBinReducedMean[i]);
      psdBinReducedMin[i]    = decibel (psdBinReducedMin[i]);
      psdBinReducedMax[i]    = decibel (psdBinReducedMax[i]);
      psdBinReducedMedian[i] = decibel (psdBinReducedMedian[i]);
    }

    /* Aggerate all reduced sample PSD */
    {
      int i = traceIdx;
      //printf ("%d\n", traceIdx);
      for (int j = 0; j < freqLen; j++)
      {
        int idx                           = i * freqLen + j;
        psdBinReducedMeanAggerated[idx]   = psdBinReducedMean[j];
        psdBinReducedMinAggerated[idx]    = psdBinReducedMin[j];
        psdBinReducedMaxAggerated[idx]    = psdBinReducedMax[j];
        psdBinReducedMedianAggerated[idx] = psdBinReducedMedian[j];
      }
    }
    free (psdBinReduced);
    free (psdBinReducedMean);
    free (psdBinReducedMin);
    free (psdBinReducedMax);
    free (psdBinReducedMedian);
    free (psdBinReducedArr);

    starttimeOfThisHour += nextTimeStampOfHoursNS;
    endtimeOfThisHour += nextTimeStampOfHoursNS;
  }

  /* Calculate PDF (power density function) */
  int PDFBins       = abs (maxdB - mindB) + 1;
  double *pdfMean   = (double *)malloc (sizeof (double) * PDFBins * freqLen);
  double *pdfMin    = (double *)malloc (sizeof (double) * PDFBins * freqLen);
  double *pdfMax    = (double *)malloc (sizeof (double) * PDFBins * freqLen);
  double *pdfMedian = (double *)malloc (sizeof (double) * PDFBins * freqLen);
  for (int i = 0; i < PDFBins * freqLen; i++)
  {
    pdfMean[i]   = 0;
    pdfMin[i]    = 0;
    pdfMax[i]    = 0;
    pdfMedian[i] = 0;
  }
  /* Histogram accumulation */
  for (int i = 0; i < freqLen; i++)
  {
    for (int j = 0; j < totalSegmentsOfOneHour; j++)
    {
      int idx = j * freqLen + i;
      pdfMean[binLocation (psdBinReducedMeanAggerated[idx], maxdB) * freqLen + i]++;
      pdfMin[binLocation (psdBinReducedMinAggerated[idx], maxdB) * freqLen + i]++;
      pdfMax[binLocation (psdBinReducedMaxAggerated[idx], maxdB) * freqLen + i]++;
      pdfMedian[binLocation (psdBinReducedMedianAggerated[idx], maxdB) * freqLen + i]++;
    }
  }
  /* Calculate the PDF, i.e., probability */
  for (int i = 0; i < PDFBins; i++)
  {
    for (int j = 0; j < freqLen; j++)
    {
      int idx = i * freqLen + j;
      pdfMean[idx] /= totalSegmentsOfOneHour;
      pdfMin[idx] /= totalSegmentsOfOneHour;
      pdfMax[idx] /= totalSegmentsOfOneHour;
      pdfMedian[idx] /= totalSegmentsOfOneHour;
    }
  }
  /* Output PDF */
  FILE *pdf_out = fopen ("pdf_out.txt", "w");
  for (int i = 0; i < PDFBins; i++)
  {
    for (int j = 0; j < freqLen; j++)
    {
      int idx = i * freqLen + j;
      fprintf (pdf_out, "%lf ", pdfMean[idx]);
    }
    fprintf (pdf_out, "\n");
  }
  fclose (pdf_out);
  free (pdfMean);
  free (pdfMin);
  free (pdfMax);
  free (pdfMedian);
  FILE *psd_reduced_out = fopen ("psd_reduced_out.txt", "w");
  for (int i = 0; i < totalSegmentsOfOneHour; i++)
  {
    for (int j = 0; j < freqLen; j++)
    {
      int idx = i * freqLen + j;
      fprintf (psd_reduced_out, "%e %e %e %e %e\n", centerPeriods[j],
               psdBinReducedMeanAggerated[idx],
               psdBinReducedMinAggerated[idx],
               psdBinReducedMaxAggerated[idx],
               psdBinReducedMedianAggerated[idx]);
    }
  }
  fclose (psd_reduced_out);
  free (psdBinReducedMeanAggerated);
  free (psdBinReducedMinAggerated);
  free (psdBinReducedMaxAggerated);
  free (psdBinReducedMedianAggerated);
  free (estimatedFreqs);

  /* Output center periods */
  FILE *center_periods_out = fopen ("center_periods_out.txt", "w");
  for (int i = 0; i < freqLen; i++)
  {
    fprintf (center_periods_out, "%e\n", centerPeriods[i]);
  }
  fclose (center_periods_out);

  /* Free center periods allocated objects */
  free (leftFreqs);
  free (rightFreqs);
  free (centerPeriods);
  /* Free response file allocated objects */
  free (freq);
  free (freqResponse);
  free (poles);
  free (zeros);
  /* Free PSD summary statistics objects */
  free (psdMin);
  free (psdMax);
  free (psdMean);
  free (psdMedian);

  /* Close input miniSEED buffer */
  free (inputmseedBuffer);

  return 0;
}
