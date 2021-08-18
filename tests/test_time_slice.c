#include <assert.h>
#include <stdio.h>

#include "libmseed.h"

#include "process_trace.h"

static nstime_t NSECS = 1000000000;

/* 1-hour long segment properties, 50% overlap*/
static int lengthOfOneHour  = 3600;
static int overlapOfOneHour = 50;

/* 900 seconds, or 15-minute long segment properties, 75% overlap */
static int lengthOfSegment  = 900;
static int overlapOfSegment = 75;

/* PDF properties */
static int mindB = -200;
static int maxdB = -50;

const char *mseedfile = "../zNACB/NACB.TW..HHZ.2021.202";

int
main ()
{
  int rv;
  char *sidpattern   = "*";
  uint8_t pubversion = 0;
  double sampleRate;

  /* Get the start time and end time of this trace */
  nstime_t starttimeOfTrace;
  nstime_t endtimeOfTrace;
  int totalSegmentsOfOneHour = 0;
  rv                         = getTraceProperties (mseedfile, &starttimeOfTrace, &endtimeOfTrace, &sampleRate, &totalSegmentsOfOneHour);
  if (rv != 0)
  {
    fprintf (stderr, "Cannot open input miniSEED for getting start time and end time\n");
    return -1;
  }

  int nextTimeStamp;
  nstime_t nextTimeStampNS;
  nstime_t starttime;
  nstime_t endtime;
  nstime_t endtimeTick;

  /* Split trace to 1-hour long segment with 50% overlapping
   * for reducing processing time */
  nextTimeStamp                    = lengthOfOneHour - (lengthOfOneHour * overlapOfOneHour / 100);
  nextTimeStampNS                  = nextTimeStamp * NSECS;
  starttime                        = starttimeOfTrace;
  endtime                          = starttime + ((nstime_t)lengthOfOneHour * NSECS);
  MS3Selections *oneHourSelections = NULL;
  while (endtime <= endtimeOfTrace)
  {
    rv = ms3_addselect (&oneHourSelections, sidpattern, starttime, endtime, pubversion);
    starttime += nextTimeStampNS;
    endtime += nextTimeStampNS;
  }
  ms3_printselections (oneHourSelections);

  /* Split 1-hour long segment to 15-minute long segment
   * with 75% overlapping for reducing data variance */
  nextTimeStamp                    = lengthOfSegment - (lengthOfSegment * overlapOfSegment / 100);
  nextTimeStampNS                  = nextTimeStamp * NSECS;
  MS3Selections *segmentSelections = NULL;
  MS3SelectTime *segmentTicks      = oneHourSelections->timewindows;
  int count                        = 0;
  while (segmentTicks)
  {
    starttime = segmentTicks->starttime;
    endtime   = starttime + ((nstime_t)lengthOfSegment * NSECS);
    char starttimestr[30];
    char endtimestr[30];
    int segments = 0;
    while (endtime <= segmentTicks->endtime)
    {
      //rv = ms3_addselect (&segmentSelections, sidpattern, starttime, endtime, pubversion);
      ms_nstime2timestr (starttime, starttimestr, ISOMONTHDAY, NANO);
      ms_nstime2timestr (endtime, endtimestr, ISOMONTHDAY, NANO);
      printf ("%s - %s\n", starttimestr, endtimestr);
      starttime += nextTimeStampNS;
      endtime += nextTimeStampNS;
      segments++;
    }
    ms3_printselections (segmentSelections);
    printf ("%d\n", segments);
    printf ("========================\n");

    if (segmentSelections)
      ms3_freeselections (segmentSelections);
    segmentTicks = segmentTicks->next;

    count++;
  }
  if (oneHourSelections)
    ms3_freeselections (oneHourSelections);

  printf ("%d %d\n", count, totalSegmentsOfOneHour);

  return 0;
}
