#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "cosine_taper.h"
#include "freq_response.h"

#define PI 3.14159265358979

#define OUTPUT_FILENAME "consine_window.m"
#define FREQ_TAPER_OUTPUT_FILENAME "freq_taper_window.m"

static void
range (double *array, double sampleRate, int totalSamples)
{
  double delta         = 1. / sampleRate;
  double totalDuration = delta * totalSamples;
  int i;

  for (i = 0; i < totalSamples; i++)
    array[i] = i / totalDuration;
}

int
main ()
{
  int n       = 51;
  float alpha = 0.5;
  float w[n];

  // check input data boundary

  // check alpha boundary

  int width = (int)floor (alpha * (n - 1) / 2.0);

  // Calculate consine (Tukey) window
  // Caculation reference from here:
  // https://github.com/scipy/scipy/blob/v1.6.3/scipy/signal/windows/windows.py#L795-L875
  int i;
  for (i = 0; i < width + 1; i++)
    w[i] = 0.5 * (1 + cos (PI * (-1 + 2.0 * i / alpha / (n - 1))));
  for (i = width + 1; i < n - width - 1; i++)
    w[i] = 1.0;
  for (i = n - width - 1; i < n; i++)
    w[i] = 0.5 * (1 + cos (PI * (-2.0 / alpha + 1 + 2.0 * i / alpha / (n - 1))));

  FILE *fid = fopen (OUTPUT_FILENAME, "w");
  fprintf (fid, "%% %s: auto-generated file\n\n", OUTPUT_FILENAME);
  fprintf (fid, "clear all;\n");
  fprintf (fid, "close all;\n\n");
  fprintf (fid, "n=%u;\n", n);
  for (i = 0; i < n; i++)
  {
    fprintf (fid, "w(%4u) = %12.4e;\n", i + 1, w[i]);
  }
  fprintf (fid, "nfft=2048;\n");
  fprintf (fid, "W=20*log10(abs(fftshift(fft(w/sum(w),nfft))));\n");
  fprintf (fid, "f=[0:(nfft-1)]/nfft-0.5;\n");
  fprintf (fid, "t=0:(n-1);\n");
  fprintf (fid, "figure;\n");
  fprintf (fid, "subplot(2,1,1);\n");
  fprintf (fid, "  plot(t,w,'Color',[0 0.25 0.5],'LineWidth',2);\n");
  fprintf (fid, "  grid on;\n");
  fprintf (fid, "  xlabel('sample index');\n");
  fprintf (fid, "  ylabel('window');\n");
  fprintf (fid, "  axis([0 n-1 -0.1 1.1]);\n");
  fprintf (fid, "subplot(2,1,2);\n");
  fprintf (fid, "  plot(f,W,'Color',[0 0.5 0.25],'LineWidth',2);\n");
  fprintf (fid, "  grid on;\n");
  fprintf (fid, "  xlabel('normalized frequency');\n");
  fprintf (fid, "  ylabel('PSD [dB]');\n");
  fprintf (fid, "  axis([-0.5 0.5 -140 20]);\n");
  fprintf (fid, "title(['Cosine window']);\n");
  fclose (fid);

  int taper_length = 18000;
  double freq[taper_length];
  double taper[taper_length];
  double sampling_rate = 20.0f;

  range (freq, sampling_rate, taper_length);
  sacCosineTaper (freq, taper_length, 1., 2., 8., 9., sampling_rate, taper);

  FILE *fid1 = fopen (FREQ_TAPER_OUTPUT_FILENAME, "w");
  fprintf (fid1, "%% %s: auto-generated file\n\n", FREQ_TAPER_OUTPUT_FILENAME);
  fprintf (fid1, "clear all;\n");
  fprintf (fid1, "close all;\n\n");
  fprintf (fid1, "n=%u;\n", taper_length);
  for (i = 0; i < taper_length; i++)
  {
    fprintf (fid1, "w(%4u) = %12.4e;\n", i + 1, taper[i]);
  }
  fprintf (fid1, "t=0:(n-1);\n");
  fprintf (fid1, "figure;\n");
  fprintf (fid1, "subplot(1,1,1);\n");
  fprintf (fid1, "  plot(t,w,'Color',[0 0.25 0.5],'LineWidth',2);\n");
  fprintf (fid1, "  grid on;\n");
  fprintf (fid1, "  xlabel('sample index');\n");
  fprintf (fid1, "  ylabel('window');\n");
  fprintf (fid1, "  axis([0 n-1 -0.1 1.1]);\n");
  fclose (fid1);
  free (taper);

  return 0;
}
