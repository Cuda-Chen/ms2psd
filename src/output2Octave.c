#include <stdio.h>
#include <stdlib.h>

#include "output2Octave.h"

void
output2Octave (const char *outputfile, int nfft, float *psd)
{
  FILE *fid = fopen (outputfile, "w");
  fprintf (fid, "%% %s : auto-generated file\n", outputfile);
  fprintf (fid, "clear all;\n");
  fprintf (fid, "close all;\n\n");
  fprintf (fid, "nfft = %u;\n", nfft);
  fprintf (fid, "f    = [0:(nfft-1)]/nfft - 0.5;\n");
  fprintf (fid, "psd  = zeros(1,nfft);\n");

  int i;
  for (i = 0; i < nfft; i++)
    fprintf (fid, "psd(%6u) = %12.4e;\n", i + 1, psd[i]);

  fprintf (fid, "figure;\n");
  fprintf (fid, "plot(f, psd, '-', 'LineWidth',1.5);\n");
  fprintf (fid, "xlabel('Normalized Frequency [f/F_s]');\n");
  fprintf (fid, "ylabel('Power Spectral Density [dB]');\n");
  fprintf (fid, "grid on;\n");
  fprintf (fid, "axis([-0.5 0.5 ymin ymax]);\n");

  fclose (fid);
}
