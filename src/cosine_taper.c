#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define PI 3.14159265358979

void
cosineTaper (float *data, int n, float alpha, float *tapered)
{
  float *w = (float *)malloc (sizeof (float) * n);

  // check input data boundary

  // check alpha boundary

  int width = (int)floor (alpha * (n - 1) / 2.0);
  int i;

  // Calculate consine (Tukey) window
  // Caculation reference from here:
  // https://github.com/scipy/scipy/blob/v1.6.3/scipy/signal/windows/windows.py#L795-L875
  for (i = 0; i < width + 1; i++)
    w[i] = 0.5 * (1 + cos (PI * (-1 + 2.0 * i / alpha / (n - 1))));
  for (i = width + 1; i < n - width - 1; i++)
    w[i] = 1.0;
  for (i = n - width - 1; i < n; i++)
    w[i] = 0.5 * (1 + cos (PI * (-2.0 / alpha + 1 + 2.0 * i / alpha / (n - 1))));

  FILE *out = fopen ("taperprocess.txt", "w");
  for (i = 0; i < n; i++)
  {
    tapered[i] = data[i] * w[i];
    fprintf (out, "%f %f %f\n", data[i], w[i], tapered[i]);
  }
  fclose (out);

  free (w);
}
