#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define PI 3.14159265358979

/* Cosine taper window with given alpha value */
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

/* Return a SAC-style cosine taper window with given four corner frequencies.
 * Referenced from: https://github.com/obspy/obspy/blob/016b4a325319ab4f98479ae9490695896996eac5/obspy/signal/invsim.py#L149
 * As we are processing seismic wave, this function will also generate the other taper window
 * for negative frequency.
 * That is, the output taper window consists of two cosine taper windows.
 */
int
sacCosineTaper (double *freqs, int n, float f1, float f2, float f3, float f4, double sampling_rate, double **taper)
{
  double nyquist_freq = sampling_rate / 2.0;

  *taper = (double *)malloc (sizeof (double) * n);
  if (*taper == NULL)
  {
    fprintf (stderr, "taper window allocation failed\n");
    return -1;
  }
  int i;
  for (i = 0; i < n; i++)
    (*taper)[i] = 0.0f;

  /* Set taper window */
  for (i = 0; i < n / 2 + 1; i++)
  {
    double temp = freqs[i];

    /* Case 1 */
    if ((f1 <= temp) && (temp <= f2))
      (*taper)[i] = 0.5 * (1.0 - cos (PI * (temp - f1) / (f2 - f1)));
    /* Case 2 */
    else if ((f2 < temp) && (temp < f3))
      (*taper)[i] = 1.0f;
    /* Case 3 */
    else if ((f3 <= temp) && (temp <= f4))
      (*taper)[i] = 0.5 * (1.0 + cos (PI * (temp - f3) / (f4 - f3)));
  }
  /* Set taper window of negative frequencies */
#if 0
  f1 += nyquist_freq;
  f2 += nyquist_freq;
  f3 += nyquist_freq;
  f4 += nyquist_freq;
  for (i = n / 2; i < n; i++)
  {
    double temp = freqs[i];

    /* Case 1 */
    if ((f1 <= temp) && (temp <= f2))
      (*taper)[i] = 0.5 * (1.0 - cos (PI * (temp - f1) / (f2 - f1)));
    /* Case 2 */
    else if ((f2 < temp) && (temp < f3))
      (*taper)[i] = 1.0f;
    /* Case 3 */
    else if ((f3 <= temp) && (temp <= f4))
      (*taper)[i] = 0.5 * (1.0 + cos (PI * (temp - f3) / (f4 - f3)));
  }
#endif
  for (i = n / 2 + 1; i < n; i++)
  {
    (*taper)[i] = (*taper)[n - i];
  }

  return 0;
}
