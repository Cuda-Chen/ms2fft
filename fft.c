#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>

#include "fftw3.h"

#include "fft.h"

void
fftToFile (double *data, uint64_t dataSamples, FILE *fptr)
{
  fftw_plan fft, ifft;
  uint64_t i;
  /* allocate memory */
  fftw_complex *in  = (fftw_complex *)fftw_malloc (sizeof (fftw_complex) * dataSamples);
  fftw_complex *out = (fftw_complex *)fftw_malloc (sizeof (fftw_complex) * dataSamples);
  fftw_complex *ref = (fftw_complex *)fftw_malloc (sizeof (fftw_complex) * dataSamples);

  /* prepare input data */
  for (i = 0; i < dataSamples; i++)
  {
    in[i][0] = data[i];
    in[i][1] = 0;
  }

  /* Fourier transform and save result to `out` */
  fft = fftw_plan_dft_1d (dataSamples, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute (fft);
  fftw_destroy_plan (fft);

  /* inverse Fourier transform and save result to `ref` */
  printf ("\ninverse transform:\n");
  ifft = fftw_plan_dft_1d (dataSamples, out, ref, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute (ifft);
  /* normalize */
  for (i = 0; i < dataSamples; i++)
  {
    ref[i][0] *= 1. / dataSamples;
    ref[i][1] *= 1. / dataSamples;
  }
  for (i = 0; i < dataSamples; i++)
  {
    printf ("recover: %" PRId64 " %+9.5f %+9.5f I v.s. %+9.5f %+9.5f I\n",
            i, in[i][0], in[i][1], ref[i][0], ref[i][1]);
    fprintf (fptr, "%lf,%lf\r\n", out[i][0], out[i][1]);
  }
  fftw_destroy_plan (ifft);
}
