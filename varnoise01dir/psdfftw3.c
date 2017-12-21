#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <math.h>


/** calculate the psd and copy it to out.
 * out must pre-allocated and is n/2 long.
 * if (add==1) then add the spectrum to out
 * rather than replacing it (for averaging psds.
 * If the input buffer  is reused, then a new plan does not need
 * to be made.  If a new buffer is used, then the old plan
 * is destroyed and a new one is made.
**/
void PSD(double *in, double *out, int n, int add) {
  static double *last_in = 0;
  static fftw_plan plan;
  fftw_complex out_c[n/2];
  int n_out, i;

  if (last_in != in) {
    double *tmp_in;
    tmp_in = (double*)malloc(n*sizeof(double));

    for (i=0; i<n; i++) tmp_in[i] = in[i]; // FFTW_MEASURE overwrites 'in', so cache it.
    
    if (last_in !=0) {
      fftw_destroy_plan(plan);
    }
    last_in = in;
    fprintf(stderr, "making plan\n");
    plan = fftw_plan_dft_r2c_1d(n, in, out_c, FFTW_MEASURE);

    for (i=0; i<n; i++) in[i] = tmp_in[i];
    free(tmp_in);
  }
  
  fftw_execute(plan);

  n_out = n/2;
  

  if (add) {
    for (i=0; i<n_out; i++) {
      out[i] += (out_c[i][0]*out_c[i][0] + out_c[i][1]*out_c[i][1])/(double)(n);
    }    
  } else {
    for (i=0; i<n_out; i++) {
      out[i] = (out_c[i][0]*out_c[i][0] + out_c[i][1]*out_c[i][1])/(double)(n);
    }    
  }  
}

