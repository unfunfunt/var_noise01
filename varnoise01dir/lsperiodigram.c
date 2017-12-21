//a file for computing the ls periodigram
//
#include "lsperiodigram.h"

//reduces sampling rate by given factor using averaging
//operates over the good data given by the chunkstruct
void decimateChunks(struct ChunkStruct* C, int factor, double* data_in, double* data_out, int len){

	memset(data_out, 0, len/factor*sizeof(double));

	//loop through good data chunks
	for(int i = 0; i<C->n; i++){
		//loop through output data array
		for(int j = C->start[i]/factor; j<C->end[i]/factor; j++){
			//find front and back edges of the bin
			int lower = (j-1)*factor;
			if(lower < C->start[i]){
				lower = C->start[i];
			}
			int upper = (j+1)*factor;
			if(upper > C->end[i]){
				upper = C->end[i];
			}
			//average data togeather
			data_out[j] = 0;
			for(int k = lower; k<=upper; k++){
				data_out[j] += data_in[k]/(upper - lower+1);
			}
		}
	}
}

//perform nfft on given data
void nfft(const double* t, const double* y, int n, int m, double complex* d){
	nfft_plan p;
	nfft_init_1d(&p, 2 * m, n);
	if (y)	/* spectrum */
	{
		for (int i = 0; i < n; i++)
		{
			p.x[i] = t[i];
			p.f[i] = y[i];
		}
	}
	else	/* window */
	{
		for (int i = 0; i < n; i++)
		{
			p.x[i] = t[i];
			p.f[i] = 1.0;
		}
	}

	if (p.flags & PRE_ONE_PSI)
		nfft_precompute_one_psi(&p);

	nfft_adjoint(&p);

	for (int i = 0; i < m; i++)
		d[i] = p.f_hat[i + m];

	d[m] = conj(p.f_hat[0]);
	nfft_finalize(&p);
}

inline double square(double a){
	return a*a;
}

inline double sign(double a, double b){
	return ((b >= 0) ? 1 : -1) * fabs(a);
}

//allocates the lomb-scargle struct
LS* LSalloc(int n)
{
	LS* ls = malloc(sizeof(LS));
	ls->freqs = (double*)malloc(n * sizeof(double));
	ls->Pn = (double*)malloc(n * sizeof(double));
	ls->nfreqs = n;
	return ls;
}

//frees lomb-scargle struct
void LSFree(LS* ls){
	free(ls->freqs);
	free(ls->Pn);
	free(ls);
}

//reduces sample range to -1/2, 1/2
void reduce(double* t, int npts, double oversampling)
{
	double tmax = t[npts - 1];
	double tmin = t[0];

	double range = oversampling * (tmax - tmin);

	const double eps_border = 1e-5;
	const double a = 0.5 - eps_border;

	// Reduce to [-1/2, 1/2[
	for (int i = 0; i < npts; i++)
		t[i] = 2 * a * (t[i] - tmin) / range - a;
}

//computes mean and varience of the data
void meanAndVariance(int n, const double* y, double* mean, double* variance)
{
	*mean = 0;
	double M2 = 0;

	int nn = 1;
	for (int i = 0; i < n; i++, nn++){
		double delta = y[i] - *mean;
		*mean +=  delta / nn;
		M2 += delta * (y[i] - *mean);
	}
	*variance = M2 / (n - 1);
}

//removes the mean from the data, returns the variance
double centerData(int n, double *y, double *var)
{
	double average = 0.0;
	meanAndVariance(n, y, &average, var);

	for (int i = 0; i < n; i++)
		y[i] -= average;

	return average;
}

//computes a lomb-scargle periodigram of the given data y(t)
//taken from https://www.aanda.org/articles/aa/pdf/2012/09/aa19076-12.pdf
LS* periodigram(const double* tobs, const double* yobs, int npts, double over, double hifac){

  if(npts <= 0){
    return LSalloc(0);
  }

  //Clones the data.
  double* t = malloc(npts * sizeof(double));
  memcpy(t, tobs, npts * sizeof(double));
  double* y = malloc(npts * sizeof(double));
  memcpy(y, yobs, npts * sizeof(double));
  
  // Centers the data.
  double var;
  centerData(npts, y, &var);
  
  // Determines the Nyquist frequency for
  // a uniform sampling and the frequency
  // resolution.
  double T = t[npts - 1] - t[0];
  double df = 1.0 / (over * T);
  
  // Determines the highest frequency, fmax,
  // and the corresponding index, m, in
  // the positive part of the spectrum. 
 
  int m = (int)floor(0.5 * npts * over * hifac);
  
  // Allocates space for the output frequencies
  // and values of the periodogram.
  int nfreqs = m;
  LS* ls = LSalloc(nfreqs);
                                                       
  // Reduces the times to [-1/2, 1/2).
  reduce(t, npts, over);

  // Unnormalised Fourier transform of the data.
  double complex* sp = malloc((m + 1) * sizeof(double complex));
  nfft(t, y, npts, m, sp);
 
  // Unnormalised Fourier transform of the window.
  double complex* win = malloc((2 * m + 1) * sizeof(double complex));
  nfft(t, NULL, npts, 2 * m, win);
                                                                             
  double pmax = -1.0;
  for (int j = 1; j <= m; j++)
  {
  	double complex z1 = sp[j];
  	double complex z2 = win[2 * j];
  	double hypo = cabs(z2);
  	double hc2wt = 0.5 * cimag(z2) / hypo;
  	double hs2wt = 0.5 * creal(z2) / hypo;
  	double cwt = sqrt(0.5 + hc2wt);
  	double swt = sign(sqrt(0.5 - hc2wt), hs2wt);
  	double den = 0.5*npts + hc2wt * creal(z2) + hs2wt * cimag(z2);
  	double cterm = square(cwt * creal(z1) + swt * cimag(z1)) / den;
  	double sterm = square(cwt * cimag(z1) - swt * creal(z1)) / (npts - den);
  
  	int j1 = j - 1;
  	ls->freqs[j1] = j * df;
  	ls->Pn[j1] = (cterm + sterm) / (2);
  	
  	if (ls->Pn[j1] > pmax)
  	{
  		pmax = ls->Pn[j1];
  	}
  }
  
  free(win);
  free(sp);
  free(y);
  free(t);
  return ls;
}
