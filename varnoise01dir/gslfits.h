//Gslfits.h 
//a file that does fits using the gsl

#ifndef GSL_FITS_H
#define GSL_FITS_H

#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_blas.h>

typedef struct {
	double* y_data;
	double* x_data;
	int len;
} noisefit;

void detrendChunk(double* data_in, double* t_in, int len, int order);

void allocGsl(int len, int order);

void fitNoiseSpectrum(double* a, double* b, double* c, double* d, double* cov, double* t_array, double* x_array, double* weights, int len);

void allocNoiseFit(int len);

int noise_f(const gsl_vector* x, void* data, gsl_vector* f);

int noise_df(const gsl_vector* x, void* data, gsl_matrix* j);

double computeLogChi(double a, double b, double c, double d, double* t_array, double* x_array, int len);
#endif
