#ifndef LS_PERIODIGRAM_H
#define LS_PERIODIGRAM_H

//a file which contains functions used to compute the ls periodigram

#include <complex.h>
#include <nfft3.h>
#include "processrwnoise.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
	double* freqs; // (>0) frequencies
	double* Pn;    //periodigram ordinates
	int nfreqs;    //number of frequencies
} LS;

void decimateChunks(struct ChunkStruct* C, int factor, double* data_in, double* data_out, int len);

void nfft(const double* t, const double* y, int n, int m, double complex* d);

LS* periodigram(const double* t, const double* y, int npts, double over, double hifac);

double centerData(int n, double* y, double* var);

void meanAndVarience(int n, const double* y, double* mean, double* varience);

void reduce(double* t, int npts, double oversampling);

double square(double a);

double sign(double a, double b);

void LSFree(LS* ls);

#endif
