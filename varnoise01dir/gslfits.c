#include "gslfits.h"

gsl_matrix* fit_x;
gsl_vector* fit_y;
gsl_vector* fit_c;
gsl_matrix* fit_cov;
gsl_vector* fit_w;
gsl_multifit_linear_workspace* fit_work;

gsl_matrix* noise_covar;
gsl_rng* rng;
gsl_multifit_nlinear_workspace* noise_work;

int noiseLen = 0;

//allocs the Gsl library for a fit of the given length and order
void allocGsl(int len, int order){

	static int gsl_len = 0;
	static int gsl_order = 0;

	//old fit is not the right size
	if(gsl_len != len || gsl_order != order){
		if(fit_x != NULL){
			gsl_matrix_free(fit_x);
			gsl_matrix_free(fit_cov);
			gsl_vector_free(fit_y);
			gsl_vector_free(fit_c);
			gsl_multifit_linear_free(fit_work);
			gsl_vector_free(fit_w);
		}
		//printf("initing gsl memory\n");
		fit_x = gsl_matrix_alloc(len, order+1);
		fit_y = gsl_vector_alloc(len);
		fit_w = gsl_vector_alloc(len);

		fit_c = gsl_vector_alloc(order+1);
		fit_cov = gsl_matrix_alloc(order+1, order+1);

		fit_work = gsl_multifit_linear_alloc(len, order+1);


		gsl_len = len;
		gsl_order = order;

	}
}

//removes a poly of the given order from the provided data
void detrendChunk(double* data_in, double* t_in, int len, int order){
	double chisq;	

	if(len <= 0){
	  return;
	}

	allocGsl(len, order);

	//setup matrices for inversion
	for(int i = 0; i<len; i++){
		for(int j = 0; j<= order; j++){
			gsl_matrix_set(fit_x, i, j, pow(t_in[i], j));
		}
	  
		gsl_vector_set(fit_y, i, data_in[i]);
		gsl_vector_set(fit_w, i, 1);
	}

	gsl_multifit_wlinear(fit_x, fit_w, fit_y, fit_c, fit_cov, &chisq, fit_work);

	//fit cooeficients are in c[0] + x* c[1] + x^2* c[2]...

	//detrends the chun by subtracting the fit
	for(int i = 0; i<len; i++){
		for(int j = 0; j<=order; j++){
			data_in[i] -= gsl_vector_get(fit_c, j) * pow(t_in[i], j);
		}
	}
	//cleans up
	gsl_matrix_set_zero(fit_x);
	gsl_matrix_set_zero(fit_cov);
	  
	gsl_vector_set_zero(fit_y);
	gsl_vector_set_zero(fit_c);
}

//allocs the noise fit for a chunk of given length
void allocNoiseFit(int len){
	//old fit is wrong length
	if(len != noiseLen){
		if(noise_covar != NULL){
			gsl_matrix_free(noise_covar);
			gsl_rng_free(rng);
			gsl_multifit_nlinear_free(noise_work);
		}
		noise_covar = gsl_matrix_alloc(4, 4);
		gsl_rng_env_setup();
		rng = gsl_rng_alloc(gsl_rng_default);
		gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters(); 
		noise_work = gsl_multifit_nlinear_alloc(gsl_multifit_nlinear_trust, &fdf_params, len, 4);
		noiseLen = len;
	}
}

//evalutes the noise fit model at position x with given data, sets the residual in f
int noise_f(const gsl_vector* x, void* data, gsl_vector* f){

	noisefit params = *(noisefit*)data;
	int n = params.len;
	double* xdata = params.x_data;
	double* ydata = params.y_data;
	double a = gsl_vector_get(x, 0);
	double b = gsl_vector_get(x, 1);
	double c = gsl_vector_get(x, 2);
	double d = gsl_vector_get(x, 3);

	for(int i = 0; i< n; i++){
		//evaluates fit at each sample
		double Yi = a + b * xdata[i] + c * pow(xdata[i], d);
		//sets huge weight if a<0, as this is non physical
		if(a <= 0){
			Yi += 1e30;
		}
		gsl_vector_set(f, i, Yi - ydata[i]);
	}
	return GSL_SUCCESS;
}

//evalutes the derivative of the noise fitting function
int noise_df(const gsl_vector* x, void* data, gsl_matrix* j){

	noisefit params = *(noisefit*)data;
	int n = params.len;
	double* xdata = params.x_data;
	double c = gsl_vector_get(x, 2);
	double d = gsl_vector_get(x, 3);


	//evaluates the matrix of derivatives
	for(int i = 0; i<n; i++){
		gsl_matrix_set(j, i, 0, 1.0);
		gsl_matrix_set(j, i, 1, xdata[i]);
		gsl_matrix_set(j, i, 2, pow(xdata[i], d));
		gsl_matrix_set(j, i, 3, c * log(xdata[i]) * pow(xdata[i], d));
	}

	return GSL_SUCCESS;
}

//fits a function a + b x + cx^d
void fitNoiseSpectrum(double* a, double* b, double* c, double* d, double* chisq, double* t, double* x, double* weights, int len){

	allocNoiseFit(len);

	//initial parameter guesses for a, b, c, d
	double init_vals[4] = {0.1, 0.001, 0, -2};
	
	gsl_multifit_nlinear_fdf fdf;

	noisefit data;
	data.x_data = t;
	data.y_data = x;
	data.len = len;

	fdf.f = noise_f;
	fdf.df = noise_df;
	fdf.fvv = NULL;
	fdf.n = len;
	fdf.p = 4;
	fdf.params = &data;

	gsl_vector_view xview = gsl_vector_view_array(init_vals, 4);
	gsl_vector_view wts = gsl_vector_view_array(weights, len);

	gsl_multifit_nlinear_winit(&xview.vector, &wts.vector, &fdf, noise_work);
	int info;
	gsl_multifit_nlinear_driver(20, 1e-10, 1e-10, 0, NULL, NULL, &info, noise_work);

	*a = gsl_vector_get(noise_work->x, 0);
	*b = gsl_vector_get(noise_work->x, 1);
	*c = gsl_vector_get(noise_work->x, 2);
	*d = gsl_vector_get(noise_work->x, 3);

	*chisq = computeLogChi(*a, *b, *c, *d, t, x, len);
	//clear workspace
}

//computes a goodness of fit function in log space
double computeLogChi(double a, double b, double c, double d, double* t, double* x, int len){

	double sum = 0;
	for(int i = 0; i<len; i++){
		double eval = a+b*t[i] + c*pow(t[i], d);

		sum += pow((fabs(x[i] - eval)), 2);
	}
	return log(sum)/(double)len;
}
