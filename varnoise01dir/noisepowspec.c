#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <getdata.h>
#include "processrwnoise.h"
#include "lsperiodigram.h"
#include "gslfits.h"
#include "dipole.h"

#define N_ROW 33
#define N_COL 16

#define PREALLOC 10000
#define MIN_LEN (1024)
#define MAX_LEN (65536)

#define FREQ 119.11597

#define CHUNKS_DIR "/aux/microchunks06/"
#define RWNDIRFILE "/rwn01"
#define PLANCK_DIR "/planck02/"
#define NOISE_DIR "/var_noise01/"
#define SUB_DIR "/varnoise01dir/"

#define TEST_DIR   "/data2/mathew/noisetest00"
#define JANKY_TEST 0

#define FFTLEN 1024

#define VER "01"

#define LOW_DECIMATE_FACTOR 10
#define LF_EST_WIN (120*600)
#define NUM_CHUNKS 10
#define NUM_LF_BANDS 9
#define MAX_LF_GAP 8000

#define STEP_STITCH 1

void PSD(double *in, double *out, int n, int add); // in psdfftw.c

double bandEnds[] = {0.4, 1.0, 2.0, 4.0, 10.0, 30.0};

double LFBins[] = {0.005, 0.01, 0.025, 0.05, 0.1, 0.2, 0.4, 1, 10};

DIRFILE *spider_df;
DIRFILE *rwn_df;
DIRFILE *planck_df;
DIRFILE *noisespec_df;
DIRFILE *test_df;

double* f_array;
double* p_array;
double* weights;
int numAlloced = 0;

double* fit_a_buf;
double* fit_b_buf;
double* fit_c_buf;
double* fit_d_buf;
double* chi_buf;
short* shortBuf;

int nframes;

char spider_dir[1000];
char rwn_dir[1000];
char chunks_dir[1000];
char planck_dir[1000];
char noise_dir[1000];


void USAGE() {
  fprintf(stderr, "noisepowspec [-r] [-y] [<array>] [<col>] [<row>]\n");
  fprintf(stderr, "  -r: remove reaction wheel noise (RWN" VER ")\n"); 
  fprintf(stderr, "  -y: remove YSSN (YSSN17)\n"); 
  fprintf(stderr, "  if a parameter is missed, all values are used\n"); 
  exit(0);
}

void InitChunkStruct(struct ChunkStruct *C) {
  C->n_alloc = PREALLOC;
  C->n = 0;
  C->start = (int *) malloc(PREALLOC * sizeof(int));
  C->end = (int *) malloc(PREALLOC * sizeof(int));
  C->len = (int *) malloc(PREALLOC * sizeof(int));
  C->lf = (int*) malloc(PREALLOC * sizeof(int));
}



void outputPSD(int frame, double psd[FFTLEN/2], LS* lfdata, int a, int r, int c, int numSamps, int file_offset) {
  static int serial = 0;
  static int oldFrame = 1;

  char fieldName[100];
  char rawName[100];


  //write output to file if we are at a new channel
  if (serial != a*10000 + r * 100 + c && serial != 0) {
    //write dirfile
    noisespec_df = gd_open(noise_dir, GD_RDWR | GD_CREAT | GD_PRETTY_PRINT | GD_SIE_ENCODED);
    if(JANKY_TEST){
    	noisespec_df = test_df;
    }
    printf("serial %d being compared to new %d\n", serial, a*10000 + r*100 +c);
    char names[5] = {'a', 'b', 'c', 'd', 'x'};
    double* arrays[5] = {fit_a_buf, fit_b_buf, fit_c_buf, fit_d_buf, chi_buf};
   
    for(int i = 0; i<5; i++){
      //generate the field names for each parameter
      double zero = 0.0;
      double* currArray = arrays[i];
      sprintf(fieldName, "X%dR%02dC%02d_%c_VAR_NOISE01", serial/10000, (serial%10000)/100, serial%100, toupper(names[i]));
      sprintf(rawName, "x%dr%02dc%02d_%c_var_noise01", serial/10000, (serial%10000)/100, serial%100, names[i]);
   
      if(JANKY_TEST){
        gd_add_raw(noisespec_df, rawName, GD_INT16, 20, 0);
      }

      //the calibration coefficients for each field, allows compression to 16 bits
      double cal = 10/32768.0;
      if(i == 1){
	cal = 0.1/32768.0;
      }
      if(i == 4){
	cal = 1/32768.0;
      }
      
      //truncate extremely large values
      for(int j = 0; j<nframes*20; j++){
	if(currArray[j] > 32767.0){
	  currArray[j] = 32767.0;
	}
	shortBuf[j] = (short) (currArray[j]/cal);
      }

      if(JANKY_TEST){  
        const char* lincomName = rawName;
        gd_add_lincom(noisespec_df, fieldName, 1, &lincomName, &cal, &zero, 0);
      }
      //write the data
      int nwritten = gd_putdata(noisespec_df, rawName, file_offset/20,0,nframes,0, GD_INT16, shortBuf);
      printf("writing field %s, len %d\n", rawName, nwritten);
      if(nwritten != nframes*20){
	char* error = gd_error_string(noisespec_df, NULL, 0);
        printf("not all data written in output psd, error %s\n", error);
	free(error);

        fflush(stdout);
      }

      //cleanup memory
      memset(shortBuf, 0, 20*nframes*sizeof(short));
      memset(currArray, 0, 20*nframes*sizeof(double));

    }
    oldFrame = 1;
    gd_close(noisespec_df);
  }


  if(serial != a*10000 + r * 100 + c){
    serial = a*10000 + r * 100 + c;
    printf("setting new serial %d\n", serial);
  }

  if(serial == -1){
    return;
  }

  //if we need more memory, allocate some
  if(numAlloced < FFTLEN + lfdata->nfreqs){
    if(f_array != NULL){
      free(f_array);
      free(p_array);
      free(weights);
    }
    numAlloced = ((FFTLEN + lfdata->nfreqs) * 2);
    printf("allocating %d doubles in output\n", numAlloced);
    f_array = (double*)malloc(numAlloced* sizeof(double));
    p_array = (double*)malloc(numAlloced* sizeof(double));
    weights = (double*)malloc(numAlloced*sizeof(double));
  }

  //construct the array of the power spectrum for fitting
  double fit_a, fit_b, fit_c, fit_d;
  double cov; 
  int index = 0;
  for(int i = 0; i<lfdata->nfreqs; i++){
    //only add low frequency data between 0.05 hz and 0.5 hz
    if(lfdata->freqs[i] > 0.05 && lfdata->freqs[i] < 0.5){
      f_array[index] = lfdata->freqs[i];
      //this factor fixes the power mismatch
      p_array[index] = lfdata->Pn[i] * LOW_DECIMATE_FACTOR/ (FREQ/2);
      weights[index] = 100;
      index++;
    }
  }
  for(int i = 3; i<(int)((FFTLEN * 25)/120.0); i++){
    //add high frequency data up to 25 hz, drop the first 3 bins because they are contaminated by the detrending
    f_array[index] = (double)i*120.0/FFTLEN;
    p_array[index] = psd[i]/(FREQ/2);
    weights[index] = 100;
    index ++;
  }
  
  //fit parameters of the form 
  //a + b x + cx^d
  fitNoiseSpectrum(&fit_a, &fit_b, &fit_c, &fit_d, &cov, f_array, p_array, weights, index);

  if(index > 300 && fit_a < 0.001){
    //write output file to debug because this chunk is probably bad 
    char filename[500];
    sprintf(filename, ".%stesting/fit%d.dat", SUB_DIR, frame+file_offset/20);


    FILE* testFile = fopen(filename, "w");

    printf("filename = %s, file = %p\n", filename, testFile);
    fprintf(testFile, "%lf %lf %lf\n%lf %lf %lf\n%lf %lf %lf\n%lf %lf %lf\n%lf %lf %lf\n%lf %lf %lf\n", 0.0, 0.0, 0.0, fit_a, 0.0, 0.0,fit_b, 0.0, 0.0, fit_c, 0.0, 0.0, fit_d, 0.0, 0.0, cov, 0.0, 0.0);

    printf("fit answer at frame %d found to be %lf, %lf, %lf, %lf, cov %lf, num samples = %d, len of chunk = %d\n", frame+file_offset/20, fit_a, fit_b, fit_c, fit_d, cov, index, numSamps);

    for(int i = 0; i<index; i++){
      fprintf(testFile, "%lf %lf %lf\n", f_array[i], p_array[i], weights[i]);
    }
    fclose(testFile);

    
  }

  //save data in buffer
  //only if it has enough data to have a good fit
  if(index < 300){
    fit_a = fit_a_buf[oldFrame];
    fit_b = fit_b_buf[oldFrame];
    fit_c = fit_c_buf[oldFrame];
    fit_d = fit_d_buf[oldFrame];
    cov = chi_buf[oldFrame];
    printf("rejecting samples at frame %d for being too short, len %d\n", frame+file_offset/20, index);
  }

  //write data into buffer
  int extrapLen = 20*frame - oldFrame;
  fit_a_buf[20*frame] = fit_a;
  fit_b_buf[20*frame] = fit_b;
  fit_c_buf[20*frame] = fit_c;
  fit_d_buf[20*frame] = fit_d;
  chi_buf[20*frame] = cov;
  for(int i = 0; i<extrapLen; i++){
    fit_a_buf[oldFrame+i] = fit_a_buf[oldFrame];
    fit_b_buf[oldFrame+i] = fit_b_buf[oldFrame];
    fit_c_buf[oldFrame+i] = fit_c_buf[oldFrame];
    fit_d_buf[oldFrame+i] = fit_d_buf[oldFrame];
    chi_buf[oldFrame+i] = chi_buf[oldFrame];
  }
  oldFrame = 20*frame;

}

//adds a new chunk to the chunk struct
void addChunk(struct ChunkStruct *C, int start, int len, int lf) {
  int n = C->n;

  //if we need more memory
  if (n+1 >= C->n_alloc) {
    C->n_alloc *= 1.5;
    C->start = (int *)realloc(C->start, C->n_alloc * sizeof(int));
    C->end = (int *)realloc(C->end, C->n_alloc * sizeof(int));
    C->len = (int *)realloc(C->len, C->n_alloc * sizeof(int));
    C->lf = (int*) realloc(C->lf, C->n_alloc * sizeof(int));
  }
  
  //add the chunk particulars
  C->start[n] = start;
  C->end[n] = start + len;
  C->len[n] = len;
  C->lf[n] = lf;
  if(STEP_STITCH){
    C->lf[n] = 1;
  }
  C->n++;
}

int main(int argc, char *argv[]) {
  int array_in=-1, col_in=-1, row_in = -1;
  int a0, a1, c0, c1, r0, r1;
  int array, row, col;
  int tot_len;
  
  int i_arg;
  int i_chunk;
  int i_samp;
  
  int start, end, len, lf;  
  char file_name[128];
  char field_name[128];
  char line_in[256];
  struct ChunkStruct C;
  
  double *data_in; // bolometer time stream
  double *yssn_in; // YSSN time stream
  double *rwn_in; // rwn time stream
  double *planck_in; // planck time stream
  double *stitch_in; //Stepstictch timestream
  double *low_freq_data; //buffer to hold the decimated data
  double* lf_t;
  double* lf_y;
  LS* answer = NULL;

  int nread;
  double fft_out[FFTLEN];
  double dummy_t[FFTLEN];
  double pow_spec[FFTLEN/2];
  int time_bin_chunk_start;
  int add_psd;
  int psds_added;
  int lf_left_edge;
  int lf_right_edge;

  int remove_rwn = 0;
  int remove_yssn = 0;
  
  int file_start = 0;
  int file_end = -1;

  char suffix[12];
  
  FILE *fp;
  
  /********** read command line *************/
  for (i_arg = 1; i_arg < argc; i_arg++) {
    if (argv[i_arg][0] == '-') {
      switch (argv[i_arg][1]) {
        case 'y':
          remove_yssn = 1;
          break;
        case 'r':
          remove_rwn = 1;
          break;
	case 's':
	  file_start = atoi(argv[i_arg+1]);
	  break;
	case 'e':
	  file_end = atoi(argv[i_arg+1]);
	  break;
	default:
          USAGE();
          break;
      }
    }
    if (tolower(argv[i_arg][0]) == 'x') {//specify array
      array_in = atoi(argv[i_arg]+1);
    } else if (tolower(argv[i_arg][0]) == 'c') {//specify column
      col_in = atoi(argv[i_arg]+1);
    } else if (tolower(argv[i_arg][0]) == 'r') {//specify row
      row_in = atoi(argv[i_arg]+1);
    }
  }
   
  if (array_in == 0) USAGE();
  if (array_in > 6) USAGE();
  if (col_in >= N_COL) USAGE();
  if (row_in >= N_ROW) USAGE();
    
  if (array_in < 0) {
    a0 = 1;
    a1 = 7;
  } else {
    a0 = array_in;
    a1 = array_in + 1;
  }
  if (col_in < 0) {
    c0 = 0;
    c1 = N_COL;
  } else {
    c0 = col_in;
    c1 = col_in+1;
  }
  if (row_in < 0) {
    r0 = 0;
    r1 = N_ROW;
  } else {
    r0 = row_in;
    r1 = row_in+1;
  }

  //generate paths to data products
  sprintf(spider_dir, "%s", getenv("SPIDER_DATA_ROOT"));
  sprintf(rwn_dir, "%s%s", spider_dir, RWNDIRFILE);
  sprintf(chunks_dir, "%s%s", spider_dir, CHUNKS_DIR);
  sprintf(planck_dir, "%s%s", spider_dir, PLANCK_DIR);
  sprintf(noise_dir, "%s%s", spider_dir, NOISE_DIR);

  spider_df = gd_open(spider_dir, GD_RDONLY);  
  if(JANKY_TEST){
    test_df = gd_open(TEST_DIR, GD_RDWR | GD_CREAT);
  }
  //determine data range
  nframes = gd_nframes(spider_df);
 
  if(file_end == -1){
    file_end = nframes*20;
  }else{
    nframes = (file_end - file_start)/20;
  }

  printf("start: %d, end %d\n", file_start, file_end);
  
  rwn_df = gd_open(rwn_dir, GD_RDONLY);
  
  planck_df = gd_open(planck_dir, GD_RDONLY);

  //set up arrays
  data_in = (double *) malloc(nframes * 20 * sizeof(double));
  yssn_in = (double *) malloc(nframes * 20 * sizeof(double));
  rwn_in = (double *) malloc(nframes * 20 * sizeof(double));
  planck_in = (double*) malloc(nframes * 20 * sizeof(double));
  stitch_in = (double*) malloc(nframes * 20 * sizeof(double));
  low_freq_data = (double*) malloc((nframes/LOW_DECIMATE_FACTOR + 1) * 20 * sizeof(double));
  lf_t = (double*) malloc((nframes/LOW_DECIMATE_FACTOR +1) * 20 * sizeof(double));
  lf_y = (double*) malloc((nframes/LOW_DECIMATE_FACTOR +1) * 20 * sizeof(double));
  fit_a_buf = (double*) malloc(nframes * 20 * sizeof(double));
  fit_b_buf = (double*) malloc(nframes * 20 * sizeof(double));
  fit_c_buf = (double*) malloc(nframes * 20 * sizeof(double));
  fit_d_buf = (double*) malloc(nframes * 20 * sizeof(double));
  chi_buf = (double*) malloc(nframes * 20 * sizeof(double));
  shortBuf = (short*) malloc(nframes * 20 * sizeof(short));

  for(int i = 0; i<FFTLEN; i++){
    dummy_t[i] = i;
  }

  InitChunkStruct(&C); 

  //loop through detectors, parse microchunk files
  for (array = a0; array < a1; array++) {
    for (col = c0; col < c1; col++) {
      for (row = r0; row < r1; row++) {
        C.n = 0;
        tot_len = 0;
        sprintf(file_name, "%sX%dR%02dC%02d.uc", chunks_dir, array, row, col);
        fp = fopen(file_name, "r");
        if (fp == NULL) {
          printf("could not open %s\n", file_name);
          fflush(stdout);
          continue;
        }
	//read all lines
        while (fgets(line_in, 256, fp)) {
          sscanf(line_in, "%d %d %d %d", &start, &end, &len, &lf);
          if ((len > MIN_LEN) && (len < MAX_LEN) && start > file_start && end < file_end) {
            //add chunk to chunk list
	    addChunk(&C, start-file_start, len, lf);
	    tot_len += len;
          }
        }
        fclose(fp);

	if (remove_yssn) {
      	  //get the dipole from getDipole
          //read it into data_in temporarily, then add it to yssn

          // Read in this array's yssn signal
          sprintf(field_name, "X%dR%02dC%02d_YSSN17", array, row, col);
          nread = gd_getdata(spider_df, field_name, file_start/20, 0, nframes, 0, GD_FLOAT64, yssn_in);
          if (nread < 100) {
            printf("error %s reading %s\n", gd_error_string(spider_df, NULL, 0), field_name);
            fflush(stdout);
          }
	  //read in dipole into a temporary array, add it to yssn
          getDipole(data_in, array, row, col, file_start/20, (file_end-file_start)/20, spider_df);
          for(int i = 0; i<nread; i++){
            yssn_in[i] += data_in[i];
          }
        }

	//read in raw data
        sprintf(field_name, "X%dR%02dC%02d_DCCLEAN08", array, row, col);
        nread = gd_getdata(spider_df, field_name, file_start/20, 0, nframes, 0, GD_FLOAT64, data_in);
	//this part is for testing
	if(JANKY_TEST){
	  sprintf(field_name, "x%dr%02dc%02d_noisetest00", array, row, col);
	  nread = gd_getdata(test_df, field_name, file_start/20, 0, nframes, 0, GD_FLOAT64, data_in);
	  if(nread < 100){
	    printf("problem reading field %s from test directory %s\n", field_name, TEST_DIR);
	  } 
	}

	//read in stepstitch product
	sprintf(field_name, "X%dR%02dC%02d_STEPSTITCH07", array, row, col);
	nread = gd_getdata(spider_df, field_name, file_start/20, 0, nframes, 0, GD_FLOAT64, stitch_in);

	if(nread < 100){
	  printf("problem reading in step stich timestream %s from file %s\n", field_name, spider_dir);
	}

	//remove rescanned planck to get actual noise
	sprintf(field_name, "X%dR%02dC%02d_PLANCK02", array, row, col);
	nread = gd_getdata(planck_df, field_name, file_start/20, 0, nframes, 0, GD_FLOAT64, planck_in);

	if(nread < 100) {
	  printf("problem reading in planck timestream %s, from file %s\n", field_name, planck_dir);
	  fflush(stdout);
	  memset(planck_in, 0, nframes*20*sizeof(double));
	}

	//read in rwn product
        if (remove_rwn) {
          sprintf(field_name, "X%dR%02dC%02d_RWN" VER, array, row, col);        
          nread = gd_getdata(rwn_df, field_name, file_start/20, 0, nframes, 0, GD_FLOAT64, rwn_in);
          if (nread < 100) {
            printf("error reading %s\n", field_name);
            fflush(stdout);
          }
        }
       
	if(JANKY_TEST){
	  memset(planck_in, 0, nread*sizeof(double));
	  memset(rwn_in, 0, nread*sizeof(double));
	  memset(stitch_in, 0, nread*sizeof(double));
	}

	if(STEP_STITCH){
	  for(i_samp = 0; i_samp< nread; i_samp++){
	    data_in[i_samp] = data_in[i_samp] - stitch_in[i_samp];
	  }
	}
        // optionally remove templates.
        if (remove_yssn && remove_rwn) {
          sprintf(suffix, "trend_both");
	  printf("removing yssn and rwn\n");
          for (i_samp = 0; i_samp < nread; i_samp++) {
            data_in[i_samp] = data_in[i_samp] - (yssn_in[i_samp] + rwn_in[i_samp] + planck_in[i_samp]);
          }
        } else if (remove_yssn) {
          sprintf(suffix, "trend_yssn");
          for (i_samp = 0; i_samp < nread; i_samp++) {
            data_in[i_samp] = data_in[i_samp] - (yssn_in[i_samp] +planck_in[i_samp]);
          }
        } else if (remove_rwn) {
          sprintf(suffix, "trend_rwn");
          for (i_samp = 0; i_samp < nread; i_samp++) {
            data_in[i_samp] = data_in[i_samp] - (rwn_in[i_samp] + planck_in[i_samp]);
          }
        } else {
          sprintf(suffix, "trend_raw");
	  for (i_samp = 0; i_samp < nread; i_samp++) {
	    data_in[i_samp] = data_in[i_samp] - (planck_in[i_samp]);
	  }
        }


	//decimate chunks and get the low frequencies
	decimateChunks(&C, LOW_DECIMATE_FACTOR, data_in, low_freq_data, nframes*20);

        printf("%s with %d chunks from %d to %d.  read %d. using %g%%\n", file_name, C.n, file_start, file_end, nread, 100.0*(double)tot_len/(double)nread);
        fflush(stdout);
        
        time_bin_chunk_start = C.start[0];
        
        add_psd = 0;
        psds_added = 0;
	lf_left_edge = 0;
	lf_right_edge = 0;

	//loop through chunks
	int old_lf = C.lf[0];
        int num_samps;
	for (i_chunk = 0; i_chunk < C.n; i_chunk++) {
	
	
          int n_steps;
          int step;
          int i_step;
          int i_f;
          int i_samp;
          n_steps = floor(((double)C.len[i_chunk]/(double)FFTLEN -1.0)*2.0) + 1;
          step = (C.len[i_chunk] - FFTLEN)/n_steps;
          for (i_step = 0; i_step < n_steps; i_step++) {
            for (i_samp = 0; i_samp < FFTLEN; i_samp++) {
              fft_out[i_samp] = data_in[C.start[i_chunk]+i_step*step+i_samp];
	    }
	    detrendChunk(fft_out, dummy_t, 1024, 2); 
            PSD(fft_out, pow_spec, FFTLEN, add_psd);
            add_psd = 1; 
            psds_added++;
          }

	  //this code makes sure to break macrochunks over big gaps
	  if(((C.start[i_chunk] - C.end[i_chunk -1] > MAX_LF_GAP)||(C.start[i_chunk] - time_bin_chunk_start >  LF_EST_WIN))&& i_chunk != 0 &&old_lf > 0){
	    C.lf[i_chunk] *= -1;
	  }


          if ((C.lf[i_chunk] != old_lf) && (abs(C.lf[i_chunk]) != -1* C.lf[i_chunk -1])) {

	    for (i_f = 0; i_f < FFTLEN/2; i_f++) {
              pow_spec[i_f]/=(double)psds_added;
            }
            
	    //make data and time arrays
	    lf_right_edge = i_chunk -1; 
	  
	    int q = 1;
	    while(C.lf[i_chunk-q] == C.lf[i_chunk -1]){
	      q++;
	    }
	    lf_left_edge = i_chunk - q + 1;
	 
	    if(C.lf[lf_left_edge -1] < 0 && C.lf[lf_right_edge] > 0){
	      lf_left_edge --;
	    }

	    //make the low frequency data arrays
	    num_samps = 0;
	    for(int i = lf_left_edge; i<=lf_right_edge; i++){
	      for(int j = 0; j< C.len[i]/LOW_DECIMATE_FACTOR; j++){
	        lf_t[num_samps] = (j*LOW_DECIMATE_FACTOR + C.start[i])/120.0;
	        lf_y[num_samps] = low_freq_data[j + (C.start[i]/LOW_DECIMATE_FACTOR)];
	        num_samps++;
	      }
	    }

	    //if the fill factor is greater that 1/4
	    if(num_samps < (int)(0.25*(C.end[lf_right_edge] - C.start[lf_left_edge])/LOW_DECIMATE_FACTOR)){
	      printf("num_samps = %d and total should be %d\n", num_samps, (int)(0.25*(C.end[lf_right_edge] - C.start[lf_left_edge])/LOW_DECIMATE_FACTOR)); 
	      num_samps = 0;

	    }

	    //detrend the low frequency data
	    detrendChunk(lf_y, lf_t, num_samps, 2);

	    //compute the lomb scargle periodigram
	    answer = periodigram(lf_t, lf_y, num_samps, 1, 1);//last two are the oversampling factor and the highest frequency in units of the niquist freqnency
	    memset(lf_t, 0, (nframes/LOW_DECIMATE_FACTOR +1) * 20 * sizeof(double));
	    memset(lf_y, 0, (nframes/LOW_DECIMATE_FACTOR +1) * 20 * sizeof(double));
	    //output low frequency and high frequency data
	    outputPSD(time_bin_chunk_start/20, pow_spec, answer, array, row, col, num_samps, file_start);
            if(answer != NULL){
	      LSFree(answer);
	    }
	    add_psd = 0;
            psds_added = 0;
            time_bin_chunk_start = C.start[i_chunk];
          }
	  old_lf = C.lf[i_chunk];
        }
      }
    }
  }
  //write the final buffer to disk
  outputPSD(0, NULL, NULL, 0, 0, -1, 0, file_start);
}
