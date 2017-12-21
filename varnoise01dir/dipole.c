//dipole.c
//Gets the dipole for a given detector

#include "dipole.h"

//Reads from a field and ensures you get the correct amount of data
void verifyRead(char *field, int nread, int nreq) {
    if (nreq - nread > 10000) {
      printf("fatal error reading %s (read %d instead of %d) \n", field, nread, nreq);
      fflush(stdout);
      exit(1);
    } else if (abs(nread- nreq) > 100) {
      printf("warning: short read for %s (read %d instead of %d) \n", field, nread, nreq);
      fflush(stdout);
    }
}

//Returns the dipole contribution from the spider dirfile
void getDipole(double *dipole_out, int array_in, int row_in, int col_in, int f0, int nframes, DIRFILE* spider_df) {

  char tools_dir[1000];
  char centroids_file[1000];
  char gain_file[1000];

  sprintf(tools_dir, "%s", getenv("SPIDER_TOOLS_PATH"));
  sprintf(centroids_file, "%s%s", tools_dir, CENTROIDSFILE);
  sprintf(gain_file, "%s%s", tools_dir, GAINSFILE);

  static quat_t *bore_q = NULL;
  static qp_memory_t *mem = NULL;
  static double *time = NULL;
  static double d_az[7][N_ROW][N_COL];
  static double d_el[7][N_ROW][N_COL];
  static double gain[7][N_ROW][N_COL];

  int nread;

  int nsamp;

  nsamp = 20*nframes;

  if (bore_q == NULL) {
    // **************************************************
    // read pointing offsets
    char linein[256];
    FILE *fp = fopen(centroids_file, "r");
    if (fp == NULL) {
      fprintf(stderr, "fatal: could not open centroids file %s\n", centroids_file);
      exit(0);
    }

    while (fgets(linein, 255,fp)) {
      int array, row, col;
      char str1[128], str2[128];
      double x, y;

      if (linein[0] == 'x') {
        sscanf(linein,"x%1dr%2dc%2d %s %s", &array, &row, &col, str1, str2);

        if ((array < 1) || (array > 6)) {
          fprintf(stderr, "error parsing %s\n%s\n", centroids_file, linein);
          exit(0);
        }
        if (tolower(str1[0]) == 'n') {
          x = 0;
          y = 0;
        } else {
          x = atof(str1);
          y = atof(str2);
        }
        d_az[array][row][col] = x;
        d_el[array][row][col] = y;
      }
    }
    fclose(fp);

    // **************************************************
    // read gains
    fp = fopen(gain_file, "r");
    if (fp == NULL) {
      fprintf(stderr, "fatal: could not open gains file %s\n", gain_file);
      exit(0);
    }

    double sum_gains[7][4];
    int num_gains[7][4];
    for (int i=0; i<7; i++) {
      for (int j=0; j<4; j++) {
        sum_gains[i][j] = 0;
        num_gains[i][j] = 0;
      }
    }

    while (fgets(linein, 255,fp)) {
      int array, row, col;
      char str1[128], str2[128], str3[128];
      double x, y;

      if (linein[0] == 'x') {
        sscanf(linein,"x%1dr%2dc%2d %s %s %s", &array, &row, &col, str1, str2, str3);

        if ((array < 1) || (array > 6)) {
          fprintf(stderr, "error parsing %s\n%s\n", gain_file, linein);
          exit(0);
        }
        if (tolower(str1[0]) == 'n') {
          x = 0;
          y = 1000;
        } else {
          x = atof(str1);
          y = atof(str3);
        }
        gain[array][row][col] = x;
        if (y<MAX_FITRES) {
          sum_gains[array][col/4] += x;
          num_gains[array][col/4]++;
        }
      }
    }
    fclose(fp);

    // **************************************************
    // clean up gains
    for (int array = 1; array < 7; array++) {
      for (int col = 0; col < N_COL; col++) {
        for (int row = 0; row < N_ROW; row++) {
          if (gain[array][row][col] == 0) {
            gain[array][row][col] = sum_gains[array][col/4]/(double)num_gains[array][col/4];
          }
        }
      }
    }

    // **************************************************
    // read time.
    double *time_in;

    time = (double *) malloc(nsamp * sizeof(double));
    time_in = (double *) malloc(nframes * sizeof(double));

    nread = gd_getdata(spider_df, "TIME_TIME02", f0, 0, nframes, 0, GD_FLOAT64, time_in);
    verifyRead("TIME_TIME02", nread, nframes);

    for (int i_frame = 0; i_frame < nframes-1; i_frame++) {
      for (int j = 0; j<20; j++) {
        time[i_frame*20+j] = time_in[i_frame] + (double)j/20.0 * (time_in[i_frame+1]-time_in[i_frame]);
      }
    }
    free(time_in);

    // **************************************************
    // generate pointing quaternions
    mem = qp_init_memory();
    qp_set_opt_fast_math(mem, 1);

    bore_q = (quat_t *) calloc(nsamp, sizeof(quat_t));

    double *ra;
    double *dec;
    double *phi;
    ra = (double *)malloc(DIPOLECHUNK * sizeof(double));
    dec = (double *)malloc(DIPOLECHUNK * sizeof(double));
    phi = (double *)malloc(DIPOLECHUNK * sizeof(double));

    // read in quaternions
    //printf("Caculating quaternions\n");
    for (int i = 0; i<nsamp; i+=DIPOLECHUNK) {
      //printf("."); fflush(stdout);
      int n = DIPOLECHUNK;
      if (i+n > nsamp) n = nsamp - i;

      nread = gd_getdata(spider_df, "RA_"POINT_PROD, 0, i+20*f0, 0, n, GD_FLOAT64, ra);
      verifyRead("RA_"POINT_PROD, nread, n);
      nread = gd_getdata(spider_df, "DEC_"POINT_PROD, 0, i+20*f0, 0, n, GD_FLOAT64, dec);
      verifyRead("DEC_"POINT_PROD, nread, n);
      nread = gd_getdata(spider_df, "PHI_"POINT_PROD, 0, i+20*f0, 0, n, GD_FLOAT64, phi);
      verifyRead("PHI_"POINT_PROD, nread, n);

      for (int j = 0; j<nread; j++) {
        ra[j] *= 15.0;
        if (phi[j]>180.0) phi[j] -= 180.0;
      }
      //printf("+"); fflush(stdout);

      // convert boresight pointing to quaternions
      qp_radecpa2quatn(mem, ra, dec, phi, bore_q+i, n);
    }
    free(phi);
    free(dec);
    free(ra);
    //printf("\n");    
  } // end first time product generation.


  // **************************************************
  // compute channel offset quaternion
  quat_t det_q;
  qp_det_offset(d_az[array_in][row_in][col_in], d_el[array_in][row_in][col_in], 0., det_q);

  // compute dipole in *K* (need to multiply by 1e6/cal per channel to get ADU)
  qp_bore2dipole(mem, det_q, time, bore_q, dipole_out, nsamp);
  for (int i = 0; i<nsamp; i++) {
    dipole_out[i] = dipole_out[i]*1e6/gain[array_in][row_in][col_in];
  }
}
