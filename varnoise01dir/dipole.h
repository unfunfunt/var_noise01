//dipole.h
//computes a dipole timestream for a given detector
//
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <getdata.h>
#include <math.h>

#include <qpoint.h>

#define N_ROW 33
#define N_COL 16 

#define MAX_FITRES 0.1
#define DIPOLECHUNK 10000000

#define CENTROIDSFILE "/analysis/config/bolotables/centroids_deproject02.txt"
#define GAINSFILE "/analysis/config/bolotables/cal.txt"
#define POINT_PROD "POINT06"

void verifyRead(char* field, int nread, int nreq);

void getDipole(double *dipole_out, int array_in, int row_in, int col_in, int f0, int nframes, DIRFILE* spider_df);
