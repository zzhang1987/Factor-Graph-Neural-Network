/* #include 	"dbmalloc/malloc.h" */

/* #include "smartall.h" This piece of shit didn't work for me */
/*
   These are the NR routines (modified by me) that are in the 
   file nrutil.c

   Other NR routines are specified in mynr.h
   See also macopt.h 
*/
float *vector(int,int);
float **matrix(int,int,int,int);
float **convert_matrix(float *,int,int,int,int);
double *dvector(int,int);
double **dmatrix(int,int,int,int);
int *ivector(int,int);
int **imatrix(int,int,int,int);
unsigned char *cvector(int,int);
unsigned char **cmatrix(int,int,int,int);
float **submatrix(float **,int,int,int,int,int,int);
/*
void free_vector(float *,int,int);
*/
void free_dvector(double *,int,int);
void free_cvector(unsigned char *,int,int);
void free_cmatrix(unsigned char **,int,int,int,int);
void free_ivector(int *,int,int);
void free_matrix(float **,int,int,int,int);
void free_dmatrix(double **, int,int,int,int);
void free_imatrix(int **,int,int,int,int);
void free_submatrix(float **,int,int,int,int);
void free_convert_matrix(float **,int,int,int,int);
void nrerror( const char * );
