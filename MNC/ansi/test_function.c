#include <stdio.h> 
#include "test.h"

double quadratic(double *x,   void *arg )
{
  int i,j;
  double f=0.0,g;
  gq_args *param = ( gq_args * ) arg ;
 
  printf ( "quadr. at " ) ;
  for ( i = 1 ; i <= param->n ; i++ ) printf( "%g " , x[i] ) ;
  for ( i = 1 ; i <= param->n ; i++ ) {
    g=0.0;
    for(j=1;j<=param->n;j++){
      g += param->A[i][j] * x[j];
    }
    f+=x[i]*(g/2.0 - param->b[i]);
  }
  printf(" : %g\n",f);
  return(f); 
}

double grad_quadratic(double *x, double *g, void *arg )
{
  int i,j;
  double f=0.0;
  gq_args *param = ( gq_args * ) arg ;
  printf ( "grad_q at " ) ;
  for ( i = 1 ; i <= param->n ; i++ ) printf ( "%g " , x[i] ) ;
  for ( i = 1 ; i <= param->n ; i++ ) {
    g[i] = 0.0 ;
    for ( j = 1 ; j <= param->n ; j++ ) {
      g[i] += param->A[i][j] * x[j] ;
    }
    f += x[i] * ( g[i] / 2.0 - param->b[i] ) ;
    g[i] -= param->b[i] ;
  }
  printf ( " :: %g :: " , f ) ;
  for ( i = 1 ; i <= param->n ; i++ ) printf ( "%g " , g[i] ) ;
  printf ( "\n" ) ;
  return ( f ) ;
}


void vgrad_quadratic(double *x, double *g, void *arg )
{
  int i,j;
  double f=0.0;
  gq_args *param = ( gq_args * ) arg ;
  printf ( "grad_q at " ) ;
  for ( i = 1 ; i <= param->n ; i++ ) printf ( "%g " , x[i] ) ;
  for ( i = 1 ; i <= param->n ; i++ ) {
    g[i] = 0.0 ;
    for ( j = 1 ; j <= param->n ; j++ ) {
      g[i] += param->A[i][j] * x[j] ;
    }
    f += x[i] * ( g[i] / 2.0 - param->b[i] ) ;
    g[i] -= param->b[i] ;
  }
  printf ( " :: %g :: " , f ) ;
  for ( i = 1 ; i <= param->n ; i++ ) printf ( "%g " , g[i] ) ;
  printf ( "\n" ) ;
}

void Atimesh( double *x , double *h , double *v , void *arg)
{ /* This ignores the first argument */
  int i,j;
  gq_args *param = ( gq_args * ) arg ;

  for ( i = 1 ; i <= param->n ; i++ ) {
    v[i] = 0 ;
    for ( j = 1 ; j <= param->n ; j++ ) {
      v[i] += param->A[i][j] * h[j] ;
    }
  }
}
