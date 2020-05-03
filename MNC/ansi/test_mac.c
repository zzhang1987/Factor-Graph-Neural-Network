#include "../ansi/r.h"
#include "../ansi/mynr.h"

/* test program for macopt solution of equation A x = b. */
/* 

*/
   
#include "test.h"

void main(int argc, char *argv[])
{
  gq_args param;
  double *x , tol = 0.00001 ;
  int iter , itmax = 10 , n , rich = 0 , end_on_step = 1 , type ;

  type = end_on_step * 10 + rich ; 
  printf("Solving A x = b\nDimension of A?\n");
  inputi(&(param.n));
  n=param.n; 
  param.A=dmatrix(1,n,1,n);
  param.b=dvector(1,n);
  x=dvector(1,n);
  typeindmatrix(param.A,1,n,1,n);
  printf("b vector?\n");
  typeindvector(param.b,1,n);
  printf("Initial condition x?\n");
  typeindvector(x,1,n);

  macopt ( x , param.n , type , tol , &iter , itmax ,
	  vgrad_quadratic , (void *)(&param) 
	  ) ;

  printf("Solution:\n");
  quadratic(x,&param);
}
