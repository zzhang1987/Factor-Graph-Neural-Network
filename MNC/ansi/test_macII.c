#include "../ansi/r.h"
#include "../ansi/mynr.h"
#include "../ansi/macopt.h"

/* 
   test program for macopt solution of equation A x = b. 

   this equation is solved by minimizing the function 1/2 xAx - bx
*/
   
#include "test.h"

void main(int argc, char *argv[])
{
  gq_args param;
  macopt_args a ;
  double *x ;
  int n ;
  double epsilon=0.001 ;

  /* Load up the parameters of the function that you want to optimize */
  printf("Solving A x = b\nDimension of A (eg 2)?\n");
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

  /* Check that the gradient_function is the gradient of the function  */
  /* You don't have to do this, but it is a good idea when debugging ! */
  maccheckgrad (  x , param.n , epsilon , 
		quadratic , (void *)(&param) , 
		vgrad_quadratic , (void *)(&param) , 
		0
		) ;

  /* initialize the arguments of the optimizer */
  macopt_defaults ( &a ) ; 

  /* Do an optimization */
  macoptII ( x , param.n , 
	    vgrad_quadratic , (void *)(&param) , &a
	    ) ;

  printf("Solution:\n");
  quadratic(x,&param);
}
