#include "../ansi/r.h"
#include "../ansi/mynr.h"

/* test program for cg solution of equation A x = b. */
/*  This program sets up A and b, and it knows about a magical routine 
    Atimesh which applies A to a vector.

    It calls cg_solve ( x , b , Atimesh , args_for_Atimesh ) . 

    cg_solve constructs a package to give to 

    light_speed_cg ( )

    which needs functions to evaluate  (1)  grad f
    and                                (2)  grad grad f . h

    (1) is a function called cg_grad_quadratic, which finds 1/2 xAx - bx
    and its gradient. The function Atimesh is communicated to this guy
    in a package of arguments. It also needs to know about b. 
    (2) is just Atimesh.

    This will lend itself to use by a sparse Atimesh machine. 

*/
   
#include "test.h"

void main(int argc, char *argv[])
{
  gq_args param;
  double *x , tol = 0.000001 , fret;
  int iter , itmax = 10 , n , tol_type=0 ;
  int cg = 1 ;

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

  switch(cg){
  case(1):
  default:
    cg_solve ( x , param.n, param.b , Atimesh , &param , 0.00001 ) ;
    break;
  case(0):
    light_speed_cg ( x , param.n , tol_type , tol , &iter , itmax ,
		    &fret , grad_quadratic , (void *)(&param) ,
		    Atimesh , (void *)(&param)
	 	    ) ;
    break;
  }
  printf("Solution:\n");
  quadratic(x,&param);
}
