#include "../ansi/r.h"
#include "../ansi/mynr.h"

typedef struct {
  int n;
  double **A, *b; /* Used to define f(x) = 1/2 xAx + bx */
} gq_args;        /* and g = Ax + b */

double quadratic(double *, void * );
double grad_quadratic(double *,double *, void * /* gq_args */);
void Atimesh(double *,double *,double *, void * /* gq_args */);

void checkgrad2
  (double *,   int    ,   double ,
   double (*f)(double *, void *) ,   void   * ,
   double (*g)(double *, double * , void *) ,   void   *);

void main(int argc, char *argv[])
{
  gq_args param;
  double *x , tol=0.000001 , fret;
  int iter , itmax=10 , n , tol_type=0 , CG = 0 ;

  printf("Dimension?\n");
  inputi(&(param.n));
  n=param.n;
  param.A=dmatrix(1,n,1,n);
  param.b=dvector(1,n);
  x=dvector(1,n);
  typeindmatrix(param.A,1,n,1,n);
  printf("Bias vector?\n");
  typeindvector(param.b,1,n);
  printf("Initial condition?\n");
  typeindvector(x,1,n);

  if (CG)  checkgrad2 ( x , param.n , 0.001 , quadratic , &param , grad_quadratic , &param ) ;
  light_speed_cg ( x , param.n , tol_type , tol , &iter , itmax ,
		 &fret , grad_quadratic , &param ,
		 Atimesh , &param
		 ) ;
  printf("Solution:\n");
  quadratic(x, (void *) &param);
}

/* checkgrad2.c */
/* Examines objective function and d_objective function to see if 
   they agree for a step of size epsilon 
 not quite same as in cg.c */


void checkgrad2
  (double *p,
   int    n,
   double epsilon,
   double (*func)(double *, void *),
   void   *func_arg,
   double (*dfunc)(double *,double *, void *),
   void   *dfunc_arg
)
{
  int j;
  double f1;
  double *g,*h;
  
  h=dvector(1,n);
  g=dvector(1,n);
  f1=(*func)(p,func_arg);
  (*dfunc)(p,g,dfunc_arg);

  printf("Testing gradient evaluation\n");
  printf("      analytic     1st_diffs\n");
  for ( j = 1 ; j <= n ; j ++ ) {
    p[j] += epsilon ;
    h[j] = (*func)(p,func_arg) - f1 ;
    p[j] -= epsilon ;

    printf("%2d %9.5g %9.5g\n" , j , g[j] , h[j]/epsilon );
  }
  free_dvector(h,1,n);
  free_dvector(g,1,n);
}
