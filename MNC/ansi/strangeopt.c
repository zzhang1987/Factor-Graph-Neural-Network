#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
/* #include "../ansi/r.h" */
#include "../ansi/nrutil.h"
#include "../ansi/mynr.h"
#include "../ansi/macopt.h"
#include "../ansi/strangeopt.h"

static double safedivide ( double , double ) ;

/* 
   
   strangeopt

   a special purpose optimizer for fe minimization at infinite temperature

   Does steepest descent, with step size set such that one parameter 
   changes sign on each line search.

   For fe min, strangeopt has the advantage that one expects the total 
   number of gradient evaluations to be virtually identical to the 
   number of bits that are flipped.

*/

void strangeopt
  (double *p,            /* starting vector                                */
   int    n,             /* number of dimensions                           */
   double  (*dfunc)(double *,double *, void *), 
                         /* evaluates the gradient of the optimized function,
			   and returns the value */
   void   *dfunc_arg,    /* arguments that get passed to dfunc             */
   macopt_args *a        /* structure in which optimizer arguments stored  */
   )                 
{
/*
   Method: 

   evaluate f & gradient -> g ;

   do {
     evaluate, for each i, t = p[i]/g[i]. 
     If they are all negative then terminate.

     Of the positive ones, select the smallest in a robust way, t, and make 
     it a tiny bit bigger.

     do {
        p -= t * g ;
	evaluate fnew [ & gradient -> gnew ] ; 
     }
     while ( fnew > f && 	t *= 0.5    ) ;   (hopefully just one loop)

     set f = fnew, g = gnew
   }
*/
#define MYHUGE   1e80 
#define TINYZERO 1e-20
#define MYTINY   1e-3
#define HALF     0.5
#define ONEANDBIT 1.01

  int    i , poscount , linits , linitmax=20 ;
  double *g , *gnew , *gtmp ;
  double f , fnew , t , c ;

  g    = dvector ( 1 , n ) ;
  gnew = dvector ( 1 , n ) ;

  f = (*dfunc)( p , g , dfunc_arg );

  for ( a->its = 1 ; a->its <= a->itmax ; a->its ++ ) {
    
    for ( t = MYHUGE , poscount = 0 , i = 1 ; i <= n ; i ++ ) {
      c = safedivide ( p[i] , g[i] ) ;
      if ( c > TINYZERO ) {
	poscount ++ ; 
	if ( c < t ) t = c ; 
      }
    }

    if ( poscount == 0 ) break ; 
    t *= ONEANDBIT ; 

    for ( linits = 1 ; linits <= linitmax ; t *= ( linits == 1 ) ? 
	 - HALF : HALF , linits++ ) {
      for (  i = 1 ; i <= n ; i ++ ) {
	p[i] -= t * g[i] ; 
      }
      fnew = (*dfunc)( p , gnew , dfunc_arg );
      if ( fnew <= f ) break ;
    }
    if ( ! ( fnew <= f ) )
      fprintf(stderr,"strangeopt function value went up.\n");
    f = fnew ; 
    gtmp = g ; g = gnew ; gnew = gtmp ; 

  }
  if ( poscount != 0 ) 
    fprintf(stderr,"Reached iteration limit in strangeopt; continuing.\n"); 
  free_dvector ( g    , 1 , n ) ;  
  free_dvector ( gnew , 1 , n ) ;  
  return;
}

static double safedivide ( double x , double y ) {
  if ( fabs ( y ) > TINYZERO ) return ( x / y ) ; 
  else return 0.0 ; 
}

void ddcheckgrad 
/* Examines objective function and d_objective function to see if 
   they agree for a step of size epsilon */
  (double *p,
   int    n,
   double epsilon,
   double (*dfunc)(double *,double *, void *),
   void   *dfunc_arg ,
   int    stopat          /* stop at this component. If 0, do the lot. */
)
{
  int j;
  double f1;
  double *g2 , *g , *h ;
  
  h=dvector(1,n);
  g=dvector(1,n);
  g2=dvector(1,n);

  f1=  (*dfunc)(p,g,dfunc_arg);
  if ( stopat <= 0 || stopat > n ) stopat = n ; 

  printf("Testing gradient evaluation\n");
  printf("      analytic     1st_diffs\n");
  for ( j = 1 ; j <= stopat ; j ++ ) {
    p[j] += epsilon ;
    h[j] = (*dfunc)(p,g2,dfunc_arg) - f1 ;
    p[j] -= epsilon ;

    printf("%2d %9.5g %9.5g\n" , j , g[j] , h[j]/epsilon );
    fflush(stdout) ; 
  }
  free_dvector(h,1,n);
  free_dvector(g,1,n);
  free_dvector(g2,1,n);
  printf("      --------     ---------\n");
}

void creepopt
  (double *p,            /* starting vector                                */
   int    n,             /* number of dimensions                           */
   double  (*dfunc)(double *,double *, void *), 
                         /* evaluates the gradient of the optimized function,
			   and returns the value */
   void   *dfunc_arg,    /* arguments that get passed to dfunc             */
   macopt_args *a        /* structure in which optimizer arguments stored  */
   )                 
{
/*
   creepopt

   a boring conservative optimizer

   goes downhill making a step of size such that the largest movement 
   is stepmax (or smaller)

   evaluate f & gradient -> g ;

   do {
     find biggest |g[i]|.
     set t = stepmax / |g[i]|

     special termination condition: 
     evaluate, for each i, t = p[i]/g[i]. 
     If they are all negative then terminate.

     do {
        p -= t * g ;
	evaluate fnew [ & gradient -> gnew ] ; 
     }
     while ( fnew > f && 	t *= 0.5    ) ;   (hopefully just one loop)

     set f = fnew, g = gnew
   }
*/

  int    i , poscount , linits , linitmax=20 ;
  double *g , *gnew , *gtmp ;
  double f , fnew , t , tmpd , stepmax = a->stepmax ;
  int fiddle = 0 ; 

  g    = dvector ( 1 , n ) ;
  gnew = dvector ( 1 , n ) ;

  f = (*dfunc)( p , g , dfunc_arg );

  for ( a->its = 1 ; a->its <= a->itmax ; a->its ++ ) {
    
    for ( t = TINYZERO , poscount = 0 , i = 1 ; i <= n ; i ++ ) {
      if ( ! ( ( p[i] > 0.0 ) ^ ( g[i] > 0.0 ) ) ) {
	poscount ++ ; 
      }
      if ( ( tmpd = fabs ( g[i] ) ) > t ) 
	t = tmpd ; 
    }

    t = stepmax / t ; 

    if ( poscount == 0 ) break ; 

    for ( linits = 1 ; linits <= linitmax ; t *= ( linits == 1 ) ? 
	 - HALF : HALF , linits++ ) {
      for (  i = 1 ; i <= n ; i ++ ) {
	p[i] -= t * g[i] ; 
      }
      fnew = (*dfunc)( p , gnew , dfunc_arg );
      if ( fnew <= f ) break ;
    }
    if ( linits > 1 ) {
      stepmax *= HALF ; /* shorter step next time */
      fiddle ++ ; 
      if ( stepmax < MYTINY ) fprintf(stderr,"t"); 
    } else if ( fiddle ) {
      if ( stepmax < a->stepmax ) stepmax *= 2.0 ; /* try growing it again */
      else {
	stepmax = a->stepmax ;
	fiddle = 0 ; 
      }
    }
    if ( ! ( fnew <= f ) ) {
      fprintf(stderr,"c"); /* creepopt function value went up. */
    }
    f = fnew ; 
    gtmp = g ; g = gnew ; gnew = gtmp ; 

  }
  if ( poscount != 0 ) 
    fprintf(stderr,"C\n"); 
/*    fprintf(stderr,"Reached iteration limit in creepopt; continuing.\n"); */
  free_dvector ( g    , 1 , n ) ;  
  free_dvector ( gnew , 1 , n ) ;  

  return;
}

/*
<!-- hhmts start -->
Last modified: Thu Sep 21 01:58:43 1995
<!-- hhmts end -->
*/
