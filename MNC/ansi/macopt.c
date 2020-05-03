/* #include <stdio.h>
#include <math.h> */
#include "../ansi/r.h" 
/* #include "../ansi/nrutil.h" */
#include "../ansi/mynr.h"
#include "../ansi/macopt.h"

/* 
   
   http://131.111.48.24/mackay/c/macopt.html       mackay@mrao.cam.ac.uk

   David MacKay's optimizer, based on conjugate gradient ideas, 
   but using bracketing of the zero of the inner product 

             (gradient).(line_search_direction)

   to do the line minimization. Only derivative calculations are required.
   The length of the first step in the line search (often set to "1.0"
   in other code) is adapted here so that, if 0.00001 is a better step size, 
   it soon cottons on to that and saves ~log(10000) bracketing operations.
   The result is that (with rich set to 0) the program can use 
   as few as 2 derivatives per line search. (If rich is set to 1, it does 
   an extra derivative calculation at the beginning of each line search 
   making a minimum of 3 per line search. Set rich=0 if you think 
   that the surface is locally quite quadratic.) If the program does average 
   2 derivatives per line search then it must be superior to most cg methods 
   including use of Rbackprop (which costs 2 derivatives straight off)

   A possible modification: where the function can be returned at same 
   time as the dfunction --- there is nothing clever to do with the 
   value, but it could be used as a sanity check and a convergence criterion. 

   See http://131.111.48.24/mackay/c/macopt.html for further discussion.

   NB: The value of "tol" is totally arbitrary and must be set by 
   you to a value that works well for your problem. 
   It depends completely on the typical value of the gradient / step size. 

   Tol specifies a magnitude of gradient at which a halt is called. 
   or a step size.

   This program MINIMIZES a function.
*/

void macoptII
  (double *p,            /* starting vector                                */
   int    n,             /* number of dimensions                           */
   void   (*dfunc)(double *,double *, void *), 
                         /* evaluates the gradient of the optimized function */
   void   *dfunc_arg,    /* arguments that get passed to dfunc             */
   macopt_args *a        /* structure in which optimizer arguments stored  */
   )                     /* Note, (*func)(double *,void *) is not used     */
{
  int j ;
  double gg , gam , dgg ;
  double *g , *h , *xi ;
  int end_if_small_grad = 1 - a->end_if_small_step ;
  double step , tmpd ;

  /* A total of 7 double * 1..n are used by this optimizer. 
     p           is provided when the optimizer is called 
     pt          is used by the line minimizer as the temporary vector. 
                    this could be cut out with minor rewriting, using p alone
     g, h and xi are used by the cg method as in NR - could one of these
                     be cut out?
     the line minimizer uses an extra gx and gy to evaluate two gradients. 
     */
  
  a->n = n ; 
  a->g = dvector ( 1 , n ) ; g = a->g ;    /* vectors as in NR code */
  a->h = dvector ( 1 , n ) ; h = a->h ;    /*                       */
  a->xi = dvector ( 1 , n ) ; xi = a->xi ; /*                       */
  a->pt = dvector ( 1 , n ) ; /* scratch vector for sole use of macprod */
  a->gx = dvector ( 1 , n ) ; /* scratch gradients             */
  a->gy = dvector ( 1 , n ) ; /* used by maclinmin and macprod */

  (*dfunc)( p , xi , dfunc_arg );
  macopt_restart ( a , 1 ) ; 
  for ( a->its = 1 ; a->its <= a->itmax ; a->its ++ ) {

    for ( gg = 0.0 , j = 1 ; j <= n ; j ++ ) 
      gg += g[j]*g[j];          /* find the magnitude of the old gradient */
    a->gtyp = sqrt ( gg / (double)(n) ) ; 

    if ( a->verbose > 0 ) 
      printf ( "mac_it %d of %d : gg = %6.3g tol = %6.3g: ", a->its , a->itmax , gg , a->tol ) ;

    if ( ( end_if_small_grad && ( gg <= a->tol ) ) 
	|| ( gg <= a->grad_tol_tiny ) ) {
      macopt_free ( a ) ;
      if ( a->verbose > 0 ) printf ("\n");
      return;
    }

    step = maclinminII ( p , dfunc , dfunc_arg , a ) ; 

    if ( a->restart == 0 ) {
      if ( a->verbose > 1 ) printf (" (step %9.5g)",step);
      if ( a->verbose > 0 ) printf ("\n");
      if ( ( a->end_if_small_step  && ( step <= a->tol ) ) 
	  || ( step <= a->step_tol_tiny ) ) {
	macopt_free ( a ) ;
	return;
      }
    }

    /* if we are feeling rich, evaluate the gradient at the new
       `minimum'. alternatively, linmin has already estimated this
       gradient by linear combination of the last two evaluations and
       left it in xi */
    if ( a->rich || a->restart ) { 
      (*dfunc)( p , xi , dfunc_arg ) ; 
    }
    if ( a->restart ) {
      fprintf(stderr,"Restarting macopt\n" ) ; 
      macopt_restart ( a , 0 ) ;
/* this is not quite right
   should distinguish whether there was an overrun indicating that the 
   value of lastx needs to be bigger / smaller; 
   in which case resetting lastx to default value may be a bad idea, 
   giving an endless loop of resets 
*/
    } else {
      dgg=0.0;
      for ( j = 1 ; j <= n ; j ++ ) {
	dgg += ( xi[j] + g[j] ) * xi[j] ;
      }
      gam = dgg / gg ;
      for ( tmpd = 0.0 , j = 1 ; j <= n ; j ++ ) {
	g[j] = -xi[j];                /* g stores (-) the most recent gradient */
	xi[j] = h[j] = g[j] + gam * h[j] ;
	/* h stores xi, the current line direction */
	/* check that the inner product of gradient and line search is < 0 */
	tmpd -= xi[j] * g[j] ; 
      }
      if ( tmpd > 0.0  || a->verbose > 2 ) {
	fprintf(stderr,"new line search has inner prod %9.4g\n", tmpd ) ; 
      }
      if ( tmpd > 0.0 ) { 
	if ( a->rich == 0 ) {
	  fprintf (stderr, "Setting rich to 1; " ) ; 
	  a->rich = 1 ; 
	}
	a->restart = 2 ; /* signifies that g[j] = -xi[j] is already done */
	fprintf(stderr,"Restarting macopt (2)\n" ) ; 
	macopt_restart ( a , 0 ) ;
      }
    }
  }
  fprintf(stderr,"Reached iteration limit in macopt; continuing.\n"); 
  macopt_free ( a ) ;	
  return;
} /* NB this leaves the best value of p in the p vector, but
     the function has not been evaluated there if rich=0     */

/* maclinmin.
   Method: 
       evaluate gradient at a sequence of points and calculate the inner 
       product with the line search direction. Continue until a 
       bracketing is achieved ( i.e a change in sign ). */
double maclinminII 
(
 double *p , 
 void   (*dfunc)(double *,double *, void *), /* evaluates the gradient */
 void   *arg ,
 macopt_args *a )
{
  int n = a->n ; 

  double x , y ;
  double s , t , m ;
  int    its = 1 , i ;
  double step , tmpd ; 
  double  *gx = a->gx , *gy = a->gy ;

  /* at x=0, the gradient (uphill) satisfies s < 0 */
  if ( a->verbose > 2 ) { /* check this is true: (no need to do this really
			   as it is already checked at the end of the main
			   loop of macopt) */
/*
#define TESTS 5
    x = a->lastx / a->gtyp ;
    fprintf (stderr, "inner product at:\n" ) ; 
    for ( i = -TESTS ; i <= TESTS ; i += 2 ) {
      step = x * 2.0 * (double) i / (double) TESTS ; 
      fprintf (stderr, "%9.5g %9.5g\n" , step ,
	       tmpd = macprodII ( p , gy , step , dfunc , arg , a ) ) ; 
    }
*/
    fprintf (stderr, "inner product at 0 = %9.4g\n" ,
	     tmpd = macprodII ( p , gy , 0.0 , dfunc , arg , a ) ) ; 
    if ( tmpd > 0.0 ) { 
      a->restart = 1 ; 
      return 0.0 ; 
    }
  }

  x = a->lastx / a->gtyp ;
  s = macprodII ( p , gx , x , dfunc , arg , a ) ; 
  
  if ( s < 0 )  {  /* we need to go further */
    do {
      y = x * a->linmin_g1 ;
      t = macprodII ( p , gy , y , dfunc , arg , a ) ; 
      if ( a->verbose > 1 ) 
	printf ("s = %6.3g: t = %6.3g; x = %6.3g y = %6.3g\n",s, t , x , y );
      if ( t >= 0.0 ) break ;
      x = y ; s = t ; a->gunused = gx ; gx = gy ; gy = a->gunused ; 
      its++ ;
/* replaces: for ( i = 1 ; i <= n ; i ++ ) gx[i] = gy[i] ; */
    }
    while ( its <= a->linmin_maxits ) ;
  } else if ( s > 0 ) { /* need to step back inside interval */
    do {
      y = x * a->linmin_g3 ;
      t = macprodII ( p , gy , y , dfunc , arg , a ) ; 
      if ( a->verbose > 1 ) 
	printf ("s = %6.3g: t = %6.3g; x = %6.3g y = %6.3g\n",s, t , x , y );
      if ( t <= 0.0 ) break ;
      x = y ; s = t ; a->gunused = gx ; gx = gy ; gy = a->gunused ; 
      its ++ ;
    } while ( its <= a->linmin_maxits ) ;
  } else { /* hole in one s = 0.0 */
    t = 1.0 ; y = x;
  }

  if ( its > a->linmin_maxits )  {
    fprintf (stderr, "Warning! maclinmin overran" );
/* this can happen where the function goes \_ and doesn't buck up
 again; it also happens if the initial `gradient' does not satisfy
 gradient.`gradient' > 0, so that there is no minimum in the supposed
 downhill direction.  I don't know if this actually happens... If it
 does then I guess a->rich should be 1.

 If the overrun is because too big a step was taken then
 the interpolation should be made between zero and the most 
 recent measurement. 

 If the overrun is because too small a step was taken then 
 the best place to go is the most distant point. 
 I will assume that this doesn't happen for the moment.

 Also need to check up what happens to t and s in the case of overrun.
 And gx and gy. 

 Maybe sort this out when writing a macopt that makes use of the gradient 
 at zero? 
*/
    fprintf (stderr, "- inner product at 0 = %9.4g\n" ,
	     tmpd = macprodII ( p , gy , 0.0 , dfunc , arg , a ) ) ; 
    if ( tmpd > 0 && a->rich == 0 ) {
      fprintf (stderr, "setting rich to 1\n" ) ;       a->rich = 1 ; 
    }
    if ( tmpd > 0 ) a->restart = 1 ; 
  }

 /*  Linear interpolate between the last two. 
     This assumes that x and y do bracket. */
  if ( s < 0.0 ) s = - s ;
  if ( t < 0.0 ) t = - t ;
  m = ( s + t ) ;
  s /= m ; t /= m ;
  
  m =  s * y + t * x ; 
  /* evaluate the step length, not that it necessarily means anything */
  for ( step = 0.0 , i = 1 ; i <= n ; i ++ ) {
    tmpd = m * a->xi[i] ;
    p[i] += tmpd ; /* this is the point where the parameter vector steps */
    step += fabs ( tmpd ) ; 
    a->xi[i] = s * gy[i] + t * gx[i] ;
/* send back the estimated gradient in xi (NB not like linmin) */
  }
  a->lastx = m * a->linmin_g2 *  a->gtyp ;
  
  return ( step / (double) ( n ) ) ; 
}

double macprodII 
( 
 double *p , double *gy , double y , 
 void   (*dfunc)(double *,double *, void *), 
 void   *arg , 
 macopt_args *a
) {
  double *pt = a->pt ; 
  double *xi = a->xi ; 
  /* finds pt = p + y xi and gets gy there, 
				       returning gy . xi */
  int n = a->n ; 

  int i;
  double s = 0.0 ;

  for ( i = 1 ; i <= n ; i ++ ) 
    pt[i] = p[i] + y * xi[i] ;
  
  dfunc( pt , gy , arg ) ;

  for ( i = 1 ; i <= n ; i ++ ) 
    s += gy[i] * xi[i] ;

  return s ;
}
  
void macopt_defaults ( macopt_args *a ) {
  a->verbose = 2 ; /* if verbose = 1 then there is one report for each
		      line minimization.
		      if verbose = 2 then there is an additional report for
		      each step of the line minimization.
		      if verbose = 3 then extra debugging routines kick in.
		   */

  a->tol = 0.001 ;           /* Do fiddle with this */
  a->end_if_small_step = 1 ; /* Change this to 0/1 if you prefer */
  a->itmax = 100 ;           /* You may wish to change this */
  a->rich = 1 ;              /* if this is 1, then the program runs a bit slower */
  a->stepmax = 0.5 ; 

  a->grad_tol_tiny = 1e-16 ; /* Probably not worth fiddling with */
  a->step_tol_tiny = 0.0 ;   /* Probably not worth fiddling with */
  a->linmin_maxits = 20 ;    /* Probably not worth fiddling with */
  a->lastx = 0.01 ;          /* only has a transient effect      */
  a->lastx_default = 0.01 ;  /* -- defines typical distance in parameter
				space at which the line minimum is expected;
				both these should be set. the default is 
				consulted if something goes badly wrong and 
				a reset is demanded. */

/* don't fiddle with the following, unless you really mean it */
  a->linmin_g1 = 2.0 ; 
  a->linmin_g2 = 1.25 ; 
  a->linmin_g3 = 0.5 ; 
  a->restart = 0 ; 
}

void macopt_free ( macopt_args *a ) 
{
  int n = a->n ; 
  free_dvector ( a->xi , 1 , n ) ;
  free_dvector ( a->h  , 1 , n ) ;
  free_dvector ( a->g  , 1 , n ) ;  
  free_dvector ( a->pt , 1 , n ) ;   
  free_dvector ( a->gx , 1 , n ) ;  
  free_dvector ( a->gy , 1 , n ) ;  
}

void macopt_restart ( macopt_args *a , int start ) 
/* if start == 1 then this is the start of a fresh macopt, not a restart */
{
  int j , n=a->n ; 
  double *g, *h, *xi ;
  g = a->g ;  h = a->h ;  xi = a->xi ; 

  if ( start == 0 ) a->lastx = a->lastx_default ; 
  /* it is assumed that (*dfunc)( p , xi , dfunc_arg ) ;  has happened */
  for ( j = 1 ; j <= n ; j ++ ) {
    if ( a->restart != 2 ) g[j] = -xi[j] ;
    xi[j] = h[j] = g[j] ;
  }
  a->restart = 0 ; 
}

void maccheckgrad 
/* Examines objective function and d_objective function to see if 
   they agree for a step of size epsilon */
  (double *p,
   int    n,
   double epsilon,
   double (*func)(double *, void *),
   void   *func_arg,
   void   (*dfunc)(double *,double *, void *),
   void   *dfunc_arg ,
   int    stopat          /* stop at this component. If 0, do the lot. */
)
{
  int j;
  double f1;
  double *g,*h;
  
  h=dvector(1,n);
  g=dvector(1,n);
  f1=(*func)(p,func_arg);
  (*dfunc)(p,g,dfunc_arg);
  if ( stopat <= 0 || stopat > n ) stopat = n ; 

  printf("Testing gradient evaluation\n");
  printf("      analytic     1st_diffs\n");
  for ( j = 1 ; j <= stopat ; j ++ ) {
    p[j] += epsilon ;
    h[j] = (*func)(p,func_arg) - f1 ;
    p[j] -= epsilon ;

    printf("%2d %9.5g %9.5g\n" , j , g[j] , h[j]/epsilon );
    fflush(stdout) ; 
  }
  free_dvector(h,1,n);
  free_dvector(g,1,n);
  printf("      --------     ---------\n");
}

/*
<!-- hhmts start -->
Last modified: Fri May 31 11:11:12 1996
<!-- hhmts end -->
*/
