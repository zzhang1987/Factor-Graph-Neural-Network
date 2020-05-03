#include <stdio.h>
#include <math.h>
#include "nrutil.h"
#include "mynr.h"

/* Two routines: 
   light_speed_cg       Uses grad and grad grad routines to minimise
   and
   cg_solve             Solves the equation A x = b using cg, and a routine 
                        that applies A to v. 
*/

/* This was once frprmn.c */
/* First Rewritten starting 9 June 1992 to make ANSI */
/* -> void cg */
/* And also to make them non-broken --- no external functions allowed */
/* This conjugate gradients algorithm WITHOUT line search is written 
   by David MacKay 8/10/92, inspired by John Skilling */
/* The idea is to use the standard conjugate gradients sequence of directions 
   for search but take advantage of the fact that we are able to evaluate 
   Av, IE Hessian*arbitrary vector, in time similar to the time needed 
   for a single function / gradient evaluation */
/* An additional modification is that now the function itself is assumed 
   to be returned by the gradient routine */
/* And the tolerance is expressed in absolute and in fractional terms */

#define EPS 1.0e-20
#define FREEALL free_dvector(xi,1,n);free_dvector(h,1,n);free_dvector(g,1,n);free_dvector(v,1,n);

void light_speed_cg
  (double *p, /* current parameter values */
   int    n,
   int    tol_type,
   double tol,
   int    *iter,/* will contain actual number of iterations */
   int    itmax,
   double *fret,
   double (*dfunc)(double *,double *, void *),
   void   *dfunc_arg,
   void   (*Afunc)( double * , double * , double * , void * ) , 
                 /* this function evaluates Ah at p and puts it in v */
   void   *Afunc_arg
)
{
  int j;
  double gg,gam,fp,dgg,num,denom,stepsize;
  double *g,*h,*xi,*v;
/* ,*dvector();
  void linmin(),nrerror(),free_dvector(); */

  g=dvector(1,n);
  h=dvector(1,n);
  xi=dvector(1,n);
  v=dvector(1,n);
  fp=  (*dfunc)(p,xi,dfunc_arg);
  for (j=1;j<=n;j++) {
    g[j] = -xi[j];                 /* Downhill direction is in g */
    xi[j]=h[j]=g[j];               /* and in xi and h */
  }
  for ( *iter = 1 ; *iter <= itmax ; *iter = *iter + 1 ) {
    (*Afunc)(p,h,v,Afunc_arg);     /* Sets v = Ah */
    /* stepsize= gradient.h / h.A.h */
    num = 0 ; denom = 0 ;
    for ( j = 1 ; j <= n ; j++ ) {
      num+=g[j]*h[j]; denom+=h[j]*v[j];
    }
    if ( num == 0 ) { 
      FREEALL
      return;
    }
    else if(denom<=0) {
      fprintf(stderr,"Warning: lscg: hAh negative. Terminating lscg\n");
      FREEALL
      return;
    }
    stepsize=num/denom;
    /* replacing:     linmin(p,xi,n,fret,func,func_arg,flinmintol); */
    for (j=1;j<=n;j++) {
      p[j] += stepsize*xi[j];
    }
    *fret=    (*dfunc)(p,xi,dfunc_arg);
    if ( ( tol_type /* fractional */ && 
	  2.0*fabs(*fret-fp) <=  tol*(fabs(*fret)+fabs(fp)+EPS)
	  ) || 
	( tol_type==0 /* absolute */ && 
	 2.0*fabs(*fret-fp) <= tol )) {
      FREEALL
      return;
    }
/* Add in a check that function has gone down */
    else if (*fret > fp){
      fprintf( stderr , "Warning: lscg: function value increased\n" ) ;
    }
    fp=*fret;
    dgg=gg=0.0;
    for (j=1;j<=n;j++) {
      gg += g[j]*g[j];
      dgg += (xi[j]+g[j])*xi[j];
    }
    if (gg == 0.0) {
      FREEALL
      return;
    }
    gam=dgg/gg;
    for (j=1;j<=n;j++) {
      g[j] = -xi[j];
      xi[j]=h[j]=g[j]+gam*h[j];
    }
  }
  fprintf ( stderr , "Too many iterations in LSCG, but continuing.\n" ) ;
  FREEALL
  return;
}

#undef EPS
#undef FREEALL

void cg_solve 
( double *x, 
 int n , 
 double *b , 
 void (*Afunc)( double * , double * , double * , void *) , 
 void *Afunc_arg , 
 double tol ) 
/* Solves A x = b by  minimising 1/2 x A x - b x  whose gradient is  A x - b */
/* A assumed to be symmetric, and Afunc assumed to apply A to vector */
{
  int iter ;
  cggq_args  param ;  
  double fret;

  param.n = n ; /* arguments for the grad f routine to use */
  param.b = b ;
  param.Afunc = Afunc ;
  param.Afunc_arg = Afunc_arg ;

  light_speed_cg ( x , n , 1 /* tol_type */ , tol , &iter , n+1 ,
		    &fret , cg_grad_quadratic , (void *)(&param) ,
		    Afunc , Afunc_arg /* here I had &Afunc_arg */
		    ) ;
}

double cg_grad_quadratic
( double *x, double *g, void *arg ) 
{
  int i;
  double f = 0.0;
  cggq_args *param = (cggq_args *) arg;

  param->Afunc( x , x , g , param->Afunc_arg ) ; /* the first argument is junk */
  for(i=1;i<=param->n;i++){
    f += x[i] * ( g[i] * 0.5 - param->b[i] ); 
    g[i] -= param->b[i];
  } 
  return(f);
} 


