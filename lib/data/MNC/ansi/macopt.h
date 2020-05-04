/* structure for macopt */
typedef struct {
  double tol ;    /* convergence declared when the gradient vector is smaller
		     in magnitude than this, or when the mean absolute 
		     step is less than this (see above) */
  double grad_tol_tiny ; /* if gradient is less than this, we definitely 
			    stop, even if we are not using a gradient 
			    tolerance */
  double step_tol_tiny ; /* if step is less than this, we stop, even if 
			    we are not using a step tolerance */
  int end_if_small_step ; /* defines the role of tol -- alternative is
			     end_on_small_grad */
  int its ;               /* number of its */
  int itmax ;             /* max */
  int rich ; /* whether to do the extra gradient evaluation at the beginning 
	      of each new line min */
  int verbose ; 
  double stepmax ;        /* largest step permitted (not used in macopt) */

  int linmin_maxits ;     /* in maclinmin */
  double linmin_g1 ;      /* factors for growing and shrinking the interval */
  double linmin_g2 ;
  double linmin_g3 ;
  double lastx     ;      /* keeps track of typical step length */
  double lastx_default ;  /* if maclinmin is reset, lastx is set to this */

/* These should not be touched by the user. They are handy pointers for macopt
   to use 
*/
  double gtyp ; /* stores the rms gradient for linmin */
  double *pt , *gx , *gy , *gunused ;
  double *xi , *g , *h ; 
  int n ;                 /* dimension of parameter space */
  int restart ;           /* whether to restart macopt - fresh cg directions */
} macopt_args ; 


/* lastx :--- 1.0 might make general sense, (cf N.R.)
				  but the best setting of all is to have 
				  a prior idea of the eigenvalues. If 
				  the objective function is equal to sum of N
				  terms then set this to 1/N, for example 
				  Err on the small side to be conservative. */
void macoptII
  (double *,            /* starting vector                                */
   int    ,             /* number of dimensions                           */
   void   (*dfunc)(double *,double *, void *), /* evaluates the gradient   */
   void   *, 
   macopt_args *
   )  ; 

void maccheckgrad(double *, int, double, 
	       double (*f)(double *, void *), void *,
	       void (*g)(double *,double *, void *), void * , int );	

void macopt_defaults ( macopt_args * ) ;

void macopt_free ( macopt_args * ) ;

/* the following functions could be declared static within macopt.c */

double maclinminII 
(
 double * ,
 void   (*dfunc)(double *,double *, void *), /* evaluates the gradient */
 void   * ,
 macopt_args * ) ;

double macprodII 
( 
 double * , double * , double  ,
 void   (*dfunc)(double *,double *, void *), 
 void   * , 
 macopt_args *
) ;

void macopt_restart ( macopt_args * , int ) ;
