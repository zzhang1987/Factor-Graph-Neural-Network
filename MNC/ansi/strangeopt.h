void strangeopt
  (double *,            /* starting vector                                */
   int    ,             /* number of dimensions                           */
   double  (*dfunc)(double *,double *, void *), 
                         /* evaluates the gradient of the optimized function,
			   and returns the value */
   void   *_arg,    /* arguments that get passed to dfunc             */
   macopt_args *        /* structure in which optimizer arguments stored  */
   )             ;    
void creepopt
  (double *,            /* starting vector                                */
   int    ,             /* number of dimensions                           */
   double  (*dfunc)(double *,double *, void *), 
                         /* evaluates the gradient of the optimized function,
			   and returns the value */
   void   *_arg,    /* arguments that get passed to dfunc             */
   macopt_args *        /* structure in which optimizer arguments stored  */
   )             ;    
void ddcheckgrad 
/* Examines objective function and d_objective function to see if 
   they agree for a step of size epsilon */
  (double *,
   int    ,
   double ,
   double (*dfunc)(double *,double *, void *),
   void   * ,
   int              /* stop at this component. If 0, do the lot. */
) ;
/*
<!-- hhmts start -->
Last modified: Sat Apr 15 14:01:39 1995
<!-- hhmts end -->
*/
