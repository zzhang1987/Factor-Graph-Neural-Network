/* for cg.c */
typedef struct {
  int n;
  double *p,*xi,*xt;
  double (*nfunc)(double *,void *);
  void *func_arg;
} f1dim_arg;

typedef struct {
  int n;
  double *p,*xi,*xt;
  void (*ndfunc)(double * , double * , void *);
  void *func_arg;
} df1dim_arg;

typedef struct {
  int n;
  double **m;
  double **in;
  double **lu;
  int *indx;
  double det;
} dmatrix_family;


typedef struct {
  int n;
  double *b;
  void (*Afunc)(  double * ,  double * ,double * , void * ) ; 
  void   *Afunc_arg;
} cggq_args;       

void linmin(double *, double *, int n, double *, double (*f)(double *,void *), 
	    void *, double);	
/* void dlinmin();	*/
double brent(double, double, double, double (*f)(double, void *), void *, 
	     double, double *);	
double dbrent(double, double, double, double (*f)(double), 
	      double (*g)(double),
	      double, double *);	
void frprmn(double *, int, double, double, int *, int,
	    double *, double (*f)(double *, void *), void *,
	    void (*g)(double *,double *, void *), void *);	
void checkgrad(double *, int, double, 
	       double (*f)(double *, void *), void *,
	       void (*g)(double *,double *, void *), void *);	
void mnbrak(double *,double *,double *,double *,double *,double *,
	    double (*f)(double,void *),void *);
double f1dim(double, void *);	
/* double df1dim(); */	
void dfpmin(double *, int , double , int *, double *, 
	    double (*f)(double *,void *), 
	    void *, void (*g)(double *,double *,void *), 
	    void *, double **, double);

void macopt
  (double *,
   int    ,
   int ,
   double ,
   int   * ,
   int    ,
   void   (*dfunc)(double *,double *, void *), 
   void   * );
double maclinmin 
(
 double * , double * , int  ,    
 void   (*dfunc)(double *,double *, void *), 
 void   * , double * ,  double * ,  double * );
double macprod 
( 
 double * , double * , double * , double * , double  , int  ,
 void   (*dfunc)(double *,double *, void *), 
 void   * ) ;


void light_speed_cg
  (double *,   int,   int,  double,   int    *,   int ,
   double *,   double (*dfunc)(double *,double *, void *),
   void   *,   void   (*Afunc)(double *,double *,double *, void *),void   *);
void cg_solve  
( double *,  
 int  , 
 double * , 
 void (*f)(double *,double *,double *, void *) , 
 void * , 
 double  ) ;
double cg_grad_quadratic
( double * , double * , void * ) ;

void 	ludcmp( double ** , int , int * , double * );
void 	lubksb( double ** , int , int * , double * );

void	find_eigs( double **, int , double *);
void tqli ( double * , double * , int , double ** ) ;
void tred2 ( double ** , int , double * , double * ) ;

/* for matrix.c */
/*
double  evaluate_determinant();
double	quadratic_form();
double	lumatrixproduct();
double	matrixproduct();
double  trace();
void 	find_eigs();
int 	clip_eigs();
void 	det_and_tr_from_eigs();
void 	tred2();
void 	tqli();
double	invert_dmatrix();
double	luinvert_dmatrix();
void	symmetrise_dmatrix();
double	invert_dmatrixfamily();
*/
/* for eig.c */
/*
void	report_eigs();
void	report_eigs2();
void	scale_eigs_by_pow_sqrt_lambda();
void	pos_eigs();
 */

