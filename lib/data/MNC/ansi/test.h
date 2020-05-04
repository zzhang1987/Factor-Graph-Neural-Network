double quadratic ( double * , void * ) ;
double grad_quadratic ( double *,double *, void *);
void vgrad_quadratic ( double *,double *, void *);
void Atimesh ( double *,double *,double *, void *); 
typedef struct {
  int n;
  double **A, *b; /* Used to define f(x) = 1/2 xAx - bx */
} gq_args;        /* and g = Ax - b */
