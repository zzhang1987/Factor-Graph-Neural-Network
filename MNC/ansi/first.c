/* first.c
   The first ansi program I have written.
   It will do conjugate gradients on a function 
   which has arguments
   and which has a gradient routine
   that has arguments
   too */

#include "../ansi/r.h"
#include "../ansi/mynr.h"

typedef struct {
  int n;
  double *m;
} param;

double objective(double *, void *);
void d_objective(double *, double *, void *);

void main(int argc, char *argv[])
{
  double ftol=0.001, flintol=0.0001;
  param control;
  int i,iter=1000;
  double *q,v;

  /* initialise control parameters */

  control.n=10;
  control.m=dvector(1,control.n);
  q=dvector(1,control.n);
  for(i=1;i<=control.n;i++){
    control.m[i]=(double)i;
    q[i]=100.0;
  }
  
  frprmn(q,control.n,ftol,&iter,&v,objective,&control,d_objective,&control,flintol);
}

double objective(double *q, void *f_args)
{
  param control;
  double dummy=0.0;
  int i;
  
  printf(":");
  control = *( (param *) f_args);

  for(i=1;i<=control.n;i++){
    dummy+=q[i]*q[i]*control.m[i]/2.0;
  }
  return(dummy+10.0);
}

void d_objective(double *q, double *g, void *f_args)
{
  param control;
  int i;
  
  control = *( (param *) f_args);

  for(i=1;i<=control.n;i++){
    g[i]=q[i]*control.m[i];
  }
  printf("%f %f\n",q[1],q[2]);
}




