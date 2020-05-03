/* #include <stdio.h>
#include <math.h> */
#include "../ansi/r.h" 
/* #include "../ansi/nrutil.h" */
#include "../ansi/mynr.h"

/* CONJUGATE GRADIENT ALGORITHMS */

/* Rewritten starting 9 June 1992 to make them ANSI */
/* And also to make them non-broken --- no external functions allowed 
   Temporarily reduced in size leaving full file as cg.c.all */

/* frprmn.c */
/* The NR conjugate gradients algorithm */

#define EPS 1.0e-10
#define FREEALL free_dvector(xi,1,n);free_dvector(h,1,n);free_dvector(g,1,n);
#define FRPVERBOSE 0

void frprmn
  (double *p,
   int    n,
   double ftol,
   double flinmintol,
   int    *iter,
   int    itmax,
   double *fret,
   double (*func)(double *, void *),
   void   *func_arg,
   void   (*dfunc)(double *,double *, void *),
   void   *dfunc_arg
)
/* double p[],ftol,*fret,(*func)(),flinmintol; */
{
  int j,its;
  double gg,gam,fp,dgg;
  double *g,*h,*xi;
  
  g=dvector(1,n);
  h=dvector(1,n);
  xi=dvector(1,n);
  fp=(*func)(p,func_arg);
  (*dfunc)(p,xi,dfunc_arg);
  for (j=1;j<=n;j++) {
    g[j] = -xi[j];
    xi[j]=h[j]=g[j];
  }
  for (its=1;its<=itmax;its++) {
    *iter=its;
    *fret = fp ; /* 21 05 94: pass the current value through to mnbrak */
    linmin(p,xi,n,fret,func,func_arg,flinmintol); 
    /* linmin * creates new vectors for p and xi (not necessary - cut)
              * calls mnbrak and then brent 
              * moves p to the new minimum when it is done,
              * and changes xi to the step used. (Note no use is made
	        of xi later).
	      * the initial step in the bracketing is equal to xi! 
	      * I now include a factor on that step. 
       mnbrak * evaluates func at p, unnecessarily */
    if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS)) {
      FREEALL
      return;
    }
    fp=(*func)(p,func_arg);
    (*dfunc)(p,xi,dfunc_arg);
    dgg=gg=0.0;
    for (j=1;j<=n;j++) {
      gg += g[j]*g[j];
      /*		  dgg += xi[j]*xi[j];	*/
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
  printf("Too many iterations in FRPRMN, but continuing.\n"); 
  FREEALL	
  return;
}

#undef EPS
#undef FREEALL



/* linmin */

/* 
typedef struct {
  int n;
  double *p , *xi , *xt ;
  double (*nfunc)() ;
  void *func_arg ;
} f1dim_arg ; 
*/

void linmin
  (double *p,
   double *xi,
   int    n,
   double *fret,
   double (*func)(double *,void *),
   void   *func_arg,
   double linmin_tol
)
{
  int j;
  double xx,xmin,fx,fb,fa,bx,ax;
  double f1dim(double,void *);
  static double lastxx = 0.01 ; /* 1.0 might make more sense, but hey! */
  f1dim_arg f_arg;
  
  f_arg.n=n;
  f_arg.nfunc=func;
  f_arg.func_arg=func_arg;
  f_arg.p=p;              /* NB in NR, these are copied into global. 
			     Here I just copy the pointers */
  f_arg.xi=xi;
  f_arg.xt=dvector(1,n) ; /* new addition 21 05 94 : scratch pad for f1dim
			     to use */

  fa = *fret ; /* 21 05 94: pass the current value through to mnbrak */
  ax=0.0;
  xx=lastxx ;
  bx=2.0 * xx ;
  mnbrak ( &ax , &xx , &bx , &fa , &fx , &fb , f1dim , &f_arg); 
  if ( FRPVERBOSE > 1 ) printf ( "B %6.3g %6.3g %6.3g " , ax , xx , bx  ) ; 
/* typically these numbers are 1.0 0.0 -1.6
    but they can be 0.0 1.0 2.7 also. */
/* my guess is it would be good for the latter to happen sometimes ! */
/* Keep a note of the actual stretch factor needed.
   If I set it to xx (stable state), and then if ax is bigger, set it 
   to ax/2.0 ? */

  lastxx = fabs ( xx ) ;
  if ( fabs ( ax ) > lastxx ) lastxx =  fabs ( ax ) * 0.5 ;

  *fret=brent ( ax , xx , bx , f1dim , &f_arg , linmin_tol , &xmin);
  for (j=1;j<=n;j++) {
    xi[j] *= xmin; 
    p[j] += xi[j];
  }
  free_dvector(f_arg.xt , 1 , n);
}

/* brent.c */

#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

double brent
  (double ax,
   double bx,
   double cx,
   double (*f)(double,void *),
   void *f_arg,
   double tol,
   double *xmin
)
{
  int iter;
  double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  double e=0.0;

  a=((ax < cx) ? ax : cx);
  b=((ax > cx) ? ax : cx);
  x=w=v=bx;
  fw=fv=fx=(*f)(x,f_arg);
  for (iter=1;iter<=ITMAX;iter++) {
    xm=0.5*(a+b);
    tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
    if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
      *xmin=x;
      return fx;
    }
    if (fabs(e) > tol1) {
      r=(x-w)*(fx-fv);
      q=(x-v)*(fx-fw);
      p=(x-v)*q-(x-w)*r;
      q=2.0*(q-r);
      if (q > 0.0) p = -p;
      q=fabs(q);
      etemp=e;
      e=d;    /* is this an uninitialized use of d ? */
      if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	d=CGOLD*(e=(x >= xm ? a-x : b-x));
      else {
	d=p/q;
	u=x+d;
	if (u-a < tol2 || b-u < tol2)
	  d=SIGN(tol1,xm-x);
      }
    } else {
      d=CGOLD*(e=(x >= xm ? a-x : b-x));
    }
    u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		fu=(*f)(u,f_arg);
    if (fu <= fx) {
      if (u >= x) a=x; else b=x;
      SHFT(v,w,x,u)
	SHFT(fv,fw,fx,fu)
	} else {
	  if (u < x) a=u; else b=u;
	  if (fu <= fw || w == x) {
	    v=w;
	    w=u;
	    fv=fw;
	    fw=fu;
	  } else if (fu <= fv || v == x || v == w) {
	    v=u;
	    fv=fu;
	  }
	}
  }
  fprintf(stderr,"Too many iterations in BRENT");
  *xmin=x;
  return fx;
}

#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SIGN

/* dbrent.c */

#define ITMAX 100
#define ZEPS 1.0e-10
#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define MOV3(a,b,c, d,e,f) (a)=(d);(b)=(e);(c)=(f);

double dbrent
  (double ax,
   double bx,
   double cx,
   double (*f)(double),
   double (*df)(double),
   double tol,
   double *xmin
)
{
  int iter,ok1,ok2;
  double a,b,d,d1,d2,du,dv,dw,dx,e=0.0;
  double fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;
  
  a=(ax < cx ? ax : cx);
  b=(ax > cx ? ax : cx);
  x=w=v=bx;
  fw=fv=fx=(*f)(x);
  dw=dv=dx=(*df)(x);
  for (iter=1;iter<=ITMAX;iter++) {
    xm=0.5*(a+b);
    tol1=tol*fabs(x)+ZEPS;
    tol2=2.0*tol1;
    if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
      *xmin=x;
      return fx;
    }
    if (fabs(e) > tol1) {
      d1=2.0*(b-a);
      d2=d1;
      if (dw != dx)  d1=(w-x)*dx/(dx-dw);
      if (dv != dx)  d2=(v-x)*dx/(dx-dv);
      u1=x+d1;
      u2=x+d2;
      ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
      ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
      olde=e;
      e=d;
      if (ok1 || ok2) {
	if (ok1 && ok2)
	  d=(fabs(d1) < fabs(d2) ? d1 : d2);
	else if (ok1)
	  d=d1;
	else
	  d=d2;
	if (fabs(d) <= fabs(0.5*olde)) {
	  u=x+d;
	  if (u-a < tol2 || b-u < tol2)
	    d=SIGN(tol1,xm-x);
	} else {
	  d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
	}
      } else {
	d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
      }
    } else {
      d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
    }
    if (fabs(d) >= tol1) {
      u=x+d;
      fu=(*f)(u);
    } else {
      u=x+SIGN(tol1,d);
      fu=(*f)(u);
      if (fu > fx) {
	*xmin=x;
	return fx;
      }
    }
    du=(*df)(u);
    if (fu <= fx) {
      if (u >= x) a=x; else b=x;
      MOV3(v,fv,dv, w,fw,dw)
	MOV3(w,fw,dw, x,fx,dx)
	  MOV3(x,fx,dx, u,fu,du)
    } else {
      if (u < x) a=u; else b=u;
      if (fu <= fw || w == x) {
	MOV3(v,fv,dv, w,fw,dw)
	  MOV3(w,fw,dw, u,fu,du)
      } else if (fu < fv || v == x || v == w) {
	MOV3(v,fv,dv, u,fu,du)
      }
    }
  }
  nrerror("Too many iterations in routine DBRENT");
  return (0.0) ;
}

#undef ITMAX
#undef ZEPS
#undef SIGN
#undef MOV3



/* mnbrak.c */

#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define MAX(a,b) ((a) > (b) ? (a) : (b) )

void mnbrak
(double *ax,
 double *bx,
 double *cx,
 double *fa,
 double *fb,
 double *fc,
 double (*func)(double,void *),
 void   *f_arg
 )
{
  double ulim,u,r,q,fu,dum ;
  
/*  *fa = ( *func ) ( *ax,f_arg ) ; */
/* I have cut this because in my applications, func has been evaluated just */
/* I pass this value instead */
  *fb = ( *func ) ( *bx,f_arg ) ;
  if  ( *fb > *fa ) {
    SHFT ( dum,*ax,*bx,dum )
      SHFT ( dum,*fb,*fa,dum )
  }
  *cx = ( *bx ) + GOLD * ( *bx - *ax ) ;
  *fc = ( *func ) ( *cx,f_arg ) ;
  while  ( *fb > *fc ) {
    r = ( *bx - *ax ) * ( *fb - *fc ) ;
    q = ( *bx - *cx ) * ( *fb - *fa ) ;
    u = ( *bx ) -  (  ( *bx - *cx ) *q -  ( *bx - *ax ) *r )/
       ( 2.0*SIGN ( MAX ( fabs ( q - r ),TINY ),q - r ) ) ;
    ulim = ( *bx ) + GLIMIT * ( *cx - *bx ) ;
    if  (  ( *bx - u ) * ( u - *cx ) > 0.0 ) {
      fu = ( *func ) ( u , f_arg ) ;
      if  ( fu < *fc ) {
	*ax = ( *bx ) ;
	*bx = u ;
	*fa = ( *fb ) ;
	*fb = fu ;
	return ;
      } else if  ( fu > *fb ) {
	*cx = u ;
	*fc = fu ;
	return ;
      }
      u = ( *cx ) + GOLD* ( *cx - *bx ) ;
      fu = ( *func ) ( u,f_arg ) ;
    } else if  (  ( *cx - u ) * ( u - ulim ) > 0.0 ) {
      fu = ( *func ) ( u,f_arg ) ;
      if  ( fu < *fc ) {
	SHFT ( *bx,*cx,u,*cx + GOLD* ( *cx - *bx ) )
	  SHFT ( *fb,*fc,fu, ( *func ) ( u,f_arg ) )
      }
    } else if  (  ( u - ulim ) * ( ulim - *cx ) >= 0.0 ) {
      u =ulim ;
      fu = ( *func ) ( u,f_arg ) ;
    } else {
      u = ( *cx ) + GOLD* ( *cx - *bx ) ;
      fu = ( *func ) ( u,f_arg ) ;
    }
    SHFT ( *ax,*bx,*cx,u )
      SHFT ( *fa,*fb,*fc,fu )
  }
}

#undef GOLD
#undef GLIMIT
#undef TINY
#undef MAX
#undef SIGN
#undef SHFT



double f1dim
  (double x,
   void   *f_argp
) 
{
  int j;
  double f,*xt;
  f1dim_arg f_arg;
  
  f_arg = *( (f1dim_arg *)f_argp);
  xt= f_arg.xt ; 
  for (j=1;j<=f_arg.n;j++) xt[j]=f_arg.p[j]+x*f_arg.xi[j];
  f=(*f_arg.nfunc)(xt,f_arg.func_arg);
  return f;
}

/*
double df1dim
  (double x,
   void   *f_argp
) 
{
  int j ;
  double f = 0.0 , *xt , *df ;
  df1dim_arg f_arg ;
  
  f_arg = *( (df1dim_arg *)f_argp);
  xt=dvector(1,f_arg.n);
  df=dvector(1,f_arg.n);
  for ( j = 1 ; j <= f_arg.n ; j++ ) xt[j] = f_arg.p[j] + x * f_arg.xi[j] ;
  (*f_arg.ndfunc)( xt , df , f_arg.func_arg ) ;
  for ( j = 1 ; j <= f_arg.n ; j++ ) f += df[j] * f_arg.xi[j] ;
  free_dvector(xt,1,f_arg.n);
  free_dvector(df,1,f_arg.n);
  return f;
}
*/

/* dfpmin.c */
/* Modified 30 Nov 90 so that it uses allocated space for the Hessian */  
/* 	hessin=matrix(1,n,1,n);	must be declared and initialised previously */
/* also modified to force it to do at least *iter loops, (so that a good hessian 
estimate is made) */ 
/* also modified to pass a tolerance to linmin when it is called */
/* Modified 16 3 92 from dfpmin.c so that the functions func and dfunc both 
   include an additional structure argument that contains all the other variables */
/* So that globals are no longer needed */

#define ITMAX 400
#define EPS 1.0e-10

/* void dfpmin(p,n,ftol,iter,fret,func,func_arg,dfunc,dfunc_arg,hessin,frac_lin_tol)
double p[],ftol,*fret,(*func)(),**hessin,frac_lin_tol;
void (*dfunc)();
void *func_arg,*dfunc_arg;
int n,*iter; */
void dfpmin
  (double *p,
   int    n,
   double ftol,
   int    *iter,
   double *fret,
   double (*func)(double *,void *),
   void   *func_arg,
   void   (*dfunc)(double *, double *,void *),
   void   *dfunc_arg,
   double **hessin,
   double frac_lin_tol
)
{
  int j,i,its;
  double fp,fae,fad,fac;
  double *xi,*g,*dg,*hdg;
  int	min_its;
  
  min_its= *iter; 
  
  xi=dvector(1,n);
  g=dvector(1,n);
  dg=dvector(1,n);
  hdg=dvector(1,n);
  fp=(*func)(p,func_arg);
  (*dfunc)(p,g,dfunc_arg);
  
  for (i=1;i<=n;i++) {
    xi[i]=0.0;
    for (j=1;j<=n;j++) xi[i] -= hessin[i][j]*g[j];
  }
  for (its=1;its<=ITMAX;its++) {
    *iter=its;
    linmin(p,xi,n,fret,func,func_arg,frac_lin_tol);
    if ((its>=min_its)&&(2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS))) {
      free_dvector(hdg,1,n);
      free_dvector(dg,1,n);
      free_dvector(g,1,n);
      free_dvector(xi,1,n);
      return;
    }
    fp=(*fret);
    for (i=1;i<=n;i++) dg[i]=g[i];
    *fret=(*func)(p,func_arg);
    (*dfunc)(p,g,dfunc_arg);
    for (i=1;i<=n;i++) dg[i]=g[i]-dg[i];
    for (i=1;i<=n;i++) {
      hdg[i]=0.0;
      for (j=1;j<=n;j++) hdg[i] += hessin[i][j]*dg[j];
    }
    fac=fae=0.0;
    for (i=1;i<=n;i++) {
      fac += dg[i]*xi[i];
      fae += dg[i]*hdg[i];
    }
    fac=1.0/fac;
    fad=1.0/fae;
    for (i=1;i<=n;i++) dg[i]=fac*xi[i]-fad*hdg[i];
    for (i=1;i<=n;i++)
      for (j=1;j<=n;j++)
	hessin[i][j] += fac*xi[i]*xi[j]
	  -fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];
    for (i=1;i<=n;i++) {
      xi[i]=0.0;
      for (j=1;j<=n;j++) xi[i] -= hessin[i][j]*g[j];
    }
  }
  /*	nrerror("Too many iterations in DFPMIN2"); */
  printf("Too many iterations in DFPMIN2, but continuing.\n"); 
  free_dvector(hdg,1,n);
  free_dvector(dg,1,n);
  free_dvector(g,1,n);
  free_dvector(xi,1,n);
  return;
}

#undef ITMAX
#undef EPS


/* checkgrad.c */
/* Examines objective function and d_objective function to see if 
   they agree for a step of size epsilon */


void checkgrad
  (double *p,
   int    n,
   double epsilon,
   double (*func)(double *, void *),
   void   *func_arg,
   void   (*dfunc)(double *,double *, void *),
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
    fflush(stdout) ; 
  }
  free_dvector(h,1,n);
  free_dvector(g,1,n);
}
/*
static int lnsrch(double *fret, double *xmin, 
		     double (*fn1d)( double * , void * ),
		     void   (*dfn1d)(double *, double * , void * ),
		     double n_over_grad)
{
  double f0, df0, f1, f2, f_target, y1, y2, mu, nu, a;
  double l1, l2, l_min = 0.1 , l_max = 0.5 ;
  int      i, status;
  int      done = 0;
  double offset = 0;
  double scale  = 0;

  */ /* Pick a scale size (small if it's very hilly) */ /*
  if (n_over_grad > 1.0)
    scale = n_over_grad * global->lnsrch_max_step;
  else
    scale = global->lnsrch_max_step;

  */ /* Never go further than the Newton step though */ /*
  if (scale > 1.0)
    scale = 1.0;

  */ /* Get the gradient for starting point */ /*
  status = (*dfn1d)(offset, &f0, &df0);
  df0 = df0 * scale;

  */ /* Check that +ve lambda sends us downhill */ /*
  if (df0 > 0)
    {
      scale = -scale;
      df0   = -df0;
    }

  */ /* Try the full Newton jump */ /*
  l1 = 1;
  f1 =  (*fn1d)(l1 * scale + offset, );

  f_target = f0 + global->lnsrch_alpha * df0 * l1;
  if (f1 < f_target)    {
      *fret = f1;
      *xmin = l1 * scale + offset;
      return 0;
    }
  
  */ /* Now try a quadratic approx */ /*
  l_min = global->lnsrch_min_lambda;
  l2 = - df0 / (2 * (f1 - f0 - df0));
  if (l2 < l_min)    {
      l2   = l_min;
      done = 1;
    }

  f2  = (*fn1d)(l2 * scale + offset, );

  f_target = f0 + global->lnsrch_alpha * df0 * l2;
  if ( (f2 < f_target) || (done) ) 
    {
      *fret = f2;
      *xmin = l2 * scale + offset;
      return 0;
    }
  
  */ /* OK we're in cubicville AZ */ /*
  y1 = f1 - f0 - l1 * df0;
  y2 = f2 - f0 - l2 * df0;
  
  for(i = 0; i < global->lnsrch_max_iters; i++)
    {
      mu = ( (l2 * l2 * l2 * y1 - l1 * l1 * l1 * y2) / 
	    (l2 * l2 * y1 - l1 * l1 * y2));
      a  = (l1 * l1 * l2 * l2 * (l1 - l2)) / (l2 * l2 * y1 - l1 * l1 * y2);
      
      nu = mu * mu - 3 * df0 * a;
      if (nu < 0)
	{
	  fprintf( stderr , "lnsrch: Error (4) Probable roundoff error\n");

	  *fret = f2;
	  *xmin = l2 * scale + offset;
	  return -4;
	}
      
      l_min = global->lnsrch_min_lambda * l1;
      l_max = global->lnsrch_max_lambda * l1;
      l1 = l2;
      f1 = f2;
      y1 = y2;

      l2 = (mu + sqrt(nu)) / 3.0;
      if (l2 < l_min)	{
	  l2   = l_min;
	  done = 1;
	}
      else if (l2 > l_max)	{
	  l2 = l_max;
	}

      f2 = (*fn1d)(l2 * scale + offset, );
      
      f_target = f0 + global->lnsrch_alpha * df0 * l2;
      if ( (f2 < f_target) || (done) )	{
	  *fret = f2;
	  *xmin = l2 * scale + offset;
	  return 0;
	}
  
      y2 = f2 - f0 - l2 * df0;
    }

  fprintf( stderr , "lnsrch: Error (5) Too many iterations - %12.6g at %12.6g\n",
	    f2, l2 * scale + offset);
  *xmin = l2 * scale + offset;
  *fret = f2;
  
  return -5;
}
*/

/* 
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

   See end of file for further discussion.

   NB: The value of "tol" is totally arbitrary and must be set by 
   you to a value that works well for your problem. 
   It depends completely on the typical value of the gradient / step size. 

   Tol specifies a magnitude of gradient at which a halt is called. 
   or a step size.
*/

#define EPS 1.0e-10
#define FREEALL free_dvector(xi,1,n);free_dvector(h,1,n);free_dvector(g,1,n);  free_dvector ( pt , 1 , n ) ;   free_dvector ( gx , 1 , n ) ;   free_dvector ( gy , 1 , n ) ;  


void macopt
  (double *p,            /* starting vector                                */
   int    n,             /* number of dimensions                           */
   int    type,          /* communicates what sort of opt to do            */
                         /* if type is even (eg 0) then no fresh gradient 
			    is evaluated at the start of each line search, 
			    reducing the number of gradients 
		            by one per loop; `rich' is defined to be
			    type % 2                                       

			    if type/10, then tol defines 
			    a step length criterion. 
			    else tol defines a gradient magnitude criterion */
   double tol,    /* convergence declared when the gradient vector is smaller
		     in magnitude than this, or when the mean absolute 
		     step is less than this (see above) */
   int    *its,   /* or when the number of iterations                      */
   int    itmax,  /* reaches itmax                                         */
   void   (*dfunc)(double *,double *, void *), /* evaluates the gradient   */
   void   *dfunc_arg
   )                     /* Note, (*func)(double *,void *) is not used     */
{
  int j ;
  double gg,gam,dgg;
  double *g,*h,*xi;
  int rich = type % 2 ;
  int end_if_small_step = type / 10 ; 
  int end_if_small_grad = ( end_if_small_step ) ? 0 : 1 ; 
  double step ;

  double *pt , *gx , *gy ;         /* vectors for linmin to use */
  /* A total of 7 double * 1..n are used by this optimizer. 
     p           is provided when the optimizer is called 
     pt          is used by the line minimizer as the temporary vector. 
                    this could be cut out with minor rewriting, using p alone
     g, h and xi are used by the cg method as in NR - could one of these
                     be cut out?
     the line minimizer uses an extra gx and gy to evaluate two gradients. 
     */
#define MACOPTVERBOSE 2
  
  g = dvector ( 1 , n ) ;
  h = dvector ( 1 , n ) ;
  xi = dvector ( 1 , n ) ;

  pt = dvector(1,n) ; /* scratch vector */
  gx = dvector(1,n) ; /* scratch gradient */
  gy = dvector(1,n) ; /* scratch gradient */

  (*dfunc)(p,xi,dfunc_arg);
  for (j=1;j<=n;j++) {
    g[j] = -xi[j];
    xi[j]=h[j]=g[j];
  }
  for ( *its = 1 ; *its <= itmax ; *its = *its + 1 ) {

    for ( gg = 0.0 , j = 1 ; j <= n ; j ++ ) 
      gg += g[j]*g[j];          /* find the magnitude of the old gradient */

    if ( MACOPTVERBOSE > 0 ) 
      printf ( "mac_it %d of %d : gg = %6.3g tol = %6.3g: ", *its , itmax , gg , tol ) ;

    /* do a pseudo line minimization */

    if ( end_if_small_grad && gg <= tol ) {
      FREEALL
      if ( MACOPTVERBOSE > 0 ) printf ("\n");
      return;
    }

    step = maclinmin ( p , xi , n , dfunc , dfunc_arg , pt , gx , gy ) ; 

    if ( MACOPTVERBOSE > 1 ) printf (" (step %9.5g)",step);
    if ( MACOPTVERBOSE > 0 ) printf ("\n");
    if ( end_if_small_step && step <= tol ) {
      FREEALL
      return;
    }

    /* Method: 
       evaluate gradient at a sequence of points and calculate the inner 
       product with the line search direction. Continue until a 
       bracketing is achieved ( i.e a change in sign ). */

    /* if there is a step length criterion, look at step size... */
    
    /* if we are feeling rich, evaluate the gradient at the new `minimum' */
    /* alternatively, linmin constructs this gradient by linear combination 
       of the last two evaluations and leaves it in xi */
    if ( rich ) { 
      (*dfunc)(p,xi,dfunc_arg); 
    }
    dgg=0.0;
    for (j=1;j<=n;j++) {
      dgg += (xi[j]+g[j])*xi[j];
    }
    gam=dgg/gg;
    for (j=1;j<=n;j++) {
      g[j] = -xi[j];                /* g stores (-) the most recent gradient */
      xi[j]=h[j]=g[j]+gam*h[j];     /* h stores xi, the current line direction */
    }
  }
  fprintf(stderr,"Reached iteration limit in MACOPT; continuing.\n"); 
  FREEALL	
  return;
} /* NB this leaves the best value of p in the p vector, but
     the function has not been evaluated there if rich=0     */

#undef EPS
#undef FREEALL

#define G1 2.0
#define G2 1.25
#define G3 1.5
#define MAXITS 30

double maclinmin 
(
 double *p , double *xi , int n ,    
 void   (*dfunc)(double *,double *, void *), /* evaluates the gradient */
 void   *arg ,
 double *pt , double *gx , double *gy )
{
  static double lastx = 0.01 ; /* 1.0 might make general sense, (cf N.R.)
				  but the best setting of all is to have 
				  a prior idea of the eigenvalues. If 
				  the objective function is equal to sum of N
				  terms then set this to 1/N, for example 
				  Err on the small side to be conservative. */

  double x , y ;
  double s , t , m ;
  int    its = 1 , i ;
  double step , tmpd ; 

  x = lastx ;
  s = macprod ( p , xi , pt , gx , x , n , dfunc , arg ) ; 
  
  /* at x=0, the gradient (uphill) satisfies s < 0 */

  if ( s < 0 )  {  /* we need to go further */
    do {
      y = x * G1 ;
      t = macprod ( p , xi , pt , gy , y , n , dfunc , arg ) ; 
      if ( MACOPTVERBOSE > 1 ) 
	printf ("s = %6.3g: t = %6.3g; x = %6.3g y = %6.3g\n",s, t , x , y );
      if ( t >= 0 ) break ;
      x = y ; s = t ; for ( i = 1 ; i <= n ; i ++ ) gx[i] = gy[i] ;
    }
    while ( its++ < MAXITS ) ;
  } else if ( s > 0 ) { /* need to step back inside interval */
    do {
      y = x / G3 ;
      t = macprod ( p , xi , pt , gy , y , n , dfunc , arg ) ; 
      if ( MACOPTVERBOSE > 1 ) 
	printf ("s = %6.3g: t = %6.3g; x = %6.3g y = %6.3g\n",s, t , x , y );
      if ( t <= 0 ) break ;
      x = y ; s = t ; for ( i = 1 ; i <= n ; i ++ ) gx[i] = gy[i] ;
    } while ( its++ < MAXITS ) ;
  } else { /* hole in one s = 0.0 */
    t = 1.0 ; y = x;
  }
  if ( its >= MAXITS ) {
/* this can happen where the function goes \_ and doesn't buck up again */
    fprintf (stderr, "Warning! maclinmin overran\n" );
    fprintf (stderr, "- inner product at 0 = %9.4g\n" ,
	     macprod ( p , xi , pt , gy , 0.0 , n , dfunc , arg ) ) ; 

  }

 /*  linear interpolate between the last two */
  if ( s < 0.0 ) s = - s ;
  if ( t < 0.0 ) t = - t ;
  m = ( s + t ) ;
  s /= m ; t /= m ;
  
  m =  s * y + t * x ; 
  /* evaluate the step length, not that it necessarily means anything */
  for ( step = 0.0 , i = 1 ; i <= n ; i ++ ) {
    tmpd = m * xi[i] ;
    p[i] += tmpd ;
    step += fabs ( tmpd ) ; 
    xi[i] = s * gy[i] + t * gx[i] ;
/* send back the estimated gradient in xi (NB not like linmin) */
  }
  lastx = m * G2 ;
  
  return ( step / (double) ( n ) ) ; 
}

#undef G2
#undef G1
#undef G3
#undef MAXITS

double macprod 
( 
 double *p , double *xi , double *pt , double *gx , double x , int n ,
 void   (*dfunc)(double *,double *, void *), 
 void   *arg
) {
  /* finds pt = p + x xi and gets gx there, 
				       returning gx . xi */

  int i;
  double s = 0.0 ;

  for ( i = 1 ; i <= n ; i ++ ) 
    pt[i] = p[i] + x * xi[i] ;
  
  dfunc( pt , gx , arg ) ;

  for ( i = 1 ; i <= n ; i ++ ) 
    s += gx[i] * xi[i] ;

  return s ;
}
  
  
/* 
   Explanation of the macopt code 
   ------------------------------

Over the last year or two I have realised that the conjugate gradient 
code in the Numerical Recipes book is very poorly written. 
* In their algorithm is that the initial step 
size is always `1', even though this may be a very inappropriate 
step size. If this initial step size is too big then the routine mnbrak 
makes the next step equal to `-1.6': NB, this is in the wrong direction!
And it is an even bigger step than the first one, so it is likely to 
give even worse numerical problems. 
Imagine, for the sake of discussion, that the true minimum is at a step 
size of about 0.001 on each line search. Then the Numerical Recipes 
algorithm will waste a lot of time creating silly guesses -- about 
log_2 1000 of them, in fact. So about 10 unnecessary function evaluations
are done on every line search, whereas if the initial step size were 
0.001, then only a few function evaluations would be needed to locate 
the minimum. 

* Another criticism of the NR code is that they do not use gradient 
information in the line search. They mention that it is possible to 
use gradient information, but their example code does not use it. 
In fact, with gradient information, it is easier to find the line 
minimum, because you can bracket the minimum with only two gradient 
evaluations. 

* Finally, it is not necessary to locate the line minimum as accurately
as linmin does. 

I have written my own conjugate gradient algorithm that attempts to 
improve on the NR code. My algorithm is called macopt. 
It has the following properties:
	1) an adaptive step size is used for the initial step of the line 
		search
	2) gradients are used in the line search
	3) once the minimum has been bracketed, the line search 
		terminates immediately, and a new line search commences
		using interpolation to estimate the location of the minimum 
		and the value of the gradient there. 
	4) typically, when the routine has adapted its step size, two 
		gradient evaluations per line minimization are performed.
I find that this algorithm sometimes is ten times faster than 
the NR code. The Xerion group at University of Toronto have also 
written their own conjugate gradient optimizers, and they say that 
they have adaptive conjugate gradient optimizers that only need one gradient 
evaluation per line search!

You can find macopt on wol in ~mackay/ansi/cg.c
If you also get test_mac.c and test_function.c 
and nrutil.c and r.c and r.h and mynr.h
you will be able to compile a demonstration program 
called test_mac which uses macopt to minimize a quadratic function. 

*/
/*
<!-- hhmts start -->
Last modified: Sat Apr 15 21:23:46 1995
<!-- hhmts end -->
*/
