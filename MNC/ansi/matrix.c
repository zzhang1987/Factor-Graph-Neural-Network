/* clib/matrix.c */
/* cut out lots of stuff on 950815. Because ludcmp mislaid */

#include "../ansi/r.h"
#include "../ansi/mynr.h"

void	find_eigs(
	double **m,
	int n,
		  double *l)
		 
{
	double *e;

	e=dvector(1,n);
	tred2(m,n,l,e);
        tqli(l,e,n,m);

	free_dvector(e,1,n);	
}

/* Makes diagonal decomposition of the matrix */
#define	float	double
#define	free_vector free_dvector 
#define	vector dvector 
/* #include "/home/mackay/NR/recipes/ludcmp.c" */
/* #include "/home/mackay/NR/recipes/lubksb.c" */
#include "/home/mackay/clib/tred2.c"
#include "/home/mackay/clib/tqlifast.c"
#undef	float
#undef	free_vector
#undef	vector 

/* top eig */
double top_eig( 
	double **m, double *v,
	int n , double ftol, int itmin , int itmax
	)
{
  /* multiply m by an initial vector (user supplied) until small change in ratio */
  double oldeig , eig=1.0 , *w , *truev , *temp ;
  int its ; 
  double vmag = 1.0 , wmag = 1.0 ; 
  truev = v ;

  w = dvector ( 1 , n ) ; 
  its = 0 ; 
  if ( itmin <=1 ) itmin = 2 ; 
  
  do {
    its ++ ; 
    mult ( m , v , w ) ; 
    wmag = magvector ( w ) ; 
    eig = wmag ; 
    rescale ( w , (1.0/wmag) ) ; /* reduce to unit magnitude */
    if ( its >= itmin ) {
      /* check to see if there's been much change */
    }
  } while ( its < itmax ) ; 
  /* copy final eig into v */
  free_dvector(e,1,n);	
  if ( ) ;
    return eig ; 
}
