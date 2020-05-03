/* bnd.c                                      (c) DJCM 95 08 31

   - belief net decoder  - 

   This code is (c) David J.C. MacKay 1995. It is free software 
   as defined by the free software foundation. 

   ------
   bnd.c:
   ------

   This is a set of routines to do Gallager / MS like  iterative 
 solution of 

 A x = z 

 defined by 

 A, represented in an alist_matrix. 
 a prior on x, given in a vector bias
 an initial condition

 No use is made of the true answer.
 
*/
/*
   How a bit in the A matrix is labelled.

   A one bit is labelled in two ways. 

   1) Its (m,n) pair, i.e. its row (m) and column (n).

   2) Its l and u values, which counts how many 1 bits it is along the row (l)
       and how many 1 bits it is down the column (u).

   A bit can be identified by its (m,l) value, its (n,u) value, or by its (m,n)
   value.

   e.g. the last bit in the bottom row of the (6*12) matrix below 
        has m=6 and n=10, and u=2, and l=4. 

   A = [ 0 1 0 0 0 1 0 1 0 0 0 1 ]
       [ 1 0 0 0 1 0 1 0 0 0 1 0 ]
       [ 0 0 1 1 0 0 0 1 0 0 0 1 ]
       [ 0 1 0 0 1 0 0 0 0 1 1 0 ]
       [ 1 0 0 1 0 1 0 0 1 0 0 0 ]
       [ 0 0 1 0 0 0 1 0 1 1 0 0 ]

   Notice that all rows in this matrix have equal weight lmax=4 and all 
   columns have equal weight umax=2 .

   A = [   1       1   1       1 ]
       [ 1       1   1       1   ]
       [     1 1       1       1 ]
       [   1     1         1 1   ]
       [ 1     1   1     1       ]
       [     1       1   1 1     ]

*/


#include "../ansi/r.h"
#include "../ansi/mynr.h"
#include "../ansi/cmatrix.h"
#include "./bnd.h"

void bnd_allocate ( bnd_param *p , bnd_control *c ) {
  int M = p->M ; 
  int N = p->N ;
  int m ; 

/* get z, xo to point to a vector allocated elsewhere; x is not used except in the printstate diagnostic*/
  p->t = cvector ( 1 , M ) ; 

/*  p->bias = dvector ( 1 , N ) ; */ /* cut out in mnc7 onwards; this is
                                         allocated outside -- NB this makes
					 this code incompatible with mnc6.c */
  p->q1 = dvector ( 1 , N ) ; 

  p->dpf = dmatrix ( 1 , M , 0 , p->a->biggest_num_m + 1 ) ;
  p->dpr = dmatrix ( 1 , M , 0 , p->a->biggest_num_m + 1 ) ;

  p->pc0 = dmatrix ( 1 , N , 1 , p->a->biggest_num_n ) ;
  p->pc1 = dmatrix ( 1 , N , 1 , p->a->biggest_num_n ) ;

  p->dqc = dmatrix ( 1 , N , 0 , p->a->biggest_num_n + 1 ) ;

  p->qc0 = dvector ( 0 , p->a->biggest_num_n + 1 ) ;
  p->qc1 = dvector ( 0 , p->a->biggest_num_n + 1 ) ;
  p->qt0 = dvector ( 0 , p->a->biggest_num_n + 1 ) ;
  p->qt1 = dvector ( 0 , p->a->biggest_num_n + 1 ) ;
  p->qb0 = dvector ( 0 , p->a->biggest_num_n + 1 ) ;
  p->qb1 = dvector ( 0 , p->a->biggest_num_n + 1 ) ;

  for ( m = 1 ; m <= M ; m ++ ) {
    p->dpf[m][0] = 1.0 ; 
    p->dpr[m][p->a->num_mlist[m]+1 ] = 1.0 ; 
  }
  /* vertical passes must be initialized too */
}

void bnd_defaults ( bnd_param *p , bnd_control *bndc ) {
#include "bnd_var_def.c"
}

void bnd_free ( bnd_param *p , bnd_control *c ) {
  int M = p->M ; 
  int N = p->N ;

  free_cvector ( p->t , 1, M ) ; 

/*  free_dvector ( p->bias , 1 , N ) ;  */ /* changed in mnc7 */
  free_dvector (   p->q1 , 1 , N ) ; 

  free_dmatrix ( p->dpf , 1 , M , 0 , p->a->biggest_num_m + 1 ) ;
  free_dmatrix ( p->dpr , 1 , M , 0 , p->a->biggest_num_m + 1 ) ;

  free_dmatrix ( p->pc0 , 1 , N , 1 , p->a->biggest_num_n ) ;
  free_dmatrix ( p->pc1 , 1 , N , 1 , p->a->biggest_num_n ) ;

  free_dmatrix ( p->dqc , 1 , N , 0 , p->a->biggest_num_n + 1 ) ;

  free_dvector ( p->qc0 , 0 , p->a->biggest_num_n + 1 ) ;
  free_dvector ( p->qc1 , 0 , p->a->biggest_num_n + 1 ) ;
  free_dvector ( p->qt0 , 0 , p->a->biggest_num_n + 1 ) ;
  free_dvector ( p->qt1 , 0 , p->a->biggest_num_n + 1 ) ;
  free_dvector ( p->qb0 , 0 , p->a->biggest_num_n + 1 ) ;
  free_dvector ( p->qb1 , 0 , p->a->biggest_num_n + 1 ) ;

}

void   bnd_load_dqc ( bnd_param *p , bnd_control *c ) {
  int u , n , N=p->N ; 

  for ( n = 1 ; n <= N ; n++ )  {
    for ( u = 1 ; u <= p->a->biggest_num_n ; u++ )  {
      p->dqc[n][u] = 1.0 - 2.0 * p->bias[n] ;
    }
  }

}

int bndecode ( bnd_param *p , bnd_control *c ) {

  bnd_load_dqc ( p , c ) ; /* sets all dqc entries to the prior to start off */

  for ( c->loop = 1 , p->not_decoded = 1 ; 
       ( c->loop <= c->loops ) && p->not_decoded ; 
       c->loop ++ ) {

    horizontal_pass ( p ) ;  /* turns dqc into dpr and dpf 
			     then reads out dpc -> pc0 and pc1 */
    vertical_pass ( p , c ) ;    /* turns pc0 and pc1 into qt0, qb0, etc, 
			     and works out qc0, qc1 -> dqc.
			     Also works out the latest marginal belief */
    bnd_score_state ( p ) ;  /* takes belief vector, turns it into an
			     x, and runs it through Ax =? z .
			     If equal then decoding is halted.
			     */
    if ( c->writelog ) 
      bnd_fprint_state ( p , c ) ; 

  }
  printf( "BND: viols %3d  recon %3d  its %2d\n", p->count_viol , p->count_high  , c->loop ) ; 

  return ( p->count_viol ) ; /* zero if success */
}

void bnd_fprint_state ( bnd_param *p , bnd_control *c ) {
  FILE *fp ;
  int n ;
  double *q1 = p->q1 ; 

  fp = fopen ( c->logfile , "a" ) ;
  if ( fp ) {
    for ( n = 1 ; n <= c->writelog ; n++ ) 
      fprintf ( fp ,  "%-2.0f " , 99.0 * q1[n] ) ; 
    fprintf ( fp ,  "\n" ) ; 
    for ( n = 1 ; n <= c->writelog ; n++ ) 
      fprintf ( fp ,  "%-2d " , p->x[n] ) ; 
    fprintf ( fp ,  "\n" ) ; 
    fclose ( fp ) ; 
  }
  else 
    fprintf ( stderr , "can't open %s \n" , c->logfile ) ;
}

void bnd_score_state ( bnd_param *p ) {
  double *q1 = p->q1 ; 
  int N = p->N ;
  int M = p->M ;
  int m , n ;

  for ( p->count_high = 0 , n = 1 ; n <= N ; n++ )  {
    p->xo[n] = ( q1[n] >= 0.5 ) ? 1 : 0 ;
    p->count_high += p->xo[n] ;
  }

  alist_times_cvector_mod2 ( (p->a) , p->xo , p->t ) ; 

  for ( p->count_viol = 0 , m = 1 ; m <= M ; m ++ ) {
    if ( p->t[m] != p->z[m] ) p->count_viol ++ ; 
  }

  if ( p->count_viol == 0 ) p->not_decoded = 0 ; 

}


void horizontal_pass ( bnd_param *p ) 
{
  int l , u , umax ; 
  int m , n ;
  double *dpf ;
  double *dqc ;
  double dpc ; 
  double **pdpf = p->dpf ;
  double *dpr ; 
  double **pdpr = p->dpr ;
  int *lupto  = p->a->l_up_to ;         /* keeps track of what l we are up to */
  int *nlist ; 
  alist_matrix *a = p->a ;
  unsigned char *z = p->z ; 
  double **pc0 = p->pc0 ;
  double **pc1 = p->pc1 ;
                                        /* * * * * * FORWARD PASS  * * * * * */
  for ( m = 1 ; m <= p->M ; m++ )   
    lupto[m] = 0 ;                      /* check all columns are at start    */

  for ( n = 1 ; n <= p->N ; n++ ) {     /* Run through bits x[n] left->right */
    umax = a->num_nlist[n] ;            /* number of bits in this column     */
    nlist = a->nlist[n] ;               /* list of bits in this column       */
    dqc = p->dqc[n] ;                   /* contributions to be made          */

    for ( u = 1 ; u <= umax ; u ++ ) {  /* Run down column              */
      m = nlist[u] ;                    /* what row are we on ?         */
      dpf = pdpf[m] ;                   /* -- for efficiency --         */

      l = lupto[m] ;                    /* what value of l = left-right */
      lupto[m] ++ ;                     /*   are we at on this row ?    */

      dpf[l+1] = dqc[u] * dpf[l] ;      /* the forward step             */
    }
  }
/*
  for ( m = 1 ; m <= p->M ; m++ )   {
    if ( lupto[m] != a->num_mlist[m] ) {
      fprintf ( stderr , "alist visitation problem %d %d %d\n" ,
	       m , lupto[m] , a->num_mlist[m] ) ;
      pause_for_return () ; 
    }
  }
*/
                                        /* * * * * * BACKWARD PASS * * * * * */

  for ( n = p->N ; n >= 1 ; n-- ) {     /* Run back through bits n           */
    umax = a->num_nlist[n] ;            /* number of bits in this column     */
    nlist = a->nlist[n] ;               /* list of bits in this column       */
    dqc = p->dqc[n] ;                   /* contributions to be made          */
    

    for ( u = 1 ; u <= umax ; u ++ ) {  /* Run down column              */
      m = nlist[u] ;                    /* what row are we on ?         */
      dpr = pdpr[m] ;                   /* -- for efficiency --         */
      dpf = pdpf[m] ;                   

      l = lupto[m] ;                    /* what value of l = left-right */
      lupto[m] -- ;                     /*   are we at on this row ?    */

      dpr[l] = dqc[u] * dpr[l+1] ;      /* the backward step             */
      dpc = dpf[l-1] * dpr[l+1] * 0.5 ; /* the desired difference        */
      if ( !z[m] ) {                    /* I think this is right!        */
	pc0[n][u] = 0.5 + dpc ;         /* messages to bit n */
	pc1[n][u] = 0.5 - dpc ;
      } else {
	pc1[n][u] = 0.5 + dpc ; 
	pc0[n][u] = 0.5 - dpc ;         /* if space is short then only 
					   a vector is really needed here,
					   if the vertical computation
					   done immediately */
      }
    }
  }
}


void vertical_pass ( bnd_param *p , bnd_control *c  ) 
{
  int u , umax ; 
  int n ;
  double  *dqc , *qt0 = p->qt0 , *qt1 = p->qt1 ;
  double         *qb0 = p->qb0 , *qb1 = p->qb1 ; 
  double **pc0 = p->pc0 ;  double **pc1 = p->pc1 ;
  double *qc0 = p->qc0 ;   double *qc1 = p->qc1 ;
  double sum , dif ; 
  alist_matrix *a = (p->a) ;
  int clipwarn=0 ; /* whether we have spat a warning this time */
  int fwarn=0 ; /* whether we have spat a warning this time */

  for ( n = 1 ; n <= p->N ; n++ ) {     /* Run through bits x[n] left->right */
    umax = a->num_nlist[n] ;              /* number of bits in this column   */
    dqc = p->dqc[n] ;                     /* contributions to be computed    */

    qt0[0] = 1.0 - p->bias[n] ;
    qt1[0] = p->bias[n] ;
    qb0[umax+1] = 1.0 ;
    qb1[umax+1] = 1.0 ;
                                        /* * * * * * DOWNWARD PASS * * * * * */
    for ( u = 1 ; u <= umax ; u ++ ) {  /* Run down column                   */
      qt0[u] = qt0[u-1] * pc0[n][u] ;   /*     downward step                 */
      qt1[u] = qt1[u-1] * pc1[n][u] ;   /*     downward step                 */
    }
    sum = ( qt0[umax] + qt1[umax] ) ;
    if ( sum > c->tinydiv ) {
      p->q1[n] = qt1[umax] / sum ;      /* pseudoposterior prob, bit n */
    } else { 
/*      p->q1[n] = 0.49 ;  */ /* leave it as it was */
      fprintf ( stderr , "*" ) ; fflush ( stderr ) ; 
    }

                                        /* * * * * * UPWARD PASS * * * * * * */
    for ( u = umax ; u >= 1 ; u -- ) {  /* Run up column                     */
      qb0[u] = qb0[u+1] * pc0[n][u] ;   /*    upward step                    */
      qb1[u] = qb1[u+1] * pc1[n][u] ;   /*    upward step                    */

      qc0[u] = qt0[u-1] * qb0[u+1] ;    /*    everyone else's messages       */
      qc1[u] = qt1[u-1] * qb1[u+1] ;    /*    everyone else's messages       */

      sum = qc0[u] + qc1[u] ;
      dif = qc0[u] - qc1[u] ;
      if ( sum > c->tinydiv ) {
	dqc[u] = dif / sum ;              /*   difference for the next iteration */
	if ( c->dofudge ) { /* research option - does stretching the qs help
			       the decoder? */
	  if ( fudge( &dqc[u] , c->fudge ) < 0 ) {
	    if ( !fwarn ) { fprintf ( stderr , "l" ) ;  fflush ( stderr ) ; fwarn=1; } 
	  }
	}
	if ( c->doclip ) { /* pragmatic clipping procedure to stop overflow */
	  if ( dqc[u] > c->clip ) { 
	    dqc[u] = c->clip ;
	    if ( !clipwarn ) { fprintf ( stderr , "c" ) ; fflush ( stderr ) ; clipwarn=1; }
	  }
	  else if ( dqc[u] < - c->clip ) { 
	    dqc[u] = - c->clip ; 
	    if ( !clipwarn ) { 
	      fprintf ( stderr , "c" ) ; fflush ( stderr ) ; clipwarn=1; 
	    }
	  }
	}
	
      } else { 
	dqc[u] = 0.0 ; fprintf ( stderr , "*" ) ; fflush ( stderr ) ; 
      }
    }
  }
}

 int fudge ( double *d , double f ) {
  double l=0.0 , q0 , q1 ;
  int status = 0 ; 
  /* turn d into a logit, scale it by f, then turn it back */
  q0 = 1.0 + *d ; 
  q1 = 1.0 - *d ; 
  if ( q0 > 0.0 ) { l = - log ( q0 ) ; }
  else { status -- ; }
  if ( q1 > 0.0 ) { l += log ( q1 ) ; }
  else { status -- ; }
  if ( !status ) 
    *d = - tanh ( 0.5 * f * l ) ; 
  /* else *d gets left alone */
  return status ; 
}

 void fudges ( double *d , double f ) {
  double l ;
  /* turn d into a logit, scale it by f, then turn it back */
  l = d_2_logit ( *d ) ;
  *d = logit_2_d ( f * l ) ; 
}

 double d_2_logit ( double d ) {
  double l = 0.0 , q0 , q1 ; 
  q0 = 1.0 + d ; 
  q1 = 1.0 - d ; 
  if ( q0 > 0.0 ) { l = - log ( q0 ) ; }
  else { fprintf ( stderr , "L" ) ; fflush ( stderr ) ; }
  if ( q1 > 0.0 ) { l += log ( q1 ) ; }
  else { fprintf ( stderr , "l" ) ; fflush ( stderr ) ; }
  return l ;
}

 double logit_2_d ( double l ) {
  double d ; 
  d = - tanh ( l * 0.5 ) ;
  return d ;
}

/*
<!-- hhmts start -->
Last modified: Tue Aug  6 17:19:47 1996
<!-- hhmts end -->
*/
