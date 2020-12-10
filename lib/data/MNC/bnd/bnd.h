#ifndef _BND_H_
#define _BND_H_
#include "ansi/cmatrix.h"

/* 
   bnd.h

   Belief Network Decoder for the problem A x = z
                                                                */

/* Control parameters */
typedef struct { /* bnd_control */
  int loops ;       /* Max number of loops  <20>   */
  int loop ;        /* -                    <0>    */
  int writelog ;    /* number of states to log <0>    */
  char logfile[100] ;/* for state info       <->    */
  int doclip ;      /* whether to clip extreme probablities <1>    */
  double clip ;     /* clip probs here      <0.9999999999> */
  double tinydiv ;  /* value below which division should not be done <1e-40> */
  int dofudge ;     /* whether to fudge q values up/down or not <0>    */
  double fudge ;    /* fudge scale          <1.0>  */

  int verbose  ;

} bnd_control ;

typedef struct { /* bnd_param */
                  /* The sizes of everything.     */
  int M ;         /* number of rows in the matrix A */
  int N ;         /* number of columns in A         */
  int m , n ;     /* used to keep track of which row / column we are up to 
		     during complex optimizations */
                  /* The matrix A is represented in a list structure */
  alist_matrix *a ; /* unlike in fe.* , this is a pointer */
  double *bias ;      /* the prior on the s vector */
  unsigned char *z ;  
  unsigned char *xo ;  /* decoded x */
  unsigned char *t ;  /* current value of Ax */

  double *q1 ; /* the pseudoposterior */

  double **pc0 , **pc1 ; /* probabilities contributed 
			    from each relationship pc0[1..N][1..umax] */
  double *qc0 , *qc1 ; /* probabilities contributed by each bit to each 
			    relationship qc0[1..umax] */
  double **dqc ;         /* difference qc0 - qc1 */
  double *qt0 , *qt1 ; /* product of probabilities from top
			    relationship qt0[1..N][1..umax] */
  double *qb0 , *qb1 ; /* and from bottom */
  double **dpf , **dpr ; /* forward and reverse probabilities within
			    each relationship dpf[1..M][0..lmax+1] */

  unsigned char *x ;  /* the true binary state vector x[n] */

  int count ; /* count of number of differences between 
		 the state and the true answer */
  int count_high ; /* number of high bits in the s vector */
  int count_s ;    /* number of high bits in the true s vector */
  int true_s ;   /* whether the true s is supplied */
  double v ;  /* the score of the final state relative to the truth */
  int count_viol ; /* number of checks z that aren't coming out right */
  int not_decoded ; 
  bnd_control *c ; 
} bnd_param;


/* 
         memory and so forth 
                                     */

void bnd_allocate ( bnd_param * , bnd_control * ) ;
void bnd_free ( bnd_param * , bnd_control * ) ;
void bnd_defaults ( bnd_param * , bnd_control * ) ;
void   bnd_load_dqc ( bnd_param *p , bnd_control *c ) ;
void bnd_fprint_state ( bnd_param *p , bnd_control *c ) ;

/*
        the main routine
                                     */

int bndecode ( bnd_param * , bnd_control * ) ;

void bnd_score_state ( bnd_param *p ) ;
void horizontal_pass ( bnd_param *p ) ;
void vertical_pass ( bnd_param *p  , bnd_control *c ) ;

/* recently added routines to prevent overflow  and allow tweaking */
int fudge ( double * , double  ) ;
void fudges ( double * , double  ) ;
double d_2_logit ( double  );
double logit_2_d ( double  );

/*
<!-- hhmts start -->
Last modified: Tue Aug  6 17:10:27 1996
<!-- hhmts end -->
*/
#endif
