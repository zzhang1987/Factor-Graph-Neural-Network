/* My gauss_elim routines */

typedef struct {
  unsigned char **m ; /* the matrix (this copy gets munged */
  unsigned char **mo; /* original copy -- unchanged, except by the 
		       modify_cmatrix_row function */
  unsigned char **mi; /* stores the inverse */
  unsigned char **mt; /* used to assist in undoing permutation */
  int *perm ;  /* the permutation */
  int *iperm ;  /* and its inverse */
  int i ;       /* what row we have got up to in the LU decomposition */
  int l , N ;   /* e.g. 1 , N */
} cm_inversion ;

void allocate_cm_inversion 
( unsigned char **mo,
 int l,
 int N,
 unsigned char **mi,
 cm_inversion *p
)
{
  int i,j ; 

  p->l = l ; p->N = N ; p->mo = mo ; p->mi = mi ; 
  p->m = cmatrix ( l , N , l , N ) ; 
  p->mt=(unsigned char **)malloc((unsigned) (N-l+1)*sizeof(unsigned char*));
  p->mt -= l ;
  p->perm = ivector ( l , N ) ; 
  p->iperm = ivector ( l , N ) ; 
  p->i = l ; 
  for (i=l; i<=N; i++){
    p->perm[i] = l-1 ;  /* deliberately off scale */
    p->iperm[i] = l-1 ;
    for (j=l; j<=N; j++) {
      mi[i][j] = 0 ; 
      p->m[i][j] = mo[i][j] ;
    }
    mi[i][i] = 1 ;
  }
}

void print_cm_inversion 
( cm_inversion *p
)
{
/*  printf( "mo:\n" ) ; 
  printoutcmatrix (  p->mo , p->l , p->N , p->l , p->N ) ; */
  printf( "m:\n" ) ; 
  printoutcmatrix (  p->m , p->l , p->N , p->l , p->N ) ; 
  printf( "mi:\n" ) ; 
  printoutcmatrix (  p->mi , p->l , p->N , p->l , p->N ) ; 
  printf( "i = %d\n" , p->i ) ;
  printf( "perm:\n" ) ; 
  printoutivector ( p->perm , p->l , p->N ) ; 
  printf( "iperm:\n" ) ; 
  printoutivector ( p->iperm , p->l , p->N ) ; 
  pause_for_return () ;
}

void free_cm_inversion 
( cm_inversion *p
)
{
  free_ivector ( p->perm , p->l , p->N ) ; 
  free_ivector ( p->iperm , p->l , p->N ) ; 
  free_cmatrix ( p->m , p->l  , p->N , p->l , p->N ) ;
  free( ( char * ) ( p->mt + p->l ) ) ;
}

/* invert_cmatrix attempts to invert p->m, or to continue to invert it, 
   starting from the state indicated by i,j,perm,etc. 
   if invert_cmatrix succeeds it returns 1, else 0, with p->i 
   recording the row at which non-full-rank was detected . 
*/
int invert_cmatrix
(
 cm_inversion *p
)
{
  int j , done , i2 , k ; 
  unsigned char *irow ;

  for ( ; p->i <= p->N ; p->i ++ ) {
    irow = p->m[p->i] ;
    /* search for a 1 in this row */
    for ( j = p->l , done = 0 ; j <= p->N ; j++ ) {
      if ( irow[j] &&            /* if we find a 1              and */
	  p->perm[j] < p->l ) {  /* if j is not already on the list */
	p->perm[j]  = p->i ; 
	p->iperm[p->i] = j ; 
	done = 1 ;
	break ; /* break out of j loop */
      }
    }
    if ( !done ) {
      return ( 0 ) ;
      break ;
    }
    else { /* search for 1s in subsequent rows. */
      for ( i2 = p->i + 1 ; i2 <= p->N ; i2 ++ ) {
	if ( p->m[i2][j] ) {
	  for ( k = p->l ; k <= p->N ; k ++ ) {
	    p->m[i2][k] ^= irow[k] ;
 	    p->mi[i2][k] ^= p->mi[p->i][k] ; 
	  }
	}
      }
    }
  }      
  /* at this point, if done==1 then m contains a (permuted) upper
     triangular matrix and mi contains a (permuted) lower triangular
     matrix. 
     */
/*  printf ("lu decomposition achieved\n" ) ; */
/*  print_cm_inversion ( p ) ; */
/*  printf ("now working back\n" ) ;  */
  done = invert_utriangularc ( p ) ; 
  /* finally, undo the permutation - make mi[l] = mi[perm[l]] etc. */
  undo_cm_perm ( p ) ; 
  return ( done ) ; 
}

void undo_cm_perm ( cm_inversion *p ) {
  int i , pi ;
  
  for ( i = p->N ; i >= p->l ; i-- ) {
    p->mt[i] = p->mi[i] ;                 /* copy all columns */
  } 
  for ( i = p->N ; i >= p->l ; i-- ) {    /* reallocate them */
    pi = p->perm[i] ;
    if ( ( pi >= p->l ) && ( pi <= p->N ) )
      p->mi[i] = p->mt[pi] ;
    else {
      fprintf ( stderr , "error in undo_cm_perm, perm[%d] = %d\n" , i , 
	       pi ) ; 
      exit ( 0 ) ; 
    }
  } 
}

/* invert_ltriangularc takes p->m and p->mi 
   and applies a series of row sums to both 
   such that p->m becomes the identity.
   If p->mi started out as the identity, it becomes the inverse 
   of p->m; otherwise it becomes (that inverse)*(whatever it was before) */
int invert_utriangularc 
( cm_inversion *p ) {
  int i , vj , j  , k ;
   
  for ( i = p->N ; i >= p->l ; i-- ) {
    for ( vj = i + 1 ; vj <= p->N ; vj++ ) { /* virtual j */
      j = p->iperm[vj] ;
      if (!( ( j >= p->l ) && ( j <= p->N ) ))
      {
	fprintf ( stderr , "error in invert_utriangularc, iperm[%d] = %d\n" ,
		 vj ,  j ) ; 
	exit ( 0 ) ; 
      }
      if ( p->m[i][j] ) {
	for ( k = p->l ; k <= p->N ; k ++ ) {
	  p->m[i][k] ^= p->m[vj][k] ; 
	  p->mi[i][k] ^= p->mi[vj][k] ; 
	} 
/*	printf("%d,%d -> \n",i,j);
	print_cm_inversion ( p ) ; */
      }
    }
  }
  return 1 ; /* inverse has been evaluated */
}


int modify_cmatrix_row ( cm_inversion *p ) {
  /* proceed along the current row seeking an available j */
  unsigned char *irow ;
  int j , done = 0 ;

  irow = p->m[p->i] ;
  /* search for a 1 in this row */
  for ( j = p->l ; done == 0 && j <= p->N ; j++ ) {
    if ( p->perm[j] < p->l ) {  /* if j is not already on the list */
      irow[j] = 1 ; 
      p->mo[p->i][j] ^= 1 ; /* modify original matrix */
      done = 1 ;
    }
  }
  return ( done ) ; /* if done == 0 then this 
		       routine has failed to
		       find a j; an unexpected
		       outcome */
}

/* 
   pbm format is 
   P1
   30 90
   (30 rows with 90 cols)

   which gives (under xv) a picture 
   with 90 rows and 30 cols.

   So in fact rows and columns are interchanged
*/
void cmatrix2pbm
(
 unsigned char **m,
 int l1,
 int h1,
 int l2,
 int h2,
 FILE *fp 
)
{
  int i,j;
  
  fprintf( fp , "P1\n%d %d\n" , h1 - l1 + 1 , h2 - l2 + 1 ) ;
  for (i=l1; i<=h1; i++){
    for (j=l2; j<=h2; j++)
      fprintf( fp , "%d ",m[i][j]);
    fprintf( fp , "\n");
  }                   
}

void printoutcmatrix 
(
 unsigned char **m,
 int l1,
 int h1,
 int l2,
 int h2
)
{
  int i,j;
  
  for (i=l1; i<=h1; i++){
    if ( !(i%10) ) printf ( "| " ) ; 
    else printf ( "  " ) ; 
  }
  printf("\n");
  for (i=l1; i<=h1; i++){
    for (j=l2; j<=h2; j++)
      printf("%d ",m[i][j]);
    printf("\n");
  }                   
}

void printoutcmatrix1 
(
 unsigned char **m,
 int l1,
 int h1,
 int l2,
 int h2
)
{
  int i,j;
  
  for (i=l1; i<=h1; i++){
    for (j=l2; j<=h2; j++)
      if ( m[i][j] ) printf("1") ; else printf(" ") ;
    printf("\n");
  }                   
}

void printoutivector
(
 int *m,
 int l1,
 int h1
)
{
  int i;
  
  for (i=l1; i<=h1; i++)
    printf( "%d " , m[i] );
  printf("\n");
}

void write_ivector
(
 FILE *fp,
 int *m,
 int l1,
 int h1
)
{
  int i;
  
  for (i=l1; i<=h1; i++)
    fprintf( fp , "%d " , m[i] );
  fprintf( fp , "\n");
}

void write_cvector
(
 FILE *fp ,
 unsigned char *m,
 int l1,
 int h1
)
{
  int i;
  
  for (i=l1; i<=h1; i++)
    fprintf( fp , "%d " , m[i] );
  fprintf( fp ,"\n");
