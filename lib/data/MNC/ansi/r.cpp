/* My routines */

#include "r.h"
#include "rand.h"

/* input */

void inputf( float *pointer )
{
  int junk ;
  scanf( "%f" , pointer ) ;
  junk = getchar( ) ;
}
void inputd( double *pointer )
{
  int junk ;
  scanf( "%lf" , pointer ) ;
  junk = getchar( ) ;
}
void inputc( unsigned char *pointer )
{
  int junk ;
  scanf( "%c" , pointer ) ;
  junk = getchar( ) ;
}
void inputi( int *pointer )
{
  int junk ;
  scanf( "%d" , pointer ) ;
  junk = getchar( ) ;
}

void inputrf( float *pointer )
{
  int junk ;
  scanf( "%f" , pointer ) ;
  do{}while( ( junk = getchar( ) ) != 10 ) ;
}
void inputrd( double *pointer )
{
  int junk ;
  scanf( "%lf" , pointer ) ;
  do{}while( ( junk = getchar( ) ) != 10 ) ;
}
void inputrc( unsigned char *pointer )
{
  int junk ;
  scanf( "%c" , pointer ) ;
  do{}while( ( junk = getchar( ) ) != 10 ) ;
}
void inputri( int *pointer )
{
  int junk ;
  scanf( "%d" , pointer ) ;
  do{}while( ( junk = getchar( ) ) != 10 ) ;
}

void clearscan( )
{
  int junk ;
  do{}while( ( junk = getchar( ) ) != 10 ) ;
}

void typeindvector( double *w , int l , int h )
{
  int i ;
  
  printf( "Please enter %d components \n" , h-l+1 ) ;
  for ( i = l ;i<= h ;i++ )
    scanf( "%lf" , &w[i] ) ;
  clearscan( ) ;
}

void constantdmatrix ( double **w , int l1 , int h1 , int l2, int h2 , 
		      double c )
{
  int i , j ;
  
  for ( i = l1 ; i <= h1 ; i++ )
    for ( j = l2 ; j <= h2 ; j++ )
      w[i][j] = c ;
}

void set_dvector_const( double *w , int l , int h , double c )
{
  int i ;
  
  for ( i = l ; i<= h ;i++ )
    w[i] = c ;
}

void set_dvector_c_dvector( double *w , int l , int h , double c , double *v )
{
  int i ;
  
  for ( i = l ; i<= h ;i++ )
    w[i] = c * v[i] ;
}

void set_ivector_const( int *w , int l , int h , int c )
{
  int i ;
  
  for ( i = l ; i<= h ;i++ )
    w[i] = c ;
}

void typeindmatrix ( double **b , int l1 , int h1 , int l2 , int h2 )
{
  int	i , j ;
  
  printf( "please type matrix\n" ) ;
  for ( i = l1 ; i<= h1-1 ; i++ )
    for ( j = l2 ; j<= h2 ; j++ )
      inputd( &b[i][j] ) ; /* This is all a trick to obtain clearscan at the last read-in */
  for ( j = l2 ; j<= h2-1 ; j++ )
    inputd( &b[i][j] ) ;
  inputrd( &b[i][j] ) ;
}

void typeincmatrix ( unsigned char **b , int l1 , int h1 , int l2 , int h2 )
{
  int	i , j , t ;
  
  printf( "please type matrix\n" ) ;
  for ( i = l1 ; i<= h1-1 ; i++ )
    for ( j = l2 ; j<= h2 ; j++ ) {
      inputi( &t ) ; b[i][j] = (unsigned char) t ; 
    }
  for ( j = l2 ; j<= h2-1 ; j++ ) {
    inputi( &t ) ; b[i][j] = (unsigned char) t ; 
  }
  inputri( &t ) ; b[i][j] = (unsigned char) t ; 
}

int ***imatrix3( int l1 , int h1 , int l2 , int h2 , int l3 , int h3 )
{
  int ***c ;
  int s1 , i ;
  
  s1 = h1-l1+1 ;
  
  c = ( int *** )malloc( ( unsigned ) s1*sizeof( int ** ) ) - l1 ;
  for ( i = l1 ;i<= h1 ;i++ )
    c[i] = imatrix( l2 , h2 , l3 , h3 ) ;
  return c ;
}

double ***dmatrix3( int l1 , int h1 , int l2 , int h2 , int l3 , int h3 )
{
  double ***c;
  int s1,i;
  
  s1=h1-l1+1;
  
  c=(double ***) malloc((unsigned) s1*sizeof(double **)) - l1;
  for (i=l1;i<=h1;i++)
    c[i]=dmatrix(l2,h2,l3,h3);
  return c;
}

long int ***limatrix3( int l1 , int h1 , int l2 , int h2 , int l3 , int h3 )
{
  long int ***c ;
  int s1 , i ;
  
  s1 = h1-l1+1 ;
  
  c = ( long int *** )malloc( ( unsigned ) s1*sizeof( long int ** ) ) - l1 ;
  for ( i = l1 ;i<= h1 ;i++ )
    c[i] = limatrix( l2 , h2 , l3 , h3 ) ;
  return c ;
}

long int **limatrix( int nrl , int nrh , int ncl , int nch )
{
  int i ;
  long int **m ;
  
  m=( long int ** )malloc( ( unsigned ) ( nrh-nrl+1 )*sizeof( long int* ) ) ;
  if ( !m ) nrerror( "allocation failure 1 in limatrix( )" ) ;
  m -= nrl ;
  
  for( i=nrl ;i<=nrh ;i++ ) {
    m[i]=( long int * )malloc( ( unsigned ) ( nch-ncl+1 )*sizeof( long int ) ) ;
    if ( !m[i] ) nrerror( "allocation failure 2 in limatrix( )" ) ;
    m[i] -= ncl ;
  }
  return m ;
}

int ipower( int a , int b )
{
  int     i=1;
  
  for(;b>0;b--)
    i*=a;
  return(i);
}

double readindmatrix 
(
 double **b ,
 int l1,
 int h1,
 int l2,
 int h2,
 char *file
)
{
	int	i, j;
	FILE    *fp;
	double	sumd=0.0;

	fp = fopen( file, "r" );
	if( !fp )   fprintf( stderr, "No such file: %s\n", file ), exit(0);

	printf( "reading in matrix from %s\n",file );
	for (i=l1; i<=h1; i++){
		for (j=l2; j<=h2; j++){
			fscanf(fp,"%lf ",&b[i][j]);
			sumd += b[i][j];
		}
	}
	fclose( fp );
	printf( "matrix in\n" );
	return(sumd);
}

int readinimatrix 
(
 int **b ,
 int l1,
 int h1,
 int l2,
 int h2,
 char *file
)
{
  int	status = 0;
  FILE    *fp;

  fp = fopen( file, "r" );
  if( !fp )   fprintf( stderr, "No such file: %s\n", file ), exit(0);

  printf( "reading in matrix from %s\n",file );
  status = fread_imatrix ( b , l1 , h1 , l2 , h2 , fp ) ; 

  fclose( fp );
  return status ;
}

int fread_imatrix 
(
 int **b ,
 int l1,
 int h1,
 int l2,
 int h2,
 FILE *fp 
)
{
  int	i, j , status = 0;

  for (i=l1; i<=h1; i++){
    for (j=l2; j<=h2; j++){
      if ( fscanf(fp,"%d ",&b[i][j]) == EOF ) {
	status -- ;
	break;
      }
    }
    if ( status < 0 ) break ;
  }
//  if ( status == 0 ) fprintf( stderr , "%d * %d int matrix in\n" , h1 , h2 );
//  else fprintf( stderr, 
//		  "Warning: readinimatrix failed at component %d\n",i);
  return status ;
}

void readinlumatrix 
(
 double **b,
 int *indx,
 int n,
 char *file
)
{
  int	i, j , status = 0 ;
  FILE    *fp;
  
  fp = fopen( file, "r" );
  if( !fp )   fprintf( stderr, "No such file: %s\n", file ), exit(0);
  
  printf( "reading in matrix\n" );
  for (i=1; i<=n; i++){
    for (j=1; j<=n; j++){
      if ( fscanf(fp,"%lf ",&b[i][j]) == EOF ) {
	status = -1 ;
	break ;
      }
    }
  }
  if ( status == 0 ) {
    for (i=1; i<=n; i++){
      if ( fscanf(fp,"%d ",&indx[i]) == EOF ) {
	status = -1 ;
	break ; 
      }
    }
  }
  fclose( fp );
  // if ( status == 0 ) printf( "lu matrix in\n" );
  // else fprintf( stderr, 
  //		  "Warning: readinlumatrix failed at component %d\n",i);
}

void readindvector(double *w, int lo, int hi, char *file)
{
	int i, status = 0 ;
	FILE    *fp;
	
	fp = fopen( file, "r" );
	if( !fp )   fprintf( stderr, "No such file: %s\n", file ), exit(0);

	for (i=lo;i<=hi;i++) {
	  if ( fscanf(fp,"%lf ",&w[i]) == EOF ) {
	    status = -1 ; 
	    break ;
	  }
	}
	fclose( fp );
	if ( status < 0 ) 
	  fprintf( stderr, 
		  "Warning: readindvector failed at component %d\n",i);
}

int readdvector(double *w, int lo, int hi, char *file)
{
  int i, status = 0 ;
  FILE    *fp;
  
  fp = fopen( file, "r" );
  if( !fp ) {
    fprintf( stderr, "No such file: %s\n", file ) ;
    status -- ;
  } 
  else {
    for ( i = lo ; i <= hi ; i ++ ) {
      if ( fscanf ( fp , "%lf " , &w[i] ) == EOF ) {
	status = -1 ; 
	break ;
      }
    }
    fclose( fp );
    if ( status < 0 ) 
      fprintf( stderr, 
	      "Warning: readdvector failed at component %d\n",i);
  }
  return status ; 
}

int writedvector ( double *w, int lo, int hi, char *file)
{
  int i, status = 0 ;
  FILE    *fp;
  
  fp = fopen( file, "w" );
  if( !fp ) {
    fprintf( stderr, "writedvector---No such file: %s\n", file ) ;
    status -- ;
  } 
  else {
    for ( i = lo ; i <= hi ; i ++ ) {
      fprintf ( fp , "%g " , w[i] ) ;
    }
    fclose( fp );
  }
  return status ; 
}

int writedmatrix ( 
 double **m,
 int l1,
 int h1,
 int l2,
 int h2,
 char *file
)
{
  int status = 0 ;
  FILE    *fp;
  int i,j;
	
  fp = fopen( file, "w" );
  if( !fp ) {
    fprintf( stderr, "writedmatrix---No such file: %s\n", file ) ;
    status -- ;
  } 
  else {
    for (i=l1; i<=h1; i++){
      for (j=l2; j<=h2; j++)
	fprintf ( fp , "%g " , m[i][j] ) ;

      fprintf ( fp , "\n" ) ;
    }
    fclose( fp );
  }
  return status ; 
}

void readinivector(int *w, int lo, int hi, char *file)
{
  int status = 0 ;
  FILE    *fp;
  
  fp = fopen( file, "r" );
  if( !fp )   fprintf( stderr, "No such file: %s\n", file ), exit(0);

  status = fread_ivector ( w , lo , hi , fp ) ; 

  fclose( fp );
}

int fread_ivector(int *w, int lo, int hi, FILE *fp )
{
  int i, status = 0 ;
  
  for (i=lo;i<=hi;i++) {
    if ( fscanf(fp,"%d ",&w[i]) == EOF ) {
      status = -1 ; 
      break ;
    }
  }

/*  if ( status < 0 ) 
    fprintf( stderr, 
	    "Warning: fread_ivector failed at component %d\n",i);
*/
  return status ; 
}

int fread_dvector( double *w, int lo, int hi, FILE *fp )
{
  int i, status = 0 ;
  
  for (i=lo;i<=hi;i++) {
    if ( fscanf(fp,"%lf ",&w[i]) == EOF ) {
      status = -1 ; 
      break ;
    }
  }

/*  if ( status < 0 ) 
    fprintf( stderr, 
	    "Warning: fread_dvector failed at component %d\n",i);
*/
  return status ; 
}

int fread_cvector(unsigned char *w, int lo, int hi, FILE *fp )
{
  int i, status = 0 , bit ;
  
  for (i=lo;i<=hi;i++) {
    if ( fscanf(fp,"%d", &bit ) == EOF ) {
      status = -1 ; 
      break ;
    }
    w[i] = (unsigned char) bit ;
  }

/*  if ( status < 0 ) 
    fprintf( stderr, 
	    "Warning: fread_cvector failed at component %d\n",i);
*/
  return status ; 
}

void printoutimatrix 
(
 int **m,
 int l1,
 int h1,
 int l2,
 int h2
)
{
  write_imatrix ( stdout , m , l1 , h1 , l2 , h2 ) ; 
}

void write_imatrix 
(
 FILE *fp , 
 int **m,
 int l1,
 int h1,
 int l2,
 int h2
)
{
  int i,j;
  
  for (i=l1; i<=h1; i++){
    for (j=l2; j<=h2; j++)
      fprintf( fp , "%d ",m[i][j]);
    fprintf(fp , "\n");
  }                   
}

void write_imatrix2 
(
 FILE *fp , 
 int **m,
 int l1,
 int h1,
 int l2,
 int h2
)
{
  int i,j;
  
  for (i=l1; i<=h1; i++){
    for (j=l2; j<=h2; j++)
      fprintf( fp , "%+2d ",m[i][j]);
    fprintf(fp , "\n");
  }                   
}


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
    fprintf( fp , "%d\n" , m[i] );
/*  fprintf( fp ,"\n"); */
}

void printoutcvector
(
 unsigned char *m,
 int l1,
 int h1
)
{
  int i;
  
  for (i=l1; i<=h1; i++)
    printf( "%d " , m[i] );
  printf("\n");
}

void printoutcvector1
(
 unsigned char *m,
 int l1,
 int h1
)
{
  int i;
  
  printf("|");
  for (i=l1; i<=h1; i++)
    if ( m[i] ) printf( "1" ) ; else printf( " " );
  printf("|\n");
}

void printoutdmatrix 
(
 double **m,
 int l1,
 int h1,
 int l2,
 int h2,
 int style
)
{
	int i,j;
	
	for (i=l1; i<=h1; i++){
		for (j=l2; j<=h2; j++)
        	 	pd(m[i][j],style);
          	printf("\n");
     	}                   
}

void printoutdmatrix3 
(
 double ***m,
 int l1,
 int h1,
 int l2,
 int h2,
 int l3,
 int h3,
 int style
)
{
  int i,j,k;
	
  for (i=l1; i<=h1; i++){
    for (j=l2; j<=h2; j++) {
      for (k=l3; k<=h3; k++)
	pd(m[i][j][k],style);
      printf("\n");
    }                   
    printf("\n");
  }                   
}

void pd
(
 double w,
 int p
) 	/* prints a single double with a specified number of digits */ 
{
	switch(p){
		case(0):		/* Don't print */
			break;
		case(-1):		/* Standard format */
			printf("%g ",w);
			break;
		case(1):
			printf("%1.0f ",w);
			break;
		case(100):
			printf("%1.0f ",w*10);
			break;
		case(11):
			if (w==1.0)printf("%1.0f ",w);
			else printf("  ");
			break;
		case(10): /* 0/1 <- -1/1 */
			printf("%1.0f ",(double)(w+1.0)*0.5);
			break;
		case(2):
			printf("%2.0f ",w);
			break;
		case(29):
			printf("%2f ",w);
			break;
		case(21):
			printf("%2.1f ",w);
			break;
		case(20):
			printf("%2.0f ",w);
			break;
		case(200):
			printf("%2.0f",w*100);
			break;
		case(3):
			printf("%3.1f ",w);
			break;
		case(39):
			printf("%3f ",w);
			break;
		case(30):
			printf("%3.0f ",w);
			break;
		case(300):
			printf("%3.0f ",w*1000);
			break;
		case(31):
			printf("%3.1f ",w);
			break;
		case(32):
			printf("%3.2f ",w);
			break;
		case(4):
			printf("%4.1f ",w);
			break;
		case(49):
			printf("%4f ",w);
			break;
		case(40):
			printf("%4.0f ",w);
			break;
		case(400):
			printf("%4.0f ",w*10000);
			break;
		case(41):
			printf("%4.1f ",w);
			break;
		case(43):
			printf("%4.3f ",w);
			break;
		case(42):
			printf("%4.2f ",w);
			break;
		case(52):
			printf("%5.2f ",w);
			break;
		case(53):
			printf("%5.3f ",w);
			break;
		case(54):
			printf("%5.4f ",w);
			break;
		case(50): 
			printf("%5.0f ",w);
			break;
		case(500): /* for probabilities */
			printf("%5.0f ",w*100000);
			break;
		case(60): 
			printf("%6.0f ",w);
			break;
		case(62):
			printf("%6.2f ",w);
			break;
		case(63):
			printf("%6.3f ",w);
			break;
		case(64):
			printf("%6.4f ",w);
			break;
		case(72):
			printf("%7.2f ",w);
			break;
		case(74):
			printf("%7.4g ",w);
			break;
		case(84):
			printf("%8.4g ",w);
			break;
		case(94):
			printf("%+9.4g ",w);
			break;
		case(76):
			printf("%7.6f ",w);
			break;
		case(5):
			printf("%5g ",w);
			break;
		case(6):
			printf("%6g ",w);
			break;
		case(600):
			printf("%6.0f ",w*1000000);
			break;
		case(7):
			printf("%7g ",w);
			break;
		case(700):
			printf("%7.0f ",w*10000000);
			break;
		default:
			printf("Error in pd rule \n");
			break;
	}
}

void pdv
(
 double *w, int m, int n, int p
 )	/* prints dvector without newline, with specified number of digits */
{
	int i;
	for(i=m;i<=n;i++)
		pd(w[i],p);
}
 

void 	pause_for_return()
{
	int c;

	printf( "press return to continue" );
	do{}while((c=getchar())!=10);
}

double gammln(
	      double xx
)
{
  double x,tmp,ser;
  static double cof[6]={76.18009173,-86.50532033,24.01409822,
                -1.231739516,0.120858003e-2,-0.536382e-5};
  int j;
  
  x=xx-1.0;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.0;
  for (j=0;j<=5;j++) {
    x += 1.0;
    ser += cof[j]/x;
  }
  return -tmp+log(2.50662827465*ser);
}

int find_rank ( double act , double *v , int lo , int hi ) {
  /* given an ordered list v[lo] < v[lo+1] < ... < v[hi], 
     and assuming that act > v[lo], the task is to return 
     the ranking of act, i.e. the integer place where it will push 
     in. For example, if we find act < v[lo+1], the answer is lo. 
     If act > v[hi], answer is hi. 
     Note no comparison with v[lo] should occur, we already know it's
     bigger. 
     */

  int query ;

  if ( lo > hi ) { 
    fprintf ( stderr , "Eek, lo > hi ???\n" ) ; 
    exit ( 0 ) ; 
  }
  if ( lo == hi ) { /* we are done */
    return ( lo ) ; 
  } 
  else { 
    query = lo + ( hi + 1 - lo ) / 2 ; 
    if ( query == lo ) { 
      fprintf ( stderr , "Eek, query == lo ???\n" ) ; 
      exit ( 0 ) ; 
    }
    else if ( act > v[query] ) { 
      return ( find_rank ( act , v , query , hi ) ) ;
    }
    else {
      return ( find_rank ( act , v , lo , query - 1 ) ) ;
    }
  }
} 

int dotprod_mod2 ( int *a , int *b , int lo , int hi ) {
  int i , ans = 0 ; 
  for ( i = lo ; i <= hi ; i ++ ) {
    ans ^= a[i] * b[i] ;
  }
  return ans ;
}
int idotprod ( int *a , int *b , int lo , int hi ) {
  int i , ans = 0 ; 
  for ( i = lo ; i <= hi ; i ++ ) {
    ans += a[i] * b[i] ;
  }
  return ans ;
}

unsigned char cdotprod_mod2 ( unsigned char *a , unsigned char *b , int lo , int hi ) {
  unsigned char ans = 0 ; 
  int i ; 
  for ( i = lo ; i <= hi ; i ++ ) {
    ans ^= a[i] & b[i] ;
  }
  return ans ;
}

void mult_cms 
( unsigned char **A , unsigned char **B , unsigned char **C , int l , int N ) {
  int i , j , k ; 
  unsigned char a ;

  for ( i = l ; i <= N ; i ++ ) {
    for ( k = l ; k <= N ; k ++ ) {
      a = 0 ; 
      for ( j = l ; j <= N ; j ++ ) {
	a ^= A[i][j] & B[j][k] ;
      }
      C[i][k] = a ;
    }
  }
}

void mult_cm_cv 
( unsigned char **A , unsigned char *b , unsigned char *c , int l , int N ) {
  int i , j ; 
  unsigned char a ;

  for ( i = l ; i <= N ; i ++ ) {
    a = 0 ; 
    for ( j = l ; j <= N ; j ++ ) {
      a ^= A[i][j] & b[j];
    }
    c[i] = a ;
  }
}

int cdotprod ( unsigned char *a , unsigned char *b , int lo , int hi ) {
  int i , ans = 0 ; 
  for ( i = lo ; i <= hi ; i ++ ) {
    ans += a[i] & b[i] ;
  }
  return ans ;
}

/* For useful routines for dealing with triangular char matrices and 
   with lists of matrices of the form
   U = 10
       11
   (just one off diagonal 1)
   see
   ~/code/old/mnc.c.940106
*/


#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

float ran3(
	   int *idum
)
{
	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj,mk;
	int i,ii,k;

	if (*idum < 0 || iff == 0) {
		iff=1;
		mj=MSEED-(*idum < 0 ? -*idum : *idum);
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1;i<=54;i++) {
			ii=(21*i) % 55;
			ma[ii]=mk;
			mk=mj-mk;
			if (mk < MZ) mk += MBIG;
			mj=ma[ii];
		}
		for (k=1;k<=4;k++)
			for (i=1;i<=55;i++) {
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
		inext=0;
		inextp=31;
		*idum=1;
	}
	if (++inext == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if (mj < MZ) mj += MBIG;
	ma[inext]=mj;
	return mj*FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC
