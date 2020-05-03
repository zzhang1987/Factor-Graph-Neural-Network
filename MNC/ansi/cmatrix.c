/* cmatrix.c
   
   My routines for handling funny representations of unsigned char 
   matrices and doing modulo 2 arithmetic. 

   Some routines are in r.c, but others are in here 

   See cmatrix.h also

*/

#include "./r.h"
#include "./rand2.h"
#include "./mynr.h"
#include "./cmatrix.h"

void alist_times_cvector_mod2 
( alist_matrix *a , unsigned char *x , unsigned char *y ) {
  int m ; 
  for ( m = 1 ; m <= a->M ; m++ ) {
    y[m] = cdotprod_mod2_index ( a->mlist[m] , x , 1 , a->num_mlist[m] ) ; 
  }
}

/* version of this for sparse x that runs through x
   only doing anything where x is 1
*/
void alist_times_cvector_sparse_mod2 
( alist_matrix *a , unsigned char *x , unsigned char *y ) {
  int n , m , i ; 
  int *nlist ; 

  for ( m = 1 ; m <= a->M ; m++ ) {
    y[m] = 0 ; 
  }
  for ( n = 1 ; n <= a->N ; n++ ) {
    if ( x[n] ) {
      nlist = a->nlist[n] ;
      for ( i = a->num_nlist[n] ; i >= 1 ; i -- ) {
	y[ nlist[i] ] ^= 1 ;
      }
    }
  }
}

/* version of this for sparse x that runs through x
   only doing anything where x is 1
*/
void alist_transpose_cvector_sparse_mod2 
( alist_matrix *a , unsigned char *x , unsigned char *y ) {
  int m , n , i ; 
  int *mlist ; 

  for ( n = 1 ; n <= a->N ; n++ ) {
    y[n] = 0 ; 
  }
  for ( m = 1 ; m <= a->M ; m++ ) {
    if ( x[m] ) {
      mlist = a->mlist[m] ;
      for ( i = a->num_mlist[m] ; i >= 1 ; i -- ) {
	y[ mlist[i] ] ^= 1 ;
      }
    }
  }
}

/* version of this for sparse x that runs through x
   only doing anything where x is 1
*/
void alist_times_ivector_sparse
( alist_matrix *a , int *x , int *y ) {
  int n , m , i ; 
  int *nlist ; 

  for ( m = 1 ; m <= a->M ; m++ ) {
    y[m] = 0 ; 
  }
  for ( n = 1 ; n <= a->N ; n++ ) {
    if ( x[n] ) {
      nlist = a->nlist[n] ;
      for ( i = a->num_nlist[n] ; i >= 1 ; i -- ) {
	y[ nlist[i] ] ++ ;
      }
    }
  }
}

/* version of this for sparse x that runs through x
   only doing anything where x is 1
*/
void alist_transpose_ivector_sparse
( alist_matrix *a ,  int *x , int *y ) {
  int m , n , i ; 
  int *mlist ; 

  for ( n = 1 ; n <= a->N ; n++ ) {
    y[n] = 0 ; 
  }
  for ( m = 1 ; m <= a->M ; m++ ) {
    if ( x[m] ) {
      mlist = a->mlist[m] ;
      for ( i = a->num_mlist[m] ; i >= 1 ; i -- ) {
	y[ mlist[i] ] ++ ;
      }
    }
  }
}

unsigned char cdotprod_mod2_index 
( 
 int *a , unsigned char *b , int lo , int hi ) {
  unsigned char ans = 0 ; 
  int i ; 
  for ( i = lo ; i <= hi ; i ++ ) {
    ans ^= b[a[i]] ;
  }
  return ans ; 
}


int cdotprod_index 
( 
 int *a , unsigned char *b , int lo , int hi ) {
  int ans = 0 ; 
  int i ; 
  for ( i = lo ; i <= hi ; i ++ ) {
    ans += b[a[i]] ;
  }
  return ans ;
}

void initialize_alist ( alist_matrix *a , int NN , int MM , 
			      int biggestn , int biggestm ) {
  int mm , nn ; 

  a->N = NN ; a->M = MM ; 
  a->biggest_num_m = a->biggest_num_m_alloc = biggestm ;
  a->biggest_num_n = a->biggest_num_n_alloc = biggestn ;
    
  a->num_mlist = ivector ( 1 , MM ) ; 
  a->mlist = imatrix ( 1 , MM , 1 , a->biggest_num_m ) ; 
  a->num_nlist = ivector ( 1 , NN ) ; 
  a->nlist = imatrix ( 1 , NN , 1 , a->biggest_num_n ) ; 
  a->l_up_to = ivector ( 1 , MM ) ; 
  
  for ( mm = 1 ; mm <= MM ; mm++ ) {
    a->num_mlist[mm] = 0 ; 
  }
  for ( nn = 1 ; nn <= NN ; nn++ ) {
    a->num_nlist[nn] = 0 ; 
  }
}

void free_alist ( alist_matrix *a ) {
  int NN = a->N , MM = a->M ; 
  free_ivector ( a->num_mlist , 1 , MM ) ; 
  free_imatrix ( a->mlist , 1 , MM , 1 , a->biggest_num_m_alloc ) ; 
  free_ivector ( a->num_nlist , 1 , NN ) ; 
  free_imatrix ( a->nlist , 1 , NN , 1 , a->biggest_num_n_alloc ) ; 
  free_ivector ( a->l_up_to , 1 , MM ) ; 
}

void add_to_alist ( alist_matrix *a , int nn , int mm ) {
  a->num_nlist[nn] ++ ;
  a->num_mlist[mm] ++ ;
  a->mlist[mm][ a->num_mlist[mm] ] = nn ; 
  a->nlist[nn][ a->num_nlist[nn] ] = mm ;
}

void subtract_from_alist ( alist_matrix *a , int nn , int u ) {
  int m , l , ll ;
  /* identify m and ll */
  m = a->nlist[nn][u] ; 
  ll = 0 ; 
  for ( l = 1 ; l <= a->num_mlist[m] ; l ++ ) {
    if ( a->mlist[m][l] == nn ) {
      ll = l ;
      break ; 
    }
  }
  if ( ll == 0 ) fprintf ( stderr , "problem in alist subtract %d %d\n" , nn , u ) ; 
  /* shuffle back rest of list */
  list_slide_back ( a->mlist[m] , ll , a->num_mlist[m] ) ;
  list_slide_back ( a->nlist[nn] , u , a->num_nlist[nn] ) ;

  a->num_nlist[nn] -- ;
  a->num_mlist[m] -- ;
}

void list_slide_back ( int *list , int bot , int top ) {
  if ( bot < top ) list[bot] = list[top] ;
  list[top] = 0 ; 
}

void subtract_col_from_alist ( alist_matrix *a , int nn ) {
  int u ; 
  for ( u = a->num_nlist[nn] ; u >= 1 ; u -- ) {
    subtract_from_alist ( a , nn , u ) ;
  }
}

void kill_blank_col_from_alist ( alist_matrix *a , int nn ) {
  int u , N = a->N , m , l ; 
  int count ; 

  a->same_length = 0 ; 

  if ( a->num_nlist[nn] ) fprintf ( stderr , "warning, nonempty column %d\n" , nn ) ; 
  if ( nn == a->N ) { /* no need to do anything */
  } else {            /* move last column into this one's spot */
     u = a->num_nlist[N] ;
     a->num_nlist[nn] = u ; 
     for ( ; u >= 1 ; u -- ) { 
       m = a->nlist[N][u] ;
       a->nlist[nn][u] = m ; 
       /* get item m, and change his entry 
	  mlist[m][l] from N to nn */
       count = 0 ; 
       for ( l = a->num_mlist[m] ; l >= 1 ; l -- ) {
	 if ( a->mlist[m][l] == N ) { a->mlist[m][l] = nn ; count++ ; }
       }
       if ( count != 1 ) fprintf ( stderr , "warning, expected 1 count, but at %d it's %d\n" , nn , count ) ; 
     }
  }
  a->N -- ; 
}

void write_alist ( FILE *fp , alist_matrix *a ) {
  /* this assumes that mlist and nlist have the form of a rectangular
     matrix in the file; if lists have unequal lengths, then the 
     entries should be present but are ignored
     */
  int N = a->N , M = a->M ; 

  fprintf ( fp , "%d %d\n" , N , M ) ; 
  fprintf ( fp , "%d %d\n" , a->biggest_num_n , a->biggest_num_m ) ; 
  write_ivector ( fp , a->num_nlist , 1 , N ) ; 
  write_ivector ( fp , a->num_mlist , 1 , M ) ; 
  write_imatrix ( fp , a->nlist , 1 , N , 1 , a->biggest_num_n ) ; 
  write_imatrix ( fp , a->mlist , 1 , M , 1 , a->biggest_num_m ) ; 
}  

void write_alist_transpose ( FILE *fp , alist_matrix *a ) {
  int N = a->N , M = a->M ; 

  fprintf ( fp , "%d %d\n" , M , N ) ; 
  fprintf ( fp , "%d %d\n" , a->biggest_num_m , a->biggest_num_n ) ; 
  write_ivector ( fp , a->num_mlist , 1 , M ) ; 
  write_ivector ( fp , a->num_nlist , 1 , N ) ; 
  write_imatrix ( fp , a->mlist , 1 , M , 1 , a->biggest_num_m ) ; 
  write_imatrix ( fp , a->nlist , 1 , N , 1 , a->biggest_num_n ) ; 
}  

int read_allocate_alist ( alist_matrix *a , char *file ) { 
  /* this assumes that mlist and nlist have the form of a rectangular
     matrix in the file; if lists have unequal lengths, then the 
     entries should be present but are ignored
     */
  int status = 0 ; 
  int N , M ; 
  int biggestn , biggestm ;
  FILE    *fp;
  
  fp = fopen( file, "r" );
  if( !fp )   fprintf( stderr, "No such file: %s\n", file ), exit(0);

  do {
    if ( fscanf(fp,"%d %d " , &N , &M ) == EOF ) {
      status = -1 ;       break ;
    }
    if ( fscanf(fp,"%d %d " , &biggestn , &biggestm ) == EOF ) {
      status = -1 ;       break ;
    }

    initialize_alist ( a , N , M ,  biggestn ,  biggestm ) ;

    status += fread_ivector ( a->num_nlist , 1 , N , fp ) ; 
    status += fread_ivector ( a->num_mlist , 1 , M , fp ) ; 
    status += fread_imatrix ( a->nlist , 1 , N , 1 , biggestn , fp ) ;
    status += fread_imatrix ( a->mlist , 1 , M , 1 , biggestm , fp ) ;
  } while ( 0 ) ; 

  return status ; 
}

int code0_matrix_alist ( alist_matrix *a , int N , int M , int tr ) {
  int m , n , num_left ; 
  int status = 0 ; 
  
  if ( M < N ) {
    fprintf ( stderr , "code0 matrix problem with N and M\n" ) ; 
    pause_for_return () ;
    exit (0) ;
  }
  for ( n = 1 ; n <= N ; n++ ) {
    m = n ; 
    add_to_alist ( a , n , m ) ;
  }
  for ( m = 1 ; m <= M ; m++ ) { /* exact number per row */
    n = N ;
    num_left = tr ;
    for ( ; n >= 1 ; n-- ) {
      if ( num_left > 0 && ( ranf() <=  (double) num_left / (double) n ) ) {
	add_to_alist ( a , n , m ) ;
	num_left -- ;
      } 
    }
  }
  status += finish_off_alist ( a ) ;
  return status ; 
}

/*
   square matrices only 
*/
int cmatrix_2_alist ( unsigned char **A , alist_matrix *a , int NN ) {
  int mm , nn , uu ; 
  int status = 0 ; 
  
  for ( nn = 1 ; nn <= NN ; nn++ ) {
    for ( uu = 1 , mm = 1 ; mm <= NN ; mm++ ) {
      if ( A[mm][nn] ) {
	a->nlist[nn][uu] = mm ;
	uu ++ ; 
	a->num_mlist[mm] ++ ;
	a->mlist[mm][ a->num_mlist[mm] ] = nn ; 
      }
    }
    uu -- ;
    a->num_nlist[nn] = uu ;
  }
  status += finish_off_alist ( a ) ;
  return status ; 
}

int gen_cmatrix_2_alist ( unsigned char **A , alist_matrix *a , int M , int N , int transpose ) {
  int mm , nn , uu ; 
  int status = 0 ; 
  
  for ( nn = 1 ; nn <= N ; nn++ ) {
    for ( uu = 1 , mm = 1 ; mm <= M ; mm++ ) {
      if ( ( transpose && A[nn][mm] ) || ( !transpose && A[mm][nn] ) ) {
	a->nlist[nn][uu] = mm ;
	uu ++ ; 
	a->num_mlist[mm] ++ ;
	a->mlist[mm][ a->num_mlist[mm] ] = nn ; 
      }
    }
    uu -- ;
    a->num_nlist[nn] = uu ;
  }
  status += finish_off_alist ( a ) ;
  return status ; 
}

/* 
   puts them side by side if vertical = 0, 
   or on top of each other if 1.
                                                */
int cmatrices_2_alist ( unsigned char **A ,  unsigned char **B ,
			       alist_matrix *a , int NN , int vertical ) {
  int mm , nn , mat , tn , tm , moffset ; 
  unsigned char **C ;
  int status = 0 ; 

  tn = 1 ;  moffset = 0 ; 
  for ( C = A , mat = 1 ; 
       mat <= 2 ; 
       ) {
    for ( nn = 1 ; nn <= NN ; nn++ , tn++ ) {
      for ( mm = 1 ; mm <= NN ; mm++ ) {
	if ( C[mm][nn] ) {
	  tm = mm + moffset ;
	  add_to_alist ( a , tn , tm ) ; 
	}
      }
    }
    /* end of loop */
    mat++ ; C = B ; 
    if ( vertical ) {
      tn = 1 ;
      moffset = NN ;
    }
  }
  status += finish_off_alist ( a ) ;
  return status ; 
}

/*
   puts identity matrix alongside the provided matrix 
   [ C I ]
*/
int cmatrix_andI_2_alist ( unsigned char **A ,  
			       alist_matrix *a , int NN , int vertical ) {
  int mm , nn ; 
  int status = 0 ; 

  for ( nn = 1 ; nn <= NN ; nn++ ) {
    for ( mm = 1 ; mm <= NN ; mm++ ) {
      if ( A[mm][nn] ) {
	add_to_alist ( a , nn , mm ) ; 
      }
    }
  }
  if ( !vertical ) {
    for ( mm = 1, nn=NN+1 ; mm <= NN ; mm++ , nn++ ) {
      add_to_alist ( a , nn , mm ) ;
    }
  } else {
    for ( mm = NN+1, nn=1 ; nn <= NN ; mm++ , nn++ ) {
      add_to_alist ( a , nn , mm ) ;
    }
  }

  status += finish_off_alist ( a ) ;
  return status ; 
}

int finish_off_alist ( alist_matrix *a ) {
  int nn , mm , uu ; 
  int tbm = 0  , tbn = 0  ; /* keeps track of the true biggest_m , n */
  int status = 0 ; 

  for ( mm  = 1 ; mm <= a->M ; mm++ ) {
    uu = a->num_mlist[mm] ; 
    if ( uu > tbm ) tbm = uu ; 
  }
  a->tot = 0 ; 
  for ( nn  = 1 ; nn <= a->N ; nn++ ) {
    uu = a->num_nlist[nn] ; 
    if ( uu > tbn ) tbn = uu ; 
    a->tot += uu ;
  }
  if ( tbm >= a->biggest_num_m_alloc ) {
    fprintf (stderr , "eek m list overrun %d %d\n" ,
	     tbm , a->biggest_num_m_alloc ) ;
    status -- ; 
  }
  else a->biggest_num_m = tbm ; 
  if ( tbn >= a->biggest_num_n_alloc ) {
    fprintf ( stderr , "eek n list overrun %d %d\n" , 
	     tbn , a->biggest_num_n_alloc ) ;
    status -- ; 
  }
  else a->biggest_num_n = tbn ; 
  return status ; 
}

void report_alist ( alist_matrix *a , int verbose ) {
  int m , n ; 

  printf ( "nlist numbers :\n" ) ;
  printoutivector ( a->num_nlist , 1 , a->N ) ; 
  printf ( "mlist numbers :\n" ) ;
  printoutivector ( a->num_mlist , 1 , a->M ) ; 
  if ( verbose >= 2  )   {
    printf ( "nlist :\n" ) ;
    for ( n = 1 ; n <= a->N ; n++ ) {
      printf ( "%3d : " , n  ) ;
      printoutivector ( a->nlist[n] , 1 , a->num_nlist[n] ) ; 
    }
    printf ( "\n" ) ; 
    printf ( "mlist :\n" ) ;
    for ( m = 1 ; m <= a->M ; m++ ) {
      printf ( "%3d : " , m  ) ;
      printoutivector ( a->mlist[m] , 1 , a->num_mlist[m] ) ; 
    }
    printf ( "\n" ) ; 
  }
}
 

void sparse_random_cmatrix
( unsigned char **A , int per_row , int NN ) {
  int i , j , num_left  ;

  for ( i = 1 ; i <= NN ; i++ )  {
    j = NN ;
    num_left = per_row ;
    for ( ; j >= 1 ; j-- ) {
      if ( num_left > 0 && ( ranf() <=  (double) num_left / (double) j ) ) {
	A[i][j] = 1 ;
	num_left -- ;
      }  else A[i][j] = 0 ;
    }
  }
}

void sparse_random_cmatrix2 /* this ensures that number_per_col equal too */
( unsigned char **A , int per , int NN ) {
  int i , j , num_left  ;
  int *new_num_c , *num_c ; /* number left in column */
  
  new_num_c = ivector ( 1 , NN ) ; 
  num_c = ivector ( 1 , NN ) ; 
  for ( j = 1 ; j <= NN ; j++ )  
    num_c[j] = per ; 
  for ( i = NN ; i >= 1 ; i-- ) {
    do { /* generate the next row at random until the number in the row is
	    right */
      j = NN ;
      num_left = per ;
      for ( ; j >= 1 ; j-- ) {
	if ( num_c[j] > 0 && ( ranf() <=  (double) num_c[j] / (double) i ) ) {
	  A[i][j] = 1 ;
	  num_left -- ;
	  new_num_c[j] = num_c[j] - 1 ; 
	}  else {
	  A[i][j] = 0 ;
	  new_num_c[j] = num_c[j] ;
	}
      }
      if ( num_left == 0 ) { /* move on */
	for ( j = 1 ; j <= NN ; j++ )  
	  num_c[j] = new_num_c[j] ;
      }
/*      printf ( "%d %d; " , i , num_left ) ; fflush(stdout) ;  */
    } while ( num_left != 0 ) ;
  }
  free_ivector ( new_num_c , 1 , NN ) ;
  free_ivector ( num_c , 1 , NN ) ;
}

void sparse_rectangular_cmatrix ( unsigned char **A , int per , int M , int N , int hd ){
/* makes a sparse rectangular matrix with per per column */
/* and number per row as equal as possible 
   if hd > 0 then all columns must be more than hd distance from each other.
   if a column is hd or closer to another column, then a "h" is printed
   and we restart that column */
  int i , j , num_left  ;
  int *num_r ; /* number left in row */
  int *mylist ; 
  int per_row , extras , hdfail = 0 ;
  int dist = 0  , col , m ;
  double di ;

  /* N columns, M rows */

  per_row = ( per * N ) / M ;  
  extras = per * N - per_row * M ;
  if ( extras < 0 ) fprintf ( stderr , "eek unexpected extras < 0 %d\n" , extras )  ;

  num_r = ivector ( 1 , M ) ; 
  mylist =  ivector ( -1 , per ) ; 

  printf ( "Start:" ) ; fflush(stdout) ;

  for ( j = 1 ; j <= M ; j++ )  
    num_r[j] = per_row ; 
  if ( extras > 0 ) {
    fprintf ( stderr , "Adding %d extras because inexact division\n" , extras ) ; 
    for ( j = M ; j >= 1 && ( extras > 0 ); j-- )  {
      if ( ranf() <=  (double) extras / (double) j ) {
	num_r[j] ++ ;
      }
    }
  }

  for ( i = N ; i >= 1 ; i-- ) {
    di = 1.0 / (double) i ; 
    printf ( "%d;" , i ) ; fflush(stdout) ; 
    do { /* generate the next col at random until the number in the col is
	    right 
	    if numleft has got down to -1, then we can certainly stop
	    and start again */
      j = M ;
      num_left = per ;
      for ( ; j >= 1 && (num_left > -1)  ; j-- ) {
	if ( num_r[j] > 0 && ( ranf() <=  (double) num_r[j] * di ) ) {
	  mylist[num_left] = j ; 
	  num_left -- ;
	} 
      }
      if ( num_left == 0 ) { /* check hd then move on */
	if ( hd > 0 ) {
	  hdfail = 0 ; 
	  for ( col = N ; (col > i) && ( hdfail == 0 ) ; col-- ) {
	    dist = 0 ; 
	    for (  j = 1 ; j <= per && ( dist <= hd ) ; j ++ ) {
	      m = mylist[j] ;
	      dist += A[m][col] ? 0 : 2 ; 
	    }
	    if ( dist <= hd ) { hdfail = 1 ; }
	  }
	}
	if ( hd > 0 && hdfail == 1 ) {
	  printf ( "h" ) ; fflush(stdout) ;
	  if ( i == 1 ) { /* tough luck, have to accept the failure */
	    fprintf ( stderr, "\nwarning, final hamming check failure %d\n" , dist ) ;
	    hdfail = 0 ; 
	  }
	} 
	if ( hdfail == 0 ) {
	  /* final acts before moving on */
	  for ( j = 1 ; j <= per ; j ++ ) {
	    m = mylist[j] ;
	    A[m][i] = 1 ;
	    num_r[m] -- ; 
	  }
	}

      } else {
	printf ( "." ) ; /* fflush(stdout) ;   */
      }
    } while ( num_left != 0 || hdfail != 0 ) ;
  }
  free_ivector ( mylist , -1 , per ) ;
  free_ivector ( num_r , 1 , M ) ;
}

int sparse_rectangular_cmatrix2 ( unsigned char **A , int c0 , int per , int M , int N , int hd ){
/* Two differences from 1st routine.
   The matrix is constructed by ROW, from bottom to top.
   And a start can be made a distance c0 from the bottom if desired 
*/
/*
   if hd > 0 then all rows must be more than hd distance from each other.
   if a row is hd or closer to another row, then a "h" is printed
   and we restart that row 

   This version starts from row N - c0, and compares with all previous
   rows right up to N
*/
  int i , j , num_left  ;
  int *num_r ; /* number left in col */
  int *mylist ; 
  int per_col , extras , hdfail = 0 ;
  int dist = 0 , row , m ;
  double di ;
  int num_rows_to_do = N - c0 ; 
  int status = 0 ; 

  /* N rows, M cols */

  per_col = ( per * num_rows_to_do ) / M ;  
  extras = per * num_rows_to_do - per_col * M ;
  if ( extras < 0 ) fprintf ( stderr , "eek unexpected extras < 0 %d\n" , extras )  ;

  num_r = ivector ( 1 , M ) ; 
  mylist =  ivector ( -1 , per ) ; 

  printf ( "Start:" ) ; fflush(stdout) ;

  for ( j = 1 ; j <= M ; j++ )  
    num_r[j] = per_col ; 
  if ( extras > 0 ) {
    fprintf ( stderr , "Adding %d extras because inexact division\n" , extras ) ; 
    for ( j = M ; j >= 1 && ( extras > 0 ); j-- )  
      if ( ranf() <=  (double) extras / (double) j ) 
	num_r[j] ++ ;
  }

  for ( i = num_rows_to_do ; i >= 1 ; i-- ) {
    di = 1.0 / (double) i ; 
    printf ( "%d;" , i ) ; fflush(stdout) ; 
    do { /* generate the next row at random until the number in the row is
	    right */
      j = M ;
      num_left = per ;
      for ( ; j >= 1 && (num_left > -1)  ; j-- ) {
	if ( num_r[j] > 0 && ( ranf() <=  (double) num_r[j] * di ) ) {
	  mylist[num_left] = j ; 
	  num_left -- ;
	} 
      }
      if ( num_left == 0 ) { /* check hd then move on */
	if ( hd > 0 ) {
	  hdfail = 0 ; 
	  for ( row = N ; (row > i) && ( hdfail == 0 ) ; row-- ) {
	    dist = 0 ; 
	    for (  j = 1 ; j <= per && ( dist <= hd ) ; j ++ ) {
	      m = mylist[j] ;
	      dist += A[row][m] ? 0 : 2 ; /* worst case - assume any 1
					     implies another difference
					     somewhere; correct if both
					     rows have same weight */
	    }
	    if ( dist <= hd ) { hdfail = 1 ; }
	  }
	}
	if ( hd > 0 && hdfail == 1 ) {
	  printf ( "h" ) ; fflush(stdout) ;
	  if ( i == 1 ) { /* tough luck, have to accept the failure */
	    fprintf ( stderr, "\nwarning, final hamming check failure %d\n" , dist ) ;
	    hdfail = 0 ; 
	    status -- ; 
	  }
	} 
	if ( hdfail == 0 ) {
	  /* final acts before moving on */
	  for ( j = 1 ; j <= per ; j ++ ) {
	    m = mylist[j] ;
	    A[i][m] = 1 ;
	    num_r[m] -- ; 
	  }
	}

      } else {
	printf ( "." ) ; /* fflush(stdout) ;  */
      }
    } while ( num_left != 0 || hdfail != 0 ) ;
  }
  free_ivector ( mylist , -1 , per ) ;
  free_ivector ( num_r , 1 , M ) ;
  return status ; 
}

int old_sparse_rectangular_alist ( alist_matrix *a , int c0 , int c1 , int per , int Cols , int Rows , int hd ){
/* 
   Constructs a matrix by row, from bottom to top.
   A start can be made a distance c0 from the bottom if desired 
   A stop can be called after only c1 rows from the bottom have been 
   finished also.  

   if hd > 0 then all rows must be more than hd distance from each other.
   if a row is hd or closer to another row, then a "h" is printed
   and we restart that row 

   This version starts from row Rows - c0, and compares with all previous
   rows right up to Rows
*/
  int status = 0 ; 
  int nn , u , latest ; 
  int i , j , num_left  ;
  int *num_r ; /* number left in col */
  int *mylist ; 
  int per_col , tot_per_col , extras , hdfail = 0 ;
  int m ;
  double di ;
  int num_rows_to_do = c1 - c0 ; 
  int *relatives ; 
  int relmax , rr , ss , prevrr ;
  int *num_nlist = a->num_nlist ; 
  int **nlist = a->nlist ; 

  /* N rows, M cols */

  per_col = ( per * num_rows_to_do ) / Cols ;  
  extras = per * num_rows_to_do - per_col * Cols ;
  if ( extras < 0 ) fprintf ( stderr , "eek unexpected extras < 0 %d\n" , extras )  ;

  num_r = ivector ( 1 , Cols ) ; 
  mylist =  ivector ( -1 , per ) ; 

  tot_per_col = ( per * Rows ) / Cols ;  
  relmax = per * ( tot_per_col + 1 ) ; 
  relatives = ivector ( 1 , relmax ) ; 
  
  printf ( "Start:" ) ; fflush(stdout) ;

  for ( j = 1 ; j <= Cols ; j++ )  
    num_r[j] = per_col ; 
  if ( extras > 0 ) {
    fprintf ( stderr , "Adding %d extras because inexact division\n" , extras ) ; 
    for ( j = Cols ; j >= 1 && ( extras > 0 ); j-- )  
      if ( ranf() <=  (double) extras / (double) j ) 
	num_r[j] ++ ;
  }

  for ( i = num_rows_to_do , m = Rows - c0 ; i >= 1 ; i-- , m-- ) {
    di = 1.0 / (double) i ; 
    printf ( "%d;" , i ) ; fflush(stdout) ; 
    do { 
      j = Cols ;
      num_left = per ;
      for ( ; j >= 1 && (num_left > -1)  ; j-- ) {
	if ( num_r[j] > 0 && ( ranf() <=  (double) num_r[j] * di ) ) {
	  mylist[num_left] = j ; 
	  num_left -- ;
	} 
      }
      if ( num_left == 0 ) { /* check hd then move on */
	if ( hd > 0 ) {
	  hdfail = 0 ; 
/* 
   for each j on the list
   go and find any other rows that share this j
   check if we have already met this row today. If so, then that's 
   a failure for the traditional setting of hd
*/
	  rr = 0 ; 
	  prevrr = 0 ; 
	  for ( j = 1 ; j <= per && ( hdfail == 0 ) ; j ++ ) {
	    nn = mylist[j] ;
	    for ( u = num_nlist[nn] ; u >= 1 && ( hdfail == 0 ) ; u-- ) {
	      rr ++ ;
	      latest = relatives[rr] = nlist[nn][u] ; 
	      for ( ss = prevrr ; ss >= 1 ; ss -- ) {
		if ( relatives[ss] == latest ) {
		  hdfail = 1 ; 
		  break ; 
		}
	      }
	    }
	    prevrr = rr ; 
	    if ( rr > relmax ) {
	      fprintf ( stderr, "\nwarning, relatives overran vector %d %d\n" , rr , relmax ) ;
	      status -- ; 
	    }
	  }
	}
	if ( hd > 0 && hdfail == 1 ) {
	  printf ( "h" ) ; fflush(stdout) ;
	  if ( i == 1 ) { /* tough luck, have to accept the failure */
	    fprintf ( stderr, "\nwarning, final hamming check failure at relative %d %d out of %d\n" , rr , ss , relmax ) ;
	    hdfail = 0 ; 
	  }
	} 
	if ( hdfail == 0 ) {
	  /* final acts before moving on */
	  for ( j = 1 ; j <= per ; j ++ ) {
	    nn = mylist[j] ;
	    add_to_alist ( a , nn , m ) ; 
	    num_r[nn] -- ; 
	  }
	}
      } else {
	printf ( "." ) ; /* fflush(stdout) ;  */
      }
    } while ( num_left != 0 || hdfail != 0 ) ;
  }
  free_ivector ( mylist , -1 , per ) ;
  free_ivector ( num_r , 1 , Cols ) ;
  free_ivector ( relatives , 1 , relmax ) ;
  
  return status ; 
}

int sparse_rectangular_alist ( alist_matrix *a , int c0 , int c1 , int per , int Cols , int Rows , int hd ){
/* 
   Constructs a matrix by row, from bottom to top.
   A start can be made a distance c0 from the bottom if desired 
   A stop can be called after only c1 rows from the bottom have been 
   finished also.  

   if hd > 0 then all rows must be more than hd distance from each other.
   if a row is hd or closer to another row, then a "h" is printed
   and we restart that row 

   This version starts from row Rows - c0, and compares with all previous
   rows right up to Rows
*/
  int status = 0 ; 
  int nn , u , latest ; 
  int i , j , num_left  ;
  int *num_r ; /* number left in col */
  int *pikd_r ;   /* whether this col pikd this time */
  int *cum ;   /* cumulative probability */
  int *mylist ; 
  int per_col , tot_per_col , extras , hdfail = 0 ;
  int m ;
  double di ;
  int num_rows_to_do = c1 - c0 ; 
  int *relatives ; 
  int relmax , rr , ss , prevrr ;
  int *num_nlist = a->num_nlist ; 
  int **nlist = a->nlist ; 
  double un ; /* random 0,1 */

  /* N rows, M cols */

  per_col = ( per * num_rows_to_do ) / Cols ;  
  extras = per * num_rows_to_do - per_col * Cols ;
  if ( extras < 0 ) fprintf ( stderr , "eek unexpected extras < 0 %d\n" , extras )  ;

  num_r = ivector ( 1 , Cols ) ; 
  pikd_r = ivector ( 1 , Cols ) ; 
  cum = ivector ( 0 , Cols ) ; cum[0] = 0 ; 
  mylist =  ivector ( -1 , per ) ; 

  tot_per_col = ( per * Rows ) / Cols ;  
  relmax = per * ( tot_per_col + 1 ) ; 
  relatives = ivector ( 1 , relmax ) ; 
  
  printf ( "Start:" ) ; fflush(stdout) ;

  for ( j = 1 ; j <= Cols ; j++ )  
    num_r[j] = per_col ; 
  if ( extras > 0 ) {
    fprintf ( stderr , "Adding %d extras because inexact division\n" , extras ) ; 
    for ( j = Cols ; j >= 1 && ( extras > 0 ); j-- )  
      if ( ranf() <=  (double) extras / (double) j ) 
	num_r[j] ++ ;
  }

  for ( i = num_rows_to_do , m = Rows - c0 ; i >= 1 ; i-- , m-- ) {
    di = 1.0 / (double) i ; 
    if ( !(i%1000) ) { printf ( "%d;" , i ) ; fflush(stdout) ; }

    /* compute cumulatives */
    for ( j = 1 ; j <= Cols ; j ++ ) 
      cum[j] = cum[j-1] + num_r[j] ; 

    do { /* generate a row */

      for ( j = 1 ; j <= Cols ; j ++ ) 
	pikd_r[j] = 0 ;

      for ( num_left = per ; num_left >= 1 ; num_left -- ) {
	do {
	  un = ranf() ; 
	  j = read_int_from_cum ( cum , 1 , Cols , un ) ; 
	} while ( pikd_r[j] ) ; /* while j already picked */
	pikd_r[j] ++ ; 
	if ( num_r[j] > 0 ) {
	  mylist[num_left] = j ; 
	} else { fprintf (stderr , "cum failure %d %d %d\n" ,
			  j , cum[j-1] , cum[j] ) ; 	  num_left ++ ;
	       pause_for_return(); }
      }
      if ( num_left == 0 ) { /* check hd then move on */
	if ( hd > 0 ) {
	  hdfail = 0 ; 
/* 
   for each j on the list
   go and find any other rows that share this j
   check if we have already met this row today. If so, then that's 
   a failure for the traditional setting of hd
*/
	  rr = 0 ; 
	  prevrr = 0 ; 
	  for ( j = 1 ; j <= per && ( hdfail == 0 ) ; j ++ ) {
	    nn = mylist[j] ;
	    for ( u = num_nlist[nn] ; u >= 1 && ( hdfail == 0 ) ; u-- ) {
	      rr ++ ;
	      latest = relatives[rr] = nlist[nn][u] ; 
	      for ( ss = prevrr ; ss >= 1 ; ss -- ) {
		if ( relatives[ss] == latest ) {
		  hdfail = 1 ; 
		  break ; 
		}
	      }
	    }
	    prevrr = rr ; 
	    if ( rr > relmax ) {
	      fprintf ( stderr, "\nwarning, relatives overran vector %d %d\n" , rr , relmax ) ;
	      status -- ; 
	    }
	  }
	}
	if ( hd > 0 && hdfail == 1 ) {
	  printf ( "h" ) ; fflush(stdout) ;
	  if ( i == 1 ) { /* tough luck, have to accept the failure */
	    fprintf ( stderr, "\nwarning, final hamming check failure at relative %d %d out of %d\n" , rr , ss , relmax ) ;
	    hdfail = 0 ; 
	  }
	} 
	if ( hdfail == 0 ) {
	  /* final acts before moving on */
	  for ( j = 1 ; j <= per ; j ++ ) {
	    nn = mylist[j] ;
	    add_to_alist ( a , nn , m ) ; 
	    num_r[nn] -- ; 
	  }
	}
      } else {
	printf ( "." ) ; /* fflush(stdout) ;  */
      }
    } while ( num_left != 0 || hdfail != 0 ) ;
  }
  free_ivector ( mylist , -1 , per ) ;
  free_ivector ( num_r , 1 , Cols ) ;
  free_ivector ( relatives , 1 , relmax ) ;
  
  return status ; 
}

int staircase_alist ( alist_matrix *a , int c0 , int c1 , 
				     int per , int Cols , 
				     int Rows ){
/* 
   Constructs a matrix by row, from bottom to top.
   Each row is a step of a staircase. 
*/
  int i , n , m , u ;
  int num_rows_to_do = c1 - c0 ; 
  int status = 0 ; 
  n = 1 ;
  printf ( "Stair:" ) ;
  for ( i = num_rows_to_do , m = Rows - c0 ; i >= 1 ; i-- , m-- ) {
    if ( !(i%1000) ) { printf ( "%d;" , i ) ; fflush(stdout) ; }

    for ( u = 1 ; u <= per ; u++ ) {
      add_to_alist ( a , n , m ) ; 
      if ( u < per ) { n ++ ; if ( n > Cols ) n = 1 ; }
    }
  }
  return status ;
}

int slope_alist ( alist_matrix *a , int c0 , int c1 , 
				     int per , int Cols , 
				     int Rows ){
/* 
   Constructs a matrix by row, from bottom to top.
   Each row is  
   1 1 ----------
   --- 1 1 ------
   ------- 1 1 -- if per should be 2.
*/
  int i , n , m , u ;
  int num_rows_to_do = c1 - c0 ; 
  int status = 0 ; 
  n = 1 ;
  printf ( "Slope:" ) ;
  for ( i = num_rows_to_do , m = Rows - c0 ; i >= 1 ; i-- , m-- ) {
    if ( !(i%1000) ) { printf ( "%d;" , i ) ; fflush(stdout) ; }

    for ( u = 1 ; u <= per ; u++ ) {
      add_to_alist ( a , n , m ) ; 
      n ++ ; if ( n > Cols ) n = 1 ; 
    }
  }
  return status ;
}


int uneven_sparse_rectangular_alist ( alist_matrix *a , int c0 , int c1 , 
				     double per , int Cols , 
				     int Rows , int hd ){
/* 
   Constructs a matrix by row, from bottom to top.
   A start can be made a distance c0 from the bottom if desired 
   A stop can be called after only c1 rows from the bottom have been 
   finished also.  

   if hd > 0 then max permitted overlap is hd DIFFERENT FROM WHAT WENT B4
   if a row is hd or closer to another row, then a "h" is printed
   and we restart that row 
*/
  int status = 0 ; 
  int nn , u , latest ; 
  int i , j , num_left  ;
  int *num_c ; /* number left in col */
  int *pikd_c ;   /* whether this col pikd this time */
  int *cum ;   /* cumulative probability */
  int *mylist ; 
  double per_col ;
  int per_col1 , per_col2 , tot_per_col , tot_ones , hdfail = 0 ;
  int extras_r , extras_c ;
  int m ;
  int num_rows_to_do = c1 - c0 ; 
  int *relatives ; 
  int relmax , rr , ss , prevrr ;
  int *num_nlist = a->num_nlist ; 
  int **nlist = a->nlist ; 
  double un ; /* random 0,1 */
  int per1 , per2 , this_per ; 

  /* Compute what numbers per row and column */

  if ( num_rows_to_do == 0 ) { fprintf ( stderr , "0 rows to do\n" ) ;
			       return (status);
			     }
  tot_ones = (int)( per * (double) num_rows_to_do ) ;
  per_col = ( per * (double) num_rows_to_do ) / (double) Cols ;  
  per_col1 = tot_ones / Cols ; 
  per_col2 = per_col1 + 1 ; 
  extras_c = tot_ones - per_col1 * Cols ;  /* that's the number of cols that 
					    have to be allocated per_col2 */
  per1 = tot_ones / num_rows_to_do ; 
  per2 = per_col1 + 1 ; 
  extras_r = tot_ones - per1 * num_rows_to_do ;
  /* that's the number of rows that have to be allocated per2 */

  printf ( "%d %d extras in cols and rows relative to baseline %d\n" , extras_c , extras_r , per1 ) ; 

  num_c = ivector ( 1 , Cols ) ; 
  pikd_c = ivector ( 1 , Cols ) ; 
  cum = ivector ( 0 , Cols ) ; cum[0] = 0 ; 
  mylist =  ivector ( -1 , per2 ) ; 

  tot_per_col = ( per2 * Rows ) / Cols ;  
  relmax = per2 * ( tot_per_col + 1 ) ; 
  relatives = ivector ( 1 , relmax ) ; 
  
  printf ( "Start:" ) ; fflush(stdout) ;

  /* allocate the larger numbers of 1s (per_col2) to the cols with small j */

  for ( j = 1 ; j <= Cols ; j++ )  
    num_c[j] = ( j <= extras_c ) ? per_col2 : per_col1 ; 

  /* Run up the rows (m) from bottom upwards, starting from Rows - c0 */

  for ( i = num_rows_to_do , m = Rows - c0 ; i >= 1 ; i-- , m-- ) {

    if ( !(i%1000) ) { printf ( "%d;" , i ) ; fflush(stdout) ; }

    /* compute cumulatives */
    for ( j = 1 ; j <= Cols ; j ++ ) 
      cum[j] = cum[j-1] + num_c[j] ; 
    
    this_per = ( i > num_rows_to_do - extras_r ) ? per2 : per1 ;

    /* allocate the larger numbers of 1s (per2) to the rows with large i */

    do { /* generate a row */

      for ( j = 1 ; j <= Cols ; j ++ ) 
	pikd_c[j] = 0 ;
      for ( num_left = this_per ; 
	   num_left >= 1 ;
	   num_left -- ) {
	do {
	  un = ranf() ; 
	  j = read_int_from_cum ( cum , 1 , Cols , un ) ; 
	} while ( pikd_c[j] ) ; /* while j already picked */
	pikd_c[j] ++ ; 
	if ( num_c[j] > 0 ) {
	  mylist[num_left] = j ; 
	} else { fprintf ( stderr , "ERROR cum failure %d %d %d\n" ,
			  j , cum[j-1] , cum[j] ) ; 	  num_left ++ ;
	       pause_for_return(); }
      }
      if ( num_left == 0 ) { /* check hd then move on */
	if ( hd > 0 ) {
	  hdfail = 0 ; 
/* 
   for each j on the list
   go and find any other rows that share this j
   check if we have already met this row today. If so, then that's 
   a failure for the traditional setting of hd
*/
	  rr = 0 ; 
	  prevrr = 0 ; 
	  for ( j = 1 ; j <= this_per && ( hdfail == 0 ) ; j ++ ) {
	    nn = mylist[j] ;
	    for ( u = num_nlist[nn] ; u >= 1 && ( hdfail == 0 ) ; u-- ) {
	      rr ++ ;
	      latest = relatives[rr] = nlist[nn][u] ; 
	      for ( ss = prevrr ; ss >= 1 ; ss -- ) {
		if ( relatives[ss] == latest ) {
		  hdfail = 1 ; 
		  break ; 
		}
	      }
	    }
	    prevrr = rr ; 
	    if ( rr > relmax ) {
	      fprintf ( stderr, "\nwarning, relatives overran vector %d %d\n" , rr , relmax ) ;
	      status -- ; 
	    }
	  }
	}
	if ( hd > 0 && hdfail == 1 ) {
	  printf ( "h" ) ; fflush(stdout) ;
	  if ( i == 1 ) { /* tough luck, have to accept the failure */
	    fprintf ( stderr, "\nwarning, final hamming check failure at relative %d %d out of %d\n" , rr , ss , relmax ) ;
	    hdfail = 0 ; 
	  }
	}
	if ( hdfail == 0 ) {
	  /* final acts before moving on */
	  for ( j = 1 ; j <= this_per ; j ++ ) {
	    nn = mylist[j] ;
	    add_to_alist ( a , nn , m ) ; 
	    num_c[nn] -- ; 
	  }
	}
      } else {
	printf ( "?" ) ; /* fflush(stdout) ;  THIS should not happen any more */
      }
    } while ( num_left != 0 || hdfail != 0 ) ;
  }
  free_ivector ( mylist , -1 , per2 ) ;
  free_ivector ( num_c , 1 , Cols ) ;
  free_ivector ( pikd_c , 1 , Cols ) ;
  free_ivector ( cum , 0 , Cols ) ;
  free_ivector ( relatives , 1 , relmax ) ;
  printf ( "done\n" ) ; fflush(stdout) ;
  
  return status ; 
}

int 	  read_int_from_dcum ( double *dcum , int lo , int hi , double u ) {
  /* dcum[lo-1] is expected to exist and be 0
     and dcum[hi] is expected to be 1 ;
     u is expected to be in 0,1
     */
  int i , top = lo  , bot = hi ;
  int done = 0 ; 
  /* top and bot are the extreme eligible values for i */
  /* first guess */
  i = (int) ( (double) ( hi - lo + 1 ) * u ) + lo ;
  if ( i > hi ) { i = hi ; }
  if ( i < lo ) { i = lo ; }

  do {
    if ( ( top <= bot ) ||
	( ( u < dcum[i] ) && ( u > dcum[i-1] ) )
	) { /* we are done */ 
      done = 1 ; 
      break ; 
    } else if  ( u < dcum[i] ) {
      top = i-1 ;
    } else {
      bot = i+1 ;
    }
    i = ( top + bot ) / 2 ; 
  } while ( !done ) ; 


  return i ; 
} 

int 	  read_int_from_cum ( int *cum , int lo , int hi , double u ) {
  /* cum[lo-1] is expected to exist and be 0
     and cum[hi] is expected to be positive ;
     u is expected to be in 0,1
     */
  int i , top = hi  , bot = lo ;
  int done = 0 ; 
  /* top and bot are the extreme eligible values for i */
  int ONE = cum[hi] ; 
  int compare = (int) ( u * (double)(ONE) ) ; 

  /* first guess */
  i = (int) ( (double) ( hi - lo + 1 ) * u ) + lo ;
  if ( i > hi ) { i = hi ; }
  if ( i < lo ) { i = lo ; }

  do {
    if ( ( top <= bot ) ||
	( ( compare < cum[i] ) && ( compare >= cum[i-1] ) )
	) { /* we are done */ 
      done = 1 ; 
      break ; 
    } else if  ( compare < cum[i] ) {
      top = i-1 ;
    } else {
      bot = i+1 ;
    }
    i = ( top + bot ) / 2 ; 
  } while ( !done ) ; 


  return i ; 
} 

/* invert_cmatrix is in r.c */

void invert_left_hand_end_cmatrix ( unsigned char **A , unsigned char **B , int N ) {
/* takes a rectangular matrix (eg) of height N and length longer than N 
   and inverts the left hand bit, modifying the matrix to make 
   sure it is invertible. The inverse ends up in B.
   This can be used on square matrices too.
 memory should be allocated already. 
*/
  int i ;
  cm_inversion pr ;

  allocate_cm_inversion ( A , 1 , N , B , &pr ) ; 
  while ( ( i = invert_cmatrix ( &pr ) ) == 0 )  {
    printf( "m.%d." , pr.i ) ;
    i = modify_cmatrix_row ( &pr ) ;
  }
  free_cm_inversion ( &pr ) ; 
}

/* memory should be allocated already. 
   this creates an A B pair such that A is sparse */
void sparse_invertible_cmatrix (  unsigned char **A , unsigned char **B ,
				int per , int N ) {

  sparse_random_cmatrix2 ( A , per , N  ) ; 
  invert_left_hand_end_cmatrix ( A , B , N ) ;
}

int random_cvector ( unsigned char *v , double f , int lo , int hi ) {
  int j = hi , c = 0 ;

  for (  ; j >= lo ; j -- )  {
    v[j] = ( ranu() < f ) ? 1 : 0 ; 
    c += v[j] ; 
  }
  return c ; 
}

int fixed_wt_cvector ( unsigned char *v , int n , int lo , int hi ) {
  int N = hi - lo + 1 ; 
  int j = hi ;
  double dj = (double) N ; 
  int status = 0 ; 

  for (  ; j >= lo ; j -- , dj -- )  {
    if ( ranf() <=  (double) n / dj ) {
      n -- ;
      v[j] = 1 ;
    }
    else v[j] = 0 ; 
  }
  if ( n!= 0 ) {
    fprintf ( stderr , "warning: fixed_wt_cvector failed by %d / %d\n" ,
			n , N ) ; 
    status -- ; 
  }
  return status ; 
}

