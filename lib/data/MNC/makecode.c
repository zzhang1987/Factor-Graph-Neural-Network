/*
   mkcode.c
                                       (c) DJCM 95 11 06 

   - MacKay-Neal Error Correcting Code - 
   - Reads in Alist matrix. Forcibly inverts Cn. 
     Writes G = Cn^-1 Cs. (explicit format)
     Writes A. (alist)
     Writes Cn. (alist)

   ======================================================================

   This code is (c) David J.C. MacKay 1994, 1995. It is free software 
   as defined by the free software foundation. 

*/

#include "./ansi/r.h"
#include "./ansi/mynr.h"
#include "./ansi/cmatrix.h"
#include "./radford/mod2mat.h"

typedef struct {
#include "mkcode_var_str.h"
} mkcode_control ;

static void c_defaults ( mkcode_control * ) ; 
static int    process_command ( int , char ** , mkcode_control * ) ; 
static void   print_usage ( char ** , FILE * ,  mkcode_control *  );

void   main ( int , char ** ) ;

void main ( int argc, char *argv[] )
{
  mkcode_control           c ;

  int forces , *a_row , *a_col ; 
  mod2mat *G, *Cn, *CnI, *Cs ;
  int K , N , Kprime ;
  int n , u , m , cn , tr , tc ; 
  alist_matrix A , Aout , ACn ; 
  FILE *fp ;

  c_defaults ( &c ) ; 
  if ( process_command (argc, argv, &c ) < 0 ) exit (0) ;

  fprintf(stderr,"mkcode running\n") ;
  fflush(stderr);

  /* READ IN A matrix */

  if ( read_allocate_alist ( &A , c.Ain ) < 0 ) exit (0) ; 

  N = A.M ; Kprime = A.N ; K = Kprime - N ; 
  fprintf ( stderr , "Read in (%d*%d) A; creating an (%d,%d) code\n" ,
	   A.M , A.N , N , K ) ; 

  /* Copy A into explicit matrices [Cs Cn] and into
     Aout and ACn (alist format)
     */

  Cn  = mod2mat_allocate ( N , N ) ;
  Cs  = mod2mat_allocate ( N , K ) ;
  if ( Cn == 0 || Cs == 0 ) {
    fprintf ( stderr , "memory failure 1\n" ) ;
    exit ( 0 ) ;
  }

  mod2mat_clear ( Cn ) ;
  mod2mat_clear ( Cs ) ;

  tr = A.biggest_num_m ;
  tc = A.biggest_num_n ;
  if ( c.force ) {
    tr += c.force ;
    tc += c.force ; /* convention: N then M - col then row  :( */
    initialize_alist ( &Aout , Kprime , N , tc , tr ) ;
  }
  initialize_alist ( &ACn , N , N , tc , tr ) ;

  for ( n = 1 ; n <= K ; n ++ ) {
    for ( u = 1 ; u <= A.num_nlist[n] ; u ++ ) {
      m = A.nlist[n][u] ; 
      mod2mat_set ( Cs , m-1 , n-1 , 1 ) ; 
      if ( c.force ) add_to_alist ( &Aout , n , m ) ;
    }
  }
  for ( cn = 1 ; n <= Kprime ; n ++ , cn ++ ) {
    for ( u = 1 ; u <= A.num_nlist[n] ; u ++ ) {
      m = A.nlist[n][u] ; 
      mod2mat_set ( Cn , m-1 , cn-1 , 1 ) ; 
      if ( c.force ) add_to_alist ( &Aout , n , m ) ;
      add_to_alist ( &ACn , cn  , m) ;
    }
  }
  free_alist ( &A ) ; 

  /* Done with original A */

  /* Invert Cn */

  CnI = mod2mat_allocate ( N , N ) ;
  if ( CnI == 0 ) {
    fprintf ( stderr , "memory failure 2\n" ) ;
    exit ( 0 ) ;
  }

  if ( ! c.force ) {
    if ( mod2mat_invert ( Cn , CnI ) == 0 ) {
      fprintf ( stderr , "inversion failed\n" ) ;
      exit ( 0 ) ; 
    }
  } else {
    a_row = ivector ( 0 , N ) ;
    a_col = ivector ( 0 , N ) ;
    forces = mod2mat_forcibly_invert  ( Cn , CnI , a_row , a_col ) ; 
    fprintf ( stderr , "inversion achieved with %d forces\n" , forces ) ; 
    if ( forces > 0 ) {
      /* go through list, updating A matrix */
      for ( ; forces >= 1 ; forces -- ) {
	add_to_alist ( &Aout , 
		      a_col[forces - 1] + 1 + K , a_row[forces - 1] + 1 ) ;
	add_to_alist ( &ACn  , 
		      a_col[forces - 1] + 1 , a_row[forces - 1] + 1) ;
	/* there is a risk of over-run here */
	/* I assume here that there is not already a 
	   1 in the location */
      }
    }
    free_ivector ( a_row , 0 , N ) ;
    free_ivector ( a_col , 0 , N ) ;
  }
  mod2mat_free ( Cn ) ; 

  /* file writing */
  /* write Aout */
  if ( c.force ) {
    finish_off_alist ( &Aout ) ;
    fp = fopen ( c.Aout , "w" ) ; 
    if ( !fp ) fprintf ( stderr , "can't open %s\n" , c.Aout ) ;
    else {
      write_alist ( fp , &Aout ) ; 
      fclose ( fp ) ; 
    }
    free_alist ( &Aout ) ; 
  } else {
    fprintf ( stderr , "No changes requested to A so - \n" ) ;
    fprintf ( stderr , "cp -f %s %s\n" , c.Ain , c.Aout ) ;
  }

  /* write Cn */
  finish_off_alist ( &ACn ) ;
  fp = fopen ( c.Cn , "w" ) ; 
  if ( !fp ) fprintf ( stderr , "can't open %s\n" , c.Cn ) ;
  else {
    write_alist ( fp , &ACn ) ; 
    fclose ( fp ) ; 
  }
  free_alist ( &ACn ) ; 

  /* write G */

  G   = mod2mat_allocate ( N , K ) ;
  if ( G == 0 ) {
    fprintf ( stderr , "memory failure 3\n" ) ;
    exit ( 0 ) ;
  }

  mod2mat_multiply ( CnI , Cs , G ) ; 
  fp = fopen ( c.G , "w" ) ;
  if ( !fp ) fprintf ( stderr , "can't open %s\n" , c.G ) ;
  else {
    mod2mat_write ( fp , G ) ; 
    fclose ( fp ) ; 
  }

  mod2mat_free ( CnI ) ; 
  mod2mat_free ( Cs  ) ; 
  mod2mat_free ( G   ) ; 
}

static void c_defaults ( mkcode_control *c ) 
{
#include "mkcode_var_def.c"
}

static int process_command ( int argc , char **argv , mkcode_control *c ) {

  int p_usage = 0 ;
  int status = 0 ;
  int cs , i ;

  if ( argc < 1 )     {
    p_usage = 1 ; 
    status -- ;
  }

#define ERROR1 fprintf ( stderr , "arg to `%s' missing\n" , argv[i] ) ; \
               status --
#define ERROR2 fprintf ( stderr , "args to `%s' missing\n" , argv[i] ) ; \
               status --
#define ERRORREG fprintf ( stderr , "regtype must be defined before `%s'\n" , argv[i] ) ; \
               status --
  for (i = 1 ; i < argc; i++)    {
    cs = 1 ;
    if ( strcmp (argv[i], "-V") == 0 )        {
    }
#include "mkcode_var_clr.c"
    else {
      fprintf ( stderr , "arg `%s' not recognised\n" , argv[i] ) ; 
      p_usage = 1 ;
      status -- ;
    }
    if ( cs == 0 ) {
      fprintf ( stderr , "arg at or before `%s' has incorrect format\n" , 
	       argv[i] ) ;
      p_usage = 1 ;
      status -- ;
    }
  }
  if ( p_usage ) print_usage ( argv , stderr , c ) ;
  return ( status ) ;
}
#undef ERROR1
#undef ERROR2
#undef ERRORREG

#define DNT fprintf( fp, "\n        ")
#define NLNE  fprintf( fp, "\n")

static void print_usage ( char **argv , FILE * fp ,
			 mkcode_control *c
)
{
  fprintf( fp, "Usage: %s ",argv[0]);
  fprintf( fp, " [optional arguments]");

#include "mkcode_var_usg.c"
  fprintf( fp, "\n");
  return ;
  return ;
}

#undef DNT
#undef NLNE

/*
<!-- hhmts start -->
Last modified: Wed Nov  8 14:54:18 1995
<!-- hhmts end -->
*/
