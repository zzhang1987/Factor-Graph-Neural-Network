/*
   r2z.c
                                       (c) DJCM 95 11 06 

   - Reads r from BSC and multiplies by Cn

   ======================================================================

   This code is (c) David J.C. MacKay 1994, 1995. It is free software 
   as defined by the free software foundation. 

     suggested usage:

     r2z -rsuffix 1 -rfile r -zfile z -zsuffix 1 -Cnfile Cn

*/

#include "./ansi/r.h"
#include "./ansi/mynr.h"
#include "./ansi/cmatrix.h"

typedef struct {
#include "r2z_var_str.h"
} r2z_control ;

static void c_defaults ( r2z_control * ) ; 
static int    process_command ( int , char ** , r2z_control * ) ; 
static void   print_usage ( char ** , FILE * ,  r2z_control *  );

int   main ( int , char ** ) ;

/*
        MAIN
                     */
int main ( int argc, char *argv[] )
{
  r2z_control           c ;

  int done = 0 , message=0 ; 
  char junk[1000] ; 
  FILE *fpz , *fpr ;
  int b , bit ;
  alist_matrix Cn ; 
  unsigned char *z , *r ; 

  c_defaults ( &c ) ; 
  if ( process_command (argc, argv, &c ) < 0 ) exit (0) ;

  fprintf(stderr,"r2z running\n") ;
  fflush(stderr);

  /* allocate memory */
  if ( read_allocate_alist ( &Cn , c.Cnfile ) < 0 ) exit (0) ; 
  if ( ( Cn.N != c.N )  || (  Cn.M != c.N )    ) {
    fprintf ( stderr , "mismatch %d %d %d\n" , Cn.N , c.N , Cn.M ) ;
    exit ( 0 ) ; 
  }
  z = cvector ( 1 , c.N ) ;
  r = cvector ( 1 , c.N ) ;

  if ( !c.zsuffix ) { /* just open the file once */
    fpz = fopen ( c.zfile , "w" ) ;
    if ( !fpz ) {
      fprintf ( stderr , "FATAL: couldn't open file %s\n" , c.zfile ) ; 
      exit ( 0 ) ; 
    }  
  }
  if ( !c.rsuffix ) { /* just open the file once */
    fpr = fopen ( c.rfile , "r" ) ;
    if ( !fpr ) {
      fprintf ( stderr , "FATAL: couldn't open file %s\n" , c.rfile ) ; 
      exit ( 0 ) ; 
    }  
  }
  do {  /* loop 1 */
    message ++ ; 
    if ( c.rsuffix ) {
      sprintf ( junk , "%s.%04d" , c.rfile , message ) ;
      fpr = fopen ( junk , "r" ) ;
      if ( !fpr ) {
	if ( message == 1 ) { /* unexpected failure */
	  fprintf ( stderr , " couldn't open file %s\n" , junk ) ; 
	} else {
	  fprintf ( stderr , " r2z finished after %d blocks\n" , message - 1 ) ;
	}
	done = 1 ; break ; 
      }     
    }
    /* get next N bits of r vector from file */
    for ( b = 1 ; b <= c.N ; b ++ ) {
      if ( fscanf ( fpr , "%d" , &bit ) == EOF ) {
	fprintf ( stderr , "stream ended at bit %d block %d\n" , b , message ) ;
	done = 1 ; 
	break ; 
      }
      r[b] = (unsigned char ) bit ;
    }
    if ( ( done == 1 ) && ( b != 1 ) ) {
      fprintf ( stderr , "Warning: file ended with incomplete block \n" ) ;
    } else if ( ! done )  {
      if ( c.zsuffix ) {
	sprintf ( junk , "%s.%04d" , c.zfile , message ) ;
	fpz = fopen ( junk , "w" ) ;
	if ( !fpz ) {
	  fprintf ( stderr , "FATAL: couldn't open file %s\n" , junk ) ; 
	  exit ( 0 ) ; 
	}     
      }
      /* multiply */
      alist_times_cvector_mod2 ( &Cn , r , z ) ;

      for ( b = 1 ; b <= c.N ; b ++ ) {
	fprintf ( fpz , "%d\n" , z[b] ) ; 
      }

/*       fprintf( fpz , "\n" ) ; */
      if ( c.zsuffix ) fclose ( fpz ) ; 
    }
    if ( c.rsuffix ) fclose ( fpr ) ; 

  } while ( !done ) ;
  
  if ( !c.zsuffix )     fclose ( fpz ) ; 
  if ( !c.rsuffix )     fclose ( fpr ) ; 

  free_cvector ( r , 1 , c.N ) ; 
  free_cvector ( z , 1 , c.N ) ; 
  free_alist ( &Cn ) ; 
}

static void c_defaults ( r2z_control *c ) 
{
#include "r2z_var_def.c"
}

static int process_command ( int argc , char **argv , r2z_control *c ) {

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
  for (i = 1 ; i < argc; i++)    {
    cs = 1 ;
    if ( strcmp (argv[i], "-V") == 0 )        {
    }
#include "r2z_var_clr.c"
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

#define DNT fprintf( fp, "\n        ")
#define NLNE  fprintf( fp, "\n")

static void print_usage ( char **argv , FILE * fp , r2z_control *c
)
{
  fprintf( fp, "Usage: %s ",argv[0]);
  fprintf( fp, " [optional arguments]");

#include "r2z_var_usg.c"
  fprintf( fp, "\n");
  return ;
}

#undef DNT
#undef NLNE

/*
<!-- hhmts start -->
Last modified: Wed Dec  6 11:22:03 1995
<!-- hhmts end -->
*/
