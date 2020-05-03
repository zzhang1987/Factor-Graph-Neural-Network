/*
   r2b.c
                                       (c) DJCM 95 11 06 

   - Reads r from BSC and (setting z to zero) works out the 
   appropriate b

   ======================================================================

   This code is (c) David J.C. MacKay 1994, 1995. It is free software 
   as defined by the free software foundation. 

     suggested usage:

     r2b -rsuffix 1 -rfile r -bfile b -bsuffix 1 -fn 0.05

*/

#include "./ansi/r.h"
#include "./ansi/mynr.h"
#include "./ansi/cmatrix.h"

typedef struct {
#include "r2b_var_str.h"
} r2b_control ;

static void c_defaults ( r2b_control * ) ; 
static int    process_command ( int , char ** , r2b_control * ) ; 
static void   print_usage ( char ** , FILE * ,  r2b_control *  );

void   main ( int , char ** ) ;

/*
        MAIN
                     */
void main ( int argc, char *argv[] )
{
  r2b_control           c ;

  int done = 0 , message=0 ; 
  char junk[1000] ; 
  FILE *fpb , *fpr ;
  int b , bit ;
  unsigned char *r ;
  double bias ; 

  c_defaults ( &c ) ; 
  if ( process_command (argc, argv, &c ) < 0 ) exit (0) ;

  fprintf(stderr,"r2b running\n") ;
  fflush(stderr);

  r = cvector ( 1 , c.N ) ;

  if ( !c.bsuffix ) { /* just open the file once */
    fpb = fopen ( c.bfile , "w" ) ;
    if ( !fpb ) {
      fprintf ( stderr , "FATAL: couldn't open file %s\n" , c.bfile ) ; 
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
	  fprintf ( stderr , " r2b finished after %d blocks\n" , message - 1 ) ;
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
      if ( c.bsuffix ) {
	sprintf ( junk , "%s.%04d" , c.bfile , message ) ;
	fpb = fopen ( junk , "w" ) ;
	if ( !fpb ) {
	  fprintf ( stderr , "FATAL: couldn't open file %s\n" , junk ) ; 
	  exit ( 0 ) ; 
	}     
      }

      for ( b = 1 ; b <= c.N ; b ++ ) {
	bias = ( r[b] ) ? 1.0 - c.fn : c.fn ;
	fprintf ( fpb , "%g\n" , bias ) ; 
      }

/*       fprintf( fpb , "\n" ) ; */
      if ( c.bsuffix ) fclose ( fpb ) ; 
    }
    if ( c.rsuffix ) fclose ( fpr ) ; 

  } while ( !done ) ;
  
  if ( !c.bsuffix )     fclose ( fpb ) ; 
  if ( !c.rsuffix )     fclose ( fpr ) ; 

  free_cvector ( r , 1 , c.N ) ; 
}

static void c_defaults ( r2b_control *c ) 
{
#include "r2b_var_def.c"
}

static int process_command ( int argc , char **argv , r2b_control *c ) {

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
#include "r2b_var_clr.c"
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

static void print_usage ( char **argv , FILE * fp , r2b_control *c
)
{
  fprintf( fp, "Usage: %s ",argv[0]);
  fprintf( fp, " [optional arguments]");

#include "r2b_var_usg.c"
  fprintf( fp, "\n");
  return ;
}

#undef DNT
#undef NLNE

/*
<!-- hhmts start -->
Last modified: Wed Dec  6 11:46:40 1995
<!-- hhmts end -->
*/
