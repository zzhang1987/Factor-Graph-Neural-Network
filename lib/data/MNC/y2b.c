/*
   y2b.c
                                       (c) DJCM 95 11 06 

   - Reads y from Gaussian channel and reports the normalized likelihood 

   ======================================================================

   This code is (c) David J.C. MacKay 1994, 1995. It is free software 
   as defined by the free software foundation. 

     suggested usage:

     y2b -bsuffix 1 -bfile b -yfile y -ysuffix 1 -gcx 1.0

*/

#include "./ansi/r.h"
#include "./ansi/rand2.h"
#include "./ansi/mynr.h"
#include "./ansi/cmatrix.h"

typedef struct {
#include "y2b_var_str.h"
} y2b_control ;

static void c_defaults ( y2b_control * ) ; 
static int    process_command ( int , char ** , y2b_control * ) ; 
static void   print_usage ( char ** , FILE * ,  y2b_control *  );

void   main ( int , char ** ) ;

/*
        MAIN
                     */
void main ( int argc, char *argv[] )
{
  y2b_control           c ;

  int done = 0 , message=0 ; 
  char junk[1000] ; 
  FILE *fpb , *fpy ;
  int b ;
  double y , bias ;

  c_defaults ( &c ) ; 
  if ( process_command (argc, argv, &c ) < 0 ) exit (0) ;

  fprintf(stderr,"y2b running\n") ;
  fflush(stderr);

  if ( !c.bsuffix ) { /* just open the file once */
    fpb = fopen ( c.bfile , "w" ) ;
    if ( !fpb ) {
      fprintf ( stderr , "FATAL: couldn't open file %s\n" , c.bfile ) ; 
      exit ( 0 ) ; 
    }  
  }
  if ( !c.ysuffix ) { /* just open the file once */
    fpy = fopen ( c.yfile , "r" ) ;
    if ( !fpy ) {
      fprintf ( stderr , "FATAL: couldn't open file %s\n" , c.yfile ) ; 
      exit ( 0 ) ; 
    }  
  }
  do {  /* loop 1 */
    message ++ ; 
    if ( c.ysuffix ) {
      sprintf ( junk , "%s.%04d" , c.yfile , message ) ;
      fpy = fopen ( junk , "r" ) ;
      if ( !fpy ) {
	if ( message == 1 ) { /* unexpected failure */
	  fprintf ( stderr , " couldn't open file %s\n" , junk ) ; 
	} else {
	  fprintf ( stderr , " y2b finished after %d blocks\n" , message - 1 ) ;
	}
	done = 1 ; break ; 
      }     
    }
    /* get next N bits of t vector from file */
    for ( b = 1 ; b <= c.N ; b ++ ) {
      if ( fscanf ( fpy , "%lf " , &y ) == EOF ) {
	fprintf ( stderr , "stream ended at bit %d block %d\n" , b , message ) ;
	done = 1 ; 
	break ; 
      }
      if ( c.bsuffix && ( b == 1 ) ) { /* open new file */
	sprintf ( junk , "%s.%04d" , c.bfile , message ) ;
	fpb = fopen ( junk , "w" ) ;
	if ( !fpb ) {
	  fprintf ( stderr , "FATAL: couldn't open file %s\n" , junk ) ; 
	  exit ( 0 ) ; 
	}     
      }
      bias = 1 / ( 1.0 + exp ( - 2.0 * c.gcx * y ) ) ;
      fprintf ( fpb , "%3.6g\n" , bias ) ; 
    }

    if ( ( done == 1 ) && ( b != 1 ) ) {
      fprintf ( stderr , "Warning: file ended with incomplete block \n" ) ;
    } 

    if ( c.ysuffix ) fclose ( fpy ) ; 
    if ( c.bsuffix ) fclose ( fpb ) ; 

  } while ( !done ) ;
  
  if ( !c.bsuffix )     fclose ( fpb ) ; 
  if ( !c.ysuffix )     fclose ( fpy ) ; 

}

static void c_defaults ( y2b_control *c ) 
{
#include "y2b_var_def.c"
}

static int process_command ( int argc , char **argv , y2b_control *c ) {

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
#include "y2b_var_clr.c"
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

static void print_usage ( char **argv , FILE * fp , y2b_control *c )

{
  fprintf( fp, "Usage: %s ",argv[0]);
  fprintf( fp, " [optional arguments]");

#include "y2b_var_usg.c"
  fprintf( fp, "\n");
  return ;
}

#undef DNT
#undef NLNE

/*
<!-- hhmts start -->
Last modified: Wed Dec  6 12:34:09 1995
<!-- hhmts end -->
*/
