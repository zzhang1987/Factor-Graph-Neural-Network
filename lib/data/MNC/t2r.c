/*
   t2r.c
                                       (c) DJCM 95 11 06 

   - Binary Symmetric channel
   - Reads in t (blocks) and writes r (blocks)

   ======================================================================

   This code is (c) David J.C. MacKay 1994, 1995. It is free software 
   as defined by the free software foundation. 

     suggested usage:

     t2r -tsuffix 1 -tfile t -rfile r -rsuffix 1 -fn 0.1

*/

#include "./ansi/r.h"
#include "./ansi/rand2.h"
#include "./ansi/mynr.h"
#include "./ansi/cmatrix.h"

typedef struct {
#include "t2r_var_str.h"
} t2r_control ;

static void c_defaults ( t2r_control * ) ; 
static int    process_command ( int , char ** , t2r_control * ) ; 
static void   print_usage ( char ** , FILE * ,  t2r_control *  );

void   main ( int , char ** ) ;

/*
        MAIN
                     */
void main ( int argc, char *argv[] )
{
  t2r_control           c ;

  int done = 0 , message=0 ; 
  char junk[1000] ; 
  FILE *fpr , *fpt ;
  int b , bit , noise ,  r ;
  int count = 0 , cycles = 0 ;

  c_defaults ( &c ) ; 
  if ( process_command (argc, argv, &c ) < 0 ) exit (0) ;
  ran_seed ( c.seed ) ; 

  fprintf(stderr,"t2r running\n") ;
  fflush(stderr);

  if ( !c.rsuffix ) { /* just open the file once */
    fpr = fopen ( c.rfile , "w" ) ;
    if ( !fpr ) {
      fprintf ( stderr , "FATAL: couldn't open file %s\n" , c.rfile ) ; 
      exit ( 0 ) ; 
    }  
  }
  if ( !c.tsuffix ) { /* just open the file once */
    fpt = fopen ( c.tfile , "r" ) ;
    if ( !fpt ) {
      fprintf ( stderr , "FATAL: couldn't open file %s\n" , c.tfile ) ; 
      exit ( 0 ) ; 
    }  
  }
  do {  /* loop 1 */
    message ++ ; 
    if ( c.tsuffix ) {
      sprintf ( junk , "%s.%04d" , c.tfile , message ) ;
      fpt = fopen ( junk , "r" ) ;
      if ( !fpt ) {
	if ( message == 1 ) { /* unexpected failure */
	  fprintf ( stderr , " couldn't open file %s\n" , junk ) ; 
	} else {
	  fprintf ( stderr , " t2r finished after %d blocks\n" , message - 1 ) ;
	}
	done = 1 ; break ; 
      }     
    }
    /* get next N bits of t vector from file */
    for ( b = 1 ; b <= c.N ; b ++ ) {
      if ( fscanf ( fpt , "%d " , &bit ) == EOF ) {
	fprintf ( stderr , "stream ended at bit %d block %d\n" , b , message ) ;
	done = 1 ; 
	break ; 
      }
      if ( c.rsuffix && (b==1) ) { /* time to open a new file */
	sprintf ( junk , "%s.%04d" , c.rfile , message ) ;
	fpr = fopen ( junk , "w" ) ;
	if ( !fpr ) {
	  fprintf ( stderr , "FATAL: couldn't open file %s\n" , junk ) ; 
	  exit ( 0 ) ; 
	}     
      }
      noise = ( ranu() > c.fn ) ? 0 : 1 ; 
      r = noise ^ bit ;
      count += noise ;
      cycles ++ ;
      fprintf ( fpr , "%d\n" , r ) ; 
    }
    if ( ( done == 1 ) && ( b != 1 ) ) {
      fprintf ( stderr , "Warning: file ended with incomplete block \n" ) ;
    } 
    if ( c.tsuffix ) fclose ( fpt ) ; 
    if ( c.rsuffix ) fclose ( fpr ) ; 

  } while ( !done ) ;
  
  if ( !c.tsuffix ) fclose ( fpt ) ; 
  if ( !c.rsuffix ) fclose ( fpr ) ; 
  fprintf ( stderr , "Sent %d bits and added %d errors\n" , cycles , count ) ;
}

static void c_defaults ( t2r_control *c ) 
{
#include "t2r_var_def.c"
}

static int process_command ( int argc , char **argv , t2r_control *c ) {

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
#include "t2r_var_clr.c"
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

static void print_usage ( char **argv , FILE * fp , t2r_control *c )
{
  fprintf( fp, "Usage: %s ",argv[0]);
  fprintf( fp, " [optional arguments]");

#include "t2r_var_usg.c"
  fprintf( fp, "\n");
  return ;
}

#undef DNT
#undef NLNE

/*
<!-- hhmts start -->
Last modified: Wed Dec  6 14:31:59 1995
<!-- hhmts end -->
*/
