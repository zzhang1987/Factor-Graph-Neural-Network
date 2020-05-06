/*
   s2t.c
                                       (c) DJCM 95 11 06 

   - MacKay-Neal Error Correcting Code - 
   - Reads in s (stream) and writes t=Gs (blocks)

   ======================================================================

   This code is (c) David J.C. MacKay 1994, 1995. It is free software 
   as defined by the free software foundation. 

     suggested usage:

     s2t -sfile s -tsuffix 1 -Gfile G -tfile t

     this usage reads in s and writes t vectors to t.0001 t.0002 etc.

     The option -smn 1
     causes the s vector to get prepended to the tfile.

     -n and -k need to be supplied. These should be the n and k 
     of the NMN code, i.e. the dimensions of the G matrix. 

*/

#include "./ansi/r.h"
#include "./ansi/rand2.h"
#include "./ansi/mynr.h"
#include "./ansi/cmatrix.h"
#include "./radford/mod2mat.h"

typedef struct {
#include "s2t_var_str.h"
} s2t_control ;

static void c_defaults ( s2t_control * ) ; 
static int    process_command ( int , char ** , s2t_control * ) ; 
static void   print_usage ( char ** , FILE * ,  s2t_control *  );

int   main ( int , char ** ) ;

/*
        MAIN
                     */
int main ( int argc, char *argv[] )
{
  s2t_control           c ;

  int code , done ; 
  int message ; 
  char junk[1000] ; 
  FILE *fps , *fpt , *fp ;
  mod2mat *G, *s, *t ;
  int b , bit ;

  c_defaults ( &c ) ; 
  if ( process_command (argc, argv, &c ) < 0 ) exit (0) ;
  /* fprintf(stderr,"s2t running\n") ; */
  fflush(stderr);

/* G is read in in radford format . Expect it to be 
  ( c.N , c.K ) ; */
  s = mod2mat_allocate ( c.K , 1 ) ;
  t = mod2mat_allocate ( c.N , 1 ) ;
  fp = fopen ( c.Gfile , "rb" ) ; 
  if ( !fp ) {
    fprintf ( stderr , " couldn't open file %s\n" , c.Gfile ) ; 
    exit(0);
  }	
  G = mod2mat_read ( fp , &code ) ; 
  fps = fopen ( c.sfile , "r" ) ;
  if ( !fps ) {
    fprintf ( stderr , " couldn't open file %s\n" , c.sfile ) ; 
    exit(0);
  }	
  message = 0 ; done = 0 ; 
  do {  /* loop 1 */
    message ++ ; 
    /* get next K bits of s vector from file */
    for ( b = 1 ; b <= c.K ; b ++ ) {
      if ( fscanf ( fps , "%d" , &bit ) == EOF ) {
	done = 1  ;
	break ; 
      }
      mod2mat_set ( s , b-1 , 0 , bit ) ;
    }
    if ( done == 1 ) {
      if ( b == 1 ) {
	fprintf ( stderr , " exact number of blocks %d\n" , message - 1 ) ;
	break ; /* break out of loop 1 */

      }
      bit = 0 ;
      fprintf ( stderr , " filling up %dth block with zeroes\n" , message ) ;
      for ( ; b <= c.K ; b ++ ) {
	mod2mat_set ( s , b-1 , 0 , bit ) ;
      }
    }

	
    mod2mat_multiply ( G , s , t ) ; 
    if ( message == 1 || c.tsuffix ) { /* time to open output file */
      if ( c.tsuffix ) {
	sprintf ( junk , "%s.%04d" , c.tfile , message ) ;
      } else {
	sprintf ( junk , "%s" , c.tfile ) ;
      }
      fpt = fopen ( junk , "w" ) ;
      if ( !fpt ) {
	fprintf ( stderr , "FATAL: couldn't open file %s\n" , junk ) ; 
	exit ( 0 ) ; 
      } 
    }
    if ( c.smn ) { /* systematic code - copy out s to t-stream first */
      for ( b = 0 ; b < c.K ; b ++ ) {
	fprintf ( fpt , "%d\n" , mod2mat_get ( s , b , 0 ) ) ;
      }
    }
    for ( b = 0 ; b < c.N ; b ++ ) {
      fprintf ( fpt , "%d\n" , mod2mat_get ( t , b , 0 ) ) ;
    }
    if ( c.tsuffix ) fclose ( fpt ) ; 

  } while ( !done ) ;
  
  if ( !(c.tsuffix) ) fclose ( fpt ) ; 
  fclose ( fps ) ; 

}

static void c_defaults ( s2t_control *c ) 
{
#include "s2t_var_def.c"
}

static int process_command ( int argc , char **argv , s2t_control *c ) {

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
#include "s2t_var_clr.c"
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

static void print_usage ( char **argv , FILE * fp , s2t_control *c )
{
  fprintf( fp, "Usage: %s ",argv[0]);
  fprintf( fp, " [optional arguments]");

#include "s2t_var_usg.c"
  fprintf( fp, "\n");
  return ;
}

#undef DNT
#undef NLNE

/*
<!-- hhmts start -->
Last modified: Wed Dec  6 14:29:52 1995
<!-- hhmts end -->
*/
