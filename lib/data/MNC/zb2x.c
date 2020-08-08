/*
   zb2x.c
                                       (c) DJCM 95 11 06 
				                95 11 16 change to RCS


   - MacKay-Neal Error Correcting Code - 
   - Reads in z and bias vector and solves Ax=z using belief net decoder

   ======================================================================

   This code is (c) David J.C. MacKay 1994, 1995. It is free software 
   as defined by the free software foundation. 

   What this program does: 

     reads in a matrix A in alist format.
     do {
       read in z and / or b 
       give b,A,z to the bnd
     }
     while ( no more to read AND not many failures detected )

     suggested usage:

     zb2x -bfile b -zfixed 0 -bsuffix 1 -xsuffix 0 -xfile x -xso -Afile a

     this usage reads in a load of b vectors of the form b.0001 b.0002 etc.
     writes the decoded x vectors into a single file x, recording
     the source bits only, not the noise bits.
*/

#include "./ansi/r.h"
#include "./ansi/rand2.h"
#include "./ansi/mynr.h"
#include "./ansi/cmatrix.h"
#include "bnd/bnd.h"
/* this must come last: */
#include "./zb2x.h"

static void c_defaults ( zb2x_control * ) ; 
static int    process_command ( int , char ** , zb2x_all * ) ; 
static void   print_usage ( char ** , FILE * ,  zb2x_all *  );
static int make_sense ( zb2x_control * , zb2x_all * ) ; 
static int make_space ( zb2x_control * , zb2x_vectors * ) ;
static void set_up_priors ( zb2x_vectors * , zb2x_control * ) ; 

static void zb2x_free ( zb2x_all * ) ; 
static int check_alist_MN ( alist_matrix * , zb2x_vectors *vec ) ; 
static int hook_zb2x_vec_to_bnd ( bnd_param *p , zb2x_vectors *vec ) ;
static double h2 ( double ) ;

void   main ( int , char ** ) ;

/*
        MAIN
                     */
void main ( int argc, char *argv[] )
{
  bnd_control         bndc ; 
  bnd_param           bndp ;
  zb2x_control           c ;
  zb2x_vectors         vec  ; 
  zb2x_all             all  ;
  alist_matrix         a ;

  int stillgoing = 1 , message = 0 , failures = 0 ; 
  char junk[1000] ; 
  FILE *fpb, *fpz , *fpx ;

  all.bndp = &bndp ; 
  all.bndc = &bndc ; 
  all.vec  = &vec ; 
  all.a    = &a ; 
  all.c    = &c ; 
  bndp.a   = &a ; 

  bnd_defaults  ( &bndp , & bndc ) ; 
  c_defaults ( &c ) ; 

  if ( process_command (argc, argv, &all ) < 0 ) exit (0) ;
  if ( make_sense ( &c , &all ) < 0 ) exit (0) ;

  fprintf(stderr,"zb2x running\n") ;
  fflush(stderr);

  if ( make_space ( &c , &vec ) < 0 ) exit (0) ;

  if ( read_allocate_alist ( &a , c.Afile ) < 0 ) exit (0) ; 
  if ( check_alist_MN ( &a , &vec ) < 0 ) exit (0) ; 

  hook_zb2x_vec_to_bnd ( &bndp , &vec ) ;
  bnd_allocate ( &bndp , &bndc ) ; 

  /*  set up fixed b values and fixed z values */

  set_up_priors ( &vec , &c ) ; 

  do {  /* loop 1 */
    message ++ ; 
    /* get relevant parts of bias vector from file */
    if ( c.bfromfile ) {
      if ( message == 1 || c.bsuffix ) { /* need to open a file */
	if (  c.bsuffix ) {
	  sprintf ( junk , "%s.%04d" , c.bfile , message ) ;
	} else {
	  sprintf ( junk , "%s" , c.bfile ) ;
	}
	fpb = fopen ( junk , "r" ) ;
	if ( !fpb ) {
	  if ( message == 1 ) { /* unexpected failure */
	    fprintf ( stderr , " couldn't open file %s\n" , junk ) ; 
	  } 
	  stillgoing = 0 ; break ; /* break out to end of loop 1 */
	} 
      }
      /* read in relevant bits of b */
      if ( fread_dvector ( vec.bias , vec.b_readfrom , vec.b_readto , fpb ) < 0 ) {
	if ( message == 1 ) { /* unexpected failure */
	  fprintf ( stderr , " couldn't read bias vector from %s\n" , junk ) ; 
	} 
	stillgoing = 0 ; break ; /* break out to end of loop 1 */
      }  
      if ( c.bsuffix ) fclose ( fpb ) ; 
    }
    /* get z vector from file */
    if ( c.zfromfile ) {
      if ( message == 1 || c.zsuffix ) { /* need to open a file */
	if ( c.zsuffix ) {
	  sprintf ( junk , "%s.%04d" , c.zfile , message ) ;
	} else {
	  sprintf ( junk , "%s" , c.zfile ) ;
	}
	fpz = fopen ( junk , "r" ) ;
	if ( !fpz ) {
	  if ( message == 1 ) { /* unexpected failure */
	    fprintf ( stderr , " couldn't open file %s\n" , junk ) ; 
	  }
	  stillgoing = 0 ; break ; /* break out to end of loop 1 */
	} 
      }
      /* read in relevant bits of z */
      if ( fread_cvector ( vec.z , 1 , c.N , fpz ) < 0 ) {
	if ( message == 1 ) { /* unexpected failure */
	  fprintf ( stderr , " couldn't read z vector from %s\n" , junk ) ; 
	} 
	stillgoing = 0 ; break ; /* break out to end of loop 1 */
      }  
      if ( c.zsuffix ) fclose ( fpz ) ; 
    }
    
    if ( bndecode ( &bndp , &bndc ) > 0 ) {
      fprintf ( stderr , "Decoding failure block %d\n" , message ) ;
      failures ++ ; 
    }

    if ( c.xtofile ) {
      if ( message == 1 || c.xsuffix ) { /* need to open a file */
	if ( c.xsuffix ) {
	  sprintf ( junk , "%s.%04d" , c.xfile , message ) ;
	} else {
	  sprintf ( junk , "%s" , c.xfile ) ;
	}
	fpx = fopen ( junk , "w" ) ;

	if ( !fpx ) {
	  fprintf ( stderr , "FATAL: couldn't open file %s\n" , junk ) ; 
	  exit ( 0 ) ; 
	} 
      }
      write_cvector ( fpx , vec.x , vec.xfrom , vec.xto ) ; /* writes 1 0 0 1 0\n */
      if ( c.xsuffix ) fclose ( fpx ) ; 
    }
    

  } while ( stillgoing ) ;

  fprintf ( stderr , " finished after %d blocks; %d block errors\n" , message - 1 , failures ) ;

  if ( c.xtofile && (!(c.xsuffix)) ) fclose ( fpx ) ; 
  if ( c.zfromfile && (!(c.zsuffix)) ) fclose ( fpz ) ; 
  if ( c.bfromfile && (!(c.bsuffix)) ) fclose ( fpb ) ; 
  
  bnd_free ( &bndp , &bndc ) ;
  zb2x_free ( &all ) ; 

}

static int make_sense ( zb2x_control *c , zb2x_all *all ) 
{ /*  correct silly control parameters */
/* 
   first the data creation 
                             */
  int status = 0 ; 

  /* there is potential confusion 
     here:
     the "n,k" of the error correcting code
     c.f. the matrix is A is M * N
     where M = n
     and   N = "k + n" 
      */
  all->vec->M = c->N ; 
  all->vec->N = c->N + c->K ; 
  all->vec->NS = c->K ; 
  all->vec->NN = c->N ; 
  all->vec->nsfrom = 1 ; 
  all->vec->nsto   = c->K ; 
  all->vec->nnfrom = c->K + 1 ; 
  all->vec->nnto   = c->K + c->N  ; 
  if ( c->xsourceonly ) {
    all->vec->xfrom = all->vec->nsfrom ;
    all->vec->xto   = all->vec->nsto ;
  } else { /* write whole x vector */
    all->vec->xfrom = all->vec->nsfrom ;
    all->vec->xto   = all->vec->nnto ;
  }

  if ( c->bnfromfile || c->bsfromfile ) { 
    c->bfromfile = 1 ; 
    if ( c->bsfromfile && c->bnfromfile ) {
      all->vec->b_readfrom = 1 ; 
      all->vec->b_readto = all->vec->N ; 
    } else if ( c->bsfromfile ) { 
      all->vec->b_readfrom = all->vec->nsfrom ; 
      all->vec->b_readto = all->vec->nsto ; 
    } else if ( c->bnfromfile ) { 
      all->vec->b_readfrom = all->vec->nnfrom ; 
      all->vec->b_readto = all->vec->nnto ; 
    } else { fprintf ( stderr , "Huh?\n" ) ; }
    fprintf ( stderr , "bias[%d:%d] will be read from %s\n" , 
	     all->vec->b_readfrom , all->vec->b_readto , c->bfile ) ; 
  } else {
    c->bfromfile = 0 ; 
  }

  if ( c->xtofile ) {
    if ( c->xsourceonly ) {
      all->vec->xfrom = all->vec->nsfrom ;
      all->vec->xto   = all->vec->nsto ; 
    } else {
      all->vec->xfrom = 1 ; 
      all->vec->xto   = all->vec->N ; 
    }
  }

  return status ; 
}

static int make_space ( zb2x_control *c , zb2x_vectors *vec ) {
/* see also zb2x_free */
  int status = 0 ; 

  vec->x = cvector ( 1 , vec->N ) ; 
  vec->z = cvector ( 1 , vec->M ) ;

  vec->bias = dvector ( 1 , vec->N ) ;

  return status ; 
}

static int hook_zb2x_vec_to_bnd ( bnd_param *p , zb2x_vectors *vec ) {
  int status = 0 ; 
  p->M = vec->M ;
  p->N = vec->N ;
  p->z = vec->z ; 
  p->xo = vec->x ; /* where bnd gets to put its output */
  p->bias = vec->bias ;

  return status ; 
}

int check_alist_MN ( alist_matrix *a , zb2x_vectors *v ) {
  int status = 0 ; 
  if (     v->N != a->N 
      ||   v->M != a->M      ) {
    fprintf ( stderr , "eek %d %d %d %d\n" , v->N , a->N , v->M , a->M ) ; 
    status -- ; 
  }
  return status ; 
}

static void set_up_priors ( zb2x_vectors *v , zb2x_control *c ) {
  int n ;
  double *b = v->bias ; 
  unsigned char *z = v->z ; 

  if ( c->bsfromfile == 0 ) {
    for ( n = v->nsfrom ; n <= v->nsto ; n++ ) {
      b[n] = c->bfs ;
    }
  }
  if ( c->bnfromfile == 0 ) {
    for ( n = v->nnfrom ; n <= v->nnto ; n++ ) {
      b[n] = c->bfn ;
    }
  }
  if ( c->zfromfile == 0 ) {
    for ( n = 1 ; n <= c->N ; n ++ ) { 
      z[n] = c->zfixed ; 
    }
  }
}
  
static void c_defaults ( zb2x_control *c ) 
{
#include "zb2x_var_def.c"
}

static int process_command ( int argc , char **argv , zb2x_all *all ) {
  bnd_control *bndc = all->bndc ; 
  zb2x_control *c = all->c ; 

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
#include "zb2x_var_clr.c"
#include "bnd/bnd_var_clr.c"
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
  if ( p_usage ) print_usage ( argv , stderr , all ) ;
  return ( status ) ;
}
#undef ERROR1
#undef ERROR2
#undef ERRORREG

#define DNT fprintf( fp, "\n        ")
#define NLNE  fprintf( fp, "\n")

static void print_usage ( char **argv , FILE * fp ,
			 zb2x_all *all 
)
{
  zb2x_control *c = all->c ; 
  bnd_control *bndc = all->bndc ; 
  
  fprintf( fp, "Usage: %s ",argv[0]);
  fprintf( fp, " [optional arguments]");

#include "zb2x_var_usg.c"
  NLNE; fprintf( fp, " Belief Net decoder:     <defaults>");
#include "bnd/bnd_var_usg.c"
  fprintf( fp, "\n");
  return ;
  return ;
}

#undef DNT
#undef NLNE

static void zb2x_free ( zb2x_all *all ) {
  zb2x_vectors *vec = all->vec ; 

  free_cvector ( vec->x , 1 , vec->N ) ; 
  free_cvector ( vec->z , 1 , vec->M ) ; 

  free_dvector ( vec->bias , 1 , vec->N ) ; 

}

double h2 ( double x ) {
  double tmp ; 
  if ( x <= 0.0 || x>= 1.0 ) tmp = 0.0 ;
  tmp = x * log ( x ) + ( 1.0 - x ) * log ( 1.0 - x ) ;
  
  return - tmp / log ( 2.0) ; 
}

