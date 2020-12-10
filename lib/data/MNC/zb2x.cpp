#include "ansi/cmatrix.h"
#include "ansi/nrutil.h"
#include "zb2x.h"
#include <iostream>

int check_alist_MN ( alist_matrix *a , zb2x_vectors *v ) {
    int status = 0 ; 
    if (     v->N != a->N 
             ||   v->M != a->M      ) {
        // fprintf ( stderr , "eek %d %d %d %d\n" , v->N , a->N , v->M , a->M ) ; 
        status -- ; 
    }
    return status ; 
}

int make_space ( zb2x_control *c , zb2x_vectors *vec ) {
/* see also zb2x_free */
    int status = 0 ; 

    vec->x = cvector ( 1 , vec->N ) ; 
    vec->z = cvector ( 1 , vec->M ) ;

    vec->bias = dvector ( 1 , vec->N ) ;
    return status ; 
}

int make_sense ( zb2x_control *c , zb2x_all *all ) 
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

        //fprintf ( stderr , "bias[%d:%d] will be read from %s\n" , 
        //          all->vec->b_readfrom , all->vec->b_readto , c->bfile ) ; 
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


int hook_zb2x_vec_to_bnd ( bnd_param *p , zb2x_vectors *vec ) {
    int status = 0 ; 
    p->M = vec->M ;
    p->N = vec->N ;
    p->z = vec->z ; 
    p->xo = vec->x ; /* where bnd gets to put its output */
    p->bias = vec->bias ;

    return status ; 
}

void set_up_priors ( zb2x_vectors *v , zb2x_control *c ) {
    int n;
    double *b = v->bias ; 
    unsigned char *z = v->z ; 

    if ( c->bsfromfile == 0 ) {
        for ( n = v->nsfrom ; n <= v->nsto ; n++ ) {
            b[n] = c->bfs ;
        }
    }
  
    if ( c->zfromfile == 0 ) {
        for ( n = 1 ; n <= c->N ; n ++ ) {
            z[n] = c->zfixed ; 
        }
    }
}
void c_defaults ( zb2x_control *c )
{
    c->verbose = 0 ;            /* -v verbose         */
    c->bnfromfile = 1 ;         /* - -                */
    c->bsfromfile = 1 ;         /* - -                */
    c->bfromfile = 1 ;          /* - -                */
    c->bfn = 0.1 ;              /* -bfn fn            */
    c->bfs = 0.5 ;              /* -bfs fs            */
    c->zfromfile = 0 ;          /* - -                */
    c->zfixed = 0 ;             /* -zfixed z          */
    c->bsuffix = 0 ;            /* -bsuffix bs        */
    c->zsuffix = 0 ;            /* -zsuffix zs        */
    c->xsuffix = 0 ;            /* -xsuffix xs        */
    c->xtofile = 1 ;            /* - -                */
    c->xsourceonly = 0 ;        /* -xso xsourceonly   */
}


void zb2x_free ( zb2x_all *all ) {
  zb2x_vectors *vec = all->vec ; 

  free_cvector ( vec->x , 1 , vec->N ) ; 
  free_cvector ( vec->z , 1 , vec->M ) ; 

  free_dvector ( vec->bias , 1 , vec->N ) ; 

}
