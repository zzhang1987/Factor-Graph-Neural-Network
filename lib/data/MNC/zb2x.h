#ifndef _ZB2X_H_
#define _ZB2X_H_

#include <ansi/cmatrix.h>
#include "bnd/bnd.h"
/* 
   zb2x.h
*/

typedef struct {
    int verbose ;     /* verbosity            <0>    */
    char bfile[100] ; /* file: bias vector    <->    */
    int bnfromfile ;  /* -                    <1>    */
    int bsfromfile ;  /* -                    <1>    */
    int bfromfile ;   /* whether some of b comes from file <1>    */
    double bfn ;      /* constant value for bias of noise bits <0.1>  */
    double bfs ;      /* constant value for bias of signal bits <0.5>  */
    int K ;           /* number of source bits (NMN) <->    */
    int N ;           /* number of noise bits (NMN) <->    */
    char zfile[100] ; /* file: z vector       <->    */
    int zfromfile ;   /* -                    <0>    */
    int zfixed ;      /* constant value for z <0>    */
    char Afile[100] ; /* file: A matrix       <->    */
    int bsuffix ;     /* whether bfile should have suffices .0001,.0002... appended <0>    */
    int zsuffix ;     /* whether zfile should have suffices .0001,.0002... appended <0>    */
    int xsuffix ;     /* whether xfile should have suffices .0001,.0002... appended <0>    */
    char xfile[100] ; /* x vector             <->    */
    int xtofile ;     /* -                    <1>    */
    int xsourceonly ; /* write decoded source bits only <0>    */
} zb2x_control ;

typedef struct {
  unsigned char *x , *z ;

  double *bias ; 

  int M ; /* number of relationship bits */
  int N ; /* length of the vector that is being reconstructed */
  int NS ; /* how many of the N bits are `signal' bits */
  int NN ; /* and noise bits */
  int nsfrom , nsto ; /* which bits are the signal */
  int nnfrom , nnto ; /* and noise bits */
  int b_readfrom , b_readto ; 
  int xfrom , xto ; 
} zb2x_vectors ; 


typedef struct {
  zb2x_control *c ;
  zb2x_vectors   *vec ; 
  alist_matrix   *a ; 
  bnd_control    *bndc ; 
  bnd_param      *bndp ;
} zb2x_all ;


int make_sense(zb2x_control *c, zb2x_all *all);
int make_space(zb2x_control *c, zb2x_vectors *vec);
int check_alist_MN(alist_matrix *a, zb2x_vectors *v);
int hook_zb2x_vec_to_bnd(bnd_param *p, zb2x_vectors *vec);
void set_up_priors(zb2x_vectors *v, zb2x_control *c);
void c_defaults(zb2x_control *c);
void zb2x_free ( zb2x_all *all );
#endif
