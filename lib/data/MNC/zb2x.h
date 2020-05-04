/* 
   zb2x.h
*/

typedef struct {
#include "zb2x_var_str.h"
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
  
