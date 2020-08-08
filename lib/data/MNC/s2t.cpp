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





#include "pybind11/pybind11.h"
#include "xtensor/xmath.hpp"
#define FORCE_IMPORT_ARRAY
#include "xtensor-python/pyarray.hpp"


template <typename T>
xt::pyarray<T> s2t(const xt::pyarray<T>& source,
                   int K, int N,
                   const std::string& Gfile,
    int smn=0){
    
}



#undef DNT
#undef NLNE

/*
<!-- hhmts start -->
Last modified: Wed Dec  6 14:29:52 1995
<!-- hhmts end -->
*/
