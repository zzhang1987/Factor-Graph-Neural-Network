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


#include "pybind11/pybind11.h"
#include "xtensor/xmath.hpp"
#include "xtensor/xadapt.hpp"
#include "xtensor/xrandom.hpp"
#include <stdio.h>
#define FORCE_IMPORT_ARRAY
#include "xtensor-python/pyarray.hpp"


template <typename T>
xt::pyarray<T> s2t(const xt::pyarray<T>& source,
                   int K, int N,
                   const std::string& Gfile,
                   bool smn=true){
    mod2mat *G, *s, *t;
    std::vector<T> res;
    int code, done;
    FILE *fp = fopen(Gfile.c_str(), "rb");
    if(!fp){
        throw pybind11::value_error("Gfile does not exists!");
    }
    if(source.shape().size() != 1){
        throw pybind11::value_error("source must be 1d array");
    }
    int tsize = source.shape()[0];
    s = mod2mat_allocate(K, 1);
    t = mod2mat_allocate(N, 1);
    G = mod2mat_read(fp, &code);
    fclose(fp);
    done = 0;
    int message = 0;
    int offset = 0;
    do {
        int b = 0;
        message ++;
        for( b = 0; b < K; b++){
            if(offset + b >= tsize){
                // std::cout << "Done" << std::endl;
                done = 1;
                break;
            }
            int bit = source[offset + b];
            mod2mat_set(s, b, 0, bit);
        }
        offset += b;
        if (done){
            if(b == 0) break;
            int bit = 0;
            for(; b < K; b++){
                mod2mat_set(s, b, 0, bit);
            }
        }
        mod2mat_multiply(G, s, t);
        if(smn){
            for(b = 0; b < K; b++){
                res.push_back(mod2mat_get(s, b, 0));
            }
        }
        for(b = 0; b < N; b++){
            res.push_back(mod2mat_get(t, b, 0));
        }
    }while(!done);
    
    int res_size = res.size();
    // std::cout << res.size() << std::endl;
    xt::xarray<T> final_res = xt::adapt(res, {res_size});
    mod2mat_free(s);
    mod2mat_free(t);
    mod2mat_free(G);
    return final_res;
}

// template <typename T>
template <typename T>
xt::pyarray<T> t2y(xt::pyarray<long int> t, T snr_db, T sigma_b, T rho){
    T gcx = std::pow(10, snr_db / 20.0);
    size_t size = t.shape()[0];
    xt::xarray<T> res = 2 * gcx * (xt::pyarray<T>(t) - 0.5) + xt::random::randn<T>({size});

    T sigma = gcx * sigma_b;
    for(int i = 0; i < size; i++){
        auto r = xt::random::rand<T>({1});
        if(r[0] < rho){
            T noise = xt::random::randn<T>({1}, 0, sigma)[0];
            res[i] += noise;
        }
    }
    return res;
}

template <typename T>
xt::pyarray<T> y2b(xt::pyarray<T>& t, T snr_db){
    T gcx = std::pow(10, snr_db / 20.0);
    return 1.0 / (1.0 + xt::exp(-2.0 * gcx * t));
}

void init_seed(int seed){
    xt::random::seed(seed);
}

PYBIND11_MODULE(MNC, m) {
    xt::import_numpy();
    m.doc() = "Test module for xtensor python bindings";
    m.def("s2t", s2t<long int>, "Non-uniform motion blur convolution.");
    m.def("t2y", t2y<double>, "t2y");
    m.def("y2b", y2b<double>, "y2b");
    m.def("init_seed", init_seed, "init_seed");
}

#undef DNT
#undef NLNE

/*
<!-- hhmts start -->
Last modified: Wed Dec  6 14:29:52 1995
<!-- hhmts end -->
*/
