#include "ansi/r.h"
#include "ansi/rand2.h"
#include "ansi/mynr.h"
#include "ansi/cmatrix.h"
#include "pybind11/detail/common.h"
#include "radford/mod2mat.h"
#include "bnd/bnd.h"
#include "zb2x.h"


#include "pybind11/pybind11.h"
#include "xtensor/xmath.hpp"
#include "xtensor/xadapt.hpp"
#include "xtensor/xrandom.hpp"
#include <stdio.h>
#include <string>
#define FORCE_IMPORT_ARRAY
#include "xtensor-python/pyarray.hpp"

#include <iostream>

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
    if(sigma_b >= 1e-20) 
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

template <typename IntT, typename T>
xt::pyarray<IntT> zb2x(xt::pyarray<T>& z, int k, int n, std::string Afile, int xso, int bndloops){
    bnd_control         bndc ; 
    bnd_param           bndp ;
    zb2x_control           c ;
    zb2x_vectors         vec  ; 
    zb2x_all             all  ;
    alist_matrix         a ;
    if(z.shape().size() != 1){
        throw pybind11::value_error("source must be 1d array");
    }
    int zsize = z.shape()[0];
    
    int stillgoing = 1 , message = 0 , failures = 0 ; 
    char junk[1000] ; 

    all.bndp = &bndp ; 
    all.bndc = &bndc ; 
    all.vec  = &vec ; 
    all.a    = &a ; 
    all.c    = &c ; 
    bndp.a   = &a ;

    bnd_defaults(&bndp, &bndc);
    c_defaults(&c);

    all.c->K = k;
    all.c->N = n;
    strcpy(all.c->Afile, Afile.c_str());
    all.c->xsourceonly = xso;
    all.c->zfromfile = 0;
    all.c->zfixed = 0;
    //strcpy(all.bndc->logfile, "decode.log");
    //all.bndc->writelog = 1;
    all.bndc->loops = bndloops;
    
    make_sense(&c, &all);
    make_space(&c, &vec);

    read_allocate_alist(&a, (char *) Afile.c_str());
    check_alist_MN(&a, &vec);

    hook_zb2x_vec_to_bnd(&bndp, &vec);
    bnd_allocate(&bndp, &bndc);

    set_up_priors(&vec, &c);

    std::vector<IntT> res;
    int offset = 0;
    if((zsize % (k + n)) != 0){
        throw pybind11::value_error("The length of z vector can not be devided by k + n");
    }
    do {
        message++;
        int readto = k + n;
        std::cout << readto << std::endl;
        for(int i = 0; i < readto; i++){
            vec.bias[i + offset + 1] = z[i];
        }
        if(bndecode(&bndp, &bndc) > 0){
            // throw pybind11::value_error("Decoding fail at block " + std::to_string(message));
            std::cerr << "Decoding fail at block " + std::to_string(message) << std::endl; 
        }
        for(int i = 0; i < k; i++){
            res.push_back(vec.x[i + 1]);
        }
        offset += (k + n);
    }while(offset + k + n <= zsize);
    xt::xarray<IntT> final_res = xt::adapt(res, {res.size()});

    bnd_free(&bndp, &bndc);
    zb2x_free(&all);
    return final_res;
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
    m.def("zb2x", zb2x<long int, double>, "zb2x");
    m.def("init_seed", init_seed, "init_seed");
}


