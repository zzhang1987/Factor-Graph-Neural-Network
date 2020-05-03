/* RAND.H - Random number generators. */
/* Copyright (c) 1992 by Radford M. Neal */

/* SET RANDOM NUMBER SEED. */

#define ran_seed(s) srandom(s)

/* GENERATE RANDOM NUMBERS. */

#define ranf() \
  ((double)random()/(1.0+(double)0x7fffffff)) /* Uniform from interval [0,1) */

#define ranu() \
  ((1.0+(double)random())/(2.0+(double)0x7fffffff))    /* Uniform from (0,1) */

#define rani(n) \
  ( (int) (ranf()*(n)) )		    /* Uniform from 0, 1, ..., (n-1) */

#define rann() \
  (cos(2.0*3.141592654*ranf()) * sqrt(-2.0*log(1.0-ranf()))) /* From standard Norml */

#define rane() \
  (-log(ranu()))		                  /* From exponential */

#define ranc() \
  (tan(3.141592654*(ranu()-0.5)))		                      /* From Cauchy */
