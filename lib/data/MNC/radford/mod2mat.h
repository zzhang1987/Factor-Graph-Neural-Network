/* MOD2MAT.H - Interface to module for handling dense mod2 matrices. */


/* This module implements operations on matrices of mod2 elements (bits,
   with addition and multiplication being done modulo 2).  The matrices
   are stored with consecutive bits of a column packed into words, and 
   the procedures are implemented where possible using bit operations 
   on these words.  This is an appropriate representation when the matrices 
   are dense (ie, 0s and 1s are about equally frequent). */


/* PACKING OF BITS INTO WORDS.  Bits are packed into 32-bit words, with
   the low-order bit coming first. */

typedef unsigned long mod2word;	/* Data type that holds packed bits */

#define mod2_wordsize 32	/* Number of bits that fit in a mod2word */

/* Extract the i'th bit of a mod2word. */

#define mod2_getbit(w,i) (((w)>>(i))&1) 

/* Make a word like w, but with the i'th bit set to 1 (if it wasn't already). */

#define mod2_setbit1(w,i) ((w)|(1<<(i))) 

/* Make a word like w, but with the i'th bit set to 0 (if it wasn't already). */

#define mod2_setbit0(w,i) ((w)&(~(1<<(i)))) 


/* STRUCTURE REPRESENTING A MATRIX.  These structures are dynamically
   allocated using mod2mat_allocate (or by other procedures that call
   mod2mat_allocate).  They should be freed with mod2mat_free when no 
   longer required. */

typedef struct 
{
  int n_rows;		/* Number of rows in the matrix */
  int n_cols;		/* Number of columns in the matrix */

  int n_words;		/* Number of words used to store a column of bits */

  mod2word **col;	/* Pointer to array of pointers to columns */

  mod2word *bits;	/* Pointer to storage block for bits in this matrix 
                           (pieces of this block are pointed to from col */
} mod2mat;


/* PROCEDURES. */

mod2mat *mod2mat_allocate (int, int);
void mod2mat_free         (mod2mat *);

void mod2mat_print     (FILE *, mod2mat *);
int  mod2mat_write     (FILE *, mod2mat *);
mod2mat *mod2mat_read  (FILE *, int *);

int  mod2mat_get (mod2mat *, int, int);
void mod2mat_set (mod2mat *, int, int, int);

void mod2mat_copy  (mod2mat *, mod2mat *);
void mod2mat_clear (mod2mat *);

int mod2mat_equal (mod2mat *, mod2mat *);

void mod2mat_transpose (mod2mat *, mod2mat *);
void mod2mat_multiply  (mod2mat *, mod2mat *, mod2mat *);
int  mod2mat_invert    (mod2mat *, mod2mat *);

int mod2mat_forcibly_invert (mod2mat *, mod2mat *, int *, int *);
