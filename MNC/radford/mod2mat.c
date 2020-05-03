/* MOD2MAT.C - Procedures for handling dense mod2 matrices. */

/* All procedures in this module display an error message on standard 
   error and terminate the program if passed an invalid argument (indicative
   of a programming error).  Errors from lack of memory or from invalid
   contents of a file result in an error code being returned to the caller,
   with no message being printed. */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "mod2mat.h"


/* ALLOCATE SPACE FOR A MOD2 MATRIX.  Returns zero if there is not enough
   memory for the matrix. */

mod2mat *mod2mat_allocate 
( int n_rows, 		/* Number of rows in matrix */
  int n_cols		/* Number of columns in matrix */
)
{
  mod2mat *m;
  int j;

  if (n_rows<=0 || n_cols<=0)
  { fprintf(stderr,"mod2mat_allocate: Invalid number of rows or columns\n");
    exit(1);
  }
  m = calloc (1, sizeof *m);

  if (m==0) 
  { return 0;
  }

  m->n_rows = n_rows;
  m->n_cols = n_cols;

  m->n_words = (n_rows+mod2_wordsize-1) / mod2_wordsize;

  m->col = calloc (m->n_cols, sizeof *m->col);
  if (m->col==0)
  { free(m);
    return 0;
  }

  m->bits = calloc(m->n_words*m->n_cols, sizeof *m->bits);
  if (m->bits==0)
  { free(m->col);
    free(m);
    return 0;
  }

  for (j = 0; j<m->n_cols; j++)
  { m->col[j] = m->bits + j*m->n_words;
  }

  return m;
}


/* FREE SPACE OCCUPIED BY A MOD2 MATRIX.  The pointer passed should no
   longer be used after mod2mat_free returns. */

void mod2mat_free
( mod2mat *m		/* Matrix to free */
)
{ free(m->bits);
  free(m->col);
  free(m);
}


/* PRINT A MOD2 MATRIX IN HUMAN-READABLE FORM.  The matrix is printed
   as "0" and "1" characters, with one line of "0"s and "1"s for each
   row of the matrix. */

void mod2mat_print     
( FILE *f,
  mod2mat *m
)
{ 
  int i, j;

  for (i = 0; i<m->n_rows; i++)
  { for (j = 0; j<m->n_cols; j++)
    { fprintf(f," %d",mod2mat_get(m,i,j));
    }
    fprintf(f,"\n");
  }
}


/* WRITE A MOD2 MATRIX TO A FILE IN MACHINE-READABLE FORM.  Returns one
   if successful, zero if an error of some sort occurred. 

   The data written to the file consists of the number of rows and the
   number of columns (in machine-readable, not human-readable, form),
   followed by the bits in each column, packed into words. */

int mod2mat_write     
( FILE *f, 
  mod2mat *m
)
{ 
  int j;

  if (fwrite (&m->n_rows, sizeof m->n_rows, 1, f) != 1) 
  { return 0;
  }
  if (fwrite (&m->n_cols, sizeof m->n_rows, 1, f) != 1) 
  { return 0;
  }

  for (j = 0; j<m->n_cols; j++)
  {
    if (fwrite (m->col[j], sizeof *m->col[j], m->n_words, f) != m->n_words) 
    { return 0;
    }
  }

  return 1;
}


/* READ A MOD2 MATRIX STORED IN MACHINE-READABLE FORM FROM A FILE.  The data 
   is expected to be in the format written by mod2mat_write.  The value
   returned is zero if some sort of error occurred, in which case an error
   code is stored in the second argument, as follows:
 
       0:  No error
       1:  Error of some sort reading file
       2:  Contents of file are invalid
       3:  Unable to allocate space

   The file is left positioned after the last of the data read.  Other data
   (such as another matrix) may follow. */

mod2mat *mod2mat_read  
( FILE *f,
  int *error_code
)
{ 
  int n_rows, n_cols;
  mod2mat *m;
  int j;

  *error_code = 0;
  
  if (fread (&n_rows, sizeof n_rows, 1, f) != 1
   || fread (&n_cols, sizeof n_cols, 1, f) != 1)
  { *error_code = ferror(f) ? 1 : 2;
    return 0;
  }
  if (n_rows<=0 || n_cols<=0) 
  { *error_code = 2;
    return 0;
  }

  m = mod2mat_allocate(n_rows,n_cols);
  if (m==0)
  { *error_code = 3;
    return 0;
  }    

  for (j = 0; j<m->n_cols; j++)
  {
    if (fread (m->col[j], sizeof *m->col[j], m->n_words, f) != m->n_words) 
    { *error_code = ferror(f) ? 1 : 2;
      return 0;
    }
  }

  return m;
}


/* GET AN ELEMENT FROM A MOD2 MATRIX. */

int mod2mat_get  
( mod2mat *m, 		/* Matrix to get element from */
  int row, 		/* Row of element (starting with zero) */
  int col		/* Column of element (starting with zero) */
)
{
  if (row<0 || row>=m->n_rows || col<0 || col>=m->n_cols)
  { fprintf(stderr,"mod2mat_get: row or column index out of bounds\n");
    exit(1);
  }

  return mod2_getbit (m->col[col][row/mod2_wordsize], row%mod2_wordsize);
}


/* SET AN ELEMENT IN A MOD2 MATRIX. */

void mod2mat_set 
( mod2mat *m, 		/* Matrix to modify element of */
  int row, 		/* Row of element to modify (starting with zero) */
  int col, 		/* Column of element to modify (starting with zero) */
  int value		/* New value of element (0 or 1) */
)
{ 
  mod2word *w;

  if (row<0 || row>=m->n_rows || col<0 || col>=m->n_cols)
  { fprintf(stderr,"mod2mat_set: row or column index out of bounds\n");
    exit(1);
  }

  w = &m->col[col][row/mod2_wordsize];

  *w = value ? mod2_setbit1(*w,row%mod2_wordsize) 
             : mod2_setbit0(*w,row%mod2_wordsize);
}


/* COPY A MOD2 MATRIX.  Copies the elements of the first matrix to the
   second matrix (which must already exist, and have the same number of
   rows and columns as the first matrix). */

void mod2mat_copy
( mod2mat *m,		/* Matrix to copy */
  mod2mat *r		/* Place to store copy of matrix */
)
{ 
  int k, j;

  if (m->n_rows!=r->n_rows || m->n_cols!=r->n_cols)
  { fprintf(stderr,"mod2mat_copy: Matrices have different dimensions\n");
    exit(1);
  }

  for (j = 0; j<m->n_cols; j++)
  { for (k = 0; k<m->n_words; k++)
    { r->col[j][k] = m->col[j][k];
    }
  }
}


/* CLEAR A MOD2 MATRIX.  Sets all the elements of the matrix to 0. */

void mod2mat_clear
( mod2mat *r
)
{
  int k, j;

  for (j = 0; j<r->n_cols; j++)
  { for (k = 0; k<r->n_words; k++)
    { r->col[j][k] = 0;
    }
  }

}


/* SEE WHETHER TWO MATRICES ARE EQUAL.  Returns one if every element of
   the first matrix is equal to the corresponding element of the second
   matrix.  The two matrices must have the same number of rows and the
   same number of columns. */

int mod2mat_equal
( mod2mat *m1,
  mod2mat *m2
)
{
  int k, j, w;
  mod2word m;

  if (m1->n_rows!=m2->n_rows || m1->n_cols!=m2->n_cols)
  { fprintf(stderr,"mod2mat_equal: Matrices have different dimensions\n");
    exit(1);
  }

  w = m1->n_words;

  /* Form a mask that has 1s in the lower bit positions corresponding to
     bits that contain information in the last word of a matrix column. */

  m = (1 << (mod2_wordsize - (w*mod2_wordsize-m1->n_rows))) - 1;
  
  for (j = 0; j<m1->n_rows; j++)
  {
    for (k = 0; k<w-1; k++)
    { if (m1->col[j][k] != m2->col[j][k]) return 0;
    }

    if ((m1->col[j][k]&m) != (m2->col[j][k]&m)) return 0;
  }

  return 1;
}


/* COMPUTE THE TRANSPOSE OF A MATRIX.  Stores the transpose of the 
   first matrix in the second matrix (which must already exist, and
   have as many rows as the first matrix has columns, and as many
   columns as the first matrix has rows). */

void mod2mat_transpose
( mod2mat *m,		/* Matrix to compute transpose of (left unchanged) */
  mod2mat *r		/* Result of transpose operation */
)
{
  mod2word w, v, *p;
  int k1, j1, i2, j2;

  if (m->n_rows!=r->n_cols || m->n_cols!=r->n_rows)
  { fprintf(stderr,
     "mod2mat_transpose: Matrices have incompatible dimensions\n");
    exit(1);
  }

  mod2mat_clear(r);

  for (j1 = 0; j1<m->n_cols; j1++)
  { 
    i2 = j1 / mod2_wordsize;
    v = 1 << (j1 % mod2_wordsize);

    p = m->col[j1];
    k1 = 0;

    for (j2 = 0; j2<r->n_cols; j2++)
    { if (k1==0)
      { w = *p++;
        k1 = mod2_wordsize;
      }
      if (w&1)
      { r->col[j2][i2] |= v;
      }
      w >>= 1;     
      k1 -= 1;
    }
  }
}


/* MULTIPLY TWO MOD2 MATRICES.  Multiplies the first matrix by the second 
   matrix, storing the result in the third matrix.  Neither of the first 
   two matrices is changed by this operation.  The three matrices must have 
   compatible numbers of rows and columns. 

   The algorithm used runs faster if the second matrix (right operand of the
   multiply) is sparse, but it is also appropriate for dense matrices.  This
   procedure could be speeded up a bit by replacing the call of mod2mat_get
   with in-line code that avoids division, but this doesn't seem worthwhile
   at the moment. */

void mod2mat_multiply 
( mod2mat *m1, 		/* Left operand of multiply */
  mod2mat *m2,		/* Right operand of multiply */
  mod2mat *r		/* Place to store result of multiply */
)
{
  int i, j, k;
  if (m1->n_cols!=m2->n_rows || m1->n_rows!=r->n_rows || m2->n_cols!=r->n_cols)
  { fprintf(stderr,"mod2mat_multiply: Matrices have incompatible dimensions\n");
    exit(1);
  }

  mod2mat_clear(r);

  for (j = 0; j<m2->n_cols; j++)
  { for (i = 0; i<m2->n_rows; i++)
    { if (mod2mat_get(m2,i,j))
      { for (k = 0; k<r->n_words; k++)
        { r->col[j][k] ^= m1->col[i][k];
        }
      }
    }
  }
}


/* INVERT A MOD2 MATRIX.  Inverts the first matrix passed, destroying its
   contents in the process (though it remains a valid matrix for storing
   into later).  The matrix to be inverted must have the same number of
   rows as columns.  The inverse of this matrix is stored in the matrix 
   passed as the second argument (which must already exist, and have the 
   same number of rows and columns as the first).

   The value returned by mod2mat_invert is one if the inversion was 
   successful, zero if the matrix turned out to be singular (in which 
   case the contents of both the original matrix and the result matrix 
   will be garbage). */

int mod2mat_invert 
( mod2mat *m, 		/* The matrix to find the inverse of (destroyed) */
  mod2mat *r		/* Place to store the inverse */
)
{
  mod2word *s, *t;
  int i, j, k, n, w, k0, b0;

  if (m->n_rows!=m->n_cols)
  { fprintf(stderr,"mod2mat_invert: Matrix to invert is not square\n");
    exit(1);
  }

  n = m->n_rows;
  w = m->n_words;

  if (r->n_rows!=n || r->n_cols!=n)
  { fprintf(stderr,
     "mod2mat_invert: Matrix to receive inverse has wrong dimensions\n");
    exit(1);
  }

  mod2mat_clear(r);
  for (i = 0; i<n; i++) 
  { mod2mat_set(r,i,i,1);
  }

  for (i = 0; i<n; i++)
  { 
    k0 = i/mod2_wordsize;
    b0 = i%mod2_wordsize;

    for (j = i; j<n; j++) 
    { if (mod2_getbit(m->col[j][k0],b0)) break;
    }

    if (j==n) return 0;

    if (j!=i)
    {
      t = m->col[i];
      m->col[i] = m->col[j];
      m->col[j] = t;

      t = r->col[i];
      r->col[i] = r->col[j];
      r->col[j] = t;
    }

    for (j = 0; j<n; j++)
    { if (j!=i && mod2_getbit(m->col[j][k0],b0))
      { s = m->col[j];
        t = m->col[i];
        for (k = k0; k<w; k++) s[k] ^= t[k];
        s = r->col[j];
        t = r->col[i];
        for (k = 0; k<w; k++) s[k] ^= t[k];
      }
    }
  }

  return 1;
}


/* FORCIBLY INVERT A MOD2 MATRIX.  Inverts the first matrix passed, destroying 
   its contents in the process (though it remains a valid matrix for storing
   into later).  The matrix to be inverted must have the same number of
   rows as columns.  The inverse of this matrix is stored in the matrix 
   passed as the second argument (which must already exist, and have the 
   same number of rows and columns as the first).

   If the matrix to be inverted is singular, the inversion proceeds anyway,
   with bits of the matrix being changed as needed to create a non-singular
   matrix.  The value returned by mod2mat_invert is the number of elements 
   that had to be changed to make inversion possible (zero, if the original 
   matrix was non-singular). 

   The row and column indexes of the elements of the original matrix 
   that were changed are stored in the arrays passed as the last two
   elements.  These arrays must have as many elements as the dimension
   of the matrix.  (This is so even if it is known that fewer elements
   than this will be changed, as these arrays are also used as temporary
   storage by this routine.) */

int mod2mat_forcibly_invert 
( mod2mat *m, 		/* The matrix to find the inverse of (destroyed) */
  mod2mat *r,		/* Place to store the inverse */
  int *a_row,		/* Place to store row indexes of altered elements */
  int *a_col		/* Place to store column indexes of altered elements */
)
{
  mod2word *s, *t;
  int i, j, k, n, w, k0, b0;
  int u, c;

  if (m->n_rows!=m->n_cols)
  { fprintf(stderr,"mod2mat_invert: Matrix to invert is not square\n");
    exit(1);
  }

  n = m->n_rows;
  w = m->n_words;

  if (r->n_rows!=n || r->n_cols!=n)
  { fprintf(stderr,
     "mod2mat_invert: Matrix to receive inverse has wrong dimensions\n");
    exit(1);
  }

  mod2mat_clear(r);
  for (i = 0; i<n; i++) 
  { mod2mat_set(r,i,i,1);
  }

  for (i = 0; i<n; i++)
  { a_row[i] = -1;
    a_col[i] = i;
  }

  for (i = 0; i<n; i++)
  { 
    k0 = i/mod2_wordsize;
    b0 = i%mod2_wordsize;

    for (j = i; j<n; j++) 
    { if (mod2_getbit(m->col[j][k0],b0)) break;
    }

    if (j==n)
    { j = i;
      mod2mat_set(m,i,j,1);
      a_row[i] = i;
    }

    if (j!=i)
    { 
      t = m->col[i];
      m->col[i] = m->col[j];
      m->col[j] = t;

      t = r->col[i];
      r->col[i] = r->col[j];
      r->col[j] = t;

      u = a_col[i];
      a_col[i] = a_col[j];
      a_col[j] = u;
    }

    for (j = 0; j<n; j++)
    { if (j!=i && mod2_getbit(m->col[j][k0],b0))
      { s = m->col[j];
        t = m->col[i];
        for (k = k0; k<w; k++) s[k] ^= t[k];
        s = r->col[j];
        t = r->col[i];
        for (k = 0; k<w; k++) s[k] ^= t[k];
      }
    }
  }

  c = 0;
  for (i = 0; i<n; i++)
  { if (a_row[i]!=-1)
    { a_row[c] = a_row[i];
      a_col[c] = a_col[i];
      c += 1;
    }
  }

  return c;
}
