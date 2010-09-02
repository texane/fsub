#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <sys/types.h>


/* matrix */

typedef struct matrix
{
  size_t n;
  double data[1];
} matrix_t;

static void matrix_create_lower(matrix_t** a, size_t n)
{ 
  /* create a lower triangular matrix.
     the digaonal is filled with ones.
   */

  
}

static inline double* matrix_at()
{
}

static inline const double* matrix_const_at()
{
}

static inline size_t matrix_size(const matrix_t* m)
{
  return m->n;
}


/* cholesky column algorithm */

static void cdiv(matrix_t* a, size_t j)
{
  /* col j is divided by the sqrt of its diagonal entry */
  const double sqrtj = sqrt(*matrix_at(a, j, j));
  *matrix_at(a, j, j) = sqrtj;
  for (size_t i = j + 1; i < a->n; ++i)
    *matrix_at(a, i, j) /= sqrtj;
}

static void cmod(matrix_t* a, size_t j, size_t k)
{
  /* col j is modified by a multiple of col k, with k<j */
  const double ajk = *matrix_const_at(a, j, k);
  for (size_t i = j; i < a->n; ++i)
    *matrix_at(a, i, j) -= *matrix_at(a, i, k) * ajk;
}

static void cholesky(matrix_t* a, size_t l, size_t h)
{
  /* apply cholesky on columns l to h */

  size_t i, j, k;

  for (j = l; j < h; ++j)
  {
    for (k = 0; k < (j - 1); ++k)
      cmod(a, j, k);
    cdiv(a, j);
  }
}

static void cholesky_par(matrix_t* a, size_t l, size_t h)
{
  parfor(l, h)
  {
    cholesky(a, l, h);
  }
}

static void solve_row(matrix_t* a, size_t i)
{
}


#if 0
static void backsub_cholesky_col_par
(a, n)
{
  for (j = 0; j < n; ++j)
  {
    parwork(a, j);
    sync();
    cdiv(a, n, j);
  }
}
#endif


#include <sys/types.h>


/* forward substitution algorithm
 */

struct fsub_context
{
  /* ax = b; */
  const matrix_t* a;
  const vector_t* x;
  vector_t* b;

  /* size of a k block */
  size_t ksize;

  /* size of a l block */
  size_t lsize;

  /* current kkblock (not inclusive) */
  volatile unsigned long kkpos __attribute__((aligned));

} fsub_context_t;


static void sub_square_block
(fsub_context_t* fsc, size_t i, size_t j, size_t n)
{
  /* b -= aij * x, a square of size n */

  const size_t first_i = i;
  const size_t last_i = i + n;
  const size_t first_j = j;
  const size_t last_j = j + n;

  for (; i < last_i; ++i)
  {
    /* local accumulation */
    double bi = 0.f;

    for (j = first_j; j < last_j; ++j)
      bi += matrix_at(fsc->a, i, j) * vector_at(fsc->x, i);

    /* update b -= ax */
    *vector_at(fsc->b, i) -= bi;
  }
}


static void sub_tri_block(fsub_context_t* fsc, size_t i, size_t n)
{
  /* b -= aij * x, a triangle of size n */

  const size_t first_i = i;
  const size_t last_i = i + n;

  for (; i < last_i; ++i)
  {
    /* local accumulation */
    double bi = 0.f;

    size_t j;
    for (j = first_i; j < i; ++j)
      bi += matrix_at(fsc->a, i, j) * vector_at(fsc->x, i);

    /* update b -= ax */
    *vector_at(fsc->b, i) -= bi;
  }
}


static void fsub(fsub_context_t* fsc)
{
  /* assume llcount */

  const size_t kk_count = matrix_size(fsc->a) / fsc->ksize;
  const size_t ll_count = matrix_size(fsc->a) / fsc->lsize;

  size_t kkpos;

  /* step0: solve the first llblock */
  fsub_ll();

  /* sequential code, slide the current kkblock */
  for (kkpos = 0; kkpos < kkcount; ++kkpos)
  {
    
  }
}


static void init_fsub_context
(fsub_context_t* fsc, matrix_t* )
{
  
}


int main(int ac, char** av)
{
  fsub_context_t fsc;

  matrix_create_lower(&a, 3);
  matrix_create_lower(&x, 3);
  matrix_create_lower(&b, 3);

  fsub(&fsc);

  return 0;
}
