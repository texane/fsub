#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <sys/types.h>


/* matrix
 */

typedef size_t index_t;

typedef struct matrix
{
  size_t n;
  double data[1];
} matrix_t;

static inline const double* matrix_const_at
(const matrix_t* m, index_t i, index_t j)
{
  return m->data + i * m->n + j;
}

static inline double* matrix_at(matrix_t* m, index_t i, index_t j)
{
  return (double*)matrix_const_at(m, i, j);
}

static inline size_t matrix_size(const matrix_t* m)
{
  return m->n;
}

static int matrix_create_lower(matrix_t** m, size_t n)
{ 
  /* create a lower triangular matrix.
     the digaonal is filled with ones.
   */

  const size_t total_size =
    offsetof(matrix_t, data) + (n * n) * sizeof(double);

  *m = malloc(total_size);
  if (*m == NULL)
    return -1;

  (*m)->n = n;

  size_t i;
  for (i = 0; i < n; ++i)
  {
    size_t j;
    for (j = 1; j <= i; ++j)
      *matrix_at(*m, i, j - 1) = 1.f; /* rand() */
    *matrix_at(*m, i, i) = 1.f;
  }

  return 0;
}

static void matrix_destroy(matrix_t* m)
{
  free(m);
}

static void __attribute__((unused)) matrix_print(const matrix_t* m)
{
  size_t i, j;

  for (i = 0; i < m->n; ++i)
  {
    for (j = 0; j < m->n; ++j)
      printf(" %lf", *matrix_const_at(m, i, j));
    printf("\n");
  }
}


/* vector
 */

typedef matrix_t vector_t;

static inline const double* vector_const_at(const vector_t* v, index_t i)
{
  return v->data + i;
}

static inline double* vector_at(vector_t* v, index_t i)
{
  return (double*)vector_const_at(v, i);
}

static inline size_t vector_size(const vector_t* v)
{
  return v->n;
}

static int vector_create_empty(vector_t** v, size_t n)
{
  const size_t total_size = offsetof(vector_t, data) + n * sizeof(double);

  *v = malloc(total_size);
  if (*v == NULL)
    return -1;

  (*v)->n = n;

  return 0;
}

static int vector_create(vector_t** v, size_t n)
{
  if (vector_create_empty(v, n) == -1)
    return -1;

  /* fill vector */
  size_t i;
  for (i = 0; i < n; ++i)
    *vector_at(*v, i) = 1.f; /* rand() */

  return 0;
}

static void vector_destroy(vector_t* v)
{
  free(v);
}

static void vector_print(const vector_t* v)
{
  size_t i;
  for (i = 0; i < v->n; ++i)
    printf(" %d", (int)v->data[i]);
  printf("\n");
}


/* forward substitution algorithm
 */

typedef struct fsub_context
{
  /* ax = b; */
  const matrix_t* a;
  vector_t* b;

  /* size of a kxk block */
  size_t ksize;

  /* size of a lxl block */
  size_t lsize;

  /* current kkblock (not inclusive) */
  volatile unsigned long kkpos __attribute__((aligned));

} fsub_context_t;


static void __attribute__((unused)) sub_square_block
(fsub_context_t* fsc, index_t i, index_t j, size_t n)
{
  /* b -= aij * x, a square of size n */

  const index_t last_i = i + n;
  const index_t first_j = j;
  const index_t last_j = j + n;

  for (; i < last_i; ++i)
  {
    /* local accumulation */
    double sum = 0.f;

    for (j = first_j; j < last_j; ++j)
      sum += *matrix_const_at(fsc->a, i, j) * *vector_const_at(fsc->b, j);

    /* update b -= sum(axi) */
    *vector_at(fsc->b, i) -= sum;
  }
}


static void sub_tri_block(fsub_context_t* fsc, index_t i, size_t n)
{
  /* b -= aij * x, a triangle of size n */

  const index_t first_i = i;
  const index_t last_i = i + n;

  for (; i < last_i; ++i)
  {
    /* local accumulation */
    double sum = 0.f;

    index_t j;
    for (j = first_i; j < i; ++j)
      sum += *matrix_const_at(fsc->a, i, j) * *vector_const_at(fsc->b, j);

    /* update b -= sum(axi) */
    *vector_at(fsc->b, i) -= sum;
  }
}

static inline unsigned long read_atomic_ul(volatile unsigned long* ul)
{
  return __sync_fetch_and_or(ul, 0UL);
}

static inline void write_atomic_ul(volatile unsigned long* ul, unsigned long n)
{
  __sync_fetch_and_and(ul, n);
}

static inline void inc_atomic_ul(volatile unsigned long* ul)
{
  __sync_fetch_and_add(ul, 1UL);
}

static void fsub_apply(fsub_context_t* fsc)
{
  /* assume fsc->lsize >= matrix_size(fsc->a) */
  /* assume (fsc->lsize % fsc->ksize) == 0 */
  /* assume (matrix_size(fsc->a) % fsc->lsize) == 0 */

  const size_t llcount = matrix_size(fsc->a) / fsc->lsize;

  /* step0: solve the first llblock
   */

  sub_tri_block(fsc, 0, fsc->lsize);
  write_atomic_ul(&fsc->kkpos, 1UL);

  /* step1: slide (kcount-1) kkblocks along the diagonal
     every processed kkblock unlocks a band of parallelism
     fsc->kkpos maintains the current kk block position
     i is the same as fsc->kkpos in matrix coordinates
   */

  const index_t last_i = (llcount - 1) * fsc->lsize;

  index_t i;
  for (i = 0; i < last_i; i += fsc->ksize, inc_atomic_ul(&fsc->kkpos))
    sub_tri_block(fsc, i, fsc->ksize);

  /* step2: solve the last llblock
   */

  if (llcount > 1)
    sub_tri_block(fsc, i, fsc->lsize);
}


static void fsub_initialize
(fsub_context_t* fsc, const matrix_t* a, vector_t* b)
{
#define CONFIG_KSIZE (1)
#define CONFIG_LSIZE (CONFIG_KSIZE * 3)

  fsc->a = a;
  fsc->b = b;

  fsc->ksize = CONFIG_KSIZE;
  fsc->lsize = CONFIG_LSIZE;

  fsc->kkpos = 0UL;
}


int main(int ac, char** av)
{
  fsub_context_t fsc;

  matrix_t* a = NULL;
  vector_t* b = NULL;

#define CONFIG_ASIZE (3 * CONFIG_LSIZE)
  if (matrix_create_lower(&a, CONFIG_ASIZE) == -1)
    goto on_error;

  if (vector_create(&b, CONFIG_ASIZE) == -1)
    goto on_error;

  /* apply forward substitution: ax = b */
  fsub_initialize(&fsc, a, b);
  fsub_apply(&fsc);

  vector_print(b);

 on_error:
  if (a != NULL)
    matrix_destroy(a);
  if (b != NULL)
    vector_destroy(b);

  return 0;
}
