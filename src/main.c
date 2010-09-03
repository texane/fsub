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

#if CONFIG_PARALLEL
/* forward decl */
struct parwork;
#endif

typedef struct fsub_context
{
  /* ax = b; */
  const matrix_t* a;
  vector_t* b;

  /* size of a kxk block */
  size_t ksize;

  /* size of a lxl block */
  size_t lsize;

#if CONFIG_PARALLEL
  struct parwork* parwork;
#endif

} fsub_context_t;


static void sub_rect_block
(fsub_context_t* fsc, index_t i, index_t j, size_t m, size_t n)
{
  /* b -= aij * x, a rect of width=m, height=n */

  const index_t first_j = j;
  const index_t last_j = j + m;
  const index_t last_i = i + n;

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
  return *ul;
}

static inline void write_atomic_ul(volatile unsigned long* ul, unsigned long n)
{
  *ul = n;
}

static inline void inc_atomic_ul(volatile unsigned long* ul)
{
  __sync_fetch_and_add(ul, 1UL);
}

static inline void dec_atomic_ul(volatile unsigned long* ul)
{
  __sync_fetch_and_add(ul, -1);
}


#if CONFIG_PARALLEL

/* spinlock */

typedef struct spinlock
{
  volatile unsigned long value;
} spinlock_t;

static inline void spinlock_init(spinlock_t* l)

  l->value = 0;
}

static inline int spinlock_trylock(spinlock_t* l)
{
  return __sync_compare_and_swap(&l->value, 0, 1);
}

static inline void spinlock_unlock(spinlock_t* l)
{
  __sync_fetch_and_and(&l->lock, 0);
}


/* parallel work */

typedef struct parwork
{
  spinlock_t lock;

  /* for global synchronization */
  volatile unsigned long refn __attribute__((aligned));

  /* area to be processed: (i,0,asize,j) */
  volatile index_t i;
  volatile index_t j;

  /* last processed row index (non inclusive) */
  volatile index_t last_i __attribute__((aligned));

} parwork_t;

static inline void parwork_init
(parwork_t* w, fsub_context_t* fsc, index_t i, index_t j, size_t m, size_t n)
{
#define CONFIG_THREAD_COUNT 16
  spinlock_init(&w->lock);
  w->refn = CONFIG_THREAD_COUNT;
  w->fsc = fsc;
  w->i = i;
  w->j = j;
  w->m = m;
  w->n = n;
}

static inline void parwork_synchronize(parwork_t* w)
{
  /* assume the lock is owned */
  while (read_atomic_ul(&w->refn) != 1UL)
    ;
}

static inline void parwork_lock(parwork_t* w)
{
  while (!spinlock_trylock(&w->lock))
    ;
}

static inline void parwork_lock_with_sync(parwork_t* w)
{
  dec_atomic_ul(&w->refn);
  parwork_lock(&w->lock);
  inc_atomic_ul(&w->refn);
}

static inline void parwork_unlock(parwork_t* w)
{
  spinlock_unlock(&w->lock);
}

static void sub_row(fsub_context_t* fsc, index_t i, size_t m)
{
  /* b -= aij*x with j in [0,m[ */
  const index_t last_j = j + m;
  index_t j;

  /* local accumulation */
  double sum = 0.f;

  for (j = 0; j < last_j; ++j)
    sum += *matrix_const_at(fsc->a, i, j) * *vector_const_at(fsc->b, j);

  /* update b -= sum(axi) */
  *vector_at(fsc->b, i) -= sum;
}

static void* parwork_entry(void* p)
{
  parwork_t* const w = (parwork_t*)p;

  index_t i;
  index_t j;

  while (!read_atomic_ul(&w->isdone))
  {
    /* lock and extract a line to process */
    parwork_lock_with_sync(w);
    i = w->i++;
    j = w->j;
    parwork_unlock(w);

    /* process row */
    sub_row(w->fsc, i, j);

    /* update the last processed i. no need for
       atomic read since we are the only updater
    */
    if (i == parwork->last_i)
      write_atomic_ul(&parwork->last_i, (unsigned long)i + 1);
  }

  /* todo: synchronize */

  return NULL;
}

#endif /* CONFIG_PARALLEL */


static void fsub_apply(fsub_context_t* fsc)
{
  /* assume fsc->lsize >= matrix_size(fsc->a) */
  /* assume (fsc->lsize % fsc->ksize) == 0 */
  /* assume (matrix_size(fsc->a) % fsc->lsize) == 0 */

  const size_t llcount = matrix_size(fsc->a) / fsc->lsize;

  /* solve the first llblock */
  sub_tri_block(fsc, 0, fsc->lsize);

  size_t band_height = matrix_size(fsc->a) - fsc->lsize;

#if CONFIG_PARALLEL
  /* post the band to process */
  parwork_t parwork;
  parwork_init(&parwork, fsc, fsc->lsize, 0, fsc->lsize, band_height);
#else 
  /* process the band sequentially */
  sub_rect_block(fsc, fsc->lsize, 0, fsc->lsize, band_height);
#endif

  /* slide along all the kkblocks on the diagonal */
  const index_t last_i = (llcount - 1) * fsc->lsize;
  index_t i;
  for (i = fsc->lsize; i < last_i; i += fsc->ksize)
  {
    const index_t next_i = i + fsc->ksize;

#if CONFIG_PARALLEL
    /* wait until left band processed
       todo: contribute to the work while waiting
    */
    while ((index_t)read_atomic_ul(&parwork.last_i) < next_i)
      ;
#endif

    /* process the kk block sequentially */
    sub_tri_block(fsc, i, fsc->ksize);

    /* compute the sub block dimensions */
    next_i = i + fsc->ksize;
    const index_t sj = i;
    const index_t si = sj + fsc->ksize;
    const index_t sm = sj + fsc->ksize;
    const index_t sn = parwork.i - (i + fsc->ksize);

#if CONFIG_PARALLEL
    /* sync on the parallel work, so that no
       one is working while we are updating
       the parwork area.
     */

    parwork_lock(&fsc->parwork.lock);
    parwork_synchronize(&fsc->parwork);

    /* update parwork area */
    fsc->parwork.i = to_i(kkpos);
    fsc->parwork.j = next_i;

    unlock_work(&fsc->parwork);

    /* process the band sequentially */
    band_height = fsc->par_work;
    sub_rect_block(fsc, fsc->lsize, 0, fsc->ksize, band_height);

    /* update kpos */
    inc_atomic_ul(&fsc->kpos);

    /* post the new band to process */
    push_work(sfc, bi, bj, fsc->ksize, matrix_size(fsc->a) - bi);
#else
    /* process the band sequentially */
    sub_rect_block(fsc, bi, bj, fsc->ksize, matrix_size(fsc->a) - bi);
#endif /* (CONFIG_PARALLEL == 0) */
  }

#if CONFIG_PARALLEL
  /* wait until everyone is done */
  lock_work(&fsc->parwork);
  parwork_synchronize(&fsc->parwork);
  write_atomic_ul(&fsc->parwork.isdone, 1UL);
  unlock_work(&fsc->parwork);
#endif

  /* the last llblock sequentially */
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
