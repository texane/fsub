#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/time.h>

#ifndef __USE_GNU
#define __USE_GNU
#endif
#include <sched.h>
#include <pthread.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>


/* config */
#define CONFIG_THREAD_COUNT 4
#define CONFIG_KSIZE (4)
#define CONFIG_LSIZE (16 * CONFIG_KSIZE)
#define CONFIG_ASIZE (100 * CONFIG_LSIZE)
#define CONFIG_ITER_COUNT 5
#define CONFIG_TIME 1
#define CONFIG_CHECK 0


/* matrix */

typedef size_t index_t;

typedef gsl_matrix matrix_t;

static inline const double* matrix_const_at
(const matrix_t* m, index_t i, index_t j)
{
  return gsl_matrix_const_ptr(m, i, j);
}

static inline double* matrix_at(matrix_t* m, index_t i, index_t j)
{
  return (double*)matrix_const_at(m, i, j);
}

static inline size_t matrix_size(const matrix_t* m)
{
  return m->size1;
}

static inline int matrix_create_empty(matrix_t** m, size_t n)
{
  *m = gsl_matrix_alloc(n, n);
  if (*m == NULL)
    return -1;
  return 0;
}

static int matrix_create_lower(matrix_t** m, size_t n)
{ 
  /* create a lower triangular matrix.
     the digaonal is filled with ones.
   */

  if (matrix_create_empty(m, n) == -1)
    return -1;

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
  gsl_matrix_free(m);
}

static void matrix_copy(matrix_t* a, const matrix_t* b)
{
  gsl_matrix_memcpy(a, b);
}

static void __attribute__((unused)) matrix_print(const matrix_t* m)
{
  const size_t size = m->size1;

  size_t i, j;

  for (i = 0; i < size; ++i)
  {
    for (j = 0; j < size; ++j)
      printf(" %lf", *matrix_const_at(m, i, j));
    printf("\n");
  }
}


/* vector
 */

typedef gsl_vector vector_t;

static inline const double* vector_const_at(const vector_t* v, index_t i)
{
  return gsl_vector_const_ptr(v, i);
}

static inline double* vector_at(vector_t* v, index_t i)
{
  return (double*)vector_const_at(v, i);
}

static inline size_t vector_size(const vector_t* v)
{
  return v->size;
}

static inline int vector_create_empty(vector_t** v, size_t n)
{
  *v = gsl_vector_alloc(n);
  if (*v == NULL)
    return -1;

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
  gsl_vector_free(v);
}

static void __attribute__((unused)) vector_print(const vector_t* v)
{
  size_t i;
  for (i = 0; i < v->size; ++i)
    printf(" %d", (int)*vector_const_at(v, i));
  printf("\n");
}

static void vector_copy(vector_t* a, const vector_t* b)
{
  gsl_vector_memcpy(a, b);
}

#if CONFIG_CHECK
static int vector_cmp(const vector_t* a, const vector_t* b)
{
  size_t i;
  for (i = 0; i < a->size; ++i)
    if (*vector_const_at(a, i) != *vector_const_at(b, i))
      return -1;
  return 0;
}
#endif


/* atomic and spinlock
 */

typedef struct spinlock
{
  volatile unsigned long value;
} spinlock_t;

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
  __sync_fetch_and_add(ul, 1);
}

static inline void dec_atomic_ul(volatile unsigned long* ul)
{
  __sync_fetch_and_add(ul, -1);
}

static inline int cas_atomic_ul
(volatile unsigned long* a, unsigned long b, unsigned long c)
{
  return __sync_bool_compare_and_swap(a, b, c);
}

static inline void spinlock_init(spinlock_t* l)
{
  l->value = 0;
}

static inline int spinlock_trylock(spinlock_t* l)
{
  return __sync_bool_compare_and_swap(&l->value, 0, 1);
}

static inline void spinlock_unlock(spinlock_t* l)
{
  __sync_fetch_and_and(&l->value, 0);
}


/* thread pool
 */

struct fsub_context;

typedef struct thread_block
{
#define THREAD_STATE_ZERO 0
#define THREAD_STATE_DONE 1
  volatile unsigned long state;

  pthread_t thread;
  unsigned long tid;

  struct fsub_context* fsc;

} thread_block_t;

typedef struct thread_pool
{
  thread_block_t tbs[CONFIG_THREAD_COUNT];
} thread_pool_t;


/* forward substitution algorithm
 */

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


typedef struct fsub_context
{
  /* ax = b; */
  const matrix_t* a;
  vector_t* b;

  /* size of a kxk block */
  size_t ksize;

  /* size of a lxl block */
  size_t lsize;

  parwork_t parwork;
  thread_pool_t pool;

} fsub_context_t;


static inline void parwork_init(parwork_t* w)
{
  spinlock_init(&w->lock);
  w->refn = CONFIG_THREAD_COUNT;

#define INVALID_MATRIX_INDEX ((index_t)-1)
  w->i = INVALID_MATRIX_INDEX;
  w->j = INVALID_MATRIX_INDEX;
}

static inline void parwork_synchronize(parwork_t* w)
{
  /* assume the lock is owned. wait until refn
     drops to 1 meaning everyone reached it.
   */

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
  parwork_lock(w);
  inc_atomic_ul(&w->refn);
}

static inline void parwork_unlock(parwork_t* w)
{
  spinlock_unlock(&w->lock);
}

static void sub_row(fsub_context_t* fsc, index_t i, size_t m)
{
  /* b -= aij*x with j in [0,m[ */
  index_t j;

  /* local accumulation */
  double sum = 0.f;

  for (j = 0; j < m; ++j)
    sum += *matrix_const_at(fsc->a, i, j) * *vector_const_at(fsc->b, j);

  /* update b -= sum(axi) */
  *vector_at(fsc->b, i) -= sum;
}

static void parwork_next_row(fsub_context_t* fsc)
{
  /* continue processing the parallel work */

  const size_t asize = matrix_size(fsc->a);
  parwork_t* const w = &fsc->parwork;

  index_t i = INVALID_MATRIX_INDEX;
  index_t j = INVALID_MATRIX_INDEX;

  parwork_lock_with_sync(w);
  if (w->i < asize)
  {
    i = w->i++;
    j = w->j;
  }
  parwork_unlock(w);

  /* did not get work */
  if (i == INVALID_MATRIX_INDEX)
    return ;

  /* process the row */
  sub_row(fsc, i, j);

  /* update the last processed i */
  inc_atomic_ul(&w->last_i);
}

static void* parwork_entry(void* p)
{
  /* thief thread entry */

  thread_block_t* const tb = (thread_block_t*)p;
  fsub_context_t* const fsc = tb->fsc;

  while (1)
  {
    switch (read_atomic_ul(&tb->state))
    {
    case THREAD_STATE_DONE:
      return NULL;
      break;

    default:
      /* process next row in parallel area */
      parwork_next_row(fsc);
      break;
    }
  }

  return NULL;
}

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

static void fsub_apply(fsub_context_t* fsc)
{
  /* assume fsc->lsize >= matrix_size(fsc->a) */
  /* assume (fsc->lsize % fsc->ksize) == 0 */
  /* assume (matrix_size(fsc->a) % fsc->lsize) == 0 */

  const size_t llcount = matrix_size(fsc->a) / fsc->lsize;

  /* solve the first llblock */
  sub_tri_block(fsc, 0, fsc->lsize);

  /* post the band to process */
  parwork_lock(&fsc->parwork);
  fsc->parwork.i = fsc->lsize;
  fsc->parwork.j = fsc->lsize;
  fsc->parwork.last_i = fsc->lsize;
  parwork_unlock(&fsc->parwork);

  /* slide along all the kkblocks on the diagonal */
  const index_t last_i = (llcount - 1) * fsc->lsize;
  index_t i;
  for (i = fsc->lsize; i < last_i; i += fsc->ksize)
  {
    const index_t next_i = i + fsc->ksize;

    /* wait until left band processed */
    while ((index_t)read_atomic_ul(&fsc->parwork.last_i) < next_i)
      parwork_next_row(fsc);

    /* sync on the parallel work, so that no
       one is working while we are updating
       the parwork area.
     */
    parwork_lock(&fsc->parwork);

    parwork_synchronize(&fsc->parwork);

    /* process the kk block sequentially */
    sub_tri_block(fsc, i, fsc->ksize);

    /* update parwork area j (which is next_i) */
    fsc->parwork.j = next_i;

    /* updating the parwork area may have created
       a hole below the kkblock. this hole dim are:
       i + fsc->ksize, i, parwork.last_i, i + fsc->ksize
       we must save parwork.last_i before unlocking workers
    */
    const index_t saved_i = fsc->parwork.last_i;

    parwork_unlock(&fsc->parwork);

    /* process the band sequentially */
    const size_t band_height = saved_i - next_i;
    sub_rect_block(fsc, next_i, i, band_height, fsc->ksize);
  }

  /* wait until remaining left bands processed */
  const size_t asize = matrix_size(fsc->a);
  while ((index_t)read_atomic_ul(&fsc->parwork.last_i) < asize)
    parwork_next_row(fsc);

  /* the last llblock sequentially */
  if (llcount > 1)
    sub_tri_block(fsc, i, fsc->lsize);
}


static void fsub_initialize
(fsub_context_t* fsc, const matrix_t* a, vector_t* b)
{
  fsc->a = a;
  fsc->b = b;

  fsc->ksize = CONFIG_KSIZE;
  fsc->lsize = CONFIG_LSIZE;

  parwork_init(&fsc->parwork);

  size_t tid;
  for (tid = 0; tid < CONFIG_THREAD_COUNT; ++tid)
  {
    thread_block_t* const tb = &fsc->pool.tbs[tid];

    tb->state = THREAD_STATE_ZERO;
    tb->tid = tid;
    tb->fsc = fsc;

    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(tid, &cpuset);

    if (tid == 0)
    {
      pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
      continue ;
    }

    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &cpuset);
    pthread_create(&tb->thread, &attr, parwork_entry, (void*)tb);
  }
}


static void fsub_finalize(fsub_context_t* fsc)
{
  size_t tid;
  for (tid = 1; tid < CONFIG_THREAD_COUNT; ++tid)
  {
    write_atomic_ul(&fsc->pool.tbs[tid].state, THREAD_STATE_DONE);
    pthread_join(fsc->pool.tbs[tid].thread, NULL);
  }
}


static void __attribute__((unused))
fsub_apply_ref(matrix_t* a, vector_t* x, const vector_t* b)
{
  gsl_permutation* const p = gsl_permutation_alloc(vector_size(b));
  if (p == NULL)
    return ;

  int signum;
  gsl_linalg_LU_decomp(a, p, &signum);
  gsl_linalg_LU_solve(a, p, b, x);

  gsl_permutation_free(p);
}


static void fsub_apply_seq(matrix_t* a, vector_t* b)
{
  /* b -= aij * x, a triangle of size n */

  const size_t asize = matrix_size(a);

  index_t i, j;

  for (i = 0; i < asize; ++i)
  {
    /* local accumulation */
    double sum = 0.f;

    for (j = 0; j < i; ++j)
      sum += *matrix_const_at(a, i, j) * *vector_const_at(b, j);

    /* update b -= sum(axi) */
    *vector_at(b, i) -= sum;
  }
}


/* main */

int main(int ac, char** av)
{
  fsub_context_t fsc;

  matrix_t* a = NULL;
  matrix_t* const_a = NULL;
  vector_t* x = NULL;
  vector_t* const_b = NULL;
  vector_t* b = NULL;

  if (matrix_create_lower(&const_a, CONFIG_ASIZE) == -1)
    goto on_error;
  if (matrix_create_empty(&a, CONFIG_ASIZE) == -1)
    goto on_error;

  if (vector_create(&const_b, CONFIG_ASIZE) == -1)
    goto on_error;
  if (vector_create_empty(&b, CONFIG_ASIZE) == -1)
    goto on_error;

  if (vector_create_empty(&x, CONFIG_ASIZE) == -1)
    goto on_error;

  fsub_initialize(&fsc, const_a, b);

  /* apply forward substitution: ax = b */
  size_t i;
  for (i = 0; i < CONFIG_ITER_COUNT; ++i)
  {
    struct timeval tms[4];

    matrix_copy(a, const_a);
    vector_copy(b, const_b);

    /* ref first since fsub_apply modifies */
    gettimeofday(&tms[0], NULL);
#if CONFIG_CHECK
    fsub_apply_ref(a, x, const_b);
#else
    fsub_apply_seq(a, b);
#endif
    gettimeofday(&tms[1], NULL);
    timersub(&tms[1], &tms[0], &tms[2]);

    gettimeofday(&tms[0], NULL);
    fsub_apply(&fsc);
    gettimeofday(&tms[1], NULL);
    timersub(&tms[1], &tms[0], &tms[3]);

#if CONFIG_TIME
    printf("times: seq=%lu par=%lu\n",
	   tms[2].tv_sec * 1000000 + tms[2].tv_usec,
	   tms[3].tv_sec * 1000000 + tms[3].tv_usec);
#endif

#if CONFIG_CHECK
    if (vector_cmp(fsc.b, x))
    {
      vector_print(fsc.b);
      vector_print(x);
      exit(-1);
    }
#endif
  }

  fsub_finalize(&fsc);

 on_error:
  if (a != NULL)
    matrix_destroy(a);
  if (const_a != NULL)
    matrix_destroy(const_a);
  if (x != NULL)
    vector_destroy(x);
  if (b != NULL)
    vector_destroy(b);
  if (const_b != NULL)
    vector_destroy(const_b);

  return 0;
}
