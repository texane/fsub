#include <stdio.h>
#include <math.h>

#ifndef __USE_GNU
#define __USE_GNU
#endif
#include <sched.h>
#include <pthread.h>

#include "matrix.h"
#include "vector.h"

#include "atomic.h"
#include "spinlock.h"

#include "config.h"



/* thread pool
 */

struct fsub_context;

typedef struct thread_block
{
#define THREAD_STATE_ZERO 0
#define THREAD_STATE_STEAL 1
#define THREAD_STATE_DONE 2
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

typedef size_t index_t;

typedef struct parwork
{
  spinlock_t lock;

  /* for global synchronization */
  volatile unsigned long refn __attribute__((aligned));

  /* area to be processed: (i,0,asize,j) */
  volatile index_t i;
  volatile index_t j;

  /* number of processed rows */
  volatile size_t row_count __attribute__((aligned));

} parwork_t;


typedef struct fsub_context
{
  /* ax = b; */
  const matrix_t* volatile a;
  vector_t* volatile b;

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

static void parwork_next_rows(fsub_context_t* fsc)
{
  /* continue processing the parallel work */

  parwork_t* const w = &fsc->parwork;

  index_t i = INVALID_MATRIX_INDEX;
  index_t j = INVALID_MATRIX_INDEX;

  size_t row_count;

  parwork_lock_with_sync(w);
  const size_t asize = matrix_size(fsc->a);
  if (w->i < asize)
  {
    row_count = CONFIG_PAR_SIZE / w->j;
    if (row_count == 0)
      row_count = 1;
    else if ((w->i + row_count) > asize)
      row_count = asize - w->i;

    i = w->i;
    w->i += row_count;
    j = w->j;
  }
  parwork_unlock(w);

  /* did not get work */
  if (i == INVALID_MATRIX_INDEX)
    return ;

  /* process the row */
  const size_t saved_count = row_count;
  for (; row_count; --row_count, ++i)
    sub_row(fsc, i, j);

  /* update processed row_count */
  add_atomic_ul(&w->row_count, saved_count);
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

      /* wait until work avail */
    case THREAD_STATE_ZERO:
      break;

    case THREAD_STATE_STEAL:
    default:
      /* process next row in parallel area */
      parwork_next_rows(fsc);
      break;
    }
  }

  return NULL;
}

static void sub_rect_block
(fsub_context_t* fsc, index_t i, index_t j, size_t m, size_t n)
{
  /* b -= aij * x, a rect of width=m, height=n */

#if 1 /* a col order */

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

#else /* b order */

  const index_t saved_i = i;
  const index_t last_i = i + n;
  const index_t last_j = j + m;

  for (; j < last_j; ++j)
  {
    const double bval = *vector_const_at(fsc->b, j);

    size_t i;
    for (i = saved_i; i < last_i; ++i)
      *vector_at(fsc->b, i) -= *matrix_const_at(fsc->a, i, j) * bval;
  }

#endif

}

static void sub_tri_block(fsub_context_t* fsc, index_t i, size_t n)
{
  /* b -= aij * x, a triangle of size n */

#if 0 /* a col order */

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

#else /* b order */

  const index_t last_i = i + n - 1;

  for (size_t j = i; n > 1; ++j, --n)
  {
    const double bval = *vector_const_at(fsc->b, j);

    i = last_i;
    for (size_t k = n - 1; k; --k, --i)
      *vector_at(fsc->b, i) -= *matrix_const_at(fsc->a, i, j) * bval;
  }

#endif
}


/* exported */

static fsub_context_t global_fsc;

void fsub_pthread_apply(const matrix_t* a, vector_t* b)
{
  /* assume fsc->lsize >= matrix_size(fsc->a) */
  /* assume (fsc->lsize % fsc->ksize) == 0 */
  /* assume (matrix_size(fsc->a) % fsc->lsize) == 0 */

  fsub_context_t* const fsc = &global_fsc;

  /* fixme, should not rely on it */
  if (fsc->a == NULL)
  {
    fsc->a = a;
    fsc->b = b;
    __sync_synchronize();

    /* wake threads up */
    size_t tid;
    for (tid = 0; tid < CONFIG_THREAD_COUNT; ++tid)
      fsc->pool.tbs[tid].state = THREAD_STATE_STEAL;
  }

  const size_t llcount = matrix_size(fsc->a) / fsc->lsize;

  /* solve the first llblock */
  sub_tri_block(fsc, 0, fsc->lsize);

  /* post the band to process */
  parwork_lock(&fsc->parwork);
  fsc->parwork.i = fsc->lsize;
  fsc->parwork.j = fsc->lsize;
  fsc->parwork.row_count = fsc->lsize;
  parwork_unlock(&fsc->parwork);

  /* slide along all the kkblocks on the diagonal */
  const index_t last_i = (llcount - 1) * fsc->lsize;
  index_t i;
  for (i = fsc->lsize; i < last_i; i += fsc->ksize)
  {
    const index_t next_i = i + fsc->ksize;

    /* wait until left band processed. the loop is left
       with the lock held and no one working, thus safe
       to process the kk block and update parwork.
       synced_row_count is updated to reflect the new
       synchronized parwork.row_count.
    */
    while (1)
    {
      if ((size_t)read_atomic_ul(&fsc->parwork.row_count) >= next_i)
      {
	parwork_lock(&fsc->parwork);
	parwork_synchronize(&fsc->parwork);
	break ;
      }

      /* contribute to the parallel work */
      parwork_next_rows(fsc);
    }

    /* process the kk block sequentially */
    sub_tri_block(fsc, i, fsc->ksize);

    /* update parwork area j (which is next_i) */
    fsc->parwork.j = next_i;

    /* updating the parwork area may have created
       a hole below the kkblock. this hole dim are:
       i + fsc->ksize, i, parwork.row_count, i + fsc->ksize
       capture parwork.row_count before unlocking workers
    */
    const index_t row_count = fsc->parwork.row_count;

    parwork_unlock(&fsc->parwork);

    /* process the band sequentially */
    const size_t band_height = row_count - next_i;
    sub_rect_block(fsc, next_i, i, band_height, fsc->ksize);
  }

  /* wait until remaining left bands processed */
  const size_t asize = matrix_size(fsc->a);
  while ((index_t)read_atomic_ul(&fsc->parwork.row_count) < asize)
    parwork_next_rows(fsc);

  /* proces last llblock sequentially */
  if (llcount > 1)
    sub_tri_block(fsc, i, fsc->lsize);
}


void fsub_pthread_initialize(void)
{
  fsub_context_t* const fsc = &global_fsc;

  fsc->a = NULL;
  fsc->b = NULL;

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


void fsub_pthread_finalize(void)
{
  fsub_context_t* const fsc = &global_fsc;

  size_t tid;
  for (tid = 1; tid < CONFIG_THREAD_COUNT; ++tid)
  {
    write_atomic_ul(&fsc->pool.tbs[tid].state, THREAD_STATE_DONE);
    pthread_join(fsc->pool.tbs[tid].thread, NULL);
  }
}

