#include <sys/types.h>
#include "kaapi.h"
#include "vector.h"
#include "matrix.h"
#include "spinlock.h"
#include "atomic.h"
#include "config.h"


/* forward substitution context */

typedef size_t index_t;

#define INVALID_MATRIX_INDEX ((index_t)-1)
#define INVALID_MATRIX_SIZE ((size_t)-1)

typedef struct fsub_context
{
  /* ax = b; */
  const matrix_t* a;
  vector_t* b;

  /* for b concurrent updates */
  spinlock_t b_lock;

  /* size of a kxk block */
  size_t ksize;

  /* size of a lxl block */
  size_t lsize;

} fsub_context_t;

static void fsub_initialize
(fsub_context_t* fsc, const matrix_t* a, vector_t* b)
{
  fsc->a = a;
  fsc->b = b;

  fsc->ksize = CONFIG_KSIZE;
  fsc->lsize = CONFIG_LSIZE;

  spinlock_init(&fsc->b_lock);
}


/* represent a block to be processed */
typedef struct block
{
  index_t i;
  index_t j;
  size_t m; /* width */
  size_t n;
} block_t;

static inline block_t* make_block
(block_t* b, index_t i, index_t j, size_t m, size_t n)
{
  b->i = i;
  b->j = j;
  b->m = m;
  b->n = n;
  return b;
}

/* xkaapi adaptive algorithm */
typedef struct task_context
{
  /* block to process in parallel */
  spinlock_t lock;
  volatile block_t block;

  fsub_context_t* fsc;

  /* xkaapi related */
  kaapi_stealcontext_t* sc;
  kaapi_stealcontext_t* msc; /* master sc */
  kaapi_taskadaptive_result_t* ktr;

} task_context_t;

typedef block_t task_result_t;

static void process_seq_tri(task_context_t* tc, index_t i, size_t n)
{
  /* sequentially process a triangle */

  fsub_context_t* const fsc = tc->fsc;

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

static void process_seq_block_common
(task_context_t* tc, index_t i, index_t j, size_t m, size_t n, int do_lock)
{
  const matrix_t* const a = tc->fsc->a;
  vector_t* const b = tc->fsc->b;

  const index_t saved_j = j;

  const index_t last_i = i + n;
  const index_t last_j = j + m;

  for (; i < last_i; ++i)
  {
    double sum = 0.f;
    for (j = saved_j; j < last_j; ++j)
      sum += *matrix_const_at(a, i, j) * *vector_at(b, j);

    if (do_lock)
    {
      spinlock_lock(&tc->fsc->b_lock);
      *vector_at(b, i) -= sum;
      spinlock_unlock(&tc->fsc->b_lock);
    }
    else
    {
      *vector_at(b, i) -= sum;
    }
  }
}

static inline void process_seq_block
(task_context_t* tc, index_t i, index_t j, size_t m, size_t n)
{
  /* sequentially process the block */
  process_seq_block_common(tc, i, j, m, n, 0);
}

static inline void process_seq_block_with_lock
(task_context_t* tc, index_t i, index_t j, size_t m, size_t n)
{
  /* process the block. since b may be updated in
     parallel we accumulate in a local variable
   */
  process_seq_block_common(tc, i, j, m, n, 1);
}

static int extract_seq_block(task_context_t* tc, block_t* sb)
{
  /* extract the next seq block SEQ_GRAIN rows at a time */

  sb->n = CONFIG_SEQ_GRAIN;

  spinlock_lock(&tc->lock);

  if (sb->n > tc->block.n)
  {
    sb->n = tc->block.n;
    if (sb->n == 0)
      goto on_error;
  }

  sb->i = tc->block.i;
  sb->j = tc->block.j;
  sb->m = tc->block.m;

  tc->block.i += sb->n;
  tc->block.n -= sb->n;

 on_error:
  spinlock_unlock(&tc->lock);

  return (sb->n == 0) ? -1 : 0;
}

static void thief_entry(void*, kaapi_thread_t*);

static int split_par_block
(kaapi_stealcontext_t* sc, int nreq, kaapi_request_t* reqs, void* arg)
{
  int nrep = 0;

  /* victim task context */
  task_context_t* const vtc = (task_context_t*)arg;

  /* extract block to distribute. each task receive unit_size rows. */
  spinlock_lock(&vtc->lock);

  if (vtc->block.n <= CONFIG_PAR_GRAIN)
  {
    spinlock_unlock(&vtc->lock);
    goto on_error;
  }

  size_t unit_size = vtc->block.n / (nreq + 1);
  if (unit_size < CONFIG_PAR_GRAIN)
  {
    unit_size = CONFIG_PAR_GRAIN;
    nreq = (vtc->block.n / CONFIG_PAR_GRAIN) - 1;
  }

  block_t stolen_block;
  stolen_block.n = unit_size * nreq;
  stolen_block.i = vtc->block.i + vtc->block.n - stolen_block.n;
  stolen_block.j = vtc->block.j;
  stolen_block.m = vtc->block.m;

  /* update the victim block */
  vtc->block.n -= stolen_block.n;
  
  spinlock_unlock(&vtc->lock);

  for (; nreq; ++reqs)
  {
    if (!kaapi_request_ok(reqs))
      continue ;

    /* push the remote task */
    kaapi_thread_t* const thief_thread = kaapi_request_getthread(reqs);
    kaapi_task_t* const thief_task = kaapi_thread_toptask(thief_thread);

    kaapi_taskadaptive_result_t* const ktr =
      kaapi_allocate_thief_result(sc, sizeof(task_result_t), NULL);
    ((task_result_t*)ktr->data)->n = 0;

    task_context_t* const thief_tc = kaapi_thread_pushdata_align
      (thief_thread, sizeof(task_context_t), 8);
    spinlock_init(&thief_tc->lock);
    thief_tc->block.i = stolen_block.i;
    thief_tc->block.j = stolen_block.j;
    thief_tc->block.m = stolen_block.m;
    thief_tc->block.n = unit_size;
    thief_tc->fsc = vtc->fsc;
    thief_tc->sc = NULL;
    thief_tc->msc = vtc->sc;
    thief_tc->ktr = ktr;

    kaapi_task_init(thief_task, thief_entry, thief_tc);
    kaapi_thread_pushtask(thief_thread);
    kaapi_request_reply_head(vtc->sc, reqs, ktr);

    /* next stolen block range */
    stolen_block.i += unit_size;
    stolen_block.n -= unit_size;

    ++nrep;
    --nreq;
  }

 on_error:
  return nrep;
}

static int reduce_thief
(kaapi_stealcontext_t* sc, void* targ, void* tptr, size_t tsize, void* vptr)
{
  task_context_t* const vtc = vptr;
  task_result_t* const tr = tptr;

  /* conccurrent with splitter */
  spinlock_lock(&vtc->lock);

  /* retrieve work */
  vtc->block.i = tr->i;
  vtc->block.n = tr->n;

  spinlock_unlock(&vtc->lock);

  return 0;
}

static void set_par_block
(task_context_t* tc, index_t i, index_t j, size_t m, size_t n)
{
  /* set the block to process in parallel */

  spinlock_lock(&tc->lock);
  tc->block.i = i;
  tc->block.j = j;
  tc->block.m = m;
  tc->block.n = n;
  spinlock_unlock(&tc->lock);
}

static void sync_par_block(task_context_t* tc)
{
  /* wait until current parallel block is done */

  while (1)
  {
    block_t sb;

    while (extract_seq_block(tc, &sb) != -1)
      process_seq_block_with_lock(tc, sb.i, sb.j, sb.m, sb.n);

    kaapi_taskadaptive_result_t* const ktr = kaapi_get_thief_head(tc->sc);
    if (ktr == NULL)
      break ;

    kaapi_preempt_thief(tc->sc, ktr, NULL, reduce_thief, tc);
  }
}

/* thief task entry */
static void thief_entry(void* arg, kaapi_thread_t* thread)
{
  task_context_t* const tc = (task_context_t*)arg;
  block_t sb;

  tc->sc = kaapi_thread_pushstealcontext
    (thread, KAAPI_STEALCONTEXT_DEFAULT, split_par_block, arg, tc->msc);

  /* until all work is done and no more thief */
  while (1)
  {
    while (extract_seq_block(tc, &sb) != -1)
    {
      process_seq_block_with_lock(tc, sb.i, sb.j, sb.m, sb.n);

      /* update ktr and look for preemption */
      const int is_preempted = kaapi_preemptpoint
	(tc->ktr, tc->sc, NULL, NULL, (void*)&tc->block,
	 sizeof(task_result_t), NULL);
      if (is_preempted)
	goto on_finalize;
    }

    /* update ktr */
    ((task_result_t*)tc->ktr->data)->n = 0;

    /* retrieve thieves */
    kaapi_taskadaptive_result_t* const ktr = kaapi_get_thief_head(tc->sc);
    if (ktr == NULL)
      goto on_finalize;

    kaapi_preempt_thief(tc->sc, ktr, NULL, reduce_thief, tc);
  }

 on_finalize:
  kaapi_steal_finalize(tc->sc);
}


/* master task entry */
static void master_entry(void* arg, kaapi_thread_t* thread)
{
  task_context_t* const tc = (task_context_t*)arg;
  fsub_context_t* const fsc = tc->fsc;

  const size_t asize = matrix_size(tc->fsc->a);
  const size_t llcount = asize / fsc->lsize;

  tc->sc = kaapi_thread_pushstealcontext
    (thread, KAAPI_STEALCONTEXT_DEFAULT, split_par_block, arg, tc->msc);

  /* process llblock(0) */
  process_seq_tri(tc, 0, fsc->lsize);

  /* foreach kkblock */
  const index_t last_i = (llcount - 1) * fsc->lsize;
  index_t i;
  for (i = fsc->lsize; i < last_i; i += fsc->ksize)
  {
    /* post left block */
    set_par_block(tc, i, 0, fsc->lsize, fsc->ksize);
    sync_par_block(tc);

    /* process the kkblock */
    process_seq_tri(tc, i, fsc->ksize);

    /* update work with band */
    const index_t next_par_i = i + 2 * fsc->ksize;
    set_par_block(tc, next_par_i, i, fsc->ksize, asize - next_par_i);

    /* process kkblock below the current one */
    process_seq_block(tc, i + fsc->ksize, i, fsc->ksize, fsc->ksize);

    /* wait until hband processed */
    sync_par_block(tc);
  }

  /* process remaining left block */
  process_seq_block(tc, i, 0, i - fsc->ksize, asize - i);

  /* process llblock(llcount - 1) */
  if (llcount > 1)
    process_seq_tri(tc, i, fsc->lsize);

  kaapi_steal_finalize(tc->sc);
}


void fsub_xkaapi_apply(const matrix_t* a, vector_t* b)
{
  /* forward sub context */
  fsub_context_t fsc;
  fsub_initialize(&fsc, a, b);

  /* create master sequential task */
  kaapi_thread_t* thread;
  kaapi_task_t* task;
  kaapi_frame_t frame;

  task_context_t tc;
  spinlock_init(&tc.lock);
  tc.block.n = 0;
  tc.fsc = &fsc;
  tc.ktr = NULL;
  tc.msc = NULL;

  thread = kaapi_self_thread();
  kaapi_thread_save_frame(thread, &frame);
  task = kaapi_thread_toptask(thread);
  kaapi_task_init(task, master_entry, (void*)&tc);
  kaapi_thread_pushtask(thread);
  kaapi_sched_sync();
  kaapi_thread_restore_frame(thread, &frame);
}
