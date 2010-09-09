#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/time.h>

#include "config.h"

#include "fsub_gsl.h"
#include "fsub_seq.h"
#include "fsub_pthread.h"
#include "fsub_xkaapi.h"

#include "atomic.h"
#include "spinlock.h"

#include "vector.h"
#include "matrix.h"


/* main */

int main(int ac, char** av)
{
  int error = 0;

  matrix_t* a = NULL;
  matrix_t* const_a = NULL;
  vector_t* x = NULL;
  vector_t* const_b = NULL;
  vector_t* b = NULL;

#if CONFIG_FSUB_SEQ
  vector_t* bseq = NULL;
#endif

  if (matrix_create_lower(&const_a, CONFIG_ASIZE) == -1)
    goto on_error;
  if (matrix_create_empty(&a, CONFIG_ASIZE) == -1)
    goto on_error;

  if (vector_create(&const_b, CONFIG_ASIZE) == -1)
    goto on_error;
  if (vector_create_empty(&b, CONFIG_ASIZE) == -1)
    goto on_error;
#if CONFIG_FSUB_SEQ
  if (vector_create_empty(&bseq, CONFIG_ASIZE) == -1)
    goto on_error;
#endif

  if (vector_create_empty(&x, CONFIG_ASIZE) == -1)
    goto on_error;

#if CONFIG_FSUB_PTHREAD
  fsub_pthread_initialize();
#endif

  /* apply forward substitution: ax = b */
  size_t i;
  for (i = 0; i < CONFIG_ITER_COUNT; ++i)
  {
    struct timeval tms[4];

    matrix_copy(a, const_a);
    vector_copy(b, const_b);
#if CONFIG_FSUB_SEQ
    vector_copy(bseq, const_b);
#endif

#if CONFIG_FSUB_PTHREAD
    gettimeofday(&tms[0], NULL);
    fsub_pthread_apply(const_a, b);
    gettimeofday(&tms[1], NULL);
    timersub(&tms[1], &tms[0], &tms[2]);
#elif CONFIG_FSUB_XKAAPI
    gettimeofday(&tms[0], NULL);
    fsub_xkaapi_apply(const_a, b);
    gettimeofday(&tms[1], NULL);
    timersub(&tms[1], &tms[0], &tms[2]);
#endif

#if CONFIG_FSUB_GSL
    fsub_gsl_apply(a, x, const_b);
#endif

#if CONFIG_FSUB_SEQ
    gettimeofday(&tms[0], NULL);
    fsub_seq_apply(a, bseq);
    gettimeofday(&tms[1], NULL);
    timersub(&tms[1], &tms[0], &tms[3]);
#endif

    /* report times */
#if CONFIG_TIME
#if CONFIG_FSUB_PTHREAD || CONFIG_FSUB_XKAAPI
    printf("par: %lu\n", tms[2].tv_sec * 1000000 + tms[2].tv_usec);
#endif
#if CONFIG_FSUB_SEQ
    printf("seq: %lu\n", tms[3].tv_sec * 1000000 + tms[3].tv_usec);
#endif
#endif

    /* check against gsl implm */
#if CONFIG_FSUB_GSL
#if CONFIG_FSUB_SEQ
    if (vector_cmp(bseq, x))
    {
      vector_print(b);
      vector_print(x);
      printf("invalid seq\n");
      error = -1;
      goto on_error;
    }
#endif
#if CONFIG_FSUB_PTHREAD || CONFIG_FSUB_XKAAPI
    if (vector_cmp(b, x))
    {
      vector_print(b);
      vector_print(x);
      printf("invalid parallel\n");
      error = -1;
      goto on_error;
    }
#endif
#endif
  }

#if CONFIG_FSUB_PTHREAD
  fsub_pthread_finalize();
#endif

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
#if CONFIG_SEQ
  if (bseq != NULL)
    vector_destroy(bseq);
#endif

  return error;
}
