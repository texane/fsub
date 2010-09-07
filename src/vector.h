#ifndef VECTOR_H_INCLUDED
# define VECTOR_H_INCLUDED


#include <sys/types.h>
#include <gsl/gsl_vector.h>


typedef gsl_vector vector_t;


static inline const double* vector_const_at(const vector_t* v, size_t i)
{
  return gsl_vector_const_ptr(v, i);
}

static inline double* vector_at(vector_t* v, size_t i)
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

static inline void vector_destroy(vector_t* v)
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

static inline void vector_copy(vector_t* a, const vector_t* b)
{
  gsl_vector_memcpy(a, b);
}

static int __attribute__((unused))
vector_cmp(const vector_t* a, const vector_t* b)
{
  size_t i;
  for (i = 0; i < a->size; ++i)
    if (*vector_const_at(a, i) != *vector_const_at(b, i))
      return -1;
  return 0;
}


#endif /* ! VECTOR_H_INCLUDED */
