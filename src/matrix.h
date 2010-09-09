#ifndef MATRIX_H_INCLUDED
# define MATRIX_H_INCLUDED


#include <sys/types.h>
#include <gsl/gsl_matrix.h>


typedef gsl_matrix matrix_t;

static inline const double* __attribute__((unused))
matrix_const_at
(const matrix_t* m, size_t i, size_t j)
{
  return gsl_matrix_const_ptr(m, i, j);
}

static inline double* __attribute__((unused))
matrix_at(matrix_t* m, size_t i, size_t j)
{
  return (double*)matrix_const_at(m, i, j);
}

static inline size_t __attribute__((unused))
matrix_size(const matrix_t* m)
{
  return m->size1;
}

static inline int __attribute__((unused))
matrix_create_empty(matrix_t** m, size_t n)
{
  *m = gsl_matrix_alloc(n, n);
  if (*m == NULL)
    return -1;
  return 0;
}

static int __attribute__((unused))
matrix_create_lower(matrix_t** m, size_t n)
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

static void __attribute__((unused))
matrix_destroy(matrix_t* m)
{
  gsl_matrix_free(m);
}

static void __attribute__((unused))
matrix_copy(matrix_t* a, const matrix_t* b)
{
  gsl_matrix_memcpy(a, b);
}

static void __attribute__((unused))
matrix_print(const matrix_t* m)
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


#endif /* ! MATRIX_H_INCLUDED */
