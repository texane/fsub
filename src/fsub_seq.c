#include <sys/types.h>
#include "matrix.h"
#include "vector.h"

typedef size_t index_t;

void fsub_seq_apply(matrix_t* a, vector_t* b)
{
  /* b -= aij * x, a triangle of size n */

  const index_t last_i = matrix_size(a) - 1;

  /* dont do on the last element */
  size_t n = matrix_size(a);

  for (size_t j = 0; j < last_i; ++j, --n)
  {
    const double bval = *vector_const_at(b, j);

    index_t i = last_i;
    for (size_t k = n - 1; k; --k, --i)
      *vector_at(b, i) -= *matrix_const_at(a, i, j) * bval;
  }
}
