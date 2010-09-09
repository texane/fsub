#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include "matrix.h"
#include "vector.h"


void fsub_gsl_apply(matrix_t* a, vector_t* x, const vector_t* b)
{
  gsl_permutation* const p = gsl_permutation_alloc(vector_size(b));
  if (p == NULL)
    return ;

  int signum;
  gsl_linalg_LU_decomp(a, p, &signum);
  gsl_linalg_LU_solve(a, p, b, x);

  gsl_permutation_free(p);
}
