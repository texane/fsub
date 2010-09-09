#ifndef FSUB_PTHREAD_H_INCLUDED
# define FSUB_PTHREAD_H_INCLUDED

#include "matrix.h"
#include "vector.h"

void fsub_pthread_initialize(void);
void fsub_pthread_finalize(void);
void fsub_pthread_apply(const matrix_t*, vector_t*);

#endif /* ! FSUB_PTHREAD_H_INCLUDED */
