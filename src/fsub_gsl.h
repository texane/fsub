#ifndef FSUB_GSL_INCLUDED_H
#define FSUB_GSL_INCLUDED_H


struct matrix;
struct vector;

void fsub_gsl_apply(struct matrix*, struct vector*, struct vector*, const struct vector*);


#endif /* ! FSUB_GSL_INCLUDED_H */
