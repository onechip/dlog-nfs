
#include "mat_long.h"

/* Implementation of a matrix of long values.
 *
 * Written by: Chris Studholme
 * Copyright:  GPL (http://www.fsf.org/copyleft/gpl.html)
 */

NTL_START_IMPL;

NTL_matrix_impl(long,vec_long,vec_vec_long,mat_long);

NTL_eq_matrix_impl(long,vec_long,vec_vec_long,mat_long);

NTL_io_matrix_impl(long,vec_long,vec_vec_long,mat_long);

inline void clear(vec_long& a) {
  for (long i=0; i<a.length(); ++i)
    a[i]=0;
}

void clear(mat_long& A) {
  for (long i=0; i<A.NumRows(); ++i)
    clear(A[i]);
}

NTL_END_IMPL;
