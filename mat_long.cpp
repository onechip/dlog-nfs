
#include "mat_long.h"

/* Implementation of a matrix of long values.
 *
 * Written by: Chris Studholme
 * Copyright:  GPL (http://www.fsf.org/copyleft/gpl.html)
 */

NTL_START_IMPL;

void clear(mat_long& A) {
  for (long i=0; i<A.NumRows(); ++i)
    clear(A[i]);
}

NTL_END_IMPL;
