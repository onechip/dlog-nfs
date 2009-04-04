
#include "vec_short.h"

NTL_START_IMPL;

static inline
void BlockConstruct(short*, long) {}

static inline
void BlockDestroy(short*, long) {}

NTL_vector_impl_plain(short,vec_short);
NTL_eq_vector_impl(short,vec_short);
NTL_io_vector_impl(short,vec_short);

/* need to provide:
inline svec_ZZ& operator*=(svec_ZZ& x, long a) { 
  mul(x, x, a);  return x;
}
*/

NTL_END_IMPL;
