
#include "vec_short.h"

NTL_START_IMPL;

static inline
void BlockConstruct(short*, long) {}

static inline
void BlockDestroy(short*, long) {}

/* need to provide:
inline svec_ZZ& operator*=(svec_ZZ& x, long a) { 
  mul(x, x, a);  return x;
}
*/

NTL_END_IMPL;
