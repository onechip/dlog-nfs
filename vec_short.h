#ifndef NTLX_vec_short__H
#define NTLX_vec_short__H

#include <string.h>

#include <NTL/vector.h>

NTL_OPEN_NNS;

typedef Vec<short> vec_short;

inline void clear(vec_short& x) {
  memset(x.elts(),0,x.length()*sizeof(short));
}

/* need to provide:
inline svec_ZZ& operator*=(svec_ZZ& x, long a) { 
  mul(x, x, a);  return x;
}
*/

NTL_CLOSE_NNS;

#endif
