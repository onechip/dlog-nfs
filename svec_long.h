#ifndef NTLX_svec_long__H
#define NTLX_svec_long__H

#include <NTL/vec_long.h>
#include "svector.h"

NTL_svector_decl(long,vec_long,svec_long);
NTL_math_svector_decl(long,vec_long,svec_long);
NTL_eq_svector_decl(long,svec_long);
NTL_io_svector_decl(long,svec_long);

/* need to provide:
inline svec_ZZ& operator*=(svec_ZZ& x, long a) { 
  mul(x, x, a);  return x;
}
*/

#endif
