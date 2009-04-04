#include "svec_long.h"
#include "svec_ZZ.h"
#include "svec_ZZ_p.h"


/* Implementation of sparse vectors.
 *
 * Written by: Chris Studholme
 * Copyright:  GPL (http://www.fsf.org/copyleft/gpl.html)
 */

NTL_START_IMPL;

static const long ntlx_svec_long_zero=0;

inline bool IsZero(long i) {
  return (i==0);
}

inline void clear(long &i) {
  i=0;
}

inline void add(long&x, long a, long b) {
  x=a+b;
}

inline void sub(long&x, long a, long b) {
  x=a-b;
}

inline void mul(long&x, long a, long b) {
  x=a*b;
}

inline void negate(long& x, long a) {
  x=-a;
}

inline void clear(vec_long& v) {
  for (long i=v.length()-1; i>=0; --i)
    v[i]=0;
}

NTL_svector_impl(long,vec_long,svec_long,ntlx_svec_long_zero);
NTL_math_svector_impl(long,vec_long,svec_long);
NTL_eq_svector_impl(long,svec_long);
NTL_io_svector_impl(long,svec_long);

NTL_svector_impl(ZZ,vec_ZZ,svec_ZZ,ZZ::zero());
NTL_math_svector_impl(ZZ,vec_ZZ,svec_ZZ);
NTL_eq_svector_impl(ZZ,svec_ZZ);
NTL_io_svector_impl(ZZ,svec_ZZ);

NTL_svector_impl(ZZ_p,vec_ZZ_p,svec_ZZ_p,ZZ_p::zero());
NTL_math_svector_impl(ZZ_p,vec_ZZ_p,svec_ZZ_p);
NTL_eq_svector_impl(ZZ_p,svec_ZZ_p);
NTL_io_svector_impl(ZZ_p,svec_ZZ_p);

NTL_END_IMPL;
