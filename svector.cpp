
#include "svec_long.h"
#include "svec_ZZ.h"
#include "svec_ZZ_p.h"

/* Implementation of sparse vectors.
 *
 * Written by: Chris Studholme
 * Copyright:  GPL (http://www.fsf.org/copyleft/gpl.html)
 */

NTL_START_IMPL;

static const long _long_zero = 0;
template<> const long& zero_ref<long>() { return _long_zero; }
template<> const ZZ& zero_ref<ZZ>() { return ZZ::zero(); }
template<> const ZZ_p& zero_ref<ZZ_p>() { return ZZ_p::zero(); }

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

NTL_math_svector_impl(long,vec_long,svec_long);
NTL_math_svector_impl(ZZ,vec_ZZ,svec_ZZ);
NTL_math_svector_impl(ZZ_p,vec_ZZ_p,svec_ZZ_p);

NTL_END_IMPL;
