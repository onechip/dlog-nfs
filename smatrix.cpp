#include "vec_svec_long.h"
#include "vec_svec_ZZ.h"
#include "vec_svec_ZZ_p.h"

#include "smat_long.h"
#include "smat_ZZ.h"
#include "smat_ZZ_p.h"


/* Implementation of sparse matrices.
 *
 * Written by: Chris Studholme
 * Copyright:  GPL (http://www.fsf.org/copyleft/gpl.html)
 */

NTL_START_IMPL;

NTL_smatrix_impl(long,svec_long,vec_svec_long,smat_long);
NTL_conv_smatrix_impl(long,svec_long,vec_svec_long,mat_long,smat_long);
NTL_math_smatrix_impl(long,vec_long,svec_long,vec_svec_long,smat_long);
NTL_eq_smatrix_impl(long,svec_long,vec_svec_long,smat_long);
NTL_io_smatrix_impl(long,svec_long,vec_svec_long,smat_long);

NTL_smatrix_impl(ZZ,svec_ZZ,vec_svec_ZZ,smat_ZZ);
NTL_conv_smatrix_impl(ZZ,svec_ZZ,vec_svec_ZZ,mat_ZZ,smat_ZZ);
NTL_math_smatrix_impl(ZZ,vec_ZZ,svec_ZZ,vec_svec_ZZ,smat_ZZ);
NTL_eq_smatrix_impl(ZZ,svec_ZZ,vec_svec_ZZ,smat_ZZ);
NTL_io_smatrix_impl(ZZ,svec_ZZ,vec_svec_ZZ,smat_ZZ);

NTL_smatrix_impl(ZZ_p,svec_ZZ_p,vec_svec_ZZ_p,smat_ZZ_p);
NTL_conv_smatrix_impl(ZZ_p,svec_ZZ_p,vec_svec_ZZ_p,mat_ZZ_p,smat_ZZ_p);
NTL_math_smatrix_impl(ZZ_p,vec_ZZ_p,svec_ZZ_p,vec_svec_ZZ_p,smat_ZZ_p);
NTL_eq_smatrix_impl(ZZ_p,svec_ZZ_p,vec_svec_ZZ_p,smat_ZZ_p);
NTL_io_smatrix_impl(ZZ_p,svec_ZZ_p,vec_svec_ZZ_p,smat_ZZ_p);

// scalar multiplication (svec_long -> svec_ZZ_p)
inline void mul(svec_ZZ_p& result, const svec_long& a, const ZZ_p& b) {
  result.SetLength(a.length());
  clear(result);
  long an = a.nvalues();
  const long* ai = a.indices();
  const svec_long::value_type* av = a.values();
  for (long i=0; i<an; ++i)
    if (av[i]!=0)
      result[ai[i]] = av[i]*b;
}

// scalar multiplication (svec_ZZ -> svec_ZZ_p)
inline void mul(svec_ZZ_p& result, const svec_ZZ& a, const ZZ_p& b) {
  result.SetLength(a.length());
  clear(result);
  long an = a.nvalues();
  const long* ai = a.indices();
  const svec_ZZ::value_type* av = a.values();
  for (long i=0; i<an; ++i)
    if (av[i]!=0) 
      result[ai[i]] = to_ZZ_p(av[i])*b;
}

// inner product (svec_long and svec_ZZ_p)
inline void InnerProduct(ZZ_p& result, const svec_long& a, const vec_ZZ_p& b) {
  if (a.length()!=b.length()) {
    std::cerr<<"InnerProduct() length mismatch"<<std::endl;
    exit(1);
  }
  clear(result);
  long an = a.nvalues();
  const long* ai = a.indices();
  const svec_long::value_type* av = a.values();
  for (long i=0; i<an; ++i)
    result += av[i]*b[ai[i]];
}
// inner product (svec_ZZ and svec_ZZ_p)
inline void InnerProduct(ZZ_p& result, const svec_ZZ& a, const vec_ZZ_p& b) {
  if (a.length()!=b.length()) {
    std::cerr<<"InnerProduct() length mismatch"<<std::endl;
    exit(1);
  }
  clear(result);
  long an = a.nvalues();
  const long* ai = a.indices();
  const svec_ZZ::value_type* av = a.values();
  for (long i=0; i<an; ++i)
    result += to_ZZ_p(av[i])*b[ai[i]];
}

// Lanczos implementations
NTL_Lanczos_impl(ZZ_p,vec_ZZ_p,svec_ZZ_p,smat_long);
NTL_Lanczos_impl(ZZ_p,vec_ZZ_p,svec_ZZ_p,smat_ZZ);
NTL_Lanczos_impl(ZZ_p,vec_ZZ_p,svec_ZZ_p,smat_ZZ_p);



NTL_END_IMPL;
