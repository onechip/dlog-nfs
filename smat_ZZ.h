#ifndef NTLX_smat_ZZ_H
#define NTLX_smat_ZZ_H

#include <NTL/mat_ZZ.h>

#include "smatrix.h"
#include "vec_svec_ZZ.h"

NTL_smatrix_decl(ZZ,svec_ZZ,vec_svec_ZZ,smat_ZZ);
NTL_conv_smatrix_decl(ZZ,svec_ZZ,vec_svec_ZZ,mat_ZZ,smat_ZZ);
NTL_math_smatrix_decl(ZZ,vec_ZZ,svec_ZZ,vec_svec_ZZ,smat_ZZ);
NTL_eq_smatrix_decl(ZZ,svec_ZZ,vec_svec_ZZ,smat_ZZ);
NTL_io_smatrix_decl(ZZ,svec_ZZ,vec_svec_ZZ,smat_ZZ);

#endif
