#ifndef NTLX_smat_long_H
#define NTLX_smat_long_H

#include "mat_long.h"

#include "smatrix.h"
#include "vec_svec_long.h"

NTL_smatrix_decl(long,svec_long,vec_svec_long,smat_long);
NTL_conv_smatrix_decl(long,svec_long,vec_svec_long,mat_long,smat_long);
NTL_math_smatrix_decl(long,vec_long,svec_long,vec_svec_long,smat_long);
NTL_eq_smatrix_decl(long,svec_long,vec_svec_long,smat_long);
NTL_io_smatrix_decl(long,svec_long,vec_svec_long,smat_long);

#endif
