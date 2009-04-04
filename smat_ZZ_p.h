#ifndef NTLX_smat_ZZ_p_H
#define NTLX_smat_ZZ_p_H

#include <NTL/mat_ZZ_p.h>

#include "smatrix.h"
#include "vec_svec_ZZ_p.h"
#include "Lanczos.h"
#include "smat_long.h"

NTL_OPEN_NNS;

NTL_smatrix_decl(ZZ_p,svec_ZZ_p,vec_svec_ZZ_p,smat_ZZ_p);
NTL_conv_smatrix_decl(ZZ_p,svec_ZZ_p,vec_svec_ZZ_p,mat_ZZ_p,smat_ZZ_p);
NTL_math_smatrix_decl(ZZ_p,vec_ZZ_p,svec_ZZ_p,vec_svec_ZZ_p,smat_ZZ_p);
NTL_eq_smatrix_decl(ZZ_p,svec_ZZ_p,vec_svec_ZZ_p,smat_ZZ_p);
NTL_io_smatrix_decl(ZZ_p,svec_ZZ_p,vec_svec_ZZ_p,smat_ZZ_p);

// methods for solving sparse linear systems
NTL_Lanczos_decl(vec_ZZ_p,smat_ZZ_p);

NTL_CLOSE_NNS;

#endif
