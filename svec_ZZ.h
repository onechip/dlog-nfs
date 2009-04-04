#ifndef NTLX_svec_ZZ__H
#define NTLX_svec_ZZ__H

#include <NTL/vec_ZZ.h>
#include <NTL/vec_long.h>
#include "svector.h"

NTL_OPEN_NNS;

NTL_svector_decl(ZZ,vec_ZZ,svec_ZZ,long,vec_long);
NTL_math_svector_decl(ZZ,vec_ZZ,svec_ZZ);
NTL_eq_svector_decl(ZZ,svec_ZZ);
NTL_io_svector_decl(ZZ,svec_ZZ);

NTL_CLOSE_NNS;

#endif
