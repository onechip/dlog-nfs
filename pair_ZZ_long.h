#ifndef NTL_pair_ZZ_long__H
#define NTL_pair_ZZ_long__H

#include <NTL/pair.h>
#include <NTL/vector.h>
#include <NTL/ZZ.h>

NTL_OPEN_NNS;

NTL_pair_decl(ZZ,long,pair_ZZ_long);
NTL_pair_io_decl(ZZ,long,pair_ZZ_long);
NTL_pair_eq_decl(ZZ,long,pair_ZZ_long);

NTL_vector_decl(pair_ZZ_long,vec_pair_ZZ_long);
NTL_io_vector_decl(pair_ZZ_long,vec_pair_ZZ_long);
NTL_eq_vector_decl(pair_ZZ_long,vec_pair_ZZ_long);

NTL_CLOSE_NNS;

#endif
