#ifndef NTLX_mat_long_H
#define NTLX_mat_long_H

#include <NTL/vec_long.h>
#include <NTL/vec_vec_long.h>
#include <NTL/matrix.h>

NTL_OPEN_NNS;

NTL_matrix_decl(long,vec_long,vec_vec_long,mat_long);
NTL_eq_matrix_decl(long,vec_long,vec_vec_long,mat_long);
NTL_io_matrix_decl(long,vec_long,vec_vec_long,mat_long);

void clear(mat_long& A);

NTL_CLOSE_NNS;

#endif
