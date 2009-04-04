#ifndef NTLX_smat_long_H
#define NTLX_smat_long_H

#include <NTL/vec_ZZ_p.h>
#include "smatrix.h"
#include "mat_long.h"
#include "vec_svec_long.h"
#include "Lanczos.h"


NTL_OPEN_NNS;

NTL_smatrix_decl(long,svec_long,vec_svec_long,smat_long);
NTL_conv_smatrix_decl(long,svec_long,vec_svec_long,mat_long,smat_long);
NTL_math_smatrix_decl(long,vec_long,svec_long,vec_svec_long,smat_long);
NTL_eq_smatrix_decl(long,svec_long,vec_svec_long,smat_long);
NTL_io_smatrix_decl(long,svec_long,vec_svec_long,smat_long);

// methods for solving sparse linear systems
NTL_Lanczos_decl(vec_ZZ_p,smat_long);

/* Structured Gaussian Elimination.
 * Finds a smaller version of system Aorig*x=yorig and returns it
 * in Anew and ynew.  cols is the list of column indices from Aorig
 * that were kept in Anew.
 */
void SGauss(smat_long& Anew, vec_ZZ_p& ynew, vec_long& cols, 
	    const smat_long& Aorig, const vec_ZZ_p& yorig);

/* Undo gaussian elimination.  Assuming Aorig and yorig were reduced with
 * SGauss(), cols is the list of columns that were kept, and x is the solution
 * to the reduced system that was found, then xorig, the solution to the
 * original system, will be computed.
 */
bool SGauss_undo(vec_ZZ_p& xorig, const vec_ZZ_p& x, 
		 const vec_long& cols,
		 const smat_long& Aorig, const vec_ZZ_p& yorig);


NTL_CLOSE_NNS;

#endif
