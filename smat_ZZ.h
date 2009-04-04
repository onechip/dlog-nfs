#ifndef NTLX_smat_ZZ_H
#define NTLX_smat_ZZ_H

#include <NTL/mat_ZZ.h>
#include <NTL/vec_ZZ_p.h>
#include "smatrix.h"
#include "vec_svec_ZZ.h"
#include "Lanczos.h"


NTL_OPEN_NNS;

NTL_smatrix_decl(ZZ,svec_ZZ,vec_svec_ZZ,smat_ZZ);
NTL_conv_smatrix_decl(ZZ,svec_ZZ,vec_svec_ZZ,mat_ZZ,smat_ZZ);
NTL_math_smatrix_decl(ZZ,vec_ZZ,svec_ZZ,vec_svec_ZZ,smat_ZZ);
NTL_eq_smatrix_decl(ZZ,svec_ZZ,vec_svec_ZZ,smat_ZZ);
NTL_io_smatrix_decl(ZZ,svec_ZZ,vec_svec_ZZ,smat_ZZ);

// methods for solving sparse linear systems
NTL_Lanczos_decl(vec_ZZ_p,smat_ZZ);

/* Structured Gaussian Elimination.
 * Finds a smaller version of system Aorig*x=yorig and returns it
 * in Anew and ynew.  cols is the list of column indices from Aorig
 * that were kept in Anew.
 */
void SGauss(smat_ZZ& Anew, vec_ZZ_p& ynew, vec_long& cols, 
	    const smat_ZZ& Aorig, const vec_ZZ_p& yorig);

/* Undo gaussian elimination.  Assuming Aorig and yorig were reduced with
 * SGauss(), cols is the list of columns that were kept, and x is the solution
 * to the reduced system that was found, then xorig, the solution to the
 * original system, will be computed.
 */
bool SGauss_undo(vec_ZZ_p& xorig, const vec_ZZ_p& x, 
		 const vec_long& cols,
		 const smat_ZZ& Aorig, const vec_ZZ_p& yorig);


NTL_CLOSE_NNS;

#endif
