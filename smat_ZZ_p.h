#ifndef NTLX_smat_ZZ_p_H
#define NTLX_smat_ZZ_p_H

#include <NTL/mat_ZZ_p.h>

#include "smatrix.h"
#include "vec_svec_ZZ_p.h"
#include "Lanczos.h"
#include "smat_long.h"


NTL_smatrix_decl(ZZ_p,svec_ZZ_p,vec_svec_ZZ_p,smat_ZZ_p);
NTL_conv_smatrix_decl(ZZ_p,svec_ZZ_p,vec_svec_ZZ_p,mat_ZZ_p,smat_ZZ_p);
NTL_math_smatrix_decl(ZZ_p,vec_ZZ_p,svec_ZZ_p,vec_svec_ZZ_p,smat_ZZ_p);
NTL_eq_smatrix_decl(ZZ_p,svec_ZZ_p,vec_svec_ZZ_p,smat_ZZ_p);
NTL_io_smatrix_decl(ZZ_p,svec_ZZ_p,vec_svec_ZZ_p,smat_ZZ_p);


// Lanczos method for solving sparse linear systems
NTL_Lanczos_decl(vec_ZZ_p,smat_ZZ_p);
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


#endif
