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


NTL_vector_impl(svec_long,vec_svec_long);
NTL_eq_vector_impl(svec_long,vec_svec_long);
NTL_io_vector_impl(svec_long,vec_svec_long);

NTL_vector_impl(svec_ZZ,vec_svec_ZZ);
NTL_eq_vector_impl(svec_ZZ,vec_svec_ZZ);
NTL_io_vector_impl(svec_ZZ,vec_svec_ZZ);

NTL_vector_impl(svec_ZZ_p,vec_svec_ZZ_p);
NTL_eq_vector_impl(svec_ZZ_p,vec_svec_ZZ_p);
NTL_io_vector_impl(svec_ZZ_p,vec_svec_ZZ_p);

inline long IsZero(long i) {
  return (i==0);
}

inline void clear(vec_long& v) {
  memset(v.elts(),0,sizeof(long)*v.length());
}

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
  const long* av = a.values();
  for (long i=0; i<an; ++i)
    if (av[i]!=0)
      result[ai[i]] = av[i]*b;
}

// inner product (svec_long and svec_ZZ_p)
inline void InnerProduct(ZZ_p& result, const svec_long& a, const vec_ZZ_p& b) {
  if (a.length()!=b.length()) {
    cerr<<"InnerProduct() length mismatch\n";
    exit(1);
  }
  clear(result);
  long an = a.nvalues();
  const long* ai = a.indices();
  const long* av = a.values();
  for (long i=0; i<an; ++i)
    result += av[i]*b[ai[i]];
}

// Lanczos implementations
NTL_Lanczos_impl(ZZ_p,vec_ZZ_p,svec_ZZ_p,smat_ZZ_p);
NTL_Lanczos_impl(ZZ_p,vec_ZZ_p,svec_ZZ_p,smat_long);

/* Structured Gaussian Elimination.
 * Finds a smaller version of system Aorig*x=yorig and returns it
 * in Anew and ynew.  cols is the list of column indices from Aorig
 * that were kept in Anew.
 */
void SGauss(smat_long& Anew, vec_ZZ_p& ynew, vec_long& cols, 
	    const smat_long& Aorig, const vec_ZZ_p& yorig) {
  cols.SetLength(Aorig.NumCols());
  for (long i=0; i<cols.length(); ++i)
    cols[i]=i;
  Anew=Aorig;
  ynew=yorig;
}

/* Undo gaussian elimination.  Assuming Aorig and yorig were reduced with
 * SGauss(), cols is the list of columns that were kept, and x is the solution
 * to the reduced system that was found, then xorig, the solution to the
 * original system, will be computed.
 */
bool SGauss_undo(vec_ZZ_p& xorig, const vec_ZZ_p& x, 
		 const vec_long& cols,
		 const smat_long& Aorig, const vec_ZZ_p& yorig) {
  // copy known results to xorig
  xorig.SetLength(Aorig.NumCols());
  bool* known = new bool[xorig.length()];
  memset(known,0,xorig.length()*sizeof(bool));
  for (long i=0; i<cols.length(); ++i) {
    xorig[cols[i]] = x[i];
    known[cols[i]] = true;
  }

  // keep track of which rows we have used
  bool* used = new bool[Aorig.NumRows()];
  memset(used,0,Aorig.NumRows()*sizeof(bool));

  // attempt to compute unknown values
  bool done;
  do {
    done=true;
    for (long i=0; i<Aorig.NumRows(); ++i) 
      if (!used[i]) {
	// see if we have exactly one unknown value with non-zero coefficient
	long unknown=-1;
	long n = Aorig[i].nvalues();
	const long* aj = Aorig[i].indices();
	const long* av = Aorig[i].values();
	used[i]=true;
	for (long j=0; j<n; ++j) 
	  if ((av[j]!=0)&&!known[aj[j]]) {
	    if (unknown<0)
	      unknown=aj[j];
	    else {
	      // second one (not useless afterall)
	      used[i]=false;
	      unknown=-1;
	      break;
	    }
	  }
	if (unknown>=0) {
	  // compute unknown value
	  ZZ_p v;
	  InnerProduct(v,Aorig[i],xorig);
	  div(xorig[unknown],v,to_ZZ_p(Aorig[i][unknown]));
	  negate(xorig[unknown],xorig[unknown]);
	  known[unknown] = true;
	  used[i]=true;
	  done=false;
	}
      }
  } while (!done);

  delete[] used;

  // make sure we know them all
  for (long i=0; i<xorig.length(); ++i)
    if (!known[i]) {
      delete[] known;
      cerr<<"SGauss_used() failed with variable "<<i<<"\n";
      return false;
    }

  delete[] known;
  return true;
}
