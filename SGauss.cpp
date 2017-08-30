
#include <NTL/vec_vec_long.h>
#include "svec_long.h"

#include <NTL/vec_ZZ_p.h>
#include "smat_long.h"
#include "smat_ZZ.h"


/* Implementation of structured gaussian elimination.
 *
 * Written by: Chris Studholme
 * Copyright:  GPL (http://www.fsf.org/copyleft/gpl.html)
 */

NTL_START_IMPL;

static const bool SGAUSS_VERBOSE = false;

inline const ZZ_p& operator/=(ZZ_p& a, const NTL::ZZ& b) {
  a /= to_ZZ_p(b);
  return a;
}
inline ZZ_p operator*(const ZZ& a, const ZZ_p& b) {
  return to_ZZ_p(a)*b;
}
inline ZZ_p operator/(const ZZ_p& a, const ZZ& b) {
  return a/to_ZZ_p(b);
}


// inner product (svec_long and svec_ZZ_p)
template <class svec_T, class vec_T>
inline void InnerProduct(ZZ_p& result, const svec_T& a, const vec_T& b) {
  if (a.length()!=b.length()) {
    std::cerr<<"InnerProduct() length mismatch"<<std::endl;
    exit(1);
  }
  clear(result);
  const long* ai = a.indices();
  const typename svec_T::value_type* av = a.values();
  for (long i=0; i<a.nvalues(); ++i)
    result += av[i]*b[ai[i]];
}


/* create vector where x[i]==1 if and only if !IsZero(y[i]) */
template<class svec_T>
void nonzero(svec_long& x, const svec_T& y) {
  x.SetLength(y.length());
  x.SetAlloc(y.nvalues());
  clear(x);
  const long* yi = y.indices();
  const typename svec_T::value_type* yv = y.values();
  for (long j=0; j<y.nvalues(); ++j)
    if (!IsZero(yv[j]))
      x[yi[j]]=1;
}


/**************** Set  ****************/

class set_long {
public:
  static const long BITS_PER_LONG = 8*sizeof(unsigned long);

  set_long();
  set_long(long min, long max);
  ~set_long();

  void clear();
  void insert(long i);

  bool contains(long i) const;

private:
  unsigned long* pos;
  long npos;  // number of long's allocated
  unsigned long* neg;
  long nneg;  // number of long's allocated
};

set_long::set_long()
  : pos(NULL), npos(0), neg(NULL), nneg(0) {
}
set_long::set_long(long min, long max)
  : pos(NULL), npos(0), neg(NULL), nneg(0) {
}
set_long::~set_long() {
  free(pos);
  free(neg);
}

void set_long::clear() {
  if (pos)
    memset(pos,0,npos*sizeof(unsigned long));
  if (neg)
    memset(neg,0,nneg*sizeof(unsigned long));
}

void set_long::insert(long i) {
  if (i>=0) {
    if (i>=npos*BITS_PER_LONG) {
      long new_size = (i/BITS_PER_LONG)+1;
      pos = (unsigned long*)realloc(pos,new_size*sizeof(unsigned long));
      memset(pos+npos,0,(new_size-npos)*sizeof(unsigned long));
      npos = new_size;
    }
    pos[i/BITS_PER_LONG] |= 1<<(i%BITS_PER_LONG);
  }
  else {
    i = -(i+1);
    if (i>=nneg*BITS_PER_LONG) {
      long new_size = (i/BITS_PER_LONG)+1;
      neg = (unsigned long*)realloc(neg,new_size*sizeof(unsigned long));
      memset(neg+nneg,0,(new_size-nneg)*sizeof(unsigned long));
      nneg = new_size;
    }
    neg[i/BITS_PER_LONG] |= 1<<(i%BITS_PER_LONG);
  }
}

bool set_long::contains(long i) const {
  if (i>=0) {
    if (i<npos*BITS_PER_LONG)
      return (pos[i/BITS_PER_LONG]>>(i%BITS_PER_LONG))&1;
  }
  else {
    i = -(i+1);
    if (i<nneg*BITS_PER_LONG)
      return (neg[i/BITS_PER_LONG]>>(i%BITS_PER_LONG))&1;
  }
  return false;
}


/**************** structured gaussian elimination  ****************/

template <class smat_T, class svec_T, class vec_Tp, class Tp>
class smat_gauss {
public:

  // types
  //typedef smat_T smat_t;
  //typedef svec_T svec_t;
  typedef typename svec_T::value_type value_t;

  // matrix and column representing system A * x = y
  smat_gauss(const smat_T& A, const vec_Tp& y);
  ~smat_gauss();

  inline long NumRows() const {
    return nrows;
  }
  inline long NumCols() const {
    return A.NumCols()-first_col;
  }

  inline const svec_T& row(long i) const {
    return i>=0 ? A[i] : E[-i-1];
  }

  inline Tp& y(long i) {
    return i>=0 ? Ay[i] : Ey[-i-1];
  }

  inline long row_weight(long i) const {
    return i>=0 ? Arw[i] : Erw[-i-1];
  }

  inline void set_row_weight(long i, long w=0) {
    if (i>=0) Arw[i]=w; else Erw[-i-1]=w;
  }

  // bubble sort AEcs to ensure column c is correctly positioned
  // note: this method doesn't handle the case where AEcw[c] increased
  void resort(long c);

  // eliminate column with weight 1
  void eliminate_col1(long c);
  
  // eliminate column with weight 2
  void eliminate_col2(long c);
  
  // eliminate row with weight 1 and the column that has the non-zero entry
  void eliminate_row1(long r);
  
  // eliminate row
  void eliminate_row(long r);

  // eliminate up to max rows, starting with ones touching column first_col
  void eliminate_rows(long max=1);

  // reduce matrix to (N+e) x N
  void reduce(long extra_rows=0);

  // create new system A' * x' = y' (cols is needed to undo)
  void make_system(smat_T& Aprime, vec_Tp& yprime, vec_long& cols) const;

  // compute solution x to A*x=y, given solution x' to A'*x'=y'
  // returns true if all values of x were recovered
  static bool undo(vec_Tp& x, const vec_Tp& xprime, 
		   const smat_T& A, const vec_Tp& y, const vec_long& cols);
  

  // compare long values
  static int compar_long(const void* a, const void* b) {
    if (*(long*)a < *(long*)b)
      return -1;
    else if (*(long*)a == *(long*)b)
      return 0;
    return 1;
  }
  

public:
  long nrows;       // current number of non-trivial rows
  long first_col;   // first non-trivial column in AEcs

  const smat_T& A;  // original matrix
  vec_long Arw;     // row weights (0 for eliminated row)

  smat_T E;         // extra rows  
  vec_long Erw;     // row weights (0 for eliminated row)

  vec_vec_long AEc; // column index of matrix (negative values for extra rows) 
  // note: row numbers in AEc are kept sorted in reverse order

  vec_long AEcw;    // column weights (0 for eliminated column)
  vec_long AEcs;    // columns indices sorted by weight
  vec_long AEcr;    // reverse pointers into AEcs to enable sort adjustment
  // invariant: c = AEcs[AEcr[c]]

  vec_Tp Ay;        // constant term
  vec_Tp Ey;        // constant term for extra rows
};


template <class smat_T, class svec_T, class vec_Tp, class Tp>
smat_gauss<smat_T,svec_T,vec_Tp,Tp>::smat_gauss(const smat_T& A1, 
					       const vec_Tp& y1)
  : nrows(0), A(A1), Ay(y1) {

  Arw.SetLength(A.NumRows());
  AEc.SetLength(A.NumCols());
  AEcw.SetLength(A.NumCols());
  
  // scan matrix to generate Arw, AEc and AEcw
  bool row1 = false;
  for (long i=A.NumRows()-1; i>=0; --i) {
    const long* ri = A[i].indices();
    const value_t* rv = A[i].values();
    for (long j=0; j<A[i].nvalues(); ++j)
      if (!IsZero(rv[j])) {
	++Arw[i];
	append(AEc[ri[j]],i);
	++AEcw[ri[j]];
      }
    if (Arw[i]>=1) {
      ++nrows;
      if (Arw[i]==1)
	row1=true;
    }
  }

  // sort columns by weight
  long* ar = new long[2*A.NumCols()];
  for (long i=0; i<A.NumCols(); ++i) {
    ar[2*i] = AEcw[i];
    ar[2*i+1] = i;
  }
  qsort(ar,A.NumCols(),2*sizeof(long),compar_long);
  AEcs.SetLength(A.NumCols());
  AEcr.SetLength(A.NumCols());
  for (long i=0; i<A.NumCols(); ++i) {
    AEcs[i] = ar[2*i+1];
    AEcr[ar[2*i+1]] = i;
  }
  delete[] ar;

  // eliminate weight 1 rows
  if (row1) {
    for (long i=0; i<A.NumRows(); ++i)
      if (Arw[i]==1)
	eliminate_row1(i);
  }

  // find first non-trivial column
  first_col=0;
  while (AEcw[AEcs[first_col]]==0)
    ++first_col;
}

template <class smat_T, class svec_T, class vec_Tp, class Tp>
smat_gauss<smat_T,svec_T,vec_Tp,Tp>::~smat_gauss() {
}

// bubble sort AEcs to ensure column c is correctly positioned
template <class smat_T, class svec_T, class vec_Tp, class Tp>
void smat_gauss<smat_T,svec_T,vec_Tp,Tp>::resort(long c) {
  // my column weight is AEcw[c] = AEcw[AEcs[AEcr[c]]]
  // prev weight is AEcw[AEcs[AEcr[c]-1]]
  while (AEcr[c]>0 && AEcw[AEcs[AEcr[c]-1]]>AEcw[c]) {
    swap(AEcs[AEcr[c]],AEcs[AEcr[c]-1]);
    ++AEcr[AEcs[AEcr[c]]];
    --AEcr[c];
  }
}

template <class smat_T, class svec_T, class vec_Tp, class Tp>
void smat_gauss<smat_T,svec_T,vec_Tp,Tp>::eliminate_row1(long i) {
  set_row_weight(i,0);
  --nrows;
  const long* ri = row(i).indices();
  const value_t* rv = row(i).values();
  for (long j=0; j<row(i).nvalues(); ++j) {
    long c = ri[j];
    if (!IsZero(rv[j]) && AEcw[c]!=0) {
      bool row1 = false;
      // remove term from remaining rows
      for (long k=0; k<AEc[c].length(); ++k) {
	long r = AEc[c][k];
	if (row_weight(r)>0) {
	  y(r) -= row(r)[c]*y(i)/rv[j];
	  set_row_weight(r,row_weight(r)-1);
	  if (row_weight(r)<=1) {
	    if (row_weight(r)==1)
	      row1 = true;
	    else
	      --nrows;
	  }
	}
      }
      // eliminate column
      AEcw[c] = 0;
      resort(c);
      // eliminate new weight 1 rows
      if (row1) {
	for (long k=0; k<AEc[c].length(); ++k)
	  if (row_weight(AEc[c][k])==1) 
	    eliminate_row1(AEc[c][k]);
      }
      return;
    }
  }
  Error("smat_gauss::eliminate_row1() FATAL ERROR!");
}

template <class smat_T, class svec_T, class vec_Tp, class Tp>
void smat_gauss<smat_T,svec_T,vec_Tp,Tp>::eliminate_col1(long c) {
  // find the row that has the non-zero entry and eliminate it
  for (long i=0; i<AEc[c].length(); ++i) {
    long r = AEc[c][i];
    if (row_weight(r)>0) {
      AEcw[c]=0;
      resort(c);
      set_row_weight(r,0);
      --nrows;
      const long* ri = row(r).indices();
      const value_t* rv = row(r).values();
      for (long j=0; j<row(r).nvalues(); ++j) 
	if (AEcw[ri[j]]>0 && !IsZero(rv[j])) {
	  --AEcw[ri[j]];
	  resort(ri[j]);
	}
      return;
    }
  }
  Error("smat_gauss::eliminate_col1() FATAL ERROR!");
}

template <class smat_T, class svec_T, class vec_Tp, class Tp>
void smat_gauss<smat_T,svec_T,vec_Tp,Tp>::eliminate_col2(long c) {
  // find the two rows
  long i=0;
  long r1=0,r2;
  while (i<AEc[c].length()) {
    r1 = AEc[c][i++];
    if (row_weight(r1)>0)
      break;
  }
  while (i<AEc[c].length()) {
    r2 = AEc[c][i++];
    if (row_weight(r2)>0) {
      AEcw[c]=0;
      resort(c);
      // subtract rows to create new (extra) row having 0 in this column
      long r = E.NumRows();
      E.SetDims(r+1,A.NumCols());
      Ey.SetLength(r+1);
      value_t g;
      g = GCD(row(r1)[c],row(r2)[c]);
      E[r] = (row(r2)[c]/g)*row(r1) - (row(r1)[c]/g)*row(r2);
      Ey[r] = (row(r2)[c]/g)*y(r1) - (row(r1)[c]/g)*y(r2);
      // fixup EAc, Erw
      Erw.SetLength(r+1);
      Erw[r] = 0;
      svec_long rnz;
      nonzero(rnz,E[r]);
      const long* ri = rnz.indices();
      const svec_long::value_type* rv = rnz.values();
      for (long k=0; k<rnz.nvalues(); ++k) {
	if (rv[k]==1 && AEcw[ri[k]]>0) {
	  append(AEc[ri[k]],-r-1);
	  ++Erw[r];
	}
      }
      if (Erw[r]>0)
	++nrows;
      // fixup AEcw
      svec_long nz;
      nonzero(nz,row(r1));
      rnz -= nz;
      nonzero(nz,row(r2));
      rnz -= nz;
      ri = rnz.indices();  // these may have changed
      rv = rnz.values();
      for (long k=0; k<rnz.nvalues(); ++k) 
	if (rv[k]<0 && AEcw[ri[k]]>0) {
	  AEcw[ri[k]] += rv[k];
	  resort(ri[k]);
	}
      // remove original rows
      set_row_weight(r1,0);
      set_row_weight(r2,0);
      nrows-=2;
      if (Erw[r]==1)
	eliminate_row1(r);
      return;
    }
  }
  Error("smat_gauss::eliminate_col2() FATAL ERROR!");
}

template <class smat_T, class svec_T, class vec_Tp, class Tp>
void smat_gauss<smat_T,svec_T,vec_Tp,Tp>::eliminate_row(long r) {
  const long* ri = row(r).indices();
  const value_t* rv = row(r).values();
  for (long j=0; j<row(r).nvalues(); ++j) {
    long c = ri[j];
    if (!IsZero(rv[j]) && AEcw[c]>0) {
      --AEcw[c];
      resort(c);
    }
  }
  set_row_weight(r,0);
  --nrows;
}

template <class smat_T, class svec_T, class vec_Tp, class Tp>
void smat_gauss<smat_T,svec_T,vec_Tp,Tp>::eliminate_rows(long max) {
  if (max<1)
    return;
  long start=first_col;
  long end=first_col+1;
  while (end<AEcs.length() && AEcw[AEcs[end]]==AEcw[AEcs[start]])
    ++end;
  // the columns [start,end) all have the same weight
  set_long s(-E.NumRows(),A.NumRows()-1);
  vec_long dups;
  for (long i=start; i<end && (dups.length()==0 || max>1); ++i) {
    for (long j=0; j<AEc[AEcs[i]].length(); ++j) {
      long r = AEc[AEcs[i]][j];
      if (row_weight(r)>0) {
	if (s.contains(r)) {
	  append(dups,r);
	  if (max==1)
	    break;
	}
	else
	  s.insert(r);
      }
    }
  }
  while (dups.length()==0 && end<AEcs.length()) {
    // extend end to find rows to delete
    // todo: should we do a whole new [start,end) block?
    for (long j=0; j<AEc[AEcs[end]].length(); ++j) {
      long r = AEc[AEcs[end]][j];
      if (row_weight(r)>0) {
	if (s.contains(r)) {
	  append(dups,r);
	  if (max==1)
	    break;
	}
      }
    }
    ++end;
  }
  // delete all rows in dups
  for (long i=0; i<dups.length(); ++i)
    if (row_weight(dups[i])>0) {
      eliminate_row(dups[i]);
      if (--max<=0)
	break;
    }
}

// reduce matrix to (N+e) x N
template <class smat_T, class svec_T, class vec_Tp, class Tp>
void smat_gauss<smat_T,svec_T,vec_Tp,Tp>::reduce(long extra) {
  do {
    while (first_col<AEcs.length() && AEcw[AEcs[first_col]]<=2) {
      if (AEcw[AEcs[first_col]]==1)
	eliminate_col1(AEcs[first_col]);
      else if (AEcw[AEcs[first_col]]==2)
	eliminate_col2(AEcs[first_col]);
      ++first_col;
    }
    if (NumRows()<=NumCols()+extra)
      break;
    if (NumRows()<NumCols()+extra+AEcw[AEcs[first_col]]-2)
      break;  // eliminating more rows won't help
    eliminate_rows(1);
  } while (true);
}

template <class smat_T, class svec_T, class vec_Tp, class Tp>
void smat_gauss<smat_T,svec_T,vec_Tp,Tp>::make_system(smat_T& newA, 
						      vec_Tp& newy,
						      vec_long& cols) const {
  vec_long cmap;
  cmap.SetLength(AEcw.length());
  cols.SetLength(NumCols());
  for (long i=0,j=0; i<AEcw.length(); ++i)
    if (AEcw[i]>0) {
      cols[j] = i;
      cmap[i] = j++;
    }

  newA.SetDims(NumRows(),NumCols());
  clear(newA);
  newy.SetLength(NumRows());
  clear(newy);
  long r=0;
  for (long i=0; i<A.NumRows(); ++i) 
    if (Arw[i]>0) {
      const long* ri = A[i].indices();
      const value_t* rv = A[i].values();
      for (long j=0; j<A[i].nvalues(); ++j)
	if (!IsZero(rv[j]) && AEcw[ri[j]]>0)
	  newA[r][cmap[ri[j]]] = rv[j];
      newy[r] = Ay[i];
      ++r;
    }
  for (long i=0; i<E.NumRows(); ++i)
    if (Erw[i]>0) {
      const long* ri = E[i].indices();
      const value_t* rv = E[i].values();
      for (long j=0; j<E[i].nvalues(); ++j)
	if (!IsZero(rv[j]) && AEcw[ri[j]]>0)
	  newA[r][cmap[ri[j]]] = rv[j];
      newy[r] = Ey[i];
      ++r;
    }
}


template <class smat_T, class svec_T, class vec_Tp, class Tp>
bool smat_gauss<smat_T,svec_T,vec_Tp,Tp>::undo(vec_Tp& x, 
					       const vec_Tp& xprime, 
					       const smat_T& A, 
					       const vec_Tp& y, 
					       const vec_long& cols) {
  // copy known results to xorig
  x.SetLength(A.NumCols());
  clear(x);
  bool* known = new bool[x.length()];
  memset(known,0,x.length()*sizeof(bool));
  long nknown = 0;
  for (long i=0; i<cols.length(); ++i) {
    x[cols[i]] = xprime[i];
    known[cols[i]] = true;
    ++nknown;
  }

  // keep track of which rows we have used
  bool* used = new bool[A.NumRows()];
  memset(used,0,A.NumRows()*sizeof(bool));

  // attempt to compute unknown values
  bool done;
  do {
    done=true;
    for (long i=0; i<A.NumRows(); ++i) 
      if (!used[i]) {
	// see if we have exactly one unknown value with non-zero coefficient
	long unknown=-1;
	long n = A[i].nvalues();
	const long* aj = A[i].indices();
	const value_t* av = A[i].values();
	used[i]=true;
	for (long j=0; j<n; ++j) 
	  if (!IsZero(av[j]) && !known[aj[j]]) {
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
	  Tp v;
	  InnerProduct(v,A[i],x);
	  v -= y[i];
	  v /= A[i][unknown];
	  negate(x[unknown],v);
	  known[unknown] = true;
	  ++nknown;
	  done=false;
	}
      }
  } while (!done);

  delete[] used;

  /*
  if (nknown<x.length()) {
    for (long i=0; i<x.length(); ++i)
      if (!known[i]) {
	std::cerr<<"SGauss failed to recover variable "<<i<<std::endl;
	break;
      }
  }
  */

  delete[] known;
  return nknown>=x.length();
}




/* Structured Gaussian Elimination.
 * Finds a smaller version of system Aorig*x=yorig and returns it
 * in Anew and ynew.  cols is the list of column indices from Aorig
 * that were kept in Anew.
 */
void SGauss(smat_long& Anew, vec_ZZ_p& ynew, vec_long& cols, 
	    const smat_long& Aorig, const vec_ZZ_p& yorig) {
  
  if (SGAUSS_VERBOSE)
    std::cerr<<"SGauss: pre  "<<Aorig.NumRows()<<'x'<<Aorig.NumCols()
	     <<std::endl;
  smat_gauss<smat_long,svec_long,vec_ZZ_p,ZZ_p> SGE(Aorig,yorig);
  if (SGAUSS_VERBOSE)
    std::cerr<<"SGauss: init "<<SGE.NumRows()<<'x'<<SGE.NumCols()<<std::endl;
  SGE.reduce();
  SGE.make_system(Anew,ynew,cols);
  if (SGAUSS_VERBOSE)
    std::cerr<<"SGauss: done "<<SGE.NumRows()<<'x'<<SGE.NumCols()<<std::endl;
}
void SGauss(smat_ZZ& Anew, vec_ZZ_p& ynew, vec_long& cols, 
	    const smat_ZZ& Aorig, const vec_ZZ_p& yorig) {
  if (SGAUSS_VERBOSE)
    std::cerr<<"SGauss: pre  "<<Aorig.NumRows()<<'x'<<Aorig.NumCols()
	     <<std::endl;
  smat_gauss<smat_ZZ,svec_ZZ,vec_ZZ_p,ZZ_p> SGE(Aorig,yorig);
  if (SGAUSS_VERBOSE)
    std::cerr<<"SGauss: init "<<SGE.NumRows()<<'x'<<SGE.NumCols()<<std::endl;
  SGE.reduce();
  SGE.make_system(Anew,ynew,cols);
  if (SGAUSS_VERBOSE)
    std::cerr<<"SGauss: done "<<SGE.NumRows()<<'x'<<SGE.NumCols()<<std::endl;
}


/* Undo gaussian elimination.  Assuming Aorig and yorig were reduced with
 * SGauss(), cols is the list of columns that were kept, and x is the solution
 * to the reduced system that was found, then xorig, the solution to the
 * original system, will be computed.
 */
bool SGauss_undo(vec_ZZ_p& xorig, const vec_ZZ_p& x, 
		 const vec_long& cols,
		 const smat_long& Aorig, const vec_ZZ_p& yorig) {
  return 
    smat_gauss<smat_long,svec_long,vec_ZZ_p,ZZ_p>::undo(xorig,x,Aorig,yorig,
							cols);
}
bool SGauss_undo(vec_ZZ_p& xorig, const vec_ZZ_p& x, 
		 const vec_long& cols,
		 const smat_ZZ& Aorig, const vec_ZZ_p& yorig) {
  return 
    smat_gauss<smat_ZZ,svec_ZZ,vec_ZZ_p,ZZ_p>::undo(xorig,x,Aorig,yorig,cols);
}

NTL_END_IMPL;
