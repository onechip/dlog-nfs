#ifndef NTL_svector__H
#define NTL_svector__H

/*
Copyright (C) 2002 Chris Studholme

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

*****************************************************

MODULE: svector

SUMMARY:

  ... see doc/svector.txt ...

*/


#include <string.h>
#include <NTL/tools.h>
#include <NTL/vec_long.h>



/* Declaration of sparse vector class.
 */
#define NTL_svector_decl(T,vec_T,svec_T) \
class svec_T { \
protected: \
  vec_long index; /* indices of allocated elements */ \
  vec_T value;    /* values of allocated elements */ \
 \
  long len;       /* maximum value of index */ \
  bool len_fixed; /* is the length of this vector fixed? */ \
  long len_max;   /* for compatibility with NTL's vector class */ \
 \
public: \
  /* Default constructor.  Initial length is zero and no memory \
   * is allocated. \
   */ \
  svec_T() { \
    len=len_max=0; \
    len_fixed=false; \
  } \
 \
 \
  /* Copy constructor.  Allocates exactly the right amount of memory \
   * to hold the contents of other and uses T's assignment operator to \
   * do the copying.  If compact=true then zero values in other are \
   * not copied.  The "fixed" status of other is not copied. \
   */ \
  svec_T(const svec_T& other, bool compact=false);  \
 \
  /* Assignment.  Performs an element-wise assignment using T's assignment \
   * operator.  If this is "fixed", other must have the same length as this. \
   */ \
  svec_T& operator=(const svec_T& other);   \
 \
  /* Initialize with a specific length and allocation.  If alloc is zero then \
   * no memory is allocated until needed (and then the default allocation is \
   * made). \
   */ \
  svec_T(NTL_NNS INIT_SIZE_TYPE, long length, long alloc=0); \
 \
  /* Initialize to empty array and swap with other.  \
   */ \
  svec_T(svec_T& other, NTL_NNS INIT_TRANS_TYPE) { \
    len=len_max=0; \
    len_fixed=false; \
    swap(other); \
  } \
 \
  /* Destructor.  Free all allocated memory. \
   */  \
  ~svec_T() {} \
 \
  /* Set current length to n.  If the length is being set smaller than it was,\
   * values at the end of the vector may disappear.  If n is being set larger,\
   * existing values are not changed and new values are cleared (zeroed). \
   */ \
  void SetLength(long n);   \
 \
  /* Current length. \
   */ \
  inline long length() const { \
    return len; \
  } \
 \
  /* Indexing operation, starting from 0.  This version returns a non-const \
   * reference to a T and will allocate memory for the vector element as \
   * needed.  If you are just reading vector elements and you want to make \
   * sure you maintain the sparsity of the vector, use the const version \
   * instead. \
   */ \
  inline T& operator[](long i) { \
    NTLX_RANGE_CHECK_CODE; \
    return RawGet(i); \
  } \
 \
  /* Indexing operation, starting from 0.  This version can be used with a  \
   * const object and returns a const reference to a T.  No memory will be \
   * allocated.  If the i'th vector element is not allocated, a const \
   * reference to zero will be returned. \
   */ \
  inline const T& operator[](long i) const { \
    NTLX_RANGE_CHECK_CODE; \
    return RawGet(i); \
  } \
 \
  /* Indexing operation, starting from 1.  See operator[] for details. \
   */ \
  inline T& operator()(long i) { \
    --i; NTLX_RANGE_CHECK_CODE; \
    return RawGet(i); \
  } \
 \
  /* Indexing operation, starting from 1.  See operator[] for details. \
   */ \
  inline const T& operator()(long i) const { \
    --i; NTLX_RANGE_CHECK_CODE; \
    return RawGet(i); \
  } \
 \
  /* Direct access to arrays, use with care. \
   */ \
  inline long nvalues() const { \
    return value.length(); \
  } \
  inline const long* indices() const { \
    return index.elts(); \
  } \
  inline T* values() { \
    return value.elts(); \
  } \
  inline const T* values() const { \
    return value.elts(); \
  } \
 \
  /* Release all allocated space and set to length 0. \
   */ \
  void kill();  \
 \
  /* Ensure enough memory is allocated for n non-zero entries. \
   */ \
  void SetAlloc(long n); \
 \
  /* Set maximum length to n.  Does not allocate memory and does not change \
   * the current length.  Here for compatibility with non-sparse vector \
   * classes. \
   */ \
  void SetMaxLength(long n);  \
 \
  /* Sets length to n and prohibits all future length changes. \
   * FixLength may only be invoked immediately after the default \
   * construction or kill. \
   * \
   * The kill operation is also subsequently prohibited, and swap is \
   * allowed on fixed length vectors of the same length. \
   * \
   * FixLength is provided mainly to implement smat_T, to enforce \
   * the restriction that all rows have the same length. \
   */ \
  void FixLength(long n); \
 \
  /* Has this vector been fixed? \
   */ \
  inline bool fixed() const { \
    return len_fixed; \
  } \
 \
  /* Maximum length ever achieved.  Here for compatibility with non-sparse  \
   * vector classes. \
   */ \
  inline long MaxLength() const { \
    return len_max; \
  } \
 \
  /* The number of objects for which space has been allocated but not  \
   * necessarily initialized.  \
   */ \
  inline long allocated() const { \
    return value.allocated(); \
  } \
 \
  /* Indexing with no range checking.  See operator[] for details. \
   */ \
  T& RawGet(long i); \
  const T& RawGet(long i) const; \
 \
  /* Returns position of a in the vector, or -1 if it is not there. \
   * The search is conducted from position 0 to MaxAlloc()-1 of the vector, \
   * and an error is raised if the object is found at position MaxLength() \
   * or higher (in which case a references an uninitialized object). \
   * Note that if NTL_CLEAN_PTR flag is set, this routine takes \
   * linear time, and otherwise, it takes constant time. \
   */ \
  long position(const T& a) const; \
 \
  /* Returns position of a in the vector, or -1 if it is not there. \
   * The search is conducted from position 0 to length()-1 of the vector. \
   * Note that if NTL_CLEAN_PTR flag is set, this routine takes \
   * linear time, and otherwise, it takes constant time. \
   */ \
  long position1(const T& a) const; \
 \
  /* Eliminate zero values from vector to speed up read-only access. \
   */ \
  void compact(); \
 \
  /* Set all elements to zero.  Does not change length or release any memory. \
   */ \
 inline void clear() { \
   index.SetLength(0); \
   value.SetLength(0); \
 } \
 \
  /* Swaps this and other by swapping internal pointers.  If either this or \
   * other is fixed, then both must be have the same length.  The "fixed" \
   * status of the vectors is not swapped. \
   */ \
  void swap(svec_T& other); \
};\
 \
long IsZero(const svec_T &x); \
 \
inline void clear(svec_T &x) { \
  x.clear(); \
} \
 \
inline void swap(svec_T &x, svec_T &y) { \
  x.swap(y); \
} \
 \
void conv(svec_T& dest, const vec_T& src); \
 \
void conv(vec_T& dest, const svec_T& src); \
 \
void VectorCopy(svec_T& x, const svec_T& a, long n); \
 \
inline svec_T VectorCopy(const svec_T& a, long n) { \
  svec_T x; VectorCopy(x, a, n); NTL_OPT_RETURN(svec_T,x); \
} \



// standard math methods
#define NTL_math_svector_decl(T,vec_T,svec_T)  \
void mul(svec_T& x, const svec_T& a, const T& b); \
void mul(svec_T& x, const T& a, const svec_T& b); \
 \
void add(svec_T& x, const svec_T& a, const svec_T& b); \
void sub(svec_T& x, const svec_T& a, const svec_T& b); \
void negate(svec_T& x, const svec_T& a); \
 \
void InnerProduct(T& x, const svec_T& a, const svec_T& b); \
 \
inline svec_T operator+(const svec_T& a, const svec_T& b) { \
  svec_T x; add(x,a,b); NTL_OPT_RETURN(svec_T, x); \
} \
inline svec_T operator-(const svec_T& a, const svec_T& b) { \
  svec_T x; sub(x,a,b); NTL_OPT_RETURN(svec_T, x); \
} \
inline svec_T operator-(const svec_T& a) { \
  svec_T x; negate(x,a); NTL_OPT_RETURN(svec_T, x); \
} \
 \
inline svec_T operator*(const svec_T& a, const T& b) { \
  svec_T x; mul(x,a,b); NTL_OPT_RETURN(svec_T, x); \
} \
 \
inline svec_T operator*(const T& a, const svec_T& b) { \
  svec_T x; mul(x,a,b); NTL_OPT_RETURN(svec_T, x); \
} \
 \
inline T operator*(const svec_T& a, const svec_T& b) { \
  T x;  InnerProduct(x,a,b);  return x; \
} \
 \
inline svec_T& operator+=(svec_T& x, const svec_T& a) { \
  add(x, x, a);  return x; \
} \
 \
inline svec_T& operator-=(svec_T& x, const svec_T& a) { \
  sub(x, x, a);  return x; \
} \
 \
/* Note: operator*=(svec_T& x, const T& a) is not provided here because \
 *       its ambiguous as to whether it's multipy on left or multiply on \
 *       right. */ \
 \
/* mixed vector/svector math methods */ \
void mul(vec_T& x, const svec_T& a, const T& b); \
void mul(vec_T& x, const T& a, const svec_T& b); \
 \
void add(vec_T& x, const vec_T& a, const svec_T& b); \
void add(vec_T& x, const svec_T& a, const vec_T& b); \
void sub(vec_T& x, const vec_T& a, const svec_T& b); \
void sub(vec_T& x, const svec_T& a, const vec_T& b); \
void negate(vec_T& x, const svec_T& a); \
 \
void InnerProduct(T& x, const vec_T& a, const svec_T& b); \
void InnerProduct(T& x, const svec_T& a, const vec_T& b); \
 \
inline vec_T operator+(const vec_T& a, const svec_T& b) { \
  vec_T x; add(x,a,b); NTL_OPT_RETURN(vec_T, x); \
} \
inline vec_T operator+(const svec_T& a, const vec_T& b) { \
  vec_T x; add(x,a,b); NTL_OPT_RETURN(vec_T, x); \
} \
inline vec_T operator-(const vec_T& a, const svec_T& b) { \
  vec_T x; sub(x,a,b); NTL_OPT_RETURN(vec_T, x); \
} \
inline vec_T operator-(const svec_T& a, const vec_T& b) { \
  vec_T x; sub(x,a,b); NTL_OPT_RETURN(vec_T, x); \
} \
 \
inline T operator*(const vec_T& a, const svec_T& b) { \
  T x;  InnerProduct(x,a,b);  return x; \
} \
inline T operator*(const svec_T& a, const vec_T& b) { \
  T x;  InnerProduct(x,a,b);  return x; \
} \
 \
inline vec_T& operator+=(vec_T& x, const svec_T& a) { \
  add(x,x,a);  return x; \
} \
 \
inline vec_T& operator-=(vec_T& x, const svec_T& a) { \
  sub(x,x,a);  return x; \
} \


#define NTL_eq_svector_decl(T,svec_T)  \
long operator==(const svec_T& l__a, const svec_T& l__b);  \
inline long operator!=(const svec_T& l__a, const svec_T& l__b) { \
  return !(l__a==l__b); \
} \


/*
long operator==(const vec_T& l__a, const svec_T& l__b);  \
inline long operator!=(const vec_T& l__a, const svec_T& l__b) { \
  return !(l__a==l__b); \
} \
long operator==(const svec_T& l__a, const vec_T& l__b);  \
inline long operator!=(const svec_T& l__a, const vec_T& l__b) { \
  return !(l__a==l__b); \
} \
*/


#define NTL_io_svector_decl(T,svec_T)  \
NTL_SNS istream& operator>>(NTL_SNS istream&, svec_T&);  \
NTL_SNS ostream& operator<<(NTL_SNS ostream&, const svec_T&); \
NTL_SNS ostream& OutputVector(NTL_SNS ostream&, const svec_T& a); \


#ifndef NTL_RANGE_CHECK
#define NTLX_RANGE_CHECK_CODE 
#else
#define NTLX_RANGE_CHECK_CODE if ((i<0)||(i>=len)) RangeError(i);
#endif


// zero_T is a variable containing the value 0 that can be 
// (read-only) referenced
#define NTL_svector_impl(T,vec_T,svec_T,zero_T) \
svec_T::svec_T(const svec_T& a, bool compact) { \
  len=a.len; \
  len_fixed=false; \
  len_max=a.len_max; \
  if (!compact) { \
    index=a.index; \
    value=a.value; \
  } \
  else { \
    long nvals=0; \
    for (long i=0; i<a.value.length(); ++i) \
      if (!IsZero(a.value.RawGet(i))) \
	++nvals; \
    index.SetLength(nvals); \
    value.SetLength(nvals); \
    if (nvals>0) { \
      long* ind = index.elts(); \
      T* val = value.elts(); \
      for (long i=0,j=0; i<a.value.length(); ++i)  \
	if (!IsZero(a.value.RawGet(i))) { \
          val[j++] = a.value.RawGet(i); \
	  ind[j++] = a.index.RawGet(i); \
	} \
    } \
  } \
} \
 \
svec_T::svec_T(INIT_SIZE_TYPE, long n, long alloc) { \
  if (alloc>0) { \
    index.SetMaxLength(alloc); \
    value.SetMaxLength(alloc); \
  } \
  len_max=len=n; \
  len_fixed=false; \
} \
 \
svec_T& svec_T::operator=(const svec_T& a) { \
  if (len_fixed&&(len!=a.len)) \
    NTL_NNS Error("svector::operator=() attempt to change length of fixed vector"); \
  len=a.len; \
  if (len>len_max) len_max=len; \
  index=a.index; \
  value=a.value; \
  return *this; \
} \
 \
void svec_T::SetLength(long n) { \
  len=n; \
  if (len>len_max) len_max=len; \
  long i = index.length(); \
  while (i>0 && index.RawGet(--i)>=len) \
    index.SetLength(i); \
  value.SetLength(index.length()); \
} \
  \
T& svec_T::RawGet(long i) { \
  long nvals = index.length(); \
  long low = 0; \
  if (nvals>0) { \
    if (i>index.RawGet(nvals-1)) \
      /* special optimization (for appends) */ \
      low=nvals; \
     \
    else { \
      /* binary search to find value or location of insert */ \
      long high = nvals; \
      while (high>low) { \
	long mid = (low+high)/2; \
	if (index.RawGet(mid)<i) \
	  low = mid+1; \
	else if (i<index.RawGet(mid)) \
	  high = mid; \
	else \
	  return value.RawGet(mid); \
      } \
    } \
  } \
 \
  /* insert new value at position low */ \
  index.SetLength(nvals+1); \
  value.SetLength(nvals+1); \
  ::clear(value.RawGet(nvals)); \
  if (low<nvals) { \
    memmove(index.elts()+low+1,index.elts()+low,(nvals-low)*sizeof(long)); \
    for (long j=nvals; j>low; --j) /* swap is faster than assignment */ \
      ::swap(value.RawGet(j),value.RawGet(j-1));  \
  } \
  index.RawGet(low) = i; \
  return value.RawGet(low); \
} \
 \
const T& svec_T::RawGet(long i) const { \
  /* binary search */ \
  long low = 0; \
  long high = index.length(); \
  while (high>low) { \
    long mid = (low+high)/2; \
    if (index.RawGet(mid)<i) \
      low = mid+1; \
    else if (i<index.RawGet(mid)) \
      high = mid; \
    else \
      return value.RawGet(mid); \
  } \
  return zero_T; \
}  \
 \
void svec_T::kill() { \
  if (len_fixed) \
    NTL_NNS Error("svector::kill() cannot kill a fixed vector"); \
  index.kill(); \
  value.kill(); \
  len=0; \
  len_max=0; \
} \
 \
void svec_T::SetAlloc(long n) { \
  if (index.length()<n) { \
    index.SetMaxLength(n); \
    value.SetMaxLength(n); \
  } \
} \
 \
void svec_T::SetMaxLength(long n) { \
  if (n>len_max) len_max=n; \
} \
 \
void svec_T::FixLength(long n) { \
  if (len) \
    NTL_NNS Error("svector::FixLength() length has already been set"); \
  len_fixed=true; \
  len=n; \
  if (len>len_max) len_max=len; \
} \
 \
long svec_T::position(const T& a) const { \
  long i = value.position(a); \
  return i<0 ? -1 : index[i]; \
} \
 \
long svec_T::position1(const T& a) const { \
  long i = value.position1(a); \
  return (i<0 || index[i]>=len) ? -1 : index[i]; \
} \
 \
void svec_T::compact() { \
  long* ind = index.elts(); \
  T* val = value.elts(); \
  /* this could be sped up by using swap(T,T) instead of T::operator=() */ \
  for (long i=value.length()-1; i>=0; --i) \
    if (IsZero(val[i])) { \
      long nvals = value.length()-1; \
      for (long j=i; j<nvals; ++j) { \
	ind[j] = ind[j+1]; \
	val[j] = val[j+1]; \
      } \
      index.SetLength(nvals); \
      value.SetLength(nvals); \
    } \
} \
 \
void svec_T::swap(svec_T& a) { \
  /* if either is fixed, they must be same length and \
   * "fixed" status is not swapped */ \
  if (len_fixed||a.len_fixed) { \
    if (len!=a.len) \
      NTL_NNS Error("svector::swap() attempt to change length of fixed vector"); \
  } \
 \
  long o_len = len; \
  len = a.len; \
  a.len = o_len; \
 \
  long o_len_max = len_max; \
  len_max = a.len_max; \
  a.len_max = o_len_max; \
 \
  ::swap(index,a.index); \
  ::swap(value,a.value); \
} \
 \
long IsZero(const svec_T &x) { \
  /* search for non-zero value */ \
  const T* v = x.values(); \
  for (long i=x.nvalues()-1; i>=0; --i) \
    if (!IsZero(v[i])) \
      return false; \
  return true; \
} \
 \
void VectorCopy(svec_T& x, const svec_T& a, long n) { \
  if (n<0) Error("VectorCopy: negative length"); \
  if (n>=(1L<<(NTL_BITS_PER_LONG-2))) Error("overflow in VectorCopy"); \
 \
  x.SetLength(n); \
  clear(x); \
 \
  long an = a.nvalues(); \
  const long* ai = a.indices(); \
  const T* av = a.values(); \
 \
  for (long i=0; (i<an)&&(ai[i]<n); ++i) \
    if (!IsZero(av[i])) \
      x[ai[i]] = av[i]; \
} \
 \
void conv(svec_T& dest, const vec_T& src) { \
  dest.SetLength(src.length()); \
  clear(dest); \
  for (long i=0; i<src.length(); ++i) \
    if (!IsZero(src[i])) \
      dest[i]=src[i]; \
} \
void conv(vec_T& dest, const svec_T& src) { \
  dest.SetLength(src.length()); \
  clear(dest); \
  long n = src.nvalues(); \
  const long* si = src.indices(); \
  const T* sv = src.values(); \
  for (long i=0; i<n; ++i) \
    dest[si[i]] = sv[i]; \
} \


// efficient processing of sparse array
#define NTL_svector_impl_enum(T,a,b,a_only,b_only,a_b) \
  long an = a.nvalues(); \
  const long* ai = a.indices(); \
  const T* av = a.values(); \
 \
  long bn = b.nvalues(); \
  const long* bi = b.indices(); \
  const T* bv = b.values(); \
 \
  long aj=0; \
  long bj=0; \
 \
  do { \
    /* b's are done, process remaining a's */ \
    if (bj>=bn) { \
      while (aj<an) { \
        a_only; \
        ++aj; \
      } \
      break; \
    } \
 \
    /* a's are done, process remaining b's */ \
    if (aj>=an) { \
      while (bj<bn) { \
        b_only; \
        ++bj; \
      } \
      break; \
    } \
 \
    /* one more a than b, process the a */ \
    if (ai[aj]<bi[bj]) { \
      a_only; \
      ++aj; \
    } \
 \
    /* one more b than a, process the b */ \
    else if (bi[bj]<ai[aj]) { \
      b_only; \
      ++bj; \
    } \
 \
    /* process a and b together */ \
    else /* (ai[aj]==bi[bj]) */ { \
      a_b; \
      ++aj; \
      ++bj; \
    } \
  } while (true) \


#define NTL_math_svector_impl(T,vec_T,svec_T)  \
void mul(svec_T& x, const svec_T& a, const T& b) { \
  if (IsZero(b)) { \
    x.SetLength(a.length()); \
    clear(x); \
    return; \
  } \
 \
  long an = a.nvalues(); \
 \
  if (&x==&a) { \
    T* xv = x.values(); \
    for (long i=0; i<an; ++i) \
      mul(xv[i],xv[i],b); \
  } \
 \
  else { \
    x.SetLength(a.length()); \
    x.SetAlloc(an); \
    clear(x); \
 \
    const long* ai = a.indices(); \
    const T* av = a.values(); \
    for (long i=0; i<an; ++i) \
      if (!IsZero(av[i])) \
        mul(x[ai[i]],av[i],b); \
  } \
} \
 \
void mul(svec_T& x, const T& b, const svec_T& a) { \
  if (IsZero(b)) { \
    x.SetLength(a.length()); \
    clear(x); \
    return; \
  } \
 \
  long an = a.nvalues(); \
 \
  if (&x==&a) { \
    T* xv = x.values(); \
    for (long i=0; i<an; ++i) \
      mul(xv[i],b,xv[i]); \
  } \
 \
  else { \
    x.SetLength(a.length()); \
    x.SetAlloc(an); \
    clear(x); \
 \
    const long* ai = a.indices(); \
    const T* av = a.values(); \
    for (long i=0; i<an; ++i) \
      if (!IsZero(av[i])) \
        mul(x[ai[i]],b,av[i]); \
  } \
} \
 \
void InnerProduct(T& x, const svec_T& a, const svec_T& b) { \
  if (a.length()!=b.length()) \
    Error("svector InnerProduct: dimension mismatch"); \
 \
  long an = a.nvalues(); \
  const long* ai = a.indices(); \
  const T* av = a.values(); \
 \
  long bn = b.nvalues(); \
  const long* bi = b.indices(); \
  const T* bv = b.values(); \
 \
  long aj=0; \
  long bj=0; \
 \
  T p; \
  clear(x); \
 \
  while ((aj<an)&&(bj<bn)) { \
    if (ai[aj]<bi[bj]) \
      ++aj; \
    else if (bi[bj]<ai[aj]) \
      ++bj; \
    else { \
      mul(p,av[aj],bv[bj]); \
      x+=p; \
    } \
  } \
} \
 \
void add(svec_T& x, const svec_T& a, const svec_T& b) { \
  if (a.length()!=b.length()) \
    Error("add() vector length mismatch"); \
 \
  if (&x==&a) { \
    /* x+=b */ \
    const long* bi = b.indices(); \
    const T* bv = b.values(); \
    long bn = b.nvalues(); \
    /* worst case complexity: O(bn*xn), we could do: O(bn+xn) */ \
    for (long i=0; i<bn; ++i) \
      if (!IsZero(bv[i])) { \
        T& s = x[bi[i]]; \
        add(s,s,bv[i]); \
      } \
  } \
 \
  else if (&x==&b) { \
    /* x+=a */ \
    const long* ai = a.indices(); \
    const T* av = a.values(); \
    long an = a.nvalues(); \
    /* worst case complexity: O(an*xn), we could do: O(an+xn) */ \
    for (long i=0; i<an; ++i) \
      if (!IsZero(av[i])) { \
	T& s = x[ai[i]]; \
        add(s,av[i],s); \
      } \
  } \
 \
  else { \
    x.SetLength(a.length()); \
    x.SetAlloc(max(a.nvalues(),b.nvalues())); \
    clear(x); \
 \
    /* worst case complexity: O(an+bn) */ \
    NTL_svector_impl_enum(T,a,b, \
			  x[ai[aj]]=av[aj], \
			  x[bi[bj]]=bv[bj], \
                          add(x[ai[aj]],av[aj],bv[bj])); \
  } \
} \
 \
void sub(svec_T& x, const svec_T& a, const svec_T& b) { \
  if (a.length()!=b.length()) \
    Error("sub() vector length mismatch"); \
 \
  if (&x==&a) { \
    /* x-=b */ \
    const long* bi = b.indices(); \
    const T* bv = b.values(); \
    long bn = b.nvalues(); \
    /* worst case complexity: O(bn*xn), we could do: O(bn+xn) */ \
    for (long i=0; i<bn; ++i) \
      if (!IsZero(bv[i])) { \
        T& s = x[bi[i]]; \
        sub(s,s,bv[i]); \
      } \
  } \
 \
  else if (&x==&b) { \
    /* x-=a */ \
    const long* ai = a.indices(); \
    const T* av = a.values(); \
    long an = a.nvalues(); \
    /* worst case complexity: O(an*xn), we could do: O(an+xn) */ \
    for (long i=0; i<an; ++i) \
      if (!IsZero(av[i])) { \
	T& s = x[ai[i]]; \
        sub(s,av[i],s); \
      } \
  } \
 \
  else { \
    x.SetLength(a.length()); \
    x.SetAlloc(max(a.nvalues(),b.nvalues())); \
    clear(x); \
 \
    /* worst case complexity: O(an+bn) */ \
    NTL_svector_impl_enum(T,a,b, \
			  x[ai[aj]]=av[aj], \
			  negate(x[bi[bj]],bv[bj]), \
                          sub(x[ai[aj]],av[aj],bv[bj])); \
  } \
} \
 \
void negate(svec_T& x, const svec_T& a) { \
  long an = a.nvalues(); \
 \
  if (&x==&a) { \
    T* xv = x.values(); \
    for (long i=0; i<an; ++i) \
      negate(xv[i],xv[i]); \
  } \
 \
  else { \
    x.SetLength(a.length()); \
    x.SetAlloc(an); \
    clear(x); \
 \
    const long* ai = a.indices(); \
    const T* av = a.values(); \
    for (long i=0; i<an; ++i) \
      if (!IsZero(av[i])) \
        negate(x[ai[i]],av[i]); \
  } \
} \
 \
void mul(vec_T& x, const svec_T& a, const T& b) { \
  x.SetLength(a.length()); \
  clear(x); \
  long an = a.nvalues(); \
  const long* ai = a.indices(); \
  const T* av = a.values(); \
  for (long i=0; i<an; ++i) \
    mul(x[ai[i]],av[i],b); \
} \
 \
void mul(vec_T& x, const T& b, const svec_T& a) { \
  x.SetLength(a.length()); \
  clear(x); \
  long an = a.nvalues(); \
  const long* ai = a.indices(); \
  const T* av = a.values(); \
  for (long i=0; i<an; ++i) \
    mul(x[ai[i]],b,av[i]); \
} \
 \
void InnerProduct(T& x, const vec_T& a, const svec_T& b) { \
  if (a.length()!=b.length()) \
    Error("InnerProduct() vector length mismatch"); \
 \
  long bn = b.nvalues(); \
  const long* bi = b.indices(); \
  const T* bv = b.values(); \
 \
  T p; \
  clear(x); \
 \
  for (long i=0; i<bn; ++i) { \
    mul(p,a[bi[i]],bv[i]); \
    x+=p; \
  } \
} \
 \
void InnerProduct(T& x, const svec_T& b, const vec_T& a) { \
  if (a.length()!=b.length()) \
    Error("svector InnerProduct: dimension mismatch"); \
 \
  long bn = b.nvalues(); \
  const long* bi = b.indices(); \
  const T* bv = b.values(); \
 \
  T p; \
  clear(x); \
 \
  for (long i=0; i<bn; ++i) { \
    mul(p,bv[i],a[bi[i]]); \
    x+=p; \
  } \
} \
 \
void add(vec_T& x, const vec_T& a, const svec_T& b) { \
  if (a.length()!=b.length()) \
    Error("svector add: dimension mismatch"); \
 \
  if (&x!=&a) \
    x=a; \
 \
  /* x+=b */ \
  const long* bi = b.indices(); \
  const T* bv = b.values(); \
  long bn = b.nvalues(); \
  for (long i=0; i<bn; ++i) { \
    T& s = x[bi[i]]; \
    add(s,s,bv[i]); \
  } \
} \
 \
void add(vec_T& x, const svec_T& a, const vec_T& b) { \
  if (a.length()!=b.length()) \
    Error("svector add: dimension mismatch"); \
 \
  if (&x!=&b) \
    x=b; \
 \
  /* x+=a */ \
  const long* ai = a.indices(); \
  const T* av = a.values(); \
  long an = a.nvalues(); \
  for (long i=0; i<an; ++i) { \
    T& s = x[ai[i]]; \
    add(s,av[i],s); \
  } \
} \
 \
void sub(vec_T& x, const vec_T& a, const svec_T& b) { \
  if (a.length()!=b.length()) \
    Error("svector sub: dimension mismatch"); \
 \
  if (&x!=&a) \
    x=a; \
 \
  /* x-=b */ \
  const long* bi = b.indices(); \
  const T* bv = b.values(); \
  long bn = b.nvalues(); \
  for (long i=0; i<bn; ++i) { \
    T& s = x[bi[i]]; \
    sub(s,s,bv[i]); \
  } \
} \
 \
void sub(vec_T& x, const svec_T& a, const vec_T& b) { \
  if (a.length()!=b.length()) \
    Error("svector sub: dimension mismatch"); \
 \
  if (&x!=&b) \
    x=b; \
 \
  /* x-=a */ \
  const long* ai = a.indices(); \
  const T* av = a.values(); \
  long an = a.nvalues(); \
  for (long i=0; i<an; ++i) { \
    T& s = x[ai[i]]; \
    sub(s,av[i],s); \
  } \
} \
 \
void negate(vec_T& x, const svec_T& a) { \
  x.SetLength(a.length()); \
  clear(x); \
  long an = a.nvalues(); \
  const long* ai = a.indices(); \
  const T* av = a.values(); \
  for (long i=0; i<an; ++i) \
    negate(x[ai[i]],av[i]); \
} \


#define NTL_eq_svector_impl(T,svec_T) \
long operator==(const svec_T& a, const svec_T& b) { \
  if (a.length()!=b.length()) \
    return false; \
 \
  NTL_svector_impl_enum(T,a,b, \
	                if (!IsZero(av[aj])) return false, \
			if (!IsZero(bv[bj])) return false, \
                        if (!(av[aj]==bv[bj])) return false); \
 \
  return true; \
} \


#define NTL_io_svector_impl(T,svec_T) \
NTL_SNS istream & operator>>(NTL_SNS istream& in, svec_T& a) { \
  if (!in || !SkipWhiteSpace(in)) NTL_NNS Error("bad vector input"); \
  long c = in.peek();  \
  if (c=='<') { \
    /* sparse stream */ \
    in.get();  \
    if (!in || !SkipWhiteSpace(in)) NTL_NNS Error("bad vector input"); \
    svec_T ibuf; \
    do { \
      long index; \
      T elem; \
      if (!(in>>index)) NTL_NNS Error("bad vector input"); \
      if (!SkipWhiteSpace(in)) NTL_NNS Error("bad vector input"); \
      c = in.peek();  \
      if (c=='>') { \
        in.get(); \
        ibuf.SetLength(index); \
        break; \
      } \
      if (c==EOF) NTL_NNS Error("bad vector input"); \
      if (!(in>>elem)) NTL_NNS Error("bad vector input"); \
      if (!SkipWhiteSpace(in)) NTL_NNS Error("bad vector input"); \
      if (!IsZero(elem)) { \
        ibuf.SetLength(index+1); \
        ibuf[index]=elem; \
      } \
    } while (true); \
    a=ibuf; \
  } \
 \
  else if (c=='[') { \
    /* non-sparse stream */ \
    in.get();  \
    if (!in || !SkipWhiteSpace(in)) NTL_NNS Error("bad vector input"); \
    long n=0; \
    svec_T ibuf;  \
    ibuf.SetLength(0);  \
 \
    c = in.peek();  \
    while (c!=']' && c!=EOF) {   \
      ibuf.SetLength(++n);   \
      T elem; \
      if (!(in>>elem)) NTL_NNS Error("bad vector input");   \
      if (!IsZero(elem)) \
        ibuf[n-1]=elem; \
      if (!SkipWhiteSpace(in)) NTL_NNS Error("bad vector input"); \
      c = in.peek();  \
    } \
    if (c==EOF) NTL_NNS Error("bad vector input"); \
    in.get(); \
    a=ibuf; \
  } \
 \
  else { \
    NTL_NNS Error("bad vector input"); \
  } \
 \
  return in; \
} \
  \
NTL_SNS ostream& operator<<(NTL_SNS ostream& out, const svec_T& a) { \
  long an = a.nvalues(); \
  const long* ai = a.indices(); \
  const T* av = a.values(); \
  out << '<'; \
  for (long i=0; i<an; ++i) \
    if (!IsZero(av[i])) \
      out<<ai[i]<<' '<<av[i]<<' '; \
  out<<a.length()<<'>'; \
  return out; \
} \
 \
NTL_SNS ostream& OutputVector(NTL_SNS ostream& out, const svec_T& a) { \
  long n = a.length();  \
  out << '['; \
  if (n>0) { \
    out<<a[0]; \
    for (long i=1; i<n; ++i) \
      out<<' '<<a[i]; \
  } \
  out<<']'; \
  return out; \
} \


#endif
