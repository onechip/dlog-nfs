#ifndef NTL_smatrix__H
#define NTL_smatrix__H

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

MODULE: smatrix

SUMMARY:

Macros are defined providing template-like classes for sparse rectangular
matrices.

The macro NTL_smatrix_decl(T,svec_T,vec_svec_T,smat_T) declares a new class
smat_T representing matrices over T, where svec_T and vec_svec_T are
classes representing "NTL vectors" over T and svec_T, respectively.  
The implementation of smat_T can be instantiated with 
NTL_smatrix_impl(T,svec_T,vec_svec_T,smat_T).

The class generated with NTL_smatrix_decl and NTL_smatrix_impl is actually
identical to a class generated using the non-sparse versions of these
macros.  For complete documentation regarding these classes, see the 
matrix.txt documentation.

If T supports I/O and/or equalility testing, then smat_T can also be made 
to support these by using NTL_io_smatrix_decl(T,svec_T,vec_svec_T,smat_T),
NTL_io_smatrix_impl(T,svec_T,vec_svec_T,smat_T),
NTL_eq_smatrix_decl(T,svec_T,vec_svec_T,smat_T), and
NTL_eq_smatrix_impl(T,svec_T,vec_svec_T,smat_T)

Also, if a (non-sparse) matrix class for T has been created (mat_T), then
conv() and transpose() methods for converting between these matricies and
the sparse equivalents can be declared and implemented using
NTL_conv_smatrix_decl(T,svec_T,vec_svec_T,mat_T,smat_T) and
NTL_conv_smatrix_impl(T,svec_T,vec_svec_T,mat_T,smat_T).

Finally, functions for basic arithmatic can be decleared and implemented
using NTL_math_smatrix_decl(T,vec_T,svec_T,vec_svec_T,smat_T) and
NTL_math_smatrix_impl(T,vec_T,svec_T,vec_svec_T,smat_T).

For example of typical use, the declaration 

   smat_T M;

creates a 0 x 0 matrix.  We can make it have 10 rows and 20 columns like this:

   M.SetDims(10, 20);

A row can be accessed as M[i], indexing from 0, or as M(i), indexing from 1.
A matrix entry can be accessed as M[i][j], indexing from 0, or as
M(i, j), indexing from 1.  

A matrix is represented as a vec_svec_T: a vector of rows, where each row is
a svec_T.  Any attempt to resize one of the rows so as to create a
non-rectangular matrix will result in a run-time error.

*/


#include "ntl5_matrix.h"


#define NTL_smatrix_decl(T,svec_T,vec_svec_T,smat_T) \
  NTL_matrix_decl(T,svec_T,vec_svec_T,smat_T) \
 \
void clear(smat_T& a); \
long IsZero(const smat_T& a); \
void transpose(smat_T& dest, const smat_T& src); \


#define NTL_conv_smatrix_decl(T,svec_T,vec_svec_T,mat_T,smat_T) \
void conv(smat_T& dest, const mat_T& src); \
 \
void conv(mat_T& dest, const smat_T& src); \
 \
void transpose(smat_T& dest, const mat_T& src); \
 \
void transpose(mat_T& dest, const smat_T& src); \


#define NTL_math_smatrix_decl(T,vec_T,svec_T,vec_svec_T,smat_T) \
void add(smat_T& X, const smat_T& A, const smat_T& B); \
void sub(smat_T& X, const smat_T& A, const smat_T& B); \
void negate(smat_T& X, const smat_T& A); \
void mul(smat_T& X, const smat_T& A, const T& b); \
void mul(smat_T& X, const T& a, const smat_T& B); \
void mul(vec_T& x, const smat_T& A, const vec_T& b); \
void mul(vec_T& x, const vec_T& a, const smat_T& B); \
inline smat_T operator+(const smat_T& a, const smat_T& b) { \
  smat_T x; add(x,a,b); NTL_OPT_RETURN(smat_T, x); \
} \
inline smat_T operator-(const smat_T& a, const smat_T& b) { \
  smat_T x; sub(x,a,b); NTL_OPT_RETURN(smat_T, x); \
} \
inline smat_T operator-(const smat_T& a) { \
  smat_T x; negate(x,a); NTL_OPT_RETURN(smat_T, x); \
} \
inline smat_T operator*(const smat_T& a, const T& b) { \
  smat_T x; mul(x,a,b); NTL_OPT_RETURN(smat_T, x); \
} \
inline smat_T operator*(const T& a, const smat_T& b) { \
  smat_T x; mul(x,a,b); NTL_OPT_RETURN(smat_T, x); \
} \
inline vec_T operator*(const smat_T& a, const vec_T& b) { \
  vec_T x; mul(x,a,b); NTL_OPT_RETURN(vec_T, x); \
} \
inline vec_T operator*(const vec_T& a, const smat_T& b) { \
  vec_T x; mul(x,a,b); NTL_OPT_RETURN(vec_T, x); \
} \
inline smat_T& operator+=(smat_T& x, const smat_T& a) { \
  add(x,x,a);  return x; \
} \
inline smat_T& operator-=(smat_T& x, const smat_T& a) { \
  sub(x,x,a);  return x; \
} \


/*  Stuff we should add eventually:
void mul(smat_T& X, const smat_T& A, const smat_T& B); 
smat_T operator*(const smat_T& A, const smat_T& B);
 */


#define NTL_eq_smatrix_decl(T,svec_T,vec_svec_T,smat_T) \
  NTL_eq_matrix_decl(T,svec_T,vec_svec_T,smat_T) \


#define NTL_io_smatrix_decl(T,svec_T,vec_svec_T,smat_T) \
  NTL_io_matrix_decl(T,svec_T,vec_svec_T,smat_T) \
NTL_SNS ostream& OutputMatrix(NTL_SNS ostream&, const smat_T& a); \



#define NTL_smatrix_impl(T,svec_T,vec_svec_T,smat_T) \
  NTL_matrix_impl(T,svec_T,vec_svec_T,smat_T) \
 \
void clear(smat_T& a) { \
  for (long i=a.NumRows()-1; i>=0; --i) \
    clear(a[i]); \
} \
 \
long IsZero(const smat_T& a) { \
  for (long i=a.NumRows()-1; i>=0; --i) \
    if (!IsZero(a[i])) \
      return 0; \
  return 1; \
} \
 \
void transpose(smat_T& dest, const smat_T& src) { \
  dest.SetDims(src.NumCols(),src.NumRows()); \
  clear(dest); \
  for (long i=0; i<src.NumRows(); ++i) { \
    long sn = src[i].nvalues(); \
    const long* sj = src[i].indices(); \
    const T* sv = src[i].values(); \
    for (long j=0; j<sn; ++j) \
      dest[sj[j]][i] = sv[j]; \
  } \
} \


#define NTL_conv_smatrix_impl(T,svec_T,vec_svec_T,mat_T,smat_T) \
void conv(smat_T& dest, const mat_T& src) { \
  dest.SetDims(src.NumRows(),src.NumCols()); \
  clear(dest); \
  for (long i=0; i<src.NumRows(); ++i) \
    conv(dest[i],src[i]); \
} \
 \
void conv(mat_T& dest, const smat_T& src) { \
  dest.SetDims(src.NumRows(),src.NumCols()); \
  clear(dest); \
  for (long i=0; i<src.NumRows(); ++i) \
    conv(dest[i],src[i]); \
} \
 \
void transpose(smat_T& dest, const mat_T& src) { \
  dest.SetDims(src.NumCols(),src.NumRows()); \
  clear(dest); \
  for (long i=0; i<src.NumRows(); ++i) \
    for (long j=0; j<src.NumCols(); ++j) \
      if (!IsZero(src[i][j])) \
        dest[j][i] = src[i][j]; \
} \
 \
void transpose(mat_T& dest, const smat_T& src) { \
  dest.SetDims(src.NumCols(),src.NumRows()); \
  clear(dest); \
  for (long i=0; i<src.NumRows(); ++i) { \
    long sn = src[i].nvalues(); \
    const long* sj = src[i].indices(); \
    const T* sv = src[i].values(); \
    for (long j=0; j<sn; ++j) \
      dest[sj[j]][i] = sv[j]; \
  } \
} \


#define NTL_math_smatrix_impl(T,vec_T,svec_T,vec_svec_T,smat_T) \
void add(smat_T& X, const smat_T& A, const smat_T& B) { \
  if ((A.NumRows()!=B.NumRows())||(A.NumCols()!=B.NumCols())) \
    Error("add() matrix dimention mismatch"); \
 \
  long nrows = A.NumRows(); \
  X.SetDims(nrows,A.NumCols()); \
  for (long i=0; i<nrows; ++i) \
    add(X[i],A[i],B[i]); \
} \
 \
void sub(smat_T& X, const smat_T& A, const smat_T& B) { \
  if ((A.NumRows()!=B.NumRows())||(A.NumCols()!=B.NumCols())) \
    Error("sub() matrix dimention mismatch"); \
 \
  long nrows = A.NumRows(); \
  X.SetDims(nrows,A.NumCols()); \
  for (long i=0; i<nrows; ++i) \
    sub(X[i],A[i],B[i]); \
} \
 \
void negate(smat_T& X, const smat_T& A) { \
  long nrows = A.NumRows(); \
  X.SetDims(nrows,A.NumCols()); \
  for (long i=0; i<nrows; ++i) \
    negate(X[i],A[i]); \
} \
 \
void mul(smat_T& X, const smat_T& A, const T& b) { \
  long nrows = A.NumRows(); \
  X.SetDims(nrows,A.NumCols()); \
  for (long i=0; i<nrows; ++i) \
    mul(X[i],A[i],b); \
} \
 \
void mul(smat_T& X, const T& a, const smat_T& B) { \
  long nrows = B.NumRows(); \
  X.SetDims(nrows,B.NumCols()); \
  for (long i=0; i<nrows; ++i) \
    mul(X[i],a,B[i]); \
} \
 \
void mul(vec_T& x, const smat_T& A, const vec_T& b) { \
  if (A.NumCols()!=b.length()) \
    Error("mul() matrix dimention mismatch"); \
 \
  long n = A.NumRows(); \
 \
  if (&x==&b) { \
    vec_T b_copy(b); \
    x.SetLength(n); \
    for (long i=0; i<n; ++i) \
      InnerProduct(x[i],A[i],b_copy); \
  } \
 \
  else { \
    x.SetLength(n); \
    for (long i=0; i<n; ++i) \
      InnerProduct(x[i],A[i],b); \
  } \
} \
 \
void mul(vec_T& x, const vec_T& a, const smat_T& B) { \
  if (a.length()!=B.NumRows()) \
    Error("mul() matrix dimention mismatch"); \
 \
  svec_T tmp; \
 \
  if (&x==&a) { \
    vec_T a_copy(a); \
    x.SetLength(B.NumCols()); \
    clear(x); \
    for (long i=0; i<B.NumRows(); ++i) \
      if (!IsZero(a[i])) { \
        mul(tmp,a_copy[i],B[i]); \
        add(x,x,tmp); \
      } \
  } \
 \
  else { \
    x.SetLength(B.NumCols()); \
    clear(x); \
    for (long i=0; i<B.NumRows(); ++i) \
      if (!IsZero(a[i])) { \
        mul(tmp,a[i],B[i]); \
        add(x,x,tmp); \
      } \
  } \
} \


#define NTL_eq_smatrix_impl(T,svec_T,vec_svec_T,smat_T) \
  NTL_eq_matrix_impl(T,svec_T,vec_svec_T,smat_T) \


#define NTL_io_smatrix_impl(T,svec_T,vec_svec_T,smat_T) \
  NTL_io_matrix_impl(T,svec_T,vec_svec_T,smat_T) \
NTL_SNS ostream& OutputMatrix(NTL_SNS ostream& out, const smat_T& a) { \
  long n = a.NumRows();  \
  out << '['; \
  if (n>0) { \
    OutputVector(out,a[0]); \
    for (long i=1; i<n; ++i) \
      OutputVector(out<<' ',a[i]); \
  } \
  out<<']'; \
  return out; \
} \



#endif
