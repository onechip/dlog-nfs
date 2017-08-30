#ifndef _NTL_Lanczos__H_
#define _NTL_Lanczos__H_



/* Template for methods implementing the Lanczos algorithms for solving
 * sparse linear systems.
 *
 * Written by: Chris Studholme
 * Copyright:  GPL (http://www.fsf.org/copyleft/gpl.html)
 */


#define NTL_Lanczos_decl(vec_T,smat_T) \
bool Lanczos(const smat_T& A, vec_T& x, const vec_T& y, \
	     long deterministic=0, long verbose=0); \


#define NTL_Lanczos_impl(T,vec_T,svec_T,smat_T) \
bool Lanczos(const smat_T& A, vec_T& x, const vec_T& y, \
	     long deterministic, long verbose) { \
 \
  long nrows = A.NumRows(); \
  long ncols = A.NumCols(); \
 \
  if (deterministic) { \
    Error("Lanczos() deterministic algorithm not implemented yet"); \
  } \
 \
  if (nrows<ncols) { \
    std::cerr<<"Lanczos() insufficient rows"<<std::endl; \
    return false; \
  } \
 \
  /* temporaries */ \
  T t1; \
  svec_T ts1; \
  vec_T tv1; \
 \
  /* D is a (nrows x nrows) diagonal matrix with random entries (squared) */ \
  vec_T D; \
  D.SetLength(nrows); \
  for (long i=0; i<nrows; ++i) { \
    do { \
      random(D[i]); \
    } while (IsZero(D[i])); \
    sqr(D[i],D[i]); \
  } \
 \
  /* w0 = w = transpose(A)*D*y */ \
  vec_T w; \
  w.SetLength(ncols); \
  for (long j=0; j<nrows; ++j) { \
    mul(t1,D[j],y[j]); \
    mul(ts1,A[j],t1); \
    add(w,w,ts1); \
  } \
  vec_T w0(w); \
 \
  /* v1 = transpose(A)*D*A*w0 */ \
  vec_T v1; \
  v1.SetLength(ncols); \
  for (long j=0; j<nrows; ++j) { \
    InnerProduct(t1,A[j],w0); \
    mul(t1,D[j],t1); \
    mul(ts1,A[j],t1); \
    add(v1,v1,ts1); \
  } \
 \
  /* dv = |v1|^2 */ \
  T dv; \
  InnerProduct(dv,v1,v1); \
 \
  /* dvw = (w0,v1) */ \
  T dvw; \
  InnerProduct(dvw,w0,v1); \
   \
  /* check if dvw is zero */ \
  if (IsZero(dvw)) { \
    std::cerr<<"Lanczos() (w0,v1)=0"<<std::endl; \
    return false; \
  } \
 \
  /* w1 = v1 - (dv/dvw)*w0 */ \
  vec_T w1; \
  div(t1,dv,dvw); \
  mul(tv1,t1,w0); \
  sub(w1,v1,tv1); \
   \
  /* b = (w0,w) */ \
  T b; \
  InnerProduct(b,w0,w); \
 \
  /* initialize solution */ \
  div(t1,b,dvw); \
  mul(x,t1,w0); \
 \
  /* variables used in loop */ \
  vec_T v2; \
  vec_T w2; \
  v2.SetLength(ncols); \
 \
  /* status display */ \
  long count=0; \
  long last_count=0; \
  long last_percent=0; \
 \
  do { \
    /* v2 = transpose(A)*D*A*w1 */ \
    clear(v2); \
    for (long j=0; j<nrows; ++j) { \
      InnerProduct(t1,A[j],w1); \
      mul(t1,D[j],t1); \
      mul(ts1,A[j],t1); \
      add(v2,v2,ts1); \
    } \
 \
    long percent = (++count)*100 / nrows; \
    if ((percent!=last_percent)&&(count-last_count>=20)) { \
      std::cout<<"Lanczos: "<<count<<" itterations ("<<percent<<"%)  \r"<<std::flush; \
      last_percent=percent; \
      last_count=count; \
    } \
     \
    /* dvw = (w1,v2) */ \
    InnerProduct(dvw,w1,v2); \
    if (IsZero(dvw)) { \
      std::cout<<"Lanczos: "<<count<<" itterations        "<<std::endl; \
      if (IsZero(w1)) \
	break; \
      std::cerr<<"Lanczos() algorithm failed!"<<std::endl; \
      return false; \
    } \
 \
    /* b = (w1,w) */ \
    InnerProduct(b,w1,w); \
 \
    /* add to solution */ \
    div(t1,b,dvw); \
    mul(tv1,t1,w1); \
    add(x,x,tv1); \
 \
    /* dv = |v2|^2 */ \
    InnerProduct(dv,v2,v2); \
 \
    /* dv0 = (v2,v1) */ \
    T dv0; \
    InnerProduct(dv0,v2,v1); \
     \
    /* dvw0 = (w0,v1) */ \
    T dvw0; \
    InnerProduct(dvw0,w0,v1); \
 \
    /* w2 = v2 - (dv/dvw)*w1 - (dv0/dvw0)*w0 */ \
    w2 = v2; \
    div(t1,dv,dvw); \
    mul(tv1,t1,w1); \
    sub(w2,w2,tv1); \
    div(t1,dv0,dvw0); \
    mul(tv1,t1,w0); \
    sub(w2,w2,tv1); \
 \
    /* shift variables */ \
    w0=w1; \
    w1=w2; \
    v1=v2; \
 \
  } while (true); \
 \
  return true; \
} \


#endif
