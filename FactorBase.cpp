#include <iostream>
#include <fstream>

#include <limits.h>
#include <stdlib.h>
#include <string.h>

#include <NTL/ZZ.h>
#include <NTL/ZZXFactoring.h>
#include <NTL/lzz_pXFactoring.h>
#include <NTL/vec_vec_long.h>
#include <NTL/vec_double.h>

#include "FactorBase.h"


/* Class representing a factor base of integers.
 *
 * Written by: Chris Studholme
 * Copyright:  GPL (http://www.fsf.org/copyleft/gpl.html)
 */

NTL_START_IMPL;

class sieve_run {
public:
  long pi,e,start,inc;
};
typedef Vec<sieve_run> vec_run;

inline void append(vec_run& plan, long pi, long e, long start, long inc) {
  plan.SetLength(plan.length()+1);
  plan[plan.length()-1].pi = pi;
  plan[plan.length()-1].e = e;
  plan[plan.length()-1].start = start;
  plan[plan.length()-1].inc = inc;
}

inline bool IsOne(long a) {
  return a==1;
}

// if uniq==true, dups are removed
static int compar_long(const void *a, const void *b) {
  if (*(long*)a < *(long*)b)
    return -1;
  else if (*(long*)a == *(long*)b)
    return 0;
  return 1;
}

static void sort(vec_long& x, bool uniq=false) {
  long* v = x.elts();
  qsort(v,x.length(),sizeof(long),compar_long);
  if (uniq)
    for (long i=x.length()-1; i>0; --i)
      if (v[i-1]==v[i]) {
	memmove(v+i,v+i+1,(x.length()-i-1)*sizeof(long));
	x.SetLength(x.length()-1);
      }
}

inline long DivRem(long& q, long a, long b) {
  q = a/b;
  long r = a%b;
  while (r<0)
    r+=b;
  return r;
}

inline long rem(long a, long b) {
  long r = a%b;
  while (r<0)
    r+=b;
  return r;
}

// compute f(a)
inline void evaluate(ZZ& image, const ZZX& f, const ZZ& a) {
  image = coeff(f,deg(f));
  for (long i=deg(f)-1; i>=0; --i) {
    image*=a;
    image+=coeff(f,i);
  }
}

// compute fg = f(g(X))
inline void compose(ZZX& fg, const ZZX& f, const ZZX& g) {
  fg = coeff(f,deg(f));
  for (long i=deg(f)-1; i>=0; --i) {
    fg*=g;
    fg+=coeff(f,i);
  }
}


// if f(start) and f(end-1) have different sign, then
// either z in [start,end) with f(z)=0
// or z in (start,end) with f(z-1) and f(z) having different sign 
// is found and appended to zeros
static void AppendZeroCrossing(vec_ZZ& zeros, const ZZX& f, 
			       const ZZ& start, const ZZ& end) {
  if (start>=end)
    return;

  ZZ s(start),e(end);
  --e;
  ZZ fs,fe;
  evaluate(fs,f,start);
  evaluate(fe,f,end);
  if (fs>fe) {
    // swap start and end
    swap(s,e);
    swap(fs,fe);
  }

  if ((fs>0)||(fe<0))
    return;

  // at this point: fs<=0<=fe
  // binary search
  ZZ m,fm;
  while (abs(s-e)>1) {
    m = (s+e)/2;
    evaluate(fm,f,m);
    if (fm<0) {
      s=m;
      fs=fm;
    }
    else if (fm>0) {
      e=m;
      fe=fm;
    }
    else {
      append(zeros,m);
      return;
    }
  }
  append(zeros,s>e?s:e);
}


// find all zeros f(z)=0, where z in [start,end)
// furthermore, z in (start,end) such that f(z-1) and f(z) have different
// sign are also considered zeros
static void FindZeros(vec_ZZ& zeros, const ZZX& f,
		      const ZZ& start, const ZZ& end) {
  zeros.SetLength(0);

  // base cases
  if (deg(f)<=1) {
    if (IsZero(f))
      Error("FindZeros() zero polynomial provided");
    else if (deg(f)==1)
      AppendZeroCrossing(zeros,f,start,end);
    return;
  }

  ZZX df;
  diff(df,f);
  vec_ZZ extrema;
  FindZeros(extrema,df,start,end);

  if (extrema.length()==0) {
    AppendZeroCrossing(zeros,f,start,end);
    return;
  }

  AppendZeroCrossing(zeros,f,start,extrema[0]);
  for (long i=1; i<extrema.length(); ++i)
    AppendZeroCrossing(zeros,f,extrema[i-1],extrema[i]);
  AppendZeroCrossing(zeros,f,extrema[extrema.length()-1],end);
}



/**************** class FactorBase ****************/

long FactorBase::CACHE_SIZE=0;

bool FactorBase::factor(long* f, ZZ& rm, const ZZ& _n) const {
  ZZ n(_n);
  if (f)
    memset(f,0,nprimes*sizeof(long));
  rm=1;

  if (n<=1) {
    if (IsOne(n))
      return true;
    rm=n;
    return false;
  }

  /* First do trial division using all of the "small" primes in the factor
   * base.
   */
  ZZ q;
  for (long i=0; i<sprime; ++i)
    if (DivRem(q,n,primes[i])==0) {
      n=q;
      if (f) ++f[i];
      // check for repetative factors
      while (DivRem(q,n,primes[i])==0) {
        n=q;
        if (f) ++f[i];
      }
      if (IsOne(n))
	return true;
    }
  
  /* If n is less than or equal to the largest prime in the factor base,
   * then it must be a prime and we can just search for it.
   * Note:  perhaps we should just attempt to compute the index directly 
   * using pi(n)?
   */
  if (n<=primes[nprimes-1]) {
    for (long i=sprime; i<nprimes; ++i)
      if (n==primes[i]) {
	if (f) ++f[i];
	return true;
      }
    Error("FactorBase::factor() FATAL ERROR!");
  }

  // continue with trial division
  for (long i=sprime; i<nprimes; ++i)
    if (DivRem(q,n,primes[i])==0) {
      n=q;
      if (f) ++f[i];
      // check for repetative factors
      while (DivRem(q,n,primes[i])==0) {
	n=q;
	if (f) ++f[i];
      }
      if (n<=primes[nprimes-1]) {
	if (IsOne(n))
	  return true;
	// it's in the list, search for it
	for (long j=i+1; j<nprimes; ++j)
	  if (n==primes[j]) {
	    if (f) ++f[j];
	    return true;
	  }
	Error("FactorBase::factor() FATAL ERROR!");
      }
    }
  rm=n;
  return false;
}

bool FactorBase::factor(long* f, long& rm, long n) const {
  if (f)
    memset(f,0,nprimes*sizeof(long));
  rm=1;

  if (n<=1) {
    if (n==1)
      return true;
    rm=n;
    return false;
  }

  /* First do trial division using all of the "small" primes in the factor
   * base.
   */
  long q;
  for (long i=0; i<sprime; ++i)
    if (DivRem(q,n,primes[i])==0) {
      n=q;
      if (f) ++f[i];
      // check for repetative factors
      while (DivRem(q,n,primes[i])==0) {
        n=q;
        if (f) ++f[i];
      }
      if (IsOne(n))
	return true;
    }
  
  /* If n is less than or equal to the largest prime in the factor base,
   * then it must be a prime and we can just search for it.
   * Note:  perhaps we should just attempt to compute the index directly 
   * using pi(n)?
   */
  if (n<=primes[nprimes-1]) {
    for (long i=sprime; i<nprimes; ++i)
      if (n==primes[i]) {
	if (f) ++f[i];
	return true;
      }
    Error("FactorBase::factor() FATAL ERROR!");
  }

  // continue with trial division
  for (long i=sprime; i<nprimes; ++i)
    if (DivRem(q,n,primes[i])==0) {
      n=q;
      if (f) ++f[i];
      // check for repetative factors
      while (DivRem(q,n,primes[i])==0) {
	n=q;
	if (f) ++f[i];
      }
      if (n<=primes[nprimes-1]) {
	if (IsOne(n))
	  return true;
	// it's in the list, search for it
	for (long j=i+1; j<nprimes; ++j)
	  if (n==primes[j]) {
	    if (f) ++f[j];
	    return true;
	  }
	Error("FactorBase::factor() FATAL ERROR!");
      }
    }
  rm=n;
  return false;
}

// just check if n is smooth
bool FactorBase::isSmooth(const ZZ& _n) const {
  ZZ n(_n);
  if (n<=1)
    return IsOne(n);

  // do trial division
  ZZ q;
  for (long i=0; i<nprimes; ++i) {
    while (DivRem(q,n,primes[i])==0)
      n=q;
    if (n<=primes[nprimes-1])
      return true;
  }
  return false;
}

// remove smooth factor from all coefficients
void FactorBase::reduce(ZZX& g, const ZZX& f) const {
  g=f;
  if (IsZero(g))
    return;
  ZZ c;
  c = LeadCoeff(g);
  if (IsOne(c))
    return;
  for (long i=0; i<deg(g); ++i) {
    if (!IsZero(coeff(g,i))) {
      c = GCD(c,coeff(g,i));
      if (IsOne(c))
	return;
    }
  }
  // c>1 is GCD of all non-zero coefficients
  if (c>primes[nprimes-1]) {
    ZZ rm;
    factor(0,rm,c);
    c /= rm;
  }
  g /= c;
}


// the zz_p modulus is expected to be set to p
// actually, start+prev_zero is a zero of f_{p^e}(X)
void sieve_r(vec_ZZ& smooth, const ZZ& start, const ZZX& f, 
	     long e=0, long prev_zero=0) {
  long p = zz_p::modulus();

  long length = smooth.length();
  vec_long zeros;

  if (e==0) {
    // reduce coefficients of f modulo p
    zz_pX fp;
    conv(fp,f);

    // factor fp
    vec_pair_zz_pX_long u;
    CanZass(u,fp);

    // note all zeros
    for (long j=0; j<u.length(); ++j)
      if (deg(u[j].a)==1)
	append(zeros,rem(rep(-ConstTerm(u[j].a))-start,p));
  }

  // e>0
  else {
    ZZ pe;
    power(pe,p,e);

    // polynomial: start+prev_zero + Y*p^e
    ZZX g;
    SetCoeff(g,0,start+prev_zero);
    SetCoeff(g,1,pe);
    
    // compose f and g
    ZZX fg;
    compose(fg,f,g);
    
    // all coefficients of fg must be divisible by pe
    fg/=pe;

    // reduce coefficients of fg modulo p
    zz_pX fp;
    conv(fp,fg);

    if (IsZero(fp)) {
      // all values [0,p) are zeros
      ZZ r;
      for (long j=0; j<p; ++j) {
	rem(r,prev_zero+pe*j,p*pe);
	if (r<length)
	  append(zeros,to_long(r));
      }
    }

    else if (deg(fp)<1) {
      // no zeros
      return;
    }

    else if (deg(fp)==1) {
      // exactly one zero
      MakeMonic(fp);
      ZZ r;
      rem(r,prev_zero-pe*rep(ConstTerm(fp)),p*pe);
      if (r<length)
	append(zeros,to_long(r));
    }

    else {
      // deg>1, does this ever happen?
      MakeMonic(fp);
      
      // factor fp
      vec_pair_zz_pX_long u;
      CanZass(u,fp);
      
      // note all zeros
      ZZ r;
      for (long j=0; j<u.length(); ++j)
	if (deg(u[j].a)==1) {
	  rem(r,prev_zero-pe*rep(ConstTerm(u[j].a)),p*pe);
	  if (r<length)
	    append(zeros,to_long(r));
	}
    }
  }

  ++e;
  ZZ pe;
  power(pe,p,e);
  long inc = pe<length ? to_long(pe) : length;
  ZZ q;
  for (long j=0; j<zeros.length(); ++j) {
    bool done=true;
    ZZ* fv = smooth.elts();
    for (long l=zeros[j]; l<length; l+=inc)
      if (fv[l]>1) {
	fv[l]/=p;
	done=false;
      }
    if (!done)
      sieve_r(smooth,start,f,e,zeros[j]); // recurse
  }
}

// sieve over a polynomial (deterministic)
void FactorBase::sieve(vec_ZZ& smooth, const ZZX& f, const ZZ& start) const {
  zz_pBak bak;
  bak.save();

  long length=smooth.length();
  ZZ* fv = smooth.elts();

  // efficient evaluation of f using successor functions
  long k = deg(f);
  ZZX g[k+1];
  ZZ pv[k+1]; // previous value

  g[k] = f;
  evaluate(pv[k],g[k],start);
  for (long i=k-1; i>=0; --i) {
    // compute successor function
    ZZX t;
    SetCoeff(t,0,1);
    SetCoeff(t,1,1);
    compose(g[i],g[i+1],t);
    g[i]-=g[i+1];
    evaluate(pv[i],g[i],start);
  }
  
  // image of f
  for (long i=0; ; ) {
    if (IsZero(fv[i]))
      abs(fv[i],pv[k]);
    else
      conv(fv[i],-1);
    if (++i>=length)
      break;
    // update function
    for (long j=k; j>0; --j)
      pv[j]+=pv[j-1];
  }

  for (long i=0; i<nprimes; ++i) {
    long p=primes[i];
    zz_p::init(p);
    // recursive function to do all the work
    sieve_r(smooth,start,f);
  }
}

std::ostream* s1debug=0;

double s1_plan=0;
double s1_exec=0;
double s1_fin=0;

// write out times
void sieve1_times() {
  std::cout<<"sieve1_plan: "<<s1_plan<<" seconds"<<std::endl;
  std::cout<<"sieve1_exec: "<<s1_exec<<" seconds"<<std::endl;
  std::cout<<"sieve1_done: "<<s1_fin<<" seconds"<<std::endl;
}

// the zz_p modulus is expected to be set to p
// actually, start+prev_zero is a zero of f_{p^e}(X)
long sieve_r(vec_run& plan, const ZZ& start, long length, const ZZX& f,
	     long pi, long pos, long e=0) {

  long p = zz_p::modulus();
  ++e;
  ZZ pe;
  power(pe,p,e);
  
  // polynomial: start+pos + Y*p^e
  ZZX g;
  SetCoeff(g,0,start+pos);
  SetCoeff(g,1,pe);
  
  // compose f and g
  ZZX fg;
  compose(fg,f,g);

  //if (s1debug)
  //  *s1debug<<' '<<fg;
  
  // all coefficients of fg must be divisible by pe
  fg/=pe;
  
  // reduce coefficients of fg modulo p
  zz_pX fp;
  conv(fp,fg);

  if (s1debug)
    *s1debug<<"e="<<e<<' '<<fp<<std::endl;
  
  if (deg(fp)>1)
    Error("FactorBase::sieve() FATAL ERROR! poly has degree > 1");

  long depth=0;

  if (deg(fp)==1) {
    // exactly one zero
    MakeMonic(fp);
    long ofs = rep(-ConstTerm(fp));
    if (pos+pe*ofs<length) {
      long inc = pe*p<length ? to_long(pe*p) : length;
      append(plan,pi,e,pos+to_long(pe*ofs),inc);
      depth = 1+sieve_r(plan,start,length,f,pi,pos+to_long(pe*ofs),e);
    }
  }
  
  else if (IsZero(fp)) {
    // all values [0,p) are zeros
    long inc = pe<length ? to_long(pe) : length;
    append(plan,pi,e,pos,inc);
    ++depth;
    for (long j=0; j<p && pos+pe*j<length; ++j) {
      long d = 1+sieve_r(plan,start,length,f,pi,pos+to_long(pe*j),e);
      if (depth<d) depth=d;
    }
  }

  return depth;
}

// the zz_p modulus is expected to be set to p
// actually, start+prev_zero is a zero of f_{p^e}(X)
long sieve_r(vec_run& plan, const ZZ& start, long length, const ZZX& f,
	     long pi, const vec_long& exclude) {

  // reduce coefficients of f modulo p
  long p = zz_p::modulus();
  zz_pX fp;
  conv(fp,f);

  long depth=0;
  
  if (deg(fp)<=0) {
    if (IsZero(fp)) {
      // all values [0,p) are zeros
      append(plan,pi,0,0,1);
      ++depth;
      long xi=0;
      for (long j=0; j<p && j<length; ++j) {
	if (xi>=exclude.length() || j!=exclude[xi]) {
	  long d = 1+sieve_r(plan,start,length,f,pi,j);
	  if (depth<d) depth=d;
	}
	else
	  ++xi;
      }
    }
  }
  
  else if (deg(fp)==1) {
    // exactly one zero
    MakeMonic(fp);
    long ofs = rem(rep(-ConstTerm(fp))-start,p);
    if (ofs<length) {
      for (long xi=0; xi<exclude.length(); ++xi)
	if (exclude[xi]==ofs)
	  return depth;
      append(plan,pi,0,ofs,p);
      depth = 1+sieve_r(plan,start,length,f,pi,ofs);
    }
  }
  
  else {
    // factor fp
    vec_pair_zz_pX_long u;
    MakeMonic(fp);
    CanZass(u,fp);
    // note all zeros
    for (long j=0; j<u.length(); ++j)
      if (deg(u[j].a)==1) {
	long ofs = rem(rep(-ConstTerm(u[j].a))-start,p);
	bool bad=false;
	for (long xi=0; xi<exclude.length(); ++xi)
	  if (exclude[xi]==ofs) {
	    bad=true;
	    break;
	  }
	if (ofs<length && !bad) {
	  append(plan,pi,0,ofs,p);
	  long d = 1+sieve_r(plan,start,length,f,pi,ofs);
	  if (depth<d) depth=d;
	}
      }
  }

  return depth;
}

// sieve over a polynomial (probabilistic)
void FactorBase::sieve(vec_short& smooth, const ZZX& f, const ZZ& start,
		       const ZZ& rem_factor, long bound_low, long bound_high,
		       const vec_pair_long_long& exclude) const {
  zz_pBak bak;
  bak.save();

  double t_start = GetTime();

  if (bound_high<=0)
    bound_high = bound();

  long length=smooth.length();

  //if (length==500)
  //s1debug = new std::fstream("sieve1-debug.txt",std::ios::out);

  // find maximum value f can take on over domain
  ZZ max;
  evaluate(max,f,start);
  abs(max,max);
  ZZ t;
  evaluate(t,f,start+length-1);
  abs(t,t);
  if (max<t) max=t;
  ZZX df;
  diff(df,f);
  vec_ZZ extrema;
  FindZeros(extrema,df,start,start+length);

  for (long i=0; i<extrema.length(); ++i) {
    // should we also check extrema[i]-1 ?? 
    evaluate(t,f,extrema[i]);
    abs(t,t);
    if (max<t) max=t;
  }
  
  // compute log multiplier and log of primes
  //double logmult = (NTL_MAX_LONG/2)/log(max);
  //double logmult = (1<<(8*sizeof(short)-2))/log(max);
  //double logmult = (SHRT_MAX-nprimes)/log(max);
  double logmult = (SHRT_MAX/2)/log(max);

  vec_vec_long logpe;
  logpe.SetLength(nprimes);

  if (s1debug)
    *s1debug<<"logmult = "<<logmult<<std::endl;

  // sort excludes
  vec_vec_long exclude_list;
  exclude_list.SetLength(nprimes);
  for (long i=0; i<exclude.length(); ++i) {
    long j = find(exclude[i].a);
    if (j>=0)
      append(exclude_list[j],rem(exclude[i].b-start,exclude[i].a));
  }

  // create plan
  vec_run plan;
  for (long i=0; i<nprimes; ++i) 
    if (bound_low<=primes[i] && primes[i]<=bound_high) {
      sort(exclude_list[i],true);
      zz_p::init(primes[i]);
      long depth = sieve_r(plan,start,length,f,i,exclude_list[i]);
      if (depth>0) {
	double logp = logmult*log((double)primes[i]);
	logpe[i].SetLength(depth);
	double cur = 0;
	for (long j=0,sum=0; j<depth; ++j) {
	  cur += logp;
	  logpe[i][j] = (long)ceil(cur) - sum;
	  sum += logpe[i][j];
	}
      }
    }
  
  //std::cout<<"plan size: "<<plan.length()<<"    "<<std::endl;
  
  double t_end = GetTime();
  s1_plan += t_end-t_start;
  t_start=t_end;

  // execute plan
  short* fv = smooth.elts();
  long stride = CACHE_SIZE/sizeof(short);
  if (stride<8) {
    for (long i=0; i<plan.length(); ++i) {
      short logp = logpe[plan[i].pi][plan[i].e];
      for (long j=plan[i].start; j<length; j+=plan[i].inc)
	if (fv[j]>=0)
	  fv[j] += logp;
    }
  }
  else {
    long end = 0;
    do {
      end += stride;
      if (end>length)
	end = length;
      for (long i=0; i<plan.length(); ++i) {
	short logp = logpe[plan[i].pi][plan[i].e];
	while (plan[i].start<end) {
	  if (fv[plan[i].start]>=0)
	    fv[plan[i].start] += logp;
	  plan[i].start += plan[i].inc;
	}
      }
    } while (end<length);
  }
  
  // add in remaining factor
  if (rem_factor>0) {
    long lr;
    if (rem_factor<max)
      lr = (long)ceil(logmult*log(rem_factor));  // round up
    else
      lr = (long)ceil(logmult*log(max));  // round up
    for (long i=0; i<length; ++i)
      if (fv[i]>=0)
        fv[i] += lr;
  }

  t_end = GetTime();
  s1_exec += t_end-t_start;
  t_start=t_end;

  // now, figure out which ones are smooth by checking if their logarithms
  // are sufficiently large

  // find all zeros of f
  vec_ZZ zeros;
  FindZeros(zeros,f,start,start+length);
  
  // make list of endpoints consisting of start, end, extrema, and zeros
  vec_ZZ endpoints;
  append(endpoints,start);
  for (long i=0,j=0; ; ) {
    if (i<zeros.length()) {
      if (j<extrema.length()) {
	// figure out which is smaller: zeros[i] or extrema[j]
	if (zeros[i]<extrema[j])
	  append(endpoints,zeros[i++]);
	else
	  append(endpoints,extrema[j++]);
      }
      else
	append(endpoints,zeros[i++]);
    }
    else if (j<extrema.length())
      append(endpoints,extrema[j++]);
    else
      break;
  }
  append(endpoints,start+length);

  // check each region between endpoints
  for (long i=1; i<endpoints.length(); ++i) {
    ZZ p1(endpoints[i-1]);
    ZZ p2(endpoints[i]);
    if (p1>=p2)
      continue;
    --p2;
    ZZ fp1,fp2;
    evaluate(fp1,f,p1);
    evaluate(fp2,f,p2);
    abs(fp1,fp1);
    abs(fp2,fp2);
    long dir=1;
    if (fp1>fp2) {
      swap(p1,p2);
      swap(fp1,fp2);
      dir=-1;
    }

    // since fp1<=fp2, we can use log(f(x-dir)) as lower bound for log(f(x))

    // look for smooth elements
    long j = to_long(p1-start);
    long lastlog = IsZero(fp1) ? 0 : (long)(logmult*log(fp1));
    do {
      if (fv[j]>0) {
	if (fv[j]<lastlog)
	  fv[j]=0;
	else {
	  evaluate(fp1,f,p1);
	  abs(fp1,fp1);
	  if (IsZero(fp1))
	    fv[j]=0;
	  else {
	    lastlog = (long)(logmult*log(fp1));
	    if (fp1<=primes[nprimes-1])  // NOT NEEDED??
	      fv[j]=1;
	    else if (fv[j]<lastlog)
	      fv[j]=0;
	    else
	      fv[j]=1;
	  }
	}
      }
      if (p1==p2)
	break;
      p1+=dir;
      j+=dir;
    } while (true);
  }

  t_end = GetTime();
  s1_fin += t_end-t_start;

  delete s1debug;
  s1debug = 0;
}

/*
// sieve over a two variable polynomial (probabilistic)
void FactorBase::sieve(vec_vec_long& smooth, const ZZX& f, 
		       const ZZ& c_start, const ZZ& d_start) const {

  // non-optimal implementation (to say the least)
  ZZ d(d_start);
  ZZ dd;
  ZZX fd;
  SetCoeff(fd,deg(f),coeff(f,deg(f)));
  for (long i=0; i<smooth.length(); ++i, ++d) {
    dd=d;
    for (long j=deg(f)-1; ; ) {
      SetCoeff(fd,j,dd*coeff(f,j));
      if (--j<0)
	break;
      dd*=d;
    }
    sieve(smooth[i],fd,c_start);
  }
}
*/

bool FactorBase::isPrime(long p) const {
  // binary search to find prime
  long low=0;
  long high=nprimes;
  while (high>low) {
    long mid = (low+high)/2;
    if (primes[mid]<p)
      low=mid+1;
    else // (p<=primes[mid])
      high=mid;
  }
  return ((low<nprimes)&&(primes[low]==p));
}

FactorBase::FactorBase() {
  primes=NULL;
  aprimes=nprimes=sprime=0;
}

FactorBase::FactorBase(long bound) {
  primes=NULL;
  setBound(bound);
}

FactorBase::~FactorBase() {
  if (primes) delete[] primes;
}

void FactorBase::setBound(long bound) {
  unsigned char* sieve = new unsigned char[bound+1];  // 0 => composite
  memset(sieve,0,(bound+1)*sizeof(unsigned char));
  
  long end = (long)ceil(sqrt((double)bound));
  for (long i=2; i<=end; ++i)
    if (sieve[i]==0)
      for (long j=i*i; j<=bound; j+=i)
	sieve[j]=1;
  
  aprimes=0;
  for (long i=2; i<=bound; ++i)
    if (sieve[i]==0)
      ++aprimes;

  if (primes) delete[] primes;
  primes = new long[aprimes];
  nprimes = 0;

  for (long i=2; i<=bound; ++i)
    if (sieve[i]==0)
      primes[nprimes++] = i;

  delete[] sieve;

  // figure out sprime
  sprime=0;
  while (primes[sprime]<sqrt((double)primes[nprimes-1]))
    ++sprime;
}



NTL_END_IMPL;

