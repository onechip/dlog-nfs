#include <stdio.h>
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

// estimate natural logarithm of a
// it will always be the case that log_est(a)<=log(a)
//inline double log_est(const ZZ& a) {
//  return (NumBits(a)-1)*M_LN2;
//}

// if f(start) and f(end) have different sign, the value of z in [start,end]
// that has f(z) closest to zero is found and appended to zeros
void AppendZeroCrossing(vec_ZZ& zeros, const ZZX& f, 
			const ZZ& start, const ZZ& end) {
  ZZ s(start),e(end);
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
  if ((-fs)<fe)
    append(zeros,s);
  else
    append(zeros,e);
}

// find all (simple) real zeros in the interval [start,end]
void FindZeros(vec_ZZ& zeros, const ZZX& f, const ZZ& start, const ZZ& end) {
  // base cases
  if (deg(f)<=1) {
    zeros.SetLength(0);
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

  zeros.SetLength(0);
  if (extrema.length()==0) {
    AppendZeroCrossing(zeros,f,start,end);
    return;
  }

  AppendZeroCrossing(zeros,f,start,extrema[0]);
  for (long i=1; i<extrema.length(); ++i)
    AppendZeroCrossing(zeros,f,extrema[i-1],extrema[i]);
  AppendZeroCrossing(zeros,f,extrema[extrema.length()-1],end);

  // check for duplicate zeros
  for (long i=zeros.length()-1; i>0; --i) 
    if (zeros[i]==zeros[i-1]) {
      // remove duplicate
      for (long j=i; j<zeros.length()-1; ++j)
	zeros[j]=zeros[j+1];
      zeros.SetLength(zeros.length()-1);
    }
}



/**************** class FactorBase ****************/

bool FactorBase::factor(long* f, ZZ& rm, const ZZ& _n) const {
  ZZ n(_n);
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
      ++f[i];
      // check for repetative factors
      while (DivRem(q,n,primes[i])==0) {
        n=q;
        ++f[i];
      }
    }

  if (IsOne(n))
    return true;
  
  /* If n is less than or equal to the largest prime in the factor base,
   * then it must be a prime and we can just search for it.
   * Note:  perhaps we should just attempt to compute the index directly 
   * using pi(n)?
   */
  if (n<=primes[nprimes-1]) {
    for (long i=sprime; i<nprimes; ++i)
      if (n==primes[i]) {
	++f[i];
	return true;
      }
    cerr<<"FactorBase::factor() FATAL ERROR!\n";
    exit(1);
  }

  // continue with trial division
  for (long i=sprime; i<nprimes; ++i)
    if (DivRem(q,n,primes[i])==0) {
      n=q;
      ++f[i];
      // check for repetative factors
      while (DivRem(q,n,primes[i])==0) {
	n=q;
	++f[i];
      }
      if (n<=primes[nprimes-1]) {
	if (IsOne(n))
	  return true;
	// it's in the list, search for it
	for (long j=i+1; j<nprimes; ++j)
	  if (n==primes[j]) {
	    ++f[j];
	    return true;
	  }
	cerr<<"FactorBase::factor() FATAL ERROR!\n";
	exit(1);
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

// the zz_p modulus is expected to be set to p
// actually, start+prev_zero is a zero of f_{p^e}(X)
void sieve_r(vec_long& smooth, const ZZ& start, const ZZX& f, long logp,
	     long e=0, long prev_zero=0) {
  long p = zz_p::modulus();

  long length = smooth.length();
  vec_long zeros;

  if (e==0) {
    // reduce coefficients of f modulo p
    zz_pX fp;
    conv(fp,f);

    if (IsZero(fp)) {
      // all values [0,p) are zeros
      zeros.SetMaxLength(p);
      for (long j=0; j<p; ++j)
	append(zeros,rem(j-start,p));
    }

    else if (deg(fp)<1) {
      // no zeros
      return;
    }

    else if (deg(fp)==1) {
      // exactly one zero
      MakeMonic(fp);
      append(zeros,rem(rep(-ConstTerm(fp))-start,p));
    }

    else {
      // factor fp
      vec_pair_zz_pX_long u;
      MakeMonic(fp);
      CanZass(u,fp);
      
      // note all zeros
      zeros.SetMaxLength(u.length());
      for (long j=0; j<u.length(); ++j)
	if (deg(u[j].a)==1)
	  append(zeros,rem(rep(-ConstTerm(u[j].a))-start,p));
    }
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
      zeros.SetMaxLength(p);
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

    else { // deg(fp)>1
      // does this ever happen?  if not, why not?
      cout<<"sieve() polynomial for "<<p<<"^"<<(e+1)<<" has degree "
	  <<deg(fp)<<"\n";
      
      // factor fp
      vec_pair_zz_pX_long u;
      MakeMonic(fp);
      CanZass(u,fp);
      
      // note all zeros
      ZZ r;
      zeros.SetMaxLength(u.length());
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
    long* fv = smooth.elts();
    for (long l=zeros[j]; l<length; l+=inc)
      if (fv[l]>=0) {
	fv[l]+=logp;
	done=false;
      }
    if (!done)
      sieve_r(smooth,start,f,logp,e,zeros[j]); // recurse
  }
}


// sieve over a polynomial (probabilistic)
void FactorBase::sieve(vec_long& smooth, const ZZX& f, const ZZ& start) const {
  zz_pBak bak;
  bak.save();

  long length=smooth.length();

  // find maximum value f can take on over domain
  ZZ max;
  clear(max);
  ZZ t;
  evaluate(t,f,start);
  abs(t,t);
  if (max<t) max=t;
  evaluate(t,f,start+length-1);
  abs(t,t);
  if (max<t) max=t;
  
  // find extrema of f
  ZZX df;
  diff(df,f);
  vec_ZZ extrema;
  FindZeros(extrema,df,start,start+length-1);
  for (long i=0; i<extrema.length(); ++i) {
    evaluate(t,f,extrema[i]);
    abs(t,t);
    if (max<t) max=t;
  }
  
  // compute log multiplier
  double logmult = (NTL_MAX_LONG/2)/log(max);

  //double t2 = GetTime();

  // tally of logarithms of factors
  for (long i=0; i<nprimes; ++i) {
    long p=primes[i];
    zz_p::init(p);
    long logp=1+(long)(logmult*log(p));  // always round up
    // recursive function to do all the work
    sieve_r(smooth,start,f,logp);
  }

  //double t3 = GetTime();

  // now, figure out which ones are smooth by making sure their logarithms
  // are sufficiently large

  // find all zeros of f
  vec_ZZ zeros;
  FindZeros(zeros,f,start,start+length-1);

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
  append(endpoints,start+length-1);

  // check each region between endpoints
  for (long i=1; i<endpoints.length(); ++i) {
    ZZ p1(endpoints[i-1]);
    ZZ p2(endpoints[i]);
    if (p1==p2)
      continue;
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

    // look for smooth elements
    long* fj = smooth.elts();
    long j = to_long(p1-start);
    long lastlog = IsZero(fp1) ? 0 : (long)(logmult*log(fp1));
    do {
      if (fj[j]>0) {
	if (fj[j]<lastlog)
	  fj[j]=0;
	else {
	  evaluate(fp1,f,p1);
	  abs(fp1,fp1);
	  if (IsZero(fp1))
	    fj[j]=0;
	  else {
	    lastlog = (long)(logmult*log(fp1));
	    if (fp1<=primes[nprimes-1])
	      fj[j]=1;
	    else if (fj[j]<lastlog)
	      fj[j]=0;
	    else
	      fj[j]=1;
	  }
	}
      }
      if (p1==p2)
	break;
      p1+=dir;
      j+=dir;
    } while (true);
  }

  /*
  double t4 = GetTime();
  cout<<"sieve() profiling: "
      <<(t3-t2)<<" "
      <<(t4-t3)<<"           \n";
  */
}

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
  
  long end = (long)ceil(sqrt(bound));
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
  while (primes[sprime]<sqrt(primes[nprimes-1]))
    ++sprime;
}
