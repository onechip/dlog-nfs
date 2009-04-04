
#include <NTL/ZZ_p.h>
#include "ZZFactoring.h"


/* Method for factoring integers.
 *
 * Written by: Chris Studholme
 * Copyright:  GPL (http://www.fsf.org/copyleft/gpl.html)
 */

NTL_START_IMPL;


long ProvePrime(const ZZ& _n) {
  ZZ n(_n);
  if (n<0)
    abs(n,n);
  if (n<=1)
    return 0;

  if (n<=1000000) {
    // n is small so use trial division to check primality
    long ln = to_long(n);
    long end = to_long(SqrRoot(n));
    PrimeSeq s;
    for (long p=s.next(); p<=end; p=s.next())
      if ((ln%p)==0)
	return 0;
    return 1;
  }

  // check small primes
  PrimeSeq s;
  for (long p=s.next(); p<1000; p=s.next())
    if (divide(n,p))
      return 0;

  // obviously, something is missing here!

  return ProbPrime(n);
}

// factors is kept sorted by p
void addFactor(vec_pair_ZZ_long& factors, const ZZ& p, long exponent=1) {
  // fast path: factors.length()==0
  if (factors.length()==0) {
    factors.SetLength(1);
    factors[0].a = p;
    factors[0].b = exponent;
    return;
  }

  // fast path: p>=factors[factors.length()-1].a
  if (p>=factors[factors.length()-1].a) {
    if (p==factors[factors.length()-1].a)
      factors[factors.length()-1].b += exponent;
    else {
      factors.SetLength(factors.length()+1);
      factors[factors.length()-1].a = p;
      factors[factors.length()-1].b = exponent;
    }
    return;
  }

  // binary search to find location to insert
  long low=0;
  long high=factors.length();
  while (high>low) {
    long mid = (low+high)/2;
    if (factors[mid].a<p)
      low=mid+1;
    else // (p<=factors[mid].a)
      high=mid;
  }
  if ((low<factors.length())&&(factors[low].a==p))
    factors[low].b += exponent;
  else {
    // insert factor
    factors.SetLength(factors.length()+1);
    for (long i=factors.length()-1; i>low; --i)
      factors[i] = factors[i-1];
    factors[low].a = p;
    factors[low].b = exponent;
  }
}

// factor n into a*b using Pollard Rho method
// pre-condition: n>1
inline void PollardRho(ZZ& a, ZZ& b, const ZZ& n, 
		       const ZZ& bnd=ZZ::zero()) {
  ZZ_pBak bak;
  bak.save();
  ZZ_p::init(n);

  ZZ d;
  ZZ_p x1;
  random(x1);
  ZZ_p x2(x1);

  ZZ end(IsZero(bnd)?5*SqrRoot(SqrRoot(n)):2*SqrRoot(bnd));
  for (; !IsZero(end); --end) {
    x1 = x1*x1 + 1;
    x2 = x2*x2 + 1;
    x2 = x2*x2 + 1;
    GCD(d,n,rep(x2-x1));
    if ((d>1)&&(d<n)) {
      a=d; b=n/d;
      return;
    }
  }
  // failure
  a=1; b=n;
}

// recursively factor n and add results to factors
// pre-condition: factors is of zero-length or 
//                a properly sorted list of factors
// pre-condition: n>1
void PollardRho(vec_pair_ZZ_long& factors, const ZZ& n,
		const ZZ& _bnd, long deterministic, long verbose) {
  bool pn = deterministic ? ProvePrime(n) : ProbPrime(n);
  if (pn) {
    addFactor(factors,n);
    return;
  }

  ZZ bnd(n<=_bnd ? ZZ::zero() : _bnd);

  ZZ a,b;
  do {
    PollardRho(a,b,n,bnd);
    if (!IsOne(a)&&!IsOne(b))
      break;
    if (!deterministic||!IsZero(bnd)) {
      addFactor(factors,n);
      return;
    }
  } while (true);

  PollardRho(factors,a,bnd,deterministic,verbose);
  PollardRho(factors,b,bnd,deterministic,verbose);
}

// do trial division by small primes
// pre-conditions: n>0
ZZ SmallPrimes(vec_pair_ZZ_long& factors, const ZZ& _n) {
  ZZ n(_n);
  ZZ q;
  PrimeSeq s;
  long p;

  //long prime_bnd = ComputePrimeBound(NumBits(n));
  long prime_bnd = SqrRoot(n)>1000 ? 1000 : to_long(SqrRoot(n));

  for (p=s.next(); (p>0)&&(p<=prime_bnd); p=s.next()) {
    if (DivRem(q,n,p)==0) {
      long e=1;
      n=q;
      while (DivRem(q,n,p)==0) {
	++e;
	n=q;
      }
      addFactor(factors,to_ZZ(p),e);
      if (IsOne(n))
	return n;
    }
  }
  return n;
}

// general purpose factoring method
void factor(vec_pair_ZZ_long& factors, const ZZ& _n,
            const ZZ& _bnd, long deterministic, long verbose) {
  // initialize factors
  factors.SetLength(0);

  ZZ n(_n);

  if (n<=1) {
    abs(n,n);
    if (n<=1)
      return;
  }

  ZZ bnd(_bnd>0 ? _bnd : ZZ::zero());

  // small primes test
  n = SmallPrimes(factors,n);

  // Pollard Rho method
  if (n>1)
    PollardRho(factors,n,bnd,deterministic,verbose);
}


NTL_END_IMPL;
