#include <iostream>
#include <string.h>

#include "DLog_IC_Base.h"


/* Base class for index calculus like algorithms and the classic index
 * calculus technique.
 *
 * Written by: Chris Studholme
 * Copyright:  GPL (http://www.fsf.org/copyleft/gpl.html)
 */

NTL_START_IMPL;

#define IC_INITIAL_CACHE_SIZE 16

// sieve beyond relations needed (complete sieving) to ensure 
// proper optimization
//#define IC_EXTRA_SIEVING


/**************** class DLog_IC_Base ****************/

long DLog_IC_Base::VERBOSE=0;     /* verbosity level */


// test if g is a generator of the field ZZ_p
// f is the factorization of p-1
bool DLog_IC_Base::isGenerator(const ZZ_p& g, 
			       const vec_pair_ZZ_long& f) {
  ZZ m,e;
  ZZ_p t;
  sub(m,ZZ_p::modulus(),1);
  for (long i=0; i<f.length(); ++i) {
    div(e,m,f[i].a);
    power(t,g,e);
    if (IsOne(t))
      return false;
  }
  return true;
}

// complete logarithm
ZZ DLog_IC_Base::log_complete(const ZZ& pw, const ZZ& log_modq, 
			      const ZZ& q) {
  // other factor of p-1
  ZZ other;
  div(other,ZZ_p::modulus()-1,q);

  // simple method for small searches
  if (other<100) {
    ZZ lg(log_modq);
    ZZ_p t;
    do {
      power(t,g,lg);
      if (pw==rep(t)) {
	rem(lg,lg,ZZ_p::modulus()-1);
	return lg;
      }
      lg+=q;
      --other;
    } while (!IsZero(other));
    // logarithm doesn't exist
    return ZZ::zero();
  }

  // need to use Pollard Rho method here!
  return ZZ::zero();
}

// add logarithm to cache
void DLog_IC_Base::add_cache(const ZZ& prime, const ZZ& log) {
  if (ncache>=acache) {
    // reallocate cache
    acache = (acache==0 ? IC_INITIAL_CACHE_SIZE : 2*acache);
    ZZ* new_prime = new ZZ[acache];
    ZZ* new_log = new ZZ[acache];
    for (long i=0; i<ncache; ++i) {
      new_prime[i] = cache_prime[i];
      new_log[i] = cache_log[i];
    }
    if (cache_prime) delete[] cache_prime;
    if (cache_log) delete[] cache_log;
    cache_prime = new_prime;
    cache_log = new_log;
  }
  // binary search to find location to insert
  long low=0;
  long high=ncache;
  while (high>low) {
    long mid = (low+high)/2;
    if (cache_prime[mid]<prime)
      low=mid+1;
    else // (prime<=cache_prime[mid])
      high=mid;
  }
  if ((low<ncache)&&(cache_prime[low]==prime)) {
    // this should never happen
    std::cerr<<"DLog_IC_Base::add_cache() prime already in cache"<<std::endl;
    exit(1);
  }
  // make room for new entry
  for (long i=ncache; i>low; --i) {
    cache_prime[i] = cache_prime[i-1];
    cache_log[i] = cache_log[i-1];
  }
  // set new entry
  cache_prime[low] = prime;
  cache_log[low] = log;
  ++ncache;
}

// lookup log_g(prime) in cache
ZZ DLog_IC_Base::log_cache(const ZZ& prime) {
  if (ncache==0)
    return ZZ::zero();
  // binary search to find prime
  long low=0;
  long high=ncache;
  while (high>low) {
    long mid = (low+high)/2;
    if (cache_prime[mid]<prime)
      low=mid+1;
    else // (prime<=cache_prime[mid])
      high=mid;
  }
  if ((low<ncache)&&(cache_prime[low]==prime))
    return cache_log[low];
  return ZZ::zero();
}

// logarithm of any power (to base g)
ZZ DLog_IC_Base::log_g(const ZZ_p& pw, bool cacheonly) {
  // the answer (eventually)
  ZZ result;

  // multiply by g until we find a smooth product
  vec_pair_ZZ_long factors;
  ZZ_p sm(pw);
  do {
    factor(factors,rep(sm),upper_bound);
    // find largest prime factor that's not in the cache
    long ci=ncache-1;
    long fi=factors.length()-1;
    ZZ pf;
    while ((ci>=0)&&(fi>=0)) {
      pf = factors[fi].a;
      while ((ci>=0)&&(cache_prime[ci]>pf))
	--ci;
      if ((ci<0)||(cache_prime[ci]!=pf))
	break;
      --ci;
      --fi;
    }
    // if it's small enough, we can move on
    if (pf<=upper_bound)
      break;
    sm *= g;
    --result;

    if ((VERBOSE)&&(result<-1000)) {
      std::cout<<"DLog_IC_Base::log_g() searching... ("<<(-result)<<")  \r"<<std::flush;
    }

  } while (true);

  if ((VERBOSE)&&(result<-1000)) {
    std::cout<<std::endl;
  }

  // add up logarithms of factors
  for (long i=factors.length()-1; i>=0; --i) {
    ZZ l;
    l = log_cache(factors[i].a);
    if (IsZero(l)) {
      if (cacheonly)
	return ZZ::zero();
      l = log_prime(factors[i].a);
      if (IsZero(l))
	return ZZ::zero();
      add_cache(factors[i].a,l);
    }
    l*=factors[i].b;
    result += l;
  }

  // normalize result
  rem(result,result,ZZ_p::modulus()-1);
  return result;
}

// log of any power relative to base
ZZ DLog_IC_Base::log(const ZZ_p& pw) {
  if (IsZero(pw)) {
    std::cerr<<"DLog_IC_Base::log() log of 0 requested"<<std::endl;
    return to_ZZ(-1);
  }

  if (IsOne(pw))
    return ZZ::zero();

  if (IsZero(log_base)) {
    std::cerr<<"DLog_IC_Base::log() algorithm not initialized"<<std::endl;
    return to_ZZ(-2);
  }

  ZZ l;
  l = log_g(pw);
  if (IsZero(l)) {
    std::cerr<<"DLog_IC_Base::log() log_g() returned 0"<<std::endl;
    return to_ZZ(-1);
  }

  // logarithm is l/log_base (mod p-1)
  if (!IsOne(log_base)) {
    /* NOTE: log_base may not be invertable.  In that case:
     *   if GCD(log_base,modulus) | l,
     *   then logarithm exists and can be computed,
     *   otherwise, logarithm does not exist
     */
    ZZ d;
    ZZ lb(log_base);
    GCD(d,l,lb);
    if (!IsOne(d)) {
      l/=d;
      lb/=d;
    }
    GCD(d,lb,ZZ_p::modulus()-1);
    if (!IsOne(d)) {
      std::cerr<<"DLog_IC_Base::log() logarithm does not exist"<<std::endl;
      return to_ZZ(-1);
    }
    InvMod(lb,lb,ZZ_p::modulus()-1);
    MulMod(l,l,lb,ZZ_p::modulus()-1);
  }
  return l;
}


NTL_END_IMPL;
