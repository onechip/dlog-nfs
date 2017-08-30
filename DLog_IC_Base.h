#ifndef _DLOG_IC_BASE_H_
#define _DLOG_IC_BASE_H_

#include <iostream>

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include "pair_ZZ_long.h"

#include "DiscreteLog.h"


NTL_OPEN_NNS;

/* Base class for methods of computing discrete logarithms based on the
 * Index Calculus method.  These include the classic Index Calculus method
 * and the Number Field Sieve method.  These methods are only applicable
 * to the field ZZ_p.
 *
 * Written by: Chris Studholme
 * Copyright:  GPL (http://www.fsf.org/copyleft/gpl.html)
 */
class DLog_IC_Base : public DiscreteLog {
public:
  static long VERBOSE;   /* verbosity level */
  

protected:
  // prime generator (mod p); internally all log are to this base
  ZZ_p g;

  // base of discrete logarithms
  ZZ_p base;  // hides DiscreteLog::base, but that's ok

  // discrete log of base (w.r.t. g)
  ZZ log_base;

  // bound on size of primes that log_prime() can handle
  ZZ upper_bound; 

  // cache of known logarithms
  ZZ* cache_prime; // a prime
  ZZ* cache_log;   // it's logarithm
  long acache;     // number of allocated entries
  long ncache;     // number of valid entries


public:
  DLog_IC_Base() : DiscreteLog() {
    log_base=0;
    cache_prime=NULL;
    cache_log=NULL;
    acache=ncache=0;
  }
  DLog_IC_Base(const ZZ_p& base) : DiscreteLog(base) {
    log_base=0;
    cache_prime=NULL;
    cache_log=NULL;
    acache=ncache=0;
    setBase(base);
  }

  ~DLog_IC_Base() {
    if (cache_prime) delete[] cache_prime;
    if (cache_log) delete[] cache_log;
  }

  // set the base to use for logs
  // factorbase bound is determined automatically
  //virtual void setBase(const ZZ_p& base);

  // returns -1 if the discrete log does not exist
  ZZ log(const GroupElement& power) {
    std::cerr<<"DLog_IC_Base::log() "
	     <<"cannot use index calculus method with an arbitrary group"
	     <<std::endl;
    return to_ZZ(-1);
  }

  // returns -1 if the discrete log does not exist
  ZZ log(const ZZ_p& power);

  // test if g is a generator of the field ZZ_p
  // f is the factorization of p-1
  static bool isGenerator(const ZZ_p& g, const vec_pair_ZZ_long& f);

protected:
  // complete log_g(pw) calculation (mod p-1)
  // log_modq = log_g(pw) mod q
  // it is assumed that GCD(q,(p-1)/q)==1
  ZZ log_complete(const ZZ& pw, const ZZ& log_modq, const ZZ& q);

  // compute log (to base g) of a prime
  virtual ZZ log_prime(const ZZ& prime)=0;

  // lookup log_g(prime) in cache
  // returns 0 if not found
  ZZ log_cache(const ZZ& prime);

  // find log (to base g) of power
  ZZ log_g(const ZZ_p& power, bool cacheonly=false);

  // add logarithm to cache
  void add_cache(const ZZ& prime, const ZZ& log);
};

NTL_CLOSE_NNS;

#endif
