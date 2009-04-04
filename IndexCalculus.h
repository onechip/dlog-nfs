#ifndef _INDEXCALCULUS_H_
#define _INDEXCALCULUS_H_

#include <iostream>

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/vec_long.h>
#include <NTL/vec_ZZ.h>

#include "smat_long.h"

#include "DiscreteLog.h"
#include "FactorBase.h"


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


/* The classic Index Calculus method for computing discrete logarithms.
 */
class IndexCalculus : public DLog_IC_Base {
public:
  // tuning parameters (these are evaluated with L_p()
  static double L_UB;    /* upper bound */

protected:
  // largest factor of group order
  ZZ q;

  // length of sieve
  long sieve_length;

  // factorbase
  FactorBase zfb;

  // "upper" factorbase
  vec_ZZ ufb;
  ZZ sufb;   // start of upper fb

public:
  IndexCalculus() : DLog_IC_Base() {
  }
  IndexCalculus(const ZZ_p& base, long bound=0, long length=0) 
    : DLog_IC_Base(base) {
    setBase(base,bound,length);
  }

  ~IndexCalculus() {
  }

  // set the base to use for logs and factorbase bound
  void setBase(const ZZ_p& base, long bound, long sieve_length=0);

  // set the base to use for logs
  // factorbase bound is determined automatically
  void setBase(const ZZ_p& base) {
    // the other setBase() will figure out the optimal parameters
    setBase(base,0);
  }

  // compute optimum parameters for a given p
  // either of bound or width that are non-zero are not computed
  static void parameters(long& bound, long& length, const ZZ& p);


private:
  // make the system of equations to solve the factorbase
  bool make_system();

  // solve A*x=b (mod q)
  bool solve_system(const smat_long& A, vec_ZZ& x, const vec_long& b, 
		    const ZZ& q);

  // find log (to base g) of a prime
  virtual ZZ log_prime(const ZZ& prime);

  // L_p[1/2;c] = exp(c*sqrt(log(p)*log(log(p))))
  inline static double L_p(const ZZ& p, double c) {
    double lm = NTL::log(p);
    double llm = ::log(lm);
    return exp(c*sqrt(lm*llm));
  }
};

NTL_CLOSE_NNS;

#endif
