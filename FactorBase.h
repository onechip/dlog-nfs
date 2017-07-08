#ifndef _FACTORBASE_H_
#define _FACTORBASE_H_

#include <NTL/ZZX.h>
#include <NTL/pair.h>
#include <NTL/vector.h>
#include <NTL/vec_vec_long.h>
#include "vec_short.h"


/* Class representing a factor base.
 *
 * Written by: Chris Studholme
 * Copyright:  GPL (http://www.fsf.org/copyleft/gpl.html)
 */

NTL_OPEN_NNS;

typedef Pair<long,long> pair_long_long;
typedef Vec<pair_long_long> vec_pair_long_long;

inline void append(vec_pair_long_long& x, long a, long b) {
  append(x,cons(a,b));
}

//void sieve1_times();

class FactorBase {
private:
  long *primes;  // array of primes starting with 2
  long aprimes;  // allocated size of primes
  long nprimes;  // number of primes
  long sprime;   // least index such that:
                 //   primes[sprime] > sqrt(primes[nprimes-1])
  
public:

  /* Set to half L2 cache size to optimize sieve.  Setting to 0 disables
   * cache optimization.
   */
  static long CACHE_SIZE;

  /* Create empty factorbase.  Call setBound() to create usable factorbase.
   */
  FactorBase();

  /* Create factorbase consisting of primes up to and including bound.
   */
  FactorBase(long bound);

  ~FactorBase();

  /* Create factorbase consisting of primes up to and including bound.
   */
  void setBound(long bound);

  /* Sieve over polynomial (deterministic).  The sieve is over the domain 
   * [start,start+smooth.length()).  Entries in smooth must be initialized
   * to zero for domain values you want checked and a non-zero value (such
   * as -1) for domain values you don't want checked.  When sieve() completes,
   * entries in smooth with contain either: -1 for values that were not
   * checked, 1 for smooth domain elements, and the non-smooth part of the
   * absolute value of the image of the polynomial in cases where the image
   * is not smooth.  In the case where f(a)=0, smooth[a-start]=0.
   */
  void sieve(vec_ZZ& smooth, const ZZX& polynomial, const ZZ& start) const;

  /* Sieve over polynomial (probabilistic).  The sieve is over the domain 
   * [start,start+smooth.length()).  Entries in smooth must be initialized
   * to zero for domain values you want checked and a negative value (such
   * as -1) for domain values you don't want checked.  When sieve() completes,
   * entries in smooth with contain either: -1 for values that were not
   * checked, 1 for smooth images, and zero for non-smooth images.
   * In the case where f(a)=0, smooth[a-start]=0.
   *
   * New parameters:
   *   rem_factor - maximum size of remaining factor (default 0)
   *   bound_low  - start sieve at this prime (default: 2)
   *   bound_high - end sieve at this prime (default: fb bound)
   *   exclude    - pairs (p,r) where p prime and r mod p are domain values
   *                  to avoid
   */
  void sieve(vec_short& smooth, const ZZX& f, const ZZ& start=ZZ::zero(),
	     const ZZ& rem_factor=ZZ::zero(),
	     long bound_low=0, long bound_high=0,
	     const vec_pair_long_long& exclude=vec_pair_long_long()) const;


  /* Sieve over a two variable polynomial (probabilistic).
   */
  //void sieve(vec_vec_long& smooth, const ZZX& f, 
  //     const ZZ& c_start=ZZ::zero(), const ZZ& d_start=ZZ::zero()) const;

  /* Reduce polynomial f by removing any smooth factor common to all
   * coefficients.  Result is g.
   */
  void reduce(ZZX& g, const ZZX& f) const;

  /* Factor n using trial division and return true if n is smooth.  f is the
   * exponents in the factorization of n and is expected to be of size length.
   * In the case that n is smooth, rem is set to 1 and the method returns true.
   * If n is not smooth, then f contains the exponents of the smooth part of
   * n, rem contains the non-smooth factor, and the method returns false.
   */
  bool factor(long* f, ZZ& rem, const ZZ& n) const;
  bool factor(long* f, long& rem, long n) const;

  /* Factor n using trial division and return true if n is smooth.  f is the
   * exponents in the factorization of n and is expected to be of size length.
   * Returns true if n is smooth, false otherwise.
   */
  inline bool factor(long* f, const ZZ& n) const {
    ZZ rem;
    return factor(f,rem,n);
  }
  inline bool factor(long* f, long n) const {
    long rem;
    return factor(f,rem,n);
  }

  // test if n is smooth using trial division
  bool isSmooth(const ZZ& n) const;

  // return true if p is a prime in this factorbase
  bool isPrime(long p) const;

  // returns largest prime in factorbase
  inline long bound() const {
    return primes[nprimes-1];
  }

  // number of primes in array of primes
  inline long length() const {
    return nprimes;
  }

  // a prime in the factorbase
  inline long operator[](long i) const {
    return primes[i];
  }

  // a prime in the factorbase
  inline long prime(long i) const {
    return primes[i];
  }

  // find prime, returns -1 if not found
  inline long find(long p) const {
    for (long i=0; i<nprimes; ++i)
      if (p==primes[i])
	return i;
    return -1;
  }

  // read-only reference to the array of primes
  const long* getPrimes() const {
    return primes;
  }

  // return a random prime in the factorbase
  long RandomPrime() const {
    return primes[RandomBnd(nprimes)];
  }

private:  
};

NTL_CLOSE_NNS;

#endif
