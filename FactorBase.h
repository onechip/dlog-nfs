#ifndef _FACTORBASE_H_
#define _FACTORBASE_H_

#include <NTL/ZZX.h>
#include <NTL/vec_vec_long.h>

/* Class representing a factor base.
 *
 * Written by: Chris Studholme
 * Copyright:  GPL (http://www.fsf.org/copyleft/gpl.html)
 */

class FactorBase {
private:
  long *primes;  // array of primes starting with 2
  long aprimes;  // allocated size of primes
  long nprimes;  // number of primes
  long sprime;   // least index such that:
                 //   primes[sprime] > sqrt(primes[nprimes-1])
  
public:
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
   */
  void sieve(vec_long& smooth, const ZZX& f, const ZZ& start=ZZ::zero()) const;

  /* Sieve over a two variable polynomial (probabilistic).
   */
  void sieve(vec_vec_long& smooth, const ZZX& f, 
	     const ZZ& c_start=ZZ::zero(), const ZZ& d_start=ZZ::zero()) const;

  /* Factor n using trial division and return true if n is smooth.  f is the
   * exponents in the factorization of n and is expected to be of size length.
   * In the case that n is smooth, rem is set to 1 and the method returns true.
   * If n is not smooth, then f contains the exponents of the smooth part of
   * n, rem contains the non-smooth factor, and the method returns false.
   */
  bool factor(long* f, ZZ& rem, const ZZ& n) const;

  /* Factor n using trial division and return true if n is smooth.  f is the
   * exponents in the factorization of n and is expected to be of size length.
   * Returns true if n is smooth, false otherwise.
   */
  inline bool factor(long* f, const ZZ& n) const {
    ZZ rem;
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

#endif
