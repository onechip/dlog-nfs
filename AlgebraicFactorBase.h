#ifndef _ALGEBRAICFACTORBASE_H_
#define _ALGEBRAICFACTORBASE_H_

#include <NTL/ZZX.h>
#include "vec_svec_long.h"
#include "FactorBase.h"

/* Class representing an algebraic factor base.
 *
 * Written by: Chris Studholme
 * Copyright:  GPL (http://www.fsf.org/copyleft/gpl.html)
 */

class AlgebraicFactorBase {
public:
  // minimal polynomial
  ZZX f;

  // rational factorbase
  FactorBase zfb;

  // indices into algebraic factorbase
  // index[i] is assoicated with zfb[i]
  vec_svec_long index;

  // number of good first degree prime ideals
  long nprimes;

  // modulus and exponent for character map calculation
  ZZ q,e;

public:
  AlgebraicFactorBase(const ZZX& f, long bound);

  ~AlgebraicFactorBase() {}

  // update index to reflect primes lying above zfb[i].
  void primes_above(long i);

  inline long length() const {
    return nprimes;
  }

  // factor c+d*alpha and return true if smooth
  bool factor(long* fact, long c, long d) const;
  
  // return norm of c+d*alpha
  ZZ norm(long c, long d) const;

  // modulus for character map calculation
  void SetModulus(const ZZ& q);

  // compute character map for c+d*alpha
  void CharacterMap(ZZ* l, long c, long d) const;
  
};


#endif
