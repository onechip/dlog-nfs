#ifndef _NUMBERFIELDSIEVE_H_
#define _NUMBERFIELDSIEVE_H_

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/vec_ZZ.h>
#include "smat_ZZ.h"

#include "IndexCalculus.h"
#include "FactorBase.h"

// helper class defined internally
class NFS_Relations;


/* Solve discrete log problems using the Number Field Sieve method.  This
 * method is only applicable in the field ZZ_p.  
 *
 * Written by: Chris Studholme
 * Copyright:  GPL (http://www.fsf.org/copyleft/gpl.html)
 */
class DLog_NFS : public DLog_IC_Base {
public:
  // tuning parameters
  static long MAX_SIEVE;  /* maximum size of a single sieve */

protected:
  // degree of the number field
  long k;

  // largest prime factor of p-1
  ZZ q;

  // width and length of sieve
  long sieve_width;
  long sieve_length;  // set after each run

  // factorbase
  FactorBase zfb;

public:
  DLog_NFS() : DLog_IC_Base() {
  }
  DLog_NFS(const ZZ_p& base, long k=0, long bound=0, long width=0) 
    : DLog_IC_Base(base) {
    setBase(base,k,bound,width);
  }

  ~DLog_NFS() {
  }

  // set the base to use for logs, along with the degree of the number field,
  // the factorbase smoothness bound, and the width of the sieve
  void setBase(const ZZ_p& base, long k, long bound=0, long sieve_width=0);

  // set the base to use for logs
  // factorbase bound is determined automatically
  virtual void setBase(const ZZ_p& base) {
    // the other setBase() will figure out the optimal parameters
    setBase(base,0);
  }

  long getSieveLength() const {
    return sieve_length;
  }

  // compute optimum parameters for a given p
  // any of k, bound, or width that are non-zero are not computed
  static void parameters(long& k, long& bound, long& width, const ZZ& p);

private:
  // compute logarithm of a (smallish) prime to base g
  // return 0 if unsuccessful
  ZZ log_prime(const ZZ& power);

  // solve linear system induced by rels
  bool ls_solve(vec_ZZ& X, const NFS_Relations& rels, const vec_ZZ& b, 
		const ZZ& q);

  // L_p[1/3;c] = exp(c*log(p)^{1/3}*log(log(p))^{2/3})
  inline static double L_p(const ZZ& p, double c) {
    double lm = ::log(p);
    double llm = ::log(lm);
    double lllm = ::log(llm);
    return exp(c*exp(llm/3)*exp(lllm*2/3));
  }


};

#endif
