
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/lzz_pXFactoring.h>

#include "AlgebraicFactorBase.h"

//#define AFB_VERBOSE

/* Class representing an algebraic factor base.
 *
 * Written by: Chris Studholme
 * Copyright:  GPL (http://www.fsf.org/copyleft/gpl.html)
 */

NTL_START_IMPL;

AlgebraicFactorBase::AlgebraicFactorBase(const ZZX& _f, long bound) {
  f = _f;
  zfb.setBound(bound);
  index.SetLength(zfb.length());
  nprimes = 0;
  for (long i=0; i<zfb.length(); ++i)
    primes_above(i);
}

// update index to reflect primes lying above zfb[i].
void AlgebraicFactorBase::primes_above(long i) {
  long p = zfb[i];
  
  zz_pBak bak;
  bak.save();
  zz_p::init(p);
  
  // reduce the coefficients of f mod p
  zz_pX fp;
  for (long j=deg(f); j>=0; --j)
    SetCoeff(fp,j,to_zz_p(coeff(f,j)));
  
  // factor fp
  vec_pair_zz_pX_long u;
  CanZass(u,fp);
  
  if ((u.length()==1)&&(u[0].b==1)) {
#ifdef AFB_VERBOSE
    cout<<"AlgebraicFactorBase::primes_above() "
	<<p<<" does not split\n";
#endif
    return;
  }
  
  if (u.length()==1) {
    if (deg(u[0].a)==1) {
      // p is totally ramified
      index[i].SetLength(1);
      index[i][0] = nprimes++;
#ifdef AFB_VERBOSE
      cout<<"AlgebraicFactorBase::primes_above() "
	  <<p<<" is totally ramified\n";
#endif
    }
    else {
#ifdef AFB_VERBOSE
      cout<<"AlgebraicFactorBase::primes_above() "
	  <<p<<" has no first degree factors\n";
#endif
    }
    return;
  }
  
  // p splits
  index[i].SetLength(p);
  long small=0;
  long large=0;
  for (long j=0; j<u.length(); ++j) {
    if (deg(u[j].a)==1) {
      if (!IsOne(LeadCoeff(u[j].a))) {
	cerr<<"AlgebraicFactorBase::primes_above() factor is not monic ("
	    <<u[j].a<<")\n";
	MakeMonic(u[j].a);
      }
      index[i][rep(ConstTerm(u[j].a))] = nprimes++;
      ++small;
    }
    else {
      ++large;
    }
  }
  if (small==deg(f)) {
#ifdef AFB_VERBOSE
    cout<<"AlgebraicFactorBase::primes_above() "<<p<<" splits completely\n";
#endif
  }
  else {
#ifdef AFB_VERBOSE
    cout<<"AlgebraicFactorBase::primes_above() "
	<<p<<" splits: "<<small<<" first degree";
    if (large>0)
      cout<<"; "<<large<<" larger\n";
    else
      cout<<"\n";
#endif
  }
}

bool AlgebraicFactorBase::factor(long* fact, long c, long d) const {
  ZZ n;
  abs(n,norm(c,d));
  
  // check that the norm is smooth
  long* zfact = new long[zfb.length()];
  if (!zfb.factor(zfact,n)) {
    // norm is not smooth
    delete[] zfact;
    return false;
  }
  
  memset(fact,0,nprimes*sizeof(long));
  
  for (long i=0; i<zfb.length(); ++i) 
    if (zfact[i]) {
      if (index[i].length()==0) {
	// found a prime that doesn't split
	// Note: this can't happen if c and d are non-zero and coprime
	delete[] zfact;
	return false;
      }
      
      if (index[i].length()==1) {
	// prime is totally ramified
	fact[index[i][0]] = zfact[i];
	continue;
      }
      
      // the prime splits into multiple distict factors
      if (d==0) {
	// this shouldn't happen
	cerr<<"AlgebraicFactorBase::factor() d is zero\n";
	delete[] zfact;
	return false;
      }
      
      long di = InvMod(d,zfb[i]);
      long j = index[i][MulMod(c,di,zfb[i])];
      if (j<0) {
	cerr<<"AlgebraicFactorBase::factor() unknown prime\n";
	delete[] zfact;
	return false;
      }
      fact[j] = zfact[i];
    }
  
  delete[] zfact;
  return true;
}

// return norm of c+d\alpha
ZZ AlgebraicFactorBase::norm(long c, long d) const {
  long k = deg(f);
  // powers of c
  ZZ ci[k+1];
  conv(ci[0],1);
  for (long i=1; i<=k; ++i)
    mul(ci[i],ci[i-1],c);

  // powers of -d
  ZZ di;
  conv(di,1);

  ZZ result;
  for (long i=k; ; --i) {
    result += ci[i]*di*coeff(f,i);
    if (i==0) break;
    di*=-d;
  }
  return result;
}

// modulus for character map calculation
void AlgebraicFactorBase::SetModulus(const ZZ& _q) {
  q=_q;

  ZZ_pBak bak;
  bak.save();
  ZZ_p::init(q);

  // reduce the coefficients of f mod q
  ZZ_pX fq;
  for (long j=deg(f); j>=0; --j)
    SetCoeff(fq,j,to_ZZ_p(coeff(f,j)));
    
  // factor fq
  vec_pair_ZZ_pX_long u;
  CanZass(u,fq);

#ifdef AFB_VERBOSE
  cout<<"AFB::SetModulus() fq = "<<u<<"\n";
#endif

  conv(e,1);
  for (long i=0; i<u.length(); ++i) {
    // make sure q is not ramified
    if (u[i].b>1) {
      cerr<<"AFB::SetModulus() q is (partially) ramified\n";
      clear(e);
      return;
    }
    ZZ p;
    power(p,q,deg(u[i].a));
    --p;
    ZZ d;
    GCD(d,p,e);
    p/=d;
    e*=p;
  }

#ifdef AFB_VERBOSE
  cout<<"AFB::SetModulus() e = "<<e<<"\n";
#endif
}

// compute character map for c+d\alpha, f(\alpha)=0
void AlgebraicFactorBase::CharacterMap(ZZ* l, long c, long d) const {
  ZZ_pBak bak1;
  bak1.save();
  ZZ_p::init(sqr(q));
    
  ZZ_pX m;
  conv(m,f);
    
  ZZ_pEBak bak2;
  bak2.save();
  ZZ_pE::init(m);
    
  ZZ_pX g;
  SetCoeff(g,0,c);
  SetCoeff(g,1,d);
    
  ZZ_pE gamma;
  conv(gamma,g);
    
  // gamma = (c+da)^e - 1
  power(gamma,gamma,e);
  --gamma;
    
  for (long i=0; i<deg(f); ++i) {
    ZZ r;
    DivRem(l[i],r,rep(coeff(rep(gamma),i)),q);
    if (!IsZero(r)) {
      cerr<<"lambda() Warning: not power of q!  \n";
    }
  }
}
  
NTL_END_IMPL;

