#include <iostream>

#include <string.h>

#include <NTL/vec_long.h>

#include "vec_svec_long.h"
#include "smat_long.h"
#include "smat_ZZ_p.h"
#include "pair_ZZ_long.h"
#include "ZZFactoring.h"

#include "IndexCalculus.h"
#include "FactorBase.h"


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

inline void clear(vec_long& x) {
  memset(x.elts(),0,x.length()*sizeof(long));
}

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
    cerr<<"DLog_IC_Base::add_cache() prime already in cache\n";
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
      cout<<"DLog_IC_Base::log_g() searching... ("<<(-result)<<")  \r";
      cout.flush();
    }

  } while (true);

  if ((VERBOSE)&&(result<-1000)) {
    cout<<"\n";
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
    cerr<<"DLog_IC_Base::log() log of 0 requested\n";
    return to_ZZ(-1);
  }

  if (IsOne(pw))
    return ZZ::zero();

  if (IsZero(log_base)) {
    cerr<<"DLog_IC_Base::log() algorithm not initialized\n";
    return to_ZZ(-2);
  }

  ZZ l;
  l = log_g(pw);
  if (IsZero(l)) {
    cerr<<"DLog_IC_Base::log() log_g() returned 0\n";
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
      cerr<<"DLog_IC_Base::log() logarithm does not exist\n";
      return to_ZZ(-1);
    }
    InvMod(lb,lb,ZZ_p::modulus()-1);
    MulMod(l,l,lb,ZZ_p::modulus()-1);
  }
  return l;
}



/**************** class IndexCalculus ****************/

// various c's for L_p[1/2;c]
double IndexCalculus::L_UB=1.0;    /* upper bound */


/* Helper class for IndexCalculus.  
 *
 * We have 3 kinds of integers here:
 *   a) small primes in the factorbase
 *   b) medium size (composite) factors for the large integers
 *   c) large integers H+i for some small i
 */
class IC_Relations {
public:
  // integer factorbase
  const FactorBase& fb;
  vec_ZZ log_fb;   // logarithms of elements of factorbase

  // start of large integers
  ZZ H;

  // relations
  vec_svec_long rels;    // exponent matrix for small primes (factorbase)
  vec_ZZ rels_a,rels_b;  // medium (non-prime) factors (a>=b)

  // number of "good" relations
  long nrels;

  // sorted array of medium (non-prime) factors
  // a is the factor, b is the number of "good" relations referencing it
  vec_pair_ZZ_long medium;

  // number of non-zero entries in medium[*].b
  long nmedium;

  // logarithms of medium factors
  vec_ZZ log_med;

  // factorization of large integers
  vec_svec_long large_fact;  // smooth part
  vec_ZZ large_med;          // non-smooth part (==0 if integer not factored)

  // columns used in creation of matrix
  vec_long cols_fb;
  vec_long cols_med;

public:
  IC_Relations(const FactorBase& zfb, const ZZ& h) : fb(zfb),H(h) {
    rels.SetMaxLength(fb.length());
    rels_a.SetMaxLength(fb.length());
    rels_b.SetMaxLength(fb.length());
    nrels=0;
    medium.SetMaxLength(fb.length());
    nmedium=0;
    large_med.SetLength(fb.length());
    large_fact.SetLength(fb.length());
    log_fb.SetLength(fb.length());
  }

  ~IC_Relations() {
  }

  // figure out if we have enough relations to solve the linear system
  bool done() {

    // find "non-good" relations that are actually "good"
    bool done;
    do {
      done=true;
      for (long i=nrels; i<rels.length(); ++i) {
	// find the medium factors
	long aj=-1;
	long bj=-1;
	if (rels_a[i]>1)
	  for (long j=0; j<medium.length(); ++j) {
	    if (rels_b[i]==medium[j].a)
	      bj=j;
	    if (rels_a[i]==medium[j].a) {
	      aj=j;
	      break;
	    }
	  }
	// count how many of the medium factors are already in use
	if ((aj<0)||(bj<0)||(aj==bj)||
	    (medium[aj].b>0)||(medium[bj].b>0)) {
	  // this is actually a good relation 
	  // (ie. it introduces at most one new medium factor)
	  if (i>nrels) {
	    swap(rels[i],rels[nrels]);
	    swap(rels_a[i],rels_a[nrels]);
	    swap(rels_b[i],rels_b[nrels]);
	  }
	  ++nrels;
	  if (aj>=0) {
	    if (++medium[aj].b==1)
	      ++nmedium;
	  }
	  if ((bj>=0)&&(bj!=aj))
	    if (++medium[bj].b==1)
	      ++nmedium;
	  done=false;
	}
      }
    } while (!done);

    return (nrels>=fb.length()+nmedium);
  }

  // factor H+a and save for future use
  void factor(long a) {
    if (a<0) {
      cerr<<"IC_Relations::factor() FATAL ERROR!\n";
      exit(1);
    }
    if (a>=large_med.length()) {
      // extend large_*
      large_med.SetLength(a+1);
      large_fact.SetLength(a+1);
    }
    long* fact = new long[fb.length()];
    ZZ rm;
    if (fb.factor(fact,rm,H+a))
      conv(large_med[a],1);
    else
      large_med[a] = rm;
    large_fact[a].SetLength(fb.length());
    for (long i=0; i<fb.length(); ++i)
      if (fact[i])
	large_fact[a][i] = fact[i];
    delete[] fact;
  }

  // returns -1 if unsuccessful
  long findMediumFactor(const ZZ& a) {
    if (a<=1)
      return -1;
    // binary search to find place to insert
    long low=0;
    long high=medium.length();
    while (high>low) {
      long mid = (low+high)/2;
      if (medium[mid].a<a)
	low=mid+1;
      else // (a<=medium[mid].a)
	high=mid;
    }
    if ((low<medium.length())&&(medium[low].a==a))
      return low;
    // this should never happen
    cerr<<"findMediumFactor() factor not found!\n";
    return -1;
  }

  void addMediumFactor(const ZZ& a) {
    if (a<=1)
      return;
    // binary search to find place to insert
    long low=0;
    long high=medium.length();
    while (high>low) {
      long mid = (low+high)/2;
      if (medium[mid].a<a)
	low=mid+1;
      else // (a<=medium[mid].a)
	high=mid;
    }
    if ((low<medium.length())&&(medium[low].a==a))
      return;
    // insert a
    medium.SetLength(medium.length()+1);
    for (long i=medium.length()-1; i>low; --i)
      medium[i] = medium[i-1];
    medium[low].a = a;
    medium[low].b = 0;
  }

  // fact = (H+a)(H+b)
  void add(const long* fact, long a, long b) {
    long pos = rels.length();
    rels.SetLength(pos+1);
    rels[pos].SetLength(fb.length());
    for (long i=0; i<fb.length(); ++i)
      if (fact[i]!=0)
	rels[pos][i]=fact[i];

    // subtract factorization of a
    if ((a>=large_med.length())||(IsZero(large_med[a]))) {
      factor(a);
      addMediumFactor(large_med[a]);
    }
    rels[pos] -= large_fact[a];

    // subtract factorization of b
    if ((b>=large_med.length())||(IsZero(large_med[b]))) {
      factor(b);
      addMediumFactor(large_med[b]);
    }
    rels[pos] -= large_fact[b];

    rels_a.SetLength(pos+1);
    rels_b.SetLength(pos+1);
    if (large_med[a]>=large_med[b]) {
      rels_a[pos] = large_med[a];
      rels_b[pos] = large_med[b];
    }
    else {
      rels_a[pos] = large_med[b];
      rels_b[pos] = large_med[a];
    }
  }

  // figure out matrix and column for equation A*x=b
  // gen is generator (logarithm of it is assumed to be 1)
  void makeMatrix(smat_long& A, vec_long& b, long gen) {
    // figure out the minimum number of rows and columns we need here
    // NOTE: a few extra rows is ok 

    // indices of rows (relations) we want
    vec_long rows;
    rows.SetLength(nrels);
    for (long i=0; i<nrels; ++i)
      rows[i] = i;

    // indices of columns we want from factorbase
    cols_fb.SetLength(fb.length()-1);
    long gen_index=-1;
    for (long i=0,j=0; i<fb.length(); ++i) {
      if (gen!=fb[i])
	cols_fb[j++] = i;
      else
	gen_index=i;
    }
    if (gen_index<0) {
      cerr<<"IndexCalculus::makeMatrix() generator not in factorbase!\n";
      return;
    }

    // note logarithm of generator
    log_fb[gen_index]=1;

    // indices of columns we want from medium factors
    cols_med.SetMaxLength(medium.length());
    for (long i=0; i<medium.length(); ++i)
      if (medium[i].b>0) {
	cols_med.SetLength(cols_med.length()+1);
	cols_med[cols_med.length()-1] = i;
      }

    // variables need in loop
    bool done;
    long cols_zero[max(cols_fb.length(),cols_med.length())];

    // reduction
    do {
      // optimistic?
      done=true;

      // count non-zero values in each column (we want zero or one)
      for (long j=0; j<cols_fb.length(); ++j)
	cols_zero[j] = -1;
      for (long i=0; i<rows.length(); ++i) {
	const svec_long& rel = rels[rows[i]];
	for (long j=0; j<cols_fb.length(); ++j) {
	  if (rel[cols_fb[j]]) {
	    if (cols_zero[j]==-1)
	      cols_zero[j] = i;
	    else
	      cols_zero[j] = -2;
	  }
	}
      }
      // look for duplicates
      for (long j=0; j<cols_fb.length(); ++j)
	if (cols_zero[j]>=0)
	  for (long i=j+1; i<cols_fb.length(); ++i)
	    if (cols_zero[i]==cols_zero[j])
	      cols_zero[i]=-1;
      // eliminate columns with zero or one non-zero values
      // (in case of one non-zero value, eliminiate the row too)
      for (long j=cols_fb.length()-1; j>=0; --j)
	if (cols_zero[j]!=-2) {
	  // delete column
	  for (long i=j; i<cols_fb.length()-1; ++i)
	    cols_fb[i] = cols_fb[i+1];
	  cols_fb.SetLength(cols_fb.length()-1);
	  if (cols_zero[j]>=0)
	    rows[cols_zero[j]]=-1; // mark row for deletion
	  done=false;
	}
      // actually delete rows
      for (long j=rows.length()-1; j>=0; --j)
	if (rows[j]==-1) {
	  for (long i=j; i<rows.length()-1; ++i)
	    rows[i] = rows[i+1];
	  rows.SetLength(rows.length()-1);
	}

      // do the same with the "med" columns
      for (long j=0; j<cols_med.length(); ++j)
	cols_zero[j] = -1;
      for (long j=0; j<cols_med.length(); ++j)
	for (long i=0; i<rows.length(); ++i) {
	  if ((rels_a[rows[i]]==medium[cols_med[j]].a)||
	      (rels_b[rows[i]]==medium[cols_med[j]].a)) {
	    if (cols_zero[j]==-1)
	      cols_zero[j] = i;
	    else
	      cols_zero[j] = -2;
	  }
	}
      // look for duplicates
      for (long j=0; j<cols_med.length(); ++j)
	if (cols_zero[j]>=0)
	  for (long i=j+1; i<cols_med.length(); ++i)
	    if (cols_zero[i]==cols_zero[j])
	      cols_zero[i]=-1;
      // eliminate columns with zero or one non-zero values
      // (in case of one non-zero value, eliminiate row too)
      for (long j=cols_med.length()-1; j>=0; --j)
	if (cols_zero[j]!=-2) {
	  // delete column
	  for (long i=j; i<cols_med.length()-1; ++i)
	    cols_med[i] = cols_med[i+1];
	  cols_med.SetLength(cols_med.length()-1);
	  if (cols_zero[j]>=0)
	    rows[cols_zero[j]]=-1;  // mark row for deletion
	  done=false;
	}
      // actually delete rows
      for (long j=rows.length()-1; j>=0; --j)
	if (rows[j]==-1) {
	  for (long i=j; i<rows.length()-1; ++i)
	    rows[i] = rows[i+1];
	  rows.SetLength(rows.length()-1);
	}
      
    } while (!done);

    if (DLog_IC_Base::VERBOSE)
      cout<<"IndexCalculus::makeMatrix() "
	  <<rows.length()<<" rows, "
	  <<cols_fb.length()<<" fb cols, "
	  <<cols_med.length()<<" med cols\n";
    
    // make matrix and vector
    A.SetDims(rows.length(),cols_fb.length()+cols_med.length());
    b.SetLength(rows.length());
    for (long i=0; i<rows.length(); ++i) {
      const svec_long& rel = rels[rows[i]];
      // set exponents for small primes (factorbase)
      for (long j=0; j<cols_fb.length(); ++j) 
	if (rel[cols_fb[j]])
	  A[i][j]=rel[cols_fb[j]];
      // set exponents for medium factors
      for (long j=0; j<cols_med.length(); ++j) {
	if (rels_b[rows[i]]==medium[cols_med[j]].a)
	  A[i][cols_fb.length()+j] = -1;
	if (rels_a[rows[i]]==medium[cols_med[j]].a) {
	  A[i][cols_fb.length()+j] -= 1;  // remember: a and b might equal
	  break;
	}
      }
      // set constant row
      if (rels[rows[i]][gen_index])
	b[i] = -rels[rows[i]][gen_index];
      else
	b[i] = 0;
    }
  }

  // compute logarithm of as many elements of factorbase and medium factors
  // as possible (others are set to zero)
  void computeLogs(const vec_ZZ& x, const ZZ& q) {
    // find indices of medium factors
    vec_long rels_aj, rels_bj;
    rels_aj.SetLength(rels.length());
    rels_bj.SetLength(rels.length());
    for (long i=0; i<rels.length(); ++i) {
      rels_aj[i]=findMediumFactor(rels_a[i]);
      rels_bj[i]=findMediumFactor(rels_b[i]);
    }

    log_med.SetLength(medium.length());
    for (long i=0; i<cols_fb.length(); ++i)
      log_fb[cols_fb[i]] = x[i];
    for (long i=0; i<cols_med.length(); ++i)
      log_med[cols_med[i]] = x[cols_fb.length()+i];
    // figure out other logarithms from relations
    bool rels_used[rels.length()];
    memset(rels_used,0,rels.length()*sizeof(bool));
    bool done;
    do {
      done=true;
      for (long i=0; i<rels.length(); ++i) {
	if (!rels_used[i]) {
	  // find a single unknown
	  long unknown=-1;
	  long rn = rels[i].nvalues();
	  const long* rj = rels[i].indices();
	  const long* rv = rels[i].values();
	  for (long j=0; j<rn; ++j)
	    if (IsZero(log_fb[rj[j]])&&(rv[j]!=0)) {
	      if (unknown==-1)
		unknown=rj[j];
	      else {
		// too many unknowns
		unknown=-2;
		break;
	      }
	    }
	  long bj=rels_bj[i];
	  if ((bj>=0)&&IsZero(log_med[bj])) {
	    if (unknown==-1)
	      unknown=fb.length();
	    else
	      unknown=-2;
	  }
	  long aj=rels_aj[i];
	  if ((aj>=0)&&(aj!=bj)&&IsZero(log_med[aj])) {
	    if (unknown==-1)
	      unknown=fb.length()+1;
	    else
	      unknown=-2;
	  }
	  if (unknown==-1) {
	    // useless relations
	    rels_used[i]=true;
	    continue;
	  }
	  if (unknown>=0) {
	    // here we can compute an unknown logarithm
	    ZZ lg;
	    if (aj>=0)
	      lg += log_med[aj];
	    if (bj>=0)
	      lg += log_med[bj];
	    for (long j=0; j<rn; ++j)
	      lg -= rv[j]*log_fb[rj[j]];
	    if (unknown<fb.length())
	      rem(log_fb[unknown],lg,q);
	    else if (unknown==fb.length())
	      rem(log_med[bj],-lg,q);
	    else
	      rem(log_med[aj],-lg,q);
	    done=false;
	  }
	}
      }
    } while (!done);
  }

  // log of H+a
  ZZ log_lg(long a, const ZZ& q) const {
    if ((a>=large_med.length())||IsZero(large_med[a]))
      return ZZ::zero();
    ZZ lg;
    long an = large_fact[a].nvalues();
    const long* aj = large_fact[a].indices();
    const long* av = large_fact[a].values();
    for (long i=0; i<an; ++i)
      if (av[i]!=0) {
	if (IsZero(log_fb[aj[i]]))
	  return ZZ::zero();
	lg += av[i]*log_fb[aj[i]];
      }
    if (!IsOne(large_med[a])) {
      for (long i=0; i<medium.length(); ++i)
	if (large_med[a]==medium[i].a) {
	  if (IsZero(log_med[i]))
	    return ZZ::zero();
	  lg += log_med[i];
	  break;
	}
    }
    rem(lg,lg,q);
    return lg;
  }
  
};


bool IndexCalculus::make_system() {

  // p = modulus
  const ZZ &p = ZZ_p::modulus();

  // H = ceil(sqrt(p))
  ZZ H;
  H=1+SqrRoot(p);
  if (VERBOSE)
    cout<<"IndexCalculus::make_system() H = "<<H<<"\n";

  // variables needed for sieve
  IC_Relations rels(zfb,H);
  long nprimes = zfb.length();
  long* fact = new long[zfb.length()];
  vec_short sieve;
  sieve.SetMaxLength(sieve_length);

  if (VERBOSE)
    cout<<"IndexCalculus::make_system() sieve_length = "<<sieve_length<<"\n";

  // sieve stats
  long sieve_count=0;
  long sieve_total=0;
  long sieve_smooth=0;
  long sieve_bad=0;

  // sieve comparison
  double sieve_time=0;

  ZZ res;  // residule
  ZZX f;   // polynomial to sieve over

  // sieving
  bool done=false;
  long last_percent=0;
  long last_count=0;
  for (long d=0; !done; ++d) {
    rem(res,H*(H+d),p);
    SetCoeff(f,0,res);
    SetCoeff(f,1,H+d);

    sieve.SetLength(d+1);
    clear(sieve);
    double t_start = GetTime();
    zfb.sieve(sieve,f);
    sieve_time += GetTime()-t_start;

    ++sieve_count;
    sieve_total+=sieve.length();
    
    for (long c=0; c<sieve.length(); ++c) {
      if (sieve[c]==1) {
	rem(res,(H+c)*(H+d),p);
	if (!zfb.factor(fact,res)) {
	  ++sieve_bad;
	  continue;
	}
	++sieve_smooth;
	rels.add(fact,c,d);
      }
    }
    done=rels.done();

    long percent = 100 - (nprimes+rels.nmedium-rels.nrels)*100/nprimes;
    if ((percent!=last_percent)||(sieve_count-last_count>=10)) {
      cout<<"Sieving: "
	  <<sieve_count<<" sieves, "
	  <<sieve_smooth<<"/"<<sieve_total<<" smooth, "
	  <<rels.nrels<<"/"<<(nprimes+rels.nmedium)<<" needed";
      if (percent<100)
	cout<<" ("<<percent<<"%)  \r";
      else
	cout<<"        \r";
      cout.flush();
      last_percent=percent;
      last_count=sieve_count;
    }
  }
  cout<<"\n";

  //sieve1_times();
  //std::cout<<"sieve1: "<<sieve_time<<" seconds"<<std::endl;

  delete[] fact;

  // reporting
  if (sieve_bad)
    cout<<"IndexCalculus::make_system() WARNING: "<<sieve_bad
	<<" non-smooth integers from sieve\n";
  if (!done) {
    cerr<<"IndexCalculus::make_system() sieve was too short\n";
    return false;
  }

#ifdef IC_EXTRA_SIEVING
  // complete the sieving process for accurate timing
  /*
  for (long c=sieve_count; c<sieve_length; ++c) {
    rem(res,(H+c)*(H+c),p);
    SetCoeff(f,0,res);
    SetCoeff(f,1,H+c);
    sieve.SetLength(sieve_length-c);
    for (long i=0; i<sieve_length-c; ++i)
      sieve[i]=0;
    zfb.sieve(sieve,f);
  }
  */
#endif
  
  // create matrix equation
  smat_long A;
  vec_long b;
  rels.makeMatrix(A, b, to_long(rep(g)));

  if (VERBOSE) {
    cout<<"IndexCalculus::make_system() A is "
	<<A.NumRows()<<"x"<<A.NumCols()<<"\n";
  }

  // solution to linear system
  vec_ZZ x;
  
  if (!solve_system(A,x,b,q)) {
    cerr<<"IndexCalculus::make_system() failed to solve system\n";
    return false;
  }

  rels.computeLogs(x,q);

  // number of incorrect logarithms (don't know how this happens?)
  long incorrect=0;

  // copy solution to cache
  ZZ lg,pr;
  for (long i=0; i<nprimes; ++i) {
    if (!IsZero(rels.log_fb[i])) {
      conv(pr,zfb[i]);
      lg = log_complete(pr,rels.log_fb[i],q);
      if (!IsZero(lg))
	add_cache(pr,lg);
      else
	++incorrect;
    }
  }
  // check medium sized primes
  for (long i=0; i<rels.log_med.length(); ++i) {
    if (ProbPrime(rels.medium[i].a)&&!IsZero(rels.log_med[i])) {
      lg = log_complete(rels.medium[i].a,rels.log_med[i],q);
      if (!IsZero(lg))
	add_cache(rels.medium[i].a,lg);
      else
	++incorrect;
    }
  }
  if (VERBOSE)
    cout<<"IndexCalculus::make_system() "<<ncache<<" cached logarithms\n";
  if (incorrect) {
    cout<<"IndexCalculus::make_system() "
	<<incorrect<<" INCORRECT LOGARITHMS\n";
    if (incorrect>ncache)
      return false;
  }

  // "upper factorbase"
  incorrect=0;
  sufb = H;
  long nufb=rels.large_med.length();
  while ((nufb>0)&&(IsZero(rels.large_med[nufb-1])))
    --nufb;
  ufb.SetLength(nufb);
  long cufb=0;
  for (long i=0; i<nufb; ++i) {
    lg = rels.log_lg(i,q);
    if (!IsZero(lg)) {
      lg = log_complete(sufb+i,lg,q);
      if (!IsZero(lg)) {
	ufb[i]=lg;
	++cufb;
      }
      else
	++incorrect;
    }
  }
  while ((ufb.length()>0)&&(IsZero(ufb[ufb.length()-1])))
    ufb.SetLength(ufb.length()-1);
  if (VERBOSE)
    cout<<"IndexCalculus::make_system() "<<cufb<<" upper logarithms\n";
  if (incorrect) {
    cout<<"IndexCalculus::make_system() "
	<<incorrect<<" INCORRECT LOGARITHMS\n";
    if (incorrect>cufb)
      return false;
  }

  return true;
}

// solve A*x=b (mod q)
bool IndexCalculus::solve_system(const smat_long& Aorig, vec_ZZ& x,
				 const vec_long& b, const ZZ& q) {
  ZZ_pBak bak;
  bak.save();
  ZZ_p::init(q);

  // constant term
  vec_ZZ_p yorig;
  yorig.SetLength(b.length());
  for (long i=0; i<b.length(); ++i)
    conv(yorig[i],b[i]);

  // the solution
  vec_ZZ_p X;

  // structured gaussian elimination
  vec_long cols;
  smat_long A;
  vec_ZZ_p y;
  SGauss(A,y,cols,Aorig,yorig);

  if (VERBOSE)
    cout<<"IndexCalculus::solve_system() reduced matrix is "
	<<A.NumRows()<<"x"<<A.NumCols()<<"\n";

  // solve system
  if (!Lanczos(A,X,y))
    return false;

  // undo structured gaussian elimination
  vec_ZZ_p Xorig;
  SGauss_undo(Xorig,X,cols,Aorig,yorig);

  // set solution
  x.SetLength(Xorig.length());
  for (long i=0; i<x.length(); ++i)
    x[i] = rep(Xorig[i]);
  
  return true;
}

void IndexCalculus::parameters(long& bound, long& length, const ZZ& p) {
  if (bound<10)
    bound = (long)(3.3*L_p(p,0.475));

  if (length<10)
    length = (long)(2.75*L_p(p,0.42));
}

void IndexCalculus::setBase(const ZZ_p& base, long bound, long length) {
  DiscreteLog::setBase(base);
  this->base=base;
  sieve_length = length;

  // factor order of group
  vec_pair_ZZ_long factors;
  factor(factors,ZZ_p::modulus()-1);

  // largest factor must have exponent 1
  if (factors[factors.length()-1].b!=1) {
    cerr<<"IndexCalculus::setBase() largest factor of order has multiplicity>1\n";
    return;
  }

  // modulus for solving the linear system
  q = factors[factors.length()-1].a;
  if (VERBOSE)
    cout<<"IndexCalculus::setBase() q = "<<q<<"\n";

  // other factor (solved using Pollard Rho method)
  if (VERBOSE)
    cout<<"IndexCalculus::setBase() other = "<<((ZZ_p::modulus()-1)/q)<<"\n";

  // figure out optimal parameters
  parameters(bound,sieve_length,ZZ_p::modulus());

  // factorbase
  zfb.setBound(bound);
  if (VERBOSE)
    cout<<"IndexCalculus::setBase() "<<zfb.length()<<" primes in factorbase\n";

  // check base
  ZZ b(rep(base));
  if ((b<NTL_MAX_LONG)&&(zfb.isPrime(to_long(b)))&&
      (isGenerator(base,factors))) {
    g=base;
    log_base=1;
  }
  else {
    bool found=false;
    for (long i=0; i<zfb.length(); ++i) {
      conv(g,zfb[i]);
      if (isGenerator(g,factors)) {
        found=true;
        break;
      }
    }
    if (!found) {
      cerr<<"IndexCalculus::setBase() failed to find generator (increase factorbase)\n";
      g=0;
      log_base=0;
      return;
    }
  }

  if (VERBOSE)
    cout<<"IndexCalculus::setBase() generator is "<<g<<"\n";

  // sieving
  if (!make_system()) {
    cerr<<"IndexCalculus::setBase() sieving failed!\n";
    log_base=0;
    return;
  }

  // upper bound on smoothness tests
  conv(upper_bound,1.0*L_p(ZZ_p::modulus(),L_UB));

  // find log of base (if necessary)
  if (IsZero(log_base)) {
    log_base = log_g(base);
    if (IsZero(log_base))
      return;
  }
  if (VERBOSE)
    cout<<"IndexCalculus::solve_system() log_base is "<<log_base<<"\n";
}

// logarithm of a (medium sized) prime
ZZ IndexCalculus::log_prime(const ZZ& pw) {

  // variables used below
  ZZ logy,logn;
  ZZ_p y,n;
  vec_short sieve;
  ZZX f;

  // iterate through values of y
  conv(y,SqrRoot(ZZ_p::modulus())/pw);
  while (true) {
    ++y;
    if (IsOne(y))
      clear(logy);
    else {
      logy = log_g(y,true);
      if (IsZero(logy))
	continue;
    }
    if (VERBOSE)
      cout<<"IndexCalculus::log_prime() y="<<y<<" \n";

    // sieve to find v
    sieve.SetLength(ufb.length());
    for (long i=0; i<ufb.length(); ++i)
      sieve[i] = IsZero(ufb[i]) ? -1 : 0;
    SetCoeff(f,0,rep(y)*pw*sufb-ZZ_p::modulus());
    SetCoeff(f,1,rep(y)*pw);
    zfb.sieve(sieve,f);
    for (long i=0; i<ufb.length(); ++i) {
      if (sieve[i]==1) {
	conv(n,pw*(sufb+i));
	n *= y;
	logn = log_g(n,true);
	if (!IsZero(logn)) {
	  if (VERBOSE)
	    cout<<"IndexCalculus::log_prime() v="<<(sufb+i)<<" \n";
	  logn -= logy;
	  logn -= ufb[i];
	  rem(logn,logn,ZZ_p::modulus()-1);
	  return logn;
	}
      }
    }
  }
}


NTL_END_IMPL;
