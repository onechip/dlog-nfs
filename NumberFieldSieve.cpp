#include <stdio.h>
#include <string.h>

#include <NTL/ZZXFactoring.h>
#include "ZZFactoring.h"
#include "smat_ZZ_p.h"

#include "NumberFieldSieve.h"
#include "AlgebraicFactorBase.h"


/* Number Field Sieve (NFS) algorithm.
 *
 * Written by: Chris Studholme
 * Copyright:  GPL (http://www.fsf.org/copyleft/gpl.html)
 */

NTL_START_IMPL;

// create a file named 'nfs.data' showing the distribution of smooth integers
//#define NFS_CREATE_DIST_FILE

long DLog_NFS::MAX_SIEVE=20000000;  /* maximum size of a single sieve */


inline long rem(long a, long b) {
  long r = a%b;
  while (r<0)
    r+=b;
  return r;
}

// compute fg = f(g(X))
inline void compose(ZZX& fg, const ZZX& f, const ZZX& g) {
  fg = coeff(f,deg(f));
  for (long i=deg(f)-1; i>=0; --i) {
    fg*=g;
    fg+=coeff(f,i);
  }
}

// find and remove a from the sorted vector v
inline void remove(vec_long& v, long a) {
  // binary search for a in v
  long low=0;
  long high=v.length();
  while (high>low) {
    long mid = (low+high)/2;
    if (v[mid]<a)
      low=mid+1;
    else if (a<v[mid])
      high=mid;
    else {
      // found it, delete it
      for (long i=mid+1; i<v.length(); ++i)
	v[i-1]=v[i];
      v.SetLength(v.length()-1);
      return;
    }
  }
}

class NFS_Relations {
public:
  // degree of number field (size of character map)
  long k;

  // rational factorbase
  const FactorBase& zfb;

  // algebraic factorbase
  const AlgebraicFactorBase& afb;

  // relations
  vec_svec_ZZ rels; // exponent matrix


public:
  NFS_Relations(long _k, const AlgebraicFactorBase& _afb,
		const FactorBase& _zfb) : zfb(_zfb),afb(_afb) {
    k=_k;
    rels.SetMaxLength(k+afb.length()+zfb.length());
  }

  ~NFS_Relations() {
  }

  // figure out if we have enough relations to solve the linear system
  bool done() {
    long extra = (long)log((double)(k+afb.length()+zfb.length()));
    return (rels.length()>k+afb.length()+zfb.length()+extra);
  }

  void add(const long* zfact, const long* afact, const ZZ* charmap) {
    long pos=rels.length();
    rels.SetLength(pos+1);
    for (long i=0; i<k; ++i)
      if (!IsZero(charmap[i]))
	rels[pos][i] = charmap[i];
    for (long i=0; i<afb.length(); ++i)
      if (afact[i]!=0)
	rels[pos][k+i] = afact[i];
    for (long i=0; i<zfb.length(); ++i)
      if (zfact[i]!=0)
	rels[pos][k+afb.length()+i] = zfact[i];
  }

  // figure out matrix for equation A*x=newb
  // b is the initial value for newb (newb will be the same but with some
  //   zero columns removed)
  // gen is generator (must be a prime in the rational factorbase)
  bool makeMatrix(smat_ZZ_p& A, vec_ZZ& newb, 
		  const vec_ZZ& b, long gen) const {

    // find generator in zfb
    long gen_index=-1;
    for (long j=0; j<zfb.length(); ++j)
      if (gen==zfb[j])
	gen_index=j;
    if (gen_index<0) {
      cerr<<"NumberFieldSieve::makeMatrix() generator not in factorbase!\n";
      return false;
    }

    // list of columns we must keep
    bool* must_keep = new bool[k+afb.length()+zfb.length()];
    memset(must_keep,0,(k+afb.length()+zfb.length())*sizeof(bool));
    must_keep[k+afb.length()+gen_index]=true;
    for (long j=0; j<b.length(); ++j)
      if (!IsZero(b[j]))
	must_keep[j]=true;

    // columns statistics (which rows are using each column)
    vec_vec_long col_use;
    col_use.SetLength(k+afb.length()+zfb.length());
    for (long i=0; i<rels.length(); ++i) {
      long rn = rels[i].nvalues();
      const long* rj = rels[i].indices();
      const ZZ* rv = rels[i].values();
      for (long j=0; j<rn; ++j)
	if (!IsZero(rv[j]))
	  append(col_use[rj[j]],i);
    }

    // status of each row
    bool* row_kept = new bool[rels.length()];
    for (long i=0; i<rels.length(); ++i)
      row_kept[i]=true;
    long row_count=rels.length();

    // main loop to eliminate unneeded rows and columns
    bool done;
    do {
      done=true;  // optimistic?

      // find columns with only one non-zero value
      for (long j=0; j<col_use.length(); ++j) 
	if ((col_use[j].length()==1)&&(!must_keep[j])) {
	  long row = col_use[j][0];
	  long rn = rels[row].nvalues();
	  const long* rj = rels[row].indices();
	  // remove row from all elements of col_use
	  for (long i=0; i<rn; ++i)
	    remove(col_use[rj[i]],row);
	  row_kept[row]=false;
	  --row_count;
	  done=false;
	}

      // count non-empty columns
      long col_count=0;
      for (long j=0; j<col_use.length(); ++j)
	if (col_use[j].length()>0)
	  ++col_count;

      if (row_count>=col_count) {
	// need to choose a row to delete here
	// find the row that matters least as far as must_keep columns go
	long best=0;
	long best_value=0;
	for (long i=0; i<rels.length(); ++i) 
	  if (row_kept[i]) {
	    long rn = rels[i].nvalues();
	    const long* rj = rels[i].indices();
	    const ZZ* rv = rels[i].values();
	    long value=row_count;
	    for (long j=0; j<rn; ++j)
	      if (must_keep[rj[j]]&&
		  (value>col_use[rj[j]].length())&&
		  (!IsZero(rv[j])))
		value = col_use[rj[j]].length();
	    if (value>best_value) {
	      best_value=value;
	      best=i;
	    }
	  }
	// delete row best
	long rn = rels[best].nvalues();
	const long* rj = rels[best].indices();
	for (long i=0; i<rn; ++i)
	  remove(col_use[rj[i]],best);
	row_kept[best]=false;
	--row_count;
	done=false;
      }
      
    } while (!done);

    // check out must_keep columns to ensure they are represented
    for (long j=0; j<col_use.length(); ++j)
      if (must_keep[j]&&(col_use[j].length()==0)) {
	cerr<<"NFS::makeMatrix() column "<<j<<" has no representation\n";
	delete[] must_keep;
	delete[] row_kept;
	return false;
      }

    delete[] must_keep;

    // list of rows we are keeping
    vec_long rows;
    rows.SetMaxLength(rels.length());
    for (long i=0; i<rels.length(); ++i)
      if (row_kept[i])
	append(rows,i);

    delete[] row_kept;

    // list of columns we are keeping
    vec_long cols;
    for (long j=0; j<col_use.length(); ++j)
      if (col_use[j].length()>0)
	append(cols,j);
    
    // make the matrix (actually, transpose of matrix)
    A.SetDims(cols.length(),1+rows.length());
    clear(A);
    // first column is one associated with generator
    for (long j=0; j<cols.length(); ++j)
      if (cols[j]==(k+afb.length()+gen_index))
	conv(A[j][0],1);
    // all the other columns
    for (long i=0; i<rows.length(); ++i) {
      long rn = rels[rows[i]].nvalues();
      const long* rj = rels[rows[i]].indices();
      const ZZ* rv = rels[rows[i]].values();
      for (long c=0,j=0; c<cols.length(); ++c) {
	while ((j<rn)&&(rj[j]<cols[c]))
	  ++j;
	if (j==rn)
	  break;
	if (rj[j]==cols[c])
	  conv(A[c][1+i],rv[j]);
      }
    }
    
    // copy kept columns from b to newb
    newb.SetLength(cols.length());
    for (long j=0; j<cols.length(); ++j)
      newb[j] = b[cols[j]];

    return true;
  }


};


/**************** class DLog_NFS ****************/

// compute optimum parameters for a given p
// any of k, bound, width that are non-zero are not computed
void DLog_NFS::parameters(long& k, long& bound, long& width, const ZZ& p) {
  // figure out optimum degree of number field
  if (k<2) {
    if (NumBits(p)<64)
      k=2;
    else if (NumBits(p)<160)
      k=3;
    else if (NumBits(p)<256)
      k=4;
    else if (NumBits(p)<512)
      k=5;
    else
      k=6;
  }

  // figure out optimum smoothness bound
  if (bound<10) {
    switch(k) {
    case 2:
      bound = (long)(0.34*L_p(p,1.06));
      break;
    case 3:
      bound = (long)(0.71*L_p(p,0.98));
      break;
    case 4:
      bound = (long)(1.19*L_p(p,0.96));
      break;
    case 5:
      bound = (long)(1.0*L_p(p,1.0));
      break;
    case 6:
      bound = (long)(1.0*L_p(p,1.0));
      break;
    default:
      bound = (long)(1.0*L_p(p,1.0));
    }
  }
  
  // figure out optimum sieve dim
  if (width<10) {
    switch(k) {
    case 2:
      width = (long)(2.4*L_p(p,1.67));
      break;
    case 3:
      width = (long)(4.1*L_p(p,1.22));
      break;
    case 4:
      width = (long)(0.15*L_p(p,1.47));
      break;
    case 5:
      width = (long)(0.2*L_p(p,2.0));
      break;
    case 6:
      width = (long)(0.2*L_p(p,2.0));
      break;
    default:
      width = (long)(0.2*L_p(p,2.0));
    }
  }
}

void DLog_NFS::setBase(const ZZ_p& base, long _k, long bound, 
		       long _sieve_width) {
  DiscreteLog::setBase(base);
  this->base=base;
  k=_k;
  sieve_width = _sieve_width;
  sieve_length = 0;

  // modulus
  ZZ p;
  p = ZZ_p::modulus();

  // factor order of group
  vec_pair_ZZ_long factors;
  factor(factors,p-1);

  // largest factor must have exponent 1
  if (factors[factors.length()-1].b!=1) {
    cerr<<"DLog_NFS::setBase() largest factor of order has multiplicity>1\n";
    return;
  }

  // modulus for solving the linear system
  q = factors[factors.length()-1].a;
  if (VERBOSE)
    cout<<"DLog_NFS::setBase() q = "<<q<<"\n";

  // other factor (solved using Pollard Rho method)
  if (VERBOSE)
    cout<<"DLog_NFS::setBase() other = "<<((p-1)/q)<<"\n";

  parameters(k,bound,sieve_width,p);

  cout<<"DLog_NFS::setBase() degree of number field is "<<k<<"\n";
  cout<<"DLog_NFS::setBase() lower smoothness bound is "<<bound<<"\n";

  // set upper bound (k'th root of p)
  conv(upper_bound,exp(NTL::log(p)/k));
  if (VERBOSE)
    cout<<"DLog_NFS::setBase() upper smoothness bound is "<<upper_bound<<"\n";

  // factorbase
  zfb.setBound(bound);
  if (VERBOSE)
    cout<<"DLog_NFS::setBase() "
	<<zfb.length()<<" primes in rational factorbase\n";

  // check base
  ZZ b(rep(base));
  if ((b<NTL_MAX_LONG)&&(zfb.isPrime(to_long(b)))&&
      (isGenerator(base,factors))) {
    g=base;
    log_base=1;
    add_cache(rep(g),to_ZZ(1));
    if (VERBOSE)
      cout<<"DLog_NFS::setBase() generator is "<<g<<"\n";
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
      cerr<<"DLog_NFS::setBase() failed to find generator (increase factorbase)\n";
      g=0;
      log_base=0;
      return;
    }
    if (VERBOSE)
      cout<<"DLog_NFS::setBase() generator is "<<g<<"\n";
    add_cache(rep(g),to_ZZ(1));

    // find log of base
    log_base = log_g(base);
  }
  if (VERBOSE)
    cout<<"DLog_NFS::setBase() log_base is "<<log_base<<"\n";
}


// solve linear system induced by rels and b (mod q)
bool DLog_NFS::ls_solve(vec_ZZ& X, const NFS_Relations& rels, const vec_ZZ& b,
			const ZZ& q) {

  ZZ_pBak bak;
  bak.save();
  ZZ_p::init(q);

  double t1 = GetTime();

  // the matrix
  smat_ZZ_p Aq;
  vec_ZZ nb;
  if (!rels.makeMatrix(Aq,nb,b,to_long(rep(g)))) {
    cerr<<"DLog_NFS::ls_solve() failed to create suitable matrix\n";
    return false;
  }

  if (VERBOSE)
    cout<<"DLog_NFS::ls_solve() matrix is "
	<<Aq.NumRows()<<"x"<<Aq.NumCols()<<"\n";
  
  vec_ZZ_p bq;
  bq.SetLength(nb.length());
  for (long i=0; i<nb.length(); ++i)
    conv(bq[i],nb[i]);

  vec_ZZ_p Xq;

  double t2 = GetTime();
  if (VERBOSE)
    cout<<"DLog_NFS::ls_solve() creating the matrix took "<<(t2-t1)<<" seconds\n";
  
  if (!Lanczos(Aq,Xq,bq))
    return false;

  // copy solution
  X.SetLength(Xq.length());
  for (long i=0; i<Xq.length(); ++i)
    X[i] = rep(Xq[i]);

  double t3 = GetTime();
  if (VERBOSE)
    cout<<"DLog_NFS::ls_solve() Lanczos took "<<(t3-t2)<<" seconds\n";
  
  return true;
}


// discrete log of a (medium sized) prime relative to g
// assumes mp!=g
ZZ DLog_NFS::log_prime(const ZZ& _mp) {
  ZZ mp(_mp);
  // special method for very small primes 
  //bool issmall=(mp<rep(g));
  //if (issmall)
  // mp*=rep(g);

  cout<<"DLog_NFS::log_prime() log of "<<mp<<" (base "<<g<<")\n";

  // modulus
  ZZ p;
  p = ZZ_p::modulus();

  // size of factorbase
  long zprimes = zfb.length();

  // find smallest h such that m = 2^h * mp satisfies m>p^{1/k}
  ZZ h;
  h=1;
  ZZ m;
  mul(m,mp,2);
  ZZ pk;
  conv(pk,exp(NTL::log(p)/k));
  //cout<<"DLog_NFS::log_prime() pk = "<<pk<<"\n";
  while (m<=pk) {
    m*=2;
    ++h;
  }
  if (VERBOSE)
    cout<<"DLog_NFS::log_prime() m = "<<m<<"\n";
  
  // find least positive c such that cp >= m^k
  ZZ mk;
  power(mk,m,k);
  ZZ cp(p);
  while (cp<mk)
    cp+=p;

  // construct minimum polynomial by writing cp in base m
  ZZX f;
  SetCoeff(f,k);  // monic
  // write p in base m
  ZZ j;
  j = cp - mk;
  for (long i=k-1; i>0; --i) {
    ZZ mi;
    power(mi,m,i);
    ZZ c;
    div(c,j,mi);  // c=j/mi;
    SetCoeff(f,i,c);
    j -= c*mi;
  }
  // constant term
  SetCoeff(f,0,j);
  //cout<<"DLog_NFS::log_prime() f' = "<<f<<"\n";

  // make sure constant term in f is smooth by subtracting multiples of m
  long fact[zfb.length()];
  long D=0;
  while (!zfb.factor(fact,abs(ConstTerm(f)))) {
    // add X-m to f
    SetCoeff(f,0,coeff(f,0)-m);
    SetCoeff(f,1,coeff(f,1)+1);
    if (++D>100) {
      cerr<<"DLog_NFS::log_prime() failed to make constant term smooth!\n";
      return ZZ::zero();
    }
  }
  if (VERBOSE)
    cout<<"DLog_NFS::log_prime() f = "<<f<<"\n";

  // check that f is irreducible
  ZZ uc;
  vec_pair_ZZX_long u;
  factor(uc,u,f);
  //cout<<"DLog_NFS::log_prime() factorization = "<<u<<"\n";
  if ((u.length()>1)||(u[0].b!=1)) {
    cerr<<"DLog_NFS::log_prime() f is not irreducible!\n";
    return ZZ::zero();
  }

  // check that q does not divide the discriminant
  ZZ disc;
  disc = discriminant(f);
  if (IsZero(disc%q)) {
    cerr<<"DLog_NFS::log_prime() q divides discriminant!\n";
    return ZZ::zero();
  }

  AlgebraicFactorBase afb(f,zfb.bound());
  long aprimes = afb.length();
  if (VERBOSE)
    cout<<"DLog_NFS::log_prime() "
	<<aprimes<<" primes in algebraic factorbase\n";

  // modulus for character map calculation
  afb.SetModulus(q);

  // helper class to manage relations
  NFS_Relations rels(k,afb,zfb);

  // factorizations (do we really have enough space on the stack for these?)
  long zfact[zprimes];
  long afact[aprimes];
  ZZ l[k];

  if (sieve_width>m) {
    cerr<<"DLog_NFS::log_prime() sieve_width too large ("<<sieve_width<<")\n";
    //conv(sieve_width,m);
    return ZZ::zero();
  }

  // make sure width is odd
  if (sieve_width%2==0)
    --sieve_width;

  cout<<"DLog_NFS::log_prime() c in range [-"
      <<(sieve_width/2)<<","<<(sieve_width/2)<<"]  "<<endl;

  // the sieve
  vec_short sieve;
  sieve.SetLength(sieve_width);

  // sieve stats
  long sieve_count=0;
  ZZ sieve_total;
  long sieve_smooth=0;
  long sieve_bad=0;

  // variable used in loop
  ZZX fz,fa,fs;
  ZZ t;

#ifdef NFS_CREATE_DIST_FILE
  FILE* dist_file = fopen("nfs.data","w");
#endif

  double t1=GetTime();

  cout<<"Sieving... \r";
  cout.flush();

  long d=0,d_last=0,d_lastsmooth=0;
  long percent_last=-1;
  bool done;
  do {
    ++d;
    done=false;

    // polynomial for rational integers
    SetCoeff(fz,0,d*m);
    SetCoeff(fz,1,1);

    // polynomial for algebraic integers
    conv(t,1);
    for (long i=deg(f); ; ) {
      SetCoeff(fa,i,t*coeff(f,i));
      if (--i<0)
	break;
      t*=-d;
    }

    // polynomial to sieve over
    mul(fs,fz,fa);

    // skip even c's if d is even
    if ((d%2)==0) {
      // fz = 2x+1
      SetCoeff(fz,0,1);
      SetCoeff(fz,1,2);
      ZZX _fs(fs);
      compose(fs,_fs,fz);
    }

    // remove any common smooth factor from coefficients
    zfb.reduce(fs,fs);

    // start of sieve
    long start=-(sieve_width/2);

    // initialize sieve to indicate which values of c we're interested in
    vec_pair_long_long excludes;
    clear(sieve);
    if (d>2) {
      long rm;
      zfb.factor(zfact,rm,d);
      for (long i=1; i<zfb.length(); ++i) {  // ignore p=2
        if (zfact[i]>0) {
          // remove prime
          long p = zfb[i];
          if ((d%2)==1) {
            // start+i = 0 mod p  =>  i = -start mod p;
            for (long i=rem(-start,p); i<sieve_width; i+=p)
              sieve[i] = -1;
            append(excludes,p,0);
          }
          else {
            // 2*(start+i)+1 = 0 mod p  =>   i = -1/2 - start mod p
            for (long i=rem(-InvMod(2,p)-start,p); i<sieve_width; i+=p)
              sieve[i] = -1;
            // 2x+1 = 0 mod p  =>  x = -1/2 mod p
            append(excludes,p,InvMod(p-2,p));
          }
        }
      }
      if (rm>1) {
        if ((d%2)==1) {
          sieve[-start] = -1;
          for (long i=0; i<sieve_width; ++i)
            if (GCD(start+i,rm)!=1)
              sieve[i] = -1;
        }
        else {
          for (long i=0; i<sieve_width; ++i)
            if (GCD(2*(start+i)+1,rm)!=1)
              sieve[i] = -1;
        }
      }
    }
    else
      sieve[-start] = -1;

    // sieve for smooth c+d*m and c+d*alpha
    zfb.sieve(sieve,fs,to_ZZ(start),ZZ::zero(),0,0,excludes);
    ++sieve_count;      

    // check results
    for (long i=0; i<sieve_width; ++i) {
      if (sieve[i]==1) {
	++sieve_total;
	long c;
	if ((d%2)==1)
	  c = start+i;
	else
	  c = 2*(start+i)+1;

	// smooth integer
	if (!zfb.factor(zfact,c+d*m)) {
	  cerr<<"DLog_NFS::log_prime() WARNING: bad smooth rational  \n";
	  ++sieve_bad;
	  continue;
	}

	// smooth algebraic integer
	if (!afb.factor(afact,c,d)) {
	  cerr<<"DLog_NFS::log_prime() WARNING: bad smooth algebraic  \n";
	  ++sieve_bad;
	  continue;
	}

	// character map
	afb.CharacterMap(l,c,d);
	
	rels.add(zfact,afact,l);
	++sieve_smooth;

#ifdef NFS_CREATE_DIST_FILE
	if (dist_file)
	  fprintf(dist_file,"%ld %ld\n",d,c);;
#endif

	d_lastsmooth=d;
  
	if (rels.done()) {
	  done=true;
	  break;
	}
      }
      else if (sieve[i]>=0) {
	// not-smooth (we count it as having been checked)
	++sieve_total;
      }
    }
      
    long percent = sieve_smooth*100/(aprimes+zprimes);
    if ((percent!=percent_last)||(d-d_last>20)) {
      cout<<"Sieving: [d="<<d<<"] "
	  <<sieve_count<<" sieves, "
	  <<sieve_smooth<<"/"<<sieve_total<<" smooth";
      if (!done)
	cout<<" ("<<percent<<"%)  \r";
      else
	cout<<"         \r";
      cout.flush();
      percent_last=percent;
      d_last = d;
    }
    
    if ((d>10)&&(d_lastsmooth<d/2)) {
      cerr<<"\nDLog_NFS::log_prime() sieving is taking too long!\n";
      return ZZ::zero();
    }

  } while (!done);
  cout<<"\n";

  // for later
  sieve_length = d;

  // reporting
  if (sieve_bad)
    cerr<<"DLog_NFS::log_prime() WARNING: "<<sieve_bad
        <<" non-smooth integers from sieve\n";

#ifdef NFS_CREATE_DIST_FILE
  if (dist_file)
    fclose(dist_file);
#endif

  // free memory used by sieve
  sieve.kill();

  double t2=GetTime();
  if (VERBOSE)
    cout<<"DLog_NFS::log_prime() sieving took "<<(t2-t1)<<" seconds\n";

  // constant vector
  vec_ZZ b;
  b.SetLength(k+aprimes+zprimes);
  conv(b[k+aprimes],-h);  // rational prime 2
  if (!afb.factor(afact,0,1)) {
    cerr<<"DLog_NFS::log_prime() alpha is not smooth\n";
    return ZZ::zero();
  }
  afb.CharacterMap(l,0,1);
  for (long i=0; i<k; ++i)
    conv(b[i],-l[i]);
  for (long i=0; i<aprimes; ++i)
    conv(b[k+i],-afact[i]);

  //cout<<"b="<<b<<"\n";

  // solve matrix equation
  vec_ZZ X;
  if (!ls_solve(X,rels,b,q)) {
    cerr<<"DLog_NFS::log_prime() solution to linear system not found\n";
    return ZZ::zero();
  }

  //cout<<"X="<<X<<"\n";

  // the logarithm
  ZZ lg;
  lg = log_complete(mp,X[0],q);
  if (!IsZero(lg)) {
    //    if (issmall)
    //  --lg;
    return lg;
  }

  cerr<<"DLog_NFS::log_prime() algorithm failed!\n";
  return ZZ::zero();
}


NTL_END_IMPL;
