#include <stdio.h>
#include <time.h>

#include "IndexCalculus.h"
#include "NumberFieldSieve.h"

using namespace NTL;

// test if g is a generator of the group ZZ_p
// f is the factorization of p-1
bool isGenerator(const ZZ_p& g, const vec_pair_ZZ_long& f) {
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


// simple test program for computing discrete logarithms using index calculus
int main(int argc, char* argv[]) {

  // L2 cache size of sieving
  FactorBase::CACHE_SIZE = 128*1024;

  // initialize random number generator
  SetSeed(to_ZZ(time(0)));

  // size of modulus (in bits)
  long k=64;

  // choose a random modulus
  ZZ p,q;
  GenGermainPrime(q,k-1);  // q has one less bit than p
  p = 2*q+1;
  ZZ_p::init(p);
  std::cout<<"Modulus: "<<p<<std::endl;

  // factorization of p-1
  vec_pair_ZZ_long f;
  f.SetLength(2);
  f[0].a = 2;  f[0].b = 1;
  f[1].a = q;  f[1].b = 1;

  // choose a base (generator)
  ZZ_p g;
  PrimeSeq s;
  do {
    conv(g,s.next());
  } while (!isGenerator(g,f));
  std::cout<<"Base:    "<<g<<std::endl;

  // set verboseness (set to 1 for more output)
  DLog_IC_Base::VERBOSE=0;

  // initialize NFS
  double start = GetTime();
  DLog_NFS dlog(g);
  std::cout<<"Setup time: "<<(GetTime()-start)<<" seconds"<<std::endl;

  for (long i=0; i<5; ++i) {

    // choose a power to find logarithm for
    ZZ_p y;
    if (i==0) {
      conv(y,2);
      if (y==g)
        ++y;
    }
    else {
      random(y);
    }

    // find discrete logarithm
    std::cout<<std::endl;
    std::cout<<"Finding logarithm of: "<<y<<std::endl;
    ZZ x;
    start = GetTime();
    x = dlog.log(y); 
    std::cout<<"Time to find logarithm: "<<(GetTime()-start)<<" seconds"
             <<std::endl;
    // verify logarithm
    ZZ_p yy;
    power(yy,g,x);
    if (y==yy)
      std::cout<<"Logarithm is: "<<x<<std::endl;
    else 
      std::cout<<"Logarithm ("<<x<<") is incorrect!"<<std::endl;
    if (x<0)
      break;
  }

  return 0;
}
