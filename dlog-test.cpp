#include <stdio.h>
#include <time.h>

#include "IndexCalculus.h"
#include "NumberFieldSieve.h"


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

  // initialize random number generator
  SetSeed(to_ZZ(time(0)));

  // size of modulus (in bits)
  long k=56;

  // choose a random modulus
  ZZ p,q;
  GenGermainPrime(q,k-1);  // q has one less bit than p
  p = 2*q+1;
  ZZ_p::init(p);
  cout<<"Modulus is: "<<p<<"\n";

  // factorization of p-1
  vec_pair_ZZ_long f;
  f.SetLength(2);
  f[0].a = 2;  f[0].b = 1;
  f[1].a = q;  f[1].b = 1;

  // choose a base (generator)
  ZZ_p g;
  do {
    random(g);
  } while (!isGenerator(g,f));
  cout<<"Base is: "<<g<<"\n";

  // set verboseness (set to 1 for more output)
  DLog_IC_Base::VERBOSE=0;

  // initialize index calculus
  double start = GetTime();
  DLog_NFS dlog(g);
  cout<<"Setup time: "<<(GetTime()-start)<<" seconds\n";

  for (long i=0; i<5; ++i) {

    // choose a power to find logarithm for
    ZZ_p y;
    random(y);
    
    // find discrete logarithm
    cout<<"Finding logarithm of: "<<y<<"\n";
    ZZ x;
    start = GetTime();
    x = dlog.log(y); 
    cout<<"Time to find logarithm: "<<(GetTime()-start)<<" seconds\n";
    // verify logarithm
    ZZ_p yy;
    power(yy,g,x);
    if (y==yy)
      cout<<"Logarithm is: "<<x<<"\n";
    else 
      cout<<"Logarithm ("<<x<<") is incorrect!\n";
    if (x<0)
      break;
  }

  return 0;
}
