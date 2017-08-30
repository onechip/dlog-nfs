#include <time.h>

#include "DiscreteLog.h"


/* Classes for discrete logarithm algorithms and square-root methods.
 *
 * Written by: Chris Studholme
 * Copyright:  GPL (http://www.fsf.org/copyleft/gpl.html)
 */

NTL_START_IMPL;


/**************** class DiscreteLog ****************/

ZZ_p_Group DiscreteLog::ZZ_p_group;

ZZ DiscreteLog::log(const ZZ_p& power) {
  if (IsZero(power))
    return to_ZZ(-1);
  ZZ_p_GE *p = new ZZ_p_GE(power);
  ZZ result;
  result = log(*p);
  delete p;
  return result;
}

bool DiscreteLog::verify(const GroupElement& power, const ZZ& log) const {
  GroupElement* b = group.newElementCopy(*base);
  b->power(log);
  bool result = (*b==power);
  delete b;
  return result;
}

bool DiscreteLog::verify(const ZZ_p& power, const ZZ& log) const {
  ZZ_p_GE *p = new ZZ_p_GE(power);
  bool result = verify(*p,log);
  delete p;
  return result;
}


/**************** class DLog_dumb ****************/

/* Approximate minimum number of loops per second expressed as a mask
 * to be used with & to do modulo arithmatic
 */ 
#define DLog_dumb_LOOPS_PER_SEC 0xff

ZZ DLog_dumb::log(const GroupElement& power) {
  if (power.isIdentity())
    return to_ZZ(0);
  if (power==*base)
    return to_ZZ(1);

  GroupElement* p = group.newElementCopy(*base);
  GroupElement* power_inv = group.newElementCopy(power);
  power_inv->invert();

  long maxloop = group.size()<NTL_MAX_LONG ? 
    to_long(group.size()) : NTL_MAX_LONG;

  clock_t end = clock() + CLOCKS_PER_SEC*limit;

  long exp = 1;
  long last=-1;
  do {
    long per = exp*100/maxloop;
    if (per!=last) {
      std::cout<<"DLog_dumb::log() "<<per<<"%  \r"<<std::flush;
      last=per;
    }
    *p *= *base;
    ++exp;
    if (*p==power) {
      delete p;
      std::cout<<std::endl;
      return to_ZZ(exp);
    }
    if (*p==*power_inv) {
      delete p;
      std::cout<<std::endl;
      return group.size() - to_ZZ(exp);
    }
    if (p->isIdentity()) {
      delete p;
      std::cout<<std::endl;
      return to_ZZ(-1);
    }
    // check, every now and again, if we are taking too long
    if ((limit>0)&&((exp&DLog_dumb_LOOPS_PER_SEC)==0)&&(clock()>end)) {
      delete p;
      std::cout<<std::endl;
      return to_ZZ(-2);
    }
  } while (true);
}



/**************** class DLog_Pollard ****************/

// private helper method
inline void DLog_Pollard::step(GroupElement& x, ZZ& a, ZZ& b,
			       const GroupElement& base, 
			       const GroupElement& power) {
  long r = x.index() % 3;
  if (r<0) r=-r;
  switch (r) {
  case 0:
    x*=base;
    b+=1;
    break;
  case 1:
    x*=power;
    a+=1;
    break;
  case 2:
    x.power(2);
    a*=2;
    b*=2;
    break;
  }
}

// size must be prime
ZZ DLog_Pollard::raw(const GroupElement& base,
		     const GroupElement& power,
		     const ZZ& size) const {
  if (power.isIdentity())
    return to_ZZ(0);
  if (power==base)
    return to_ZZ(1);
  
  GroupElement* x1 = group.newElementIdentity();
  GroupElement* x2 = group.newElementIdentity();
  ZZ a1, a2, b1, b2;

  long maxloop = 4*SqrRoot(size)<NTL_MAX_LONG ? 
    to_long(4*SqrRoot(size)) : NTL_MAX_LONG;

  clock_t end = clock() + CLOCKS_PER_SEC*limit;

  long count=0;
  long last=-1;
  do {
    long cur = count*100/maxloop;
    if (cur!=last) {
      std::cout<<"DLog_Pollard::raw() "<<cur<<"%  \r"<<std::flush;
      last=cur;
    }
    step(*x1,a1,b1,base,power);
    step(*x2,a2,b2,base,power);
    step(*x2,a2,b2,base,power);
    rem(a1,a1,size);
    rem(b1,b1,size);
    rem(a2,a2,size);
    rem(b2,b2,size);
    // check, every now and again, if we are taking too long
    if ((limit>0)&&((count&DLog_dumb_LOOPS_PER_SEC)==0)&&(clock()>end)) {
      delete x1;
      delete x2;
      std::cout<<std::endl;
      return to_ZZ(-2);
    }
    if (++count>maxloop) {
      // log doesn't exist
      delete x1;
      delete x2;
      std::cout<<std::endl;
      return to_ZZ(-1);
    }
  } while (*x1!=*x2);
  std::cout<<std::endl;

  delete x1;
  delete x2;

  ZZ m,n;
  m = a1 - a2;
  n = b2 - b1;
  
  ZZ d,s,t;
  XGCD(d,s,t,m,size);
  if (d!=1)
    return to_ZZ(-1);

  s*=n;
  rem(s,s,size);
  return s;
}

ZZ DLog_Pollard::log(const GroupElement& power) {
  ZZ size;
  size = group.size();

  ZZ a,p;
  
  for (long i=factors.length()-1; i>=0; --i) {
    const pair_ZZ_long& factor = factors[i];
    ZZ d; 
    NTL::power(d,factor.a,factor.b);
    ZZ s;
    s = size/d;
    GroupElement* bs = group.newElementCopy(*base);
    bs->power(s);
    GroupElement* pw = group.newElementCopy(power);
    pw->power(s);
    ZZ e;
    e = raw(*bs,*pw,d);
    delete pw;
    delete bs;
    if (e<0)
      return e;

    if (IsZero(p)) {
      a = e;
      p = d;
    }
    else {
      // Chinese Remainder Algorithm
      CRT(a,p,e,d);
    }
  }
  if (a<0) a+=p;
  return a;
}


NTL_END_IMPL;

