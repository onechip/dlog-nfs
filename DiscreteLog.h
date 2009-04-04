#ifndef _DISCRETELOG_H_
#define _DISCRETELOG_H_
 
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include "ZZFactoring.h"

#include "Group.h"

/* Classes for discrete logarithm algorithms.
 *
 * Written by: Chris Studholme
 * Copyright:  GPL (http://www.fsf.org/copyleft/gpl.html)
 */

NTL_OPEN_NNS;

/* Abstract interface to algorithms for finding the discrete log in
 * arbitrary groups and in the group ZZ_p.
 */
class DiscreteLog {
protected:
  static ZZ_p_Group ZZ_p_group;
  Group& group;
  GroupElement* base;

public:
  // set group to ZZ_p and leave base unset
  DiscreteLog() : group(ZZ_p_group) {
    base=NULL;
  }
  // set group to ZZ_p and set base
  DiscreteLog(const ZZ_p& base) : group(ZZ_p_group) {
    this->base = new ZZ_p_GE(base);
  }
  // set group to some general group and leave base unset
  DiscreteLog(Group& g) : group(g) {
    base=NULL;
  }
  // set group to some general group and set base
  DiscreteLog(Group& g, const GroupElement& base) : group(g) {
    this->base = group.newElementCopy(base);
  }
  virtual ~DiscreteLog() {
    if (base)
      delete base;
  }

  // set base to an element of a general group
  virtual void setBase(const GroupElement& base) {
    if (this->base) delete this->base;
    this->base = group.newElementCopy(base);
  }

  // set base to an element of ZZ_p
  virtual void setBase(const ZZ_p& base) {
    if (this->base) delete this->base;
    this->base = new ZZ_p_GE(base);
  }

  // find log of a power from a general group
  // returns -1 if the discrete log does not exist
  // may return -2 if some limit was exceeded before the log was found
  virtual ZZ log(const GroupElement& power)=0;

  // find log of a power from ZZ_p
  // returns -1 if the discrete log does not exist
  // may return -2 if some limit was exceeded before the log was found
  virtual ZZ log(const ZZ_p& power);

  // check discrete log calculation (returns true if all is ok)
  bool verify(const GroupElement& power, const ZZ& log) const;

  // check discrete log calculation (returns true if all is ok)
  bool verify(const ZZ_p& power, const ZZ& log) const;
  
};

/* Dumb way to find discrete logarithms.  Just start at exponent 1 and
 * loop through the powers of base one by one until the logarithm is
 * found.
 *
 * This method works with any general group.
 */
class DLog_dumb : public DiscreteLog {
private:
  long limit;  // run-time limit in seconds

public:
  DLog_dumb() : DiscreteLog() { 
    limit=0; 
  }
  DLog_dumb(const ZZ_p& base) : DiscreteLog(base) { 
    limit=0; 
  }
  DLog_dumb(Group& g) : DiscreteLog(g) { 
    limit=0; 
  }
  DLog_dumb(Group& g, const GroupElement& base) : DiscreteLog(g,base) { 
    limit=0;
  }
  DLog_dumb(long limit) : DiscreteLog() { 
    this->limit=limit;
  }
  DLog_dumb(const ZZ_p& base, long limit) : DiscreteLog(base) { 
    this->limit=limit;
  }
  DLog_dumb(Group& g, long limit) : DiscreteLog(g) { 
    this->limit=limit;
  }
  DLog_dumb(Group& g, const GroupElement& base, long limit) 
    : DiscreteLog(g,base) { 
    this->limit=limit;
  }

  void setLimit(long limit) {
    this->limit=limit;
  }

  // returns -1 if the discrete log does not exist (or limit was reached)
  virtual ZZ log(const GroupElement& power);
};

/* Use the Pollard Rho algorithm to find the discrete log in a general
 * group.  The group order is factored (also using the Pollard Rho method),
 * then the discrete log is found modulo each factor and these logs are
 * combined using the Chinese Remainder Theorem.
 *
 * The run time (in seconds) of this algorithm can be limited, but its
 * measurement is per factor.
 */
class DLog_Pollard : public DiscreteLog {
private:
  long limit;  // run-time limit in seconds
  vec_pair_ZZ_long factors;

public:
  DLog_Pollard() : DiscreteLog() { 
    limit=0; 
    factor(factors,group.size());
  }
  DLog_Pollard(const ZZ_p& base) : DiscreteLog(base) { 
    limit=0; 
    factor(factors,group.size());
  }
  DLog_Pollard(Group& g) : DiscreteLog(g) { 
    limit=0; 
    factor(factors,group.size());
  }
  DLog_Pollard(Group& g, const GroupElement& base) : DiscreteLog(g,base) { 
    limit=0;
    factor(factors,group.size());
  }
  DLog_Pollard(long limit) : DiscreteLog() { 
    this->limit=limit;
    factor(factors,group.size());
  }
  DLog_Pollard(const ZZ_p& base, long limit) : DiscreteLog(base) { 
    this->limit=limit;
    factor(factors,group.size());
  }
  DLog_Pollard(Group& g, long limit) : DiscreteLog(g) { 
    this->limit=limit;
    factor(factors,group.size());
  }
  DLog_Pollard(Group& g, const GroupElement& base, long limit) 
    : DiscreteLog(g,base) { 
    this->limit=limit;
    factor(factors,group.size());
  }
  ~DLog_Pollard() {
  }

  // Set run-time limit in seconds.  Set to zero for no limit.
  void setLimit(long limit) {
    this->limit=limit;
  }

  // returns -1 if the discrete log does not exist
  virtual ZZ log(const GroupElement& power);

private:
  // a single step of the algorithm
  static void step(GroupElement& x, ZZ& a, ZZ& b,
		   const GroupElement& base, 
		   const GroupElement& power);

  // raw application of algorithm
  ZZ raw(const GroupElement& base, const GroupElement& power,
	 const ZZ& size) const;
    
};

NTL_CLOSE_NNS;

#endif
