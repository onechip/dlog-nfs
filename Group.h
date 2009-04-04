#ifndef _GROUP_H_
#define _GROUP_H_
 
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>

NTL_OPEN_NNS;

/* Group and GroupElement abstract interfaces to a general group.
 * GroupElement represents an individual element of a general group
 * and supports tests for identity and equality, transformation to
 * the inverse element, combination with another element via the group
 * operation, and repeated application (power) or the group operation.
 *
 * Also, there should be a one-to-one mapping from group elements to
 * the integers.  This mapping should be such that it is possible to
 * parition the group elements into a number of similarly sized subsets.
 * n subsets can be computed by index() (mod n).
 *
 * Written by: Chris Studholme
 * Copyright:  GPL (http://www.fsf.org/copyleft/gpl.html)
 */
class GroupElement {
public:
  virtual ~GroupElement() {}

  // test if this element is the identity
  virtual bool isIdentity() const=0;

  // test two group elements for equality
  virtual bool operator==(const GroupElement& other) const=0;

  // test for inequality (default is not ==)
  inline bool operator!=(const GroupElement& other) const {
    return !(*this==other);
  }

  // multiply this group element by another
  virtual GroupElement& operator*=(const GroupElement& other)=0;

  // set this to this ^ exponent
  virtual void power(const ZZ& exponent)=0;
 
  // set this to this ^ exponent
  virtual void power(long exponent)=0;

  // invert this group element
  virtual void invert()=0;

  // completely arbitrary one-to-one mapping from G -> Z
  virtual ZZ index() const=0;

  // set this to group element associate with index
  virtual void set(const ZZ& index)=0;
};


/* Group and GroupElement abstract interfaces to a general group.
 * Group provides methods necessary for dealing with an entire group.
 * Methods are available for creating a new identity element, copy of
 * some other group element, an element from an index, and a random 
 * group element.  Also, a read-only reference to an identity element
 * is provided along with the order of the group.  This interface is
 * intended to be used with finite groups of known order.
 */
class Group {
public:
  virtual ~Group() {}

  // create new group element that is identity
  virtual GroupElement* newElementIdentity()=0;

  // create new group element associated with index
  virtual GroupElement* newElementFromIndex(const ZZ& index)=0;

  // create new group element that is a copy of some other element
  virtual GroupElement* newElementCopy(const GroupElement& other)=0;

  // create new random group element
  virtual GroupElement* newElementRandom()=0;

  // number of elements in group
  virtual ZZ size()=0;
};


/* ZZ_p_GE objects are congurency classes mod p (for p some prime).  The
 * identity is 1 and the group operation is multiplication.
 */
class ZZ_p_GE : public GroupElement {
private:
  ZZ_p value;

public:
  ZZ_p_GE() {
    value=1;
  }
  ZZ_p_GE(const ZZ_p& v) {
    value=v;
  }
  ZZ_p_GE(const ZZ_p_GE& other) {
    value=other.value;
  }
  virtual ~ZZ_p_GE() {
  }

  // test if this element is the identity
  virtual bool isIdentity() const {
    return IsOne(value);
  }

  // test two group elements for equality
  virtual bool operator==(const GroupElement& other) const {
    return value==(*(ZZ_p_GE*)&other).value;
  }

  // multiply this group element by another
  virtual GroupElement& operator*=(const GroupElement& other) {
    value*=(*(ZZ_p_GE*)&other).value;
    return *this;
  }

  // set this to this ^ exponent
  virtual void power(const ZZ& exponent) {
    NTL::power(value,value,exponent);
  }
 
  // set this to this ^ exponent
  virtual void power(long exponent) {
    NTL::power(value,value,exponent);
  }

  // invert this group element
  virtual void invert() {
    value = inv(value);
  }

  // completely arbitrary one-to-one mapping from G -> Z
  virtual ZZ index() const {
    return rep(value);
  }

  // set this to group element associate with index
  virtual void set(const ZZ& index) {
    value = to_ZZ_p(index);
  }

  // make this a random group element
  virtual void makeRandom() {
    value = random_ZZ_p();
  }
};


/* ZZ_p_Group is the group of congruency classes mod p (p some prime)
 * combined with the multiplication operation and identity 1.  The
 * order of the group is simply p-1.
 */
class ZZ_p_Group : public Group {
public:
  ZZ_p_Group() {
  }
  virtual ~ZZ_p_Group() {
  }

  // create new group element that is identity
  virtual GroupElement* newElementIdentity() {
    return new ZZ_p_GE();
  }

  // create new group element associated with index
  virtual GroupElement* newElementFromIndex(const ZZ& index) {
    return new ZZ_p_GE(to_ZZ_p(index));
  }

  // create new group element that is a copy of some other element
  virtual GroupElement* newElementCopy(const GroupElement& other) {
    return new ZZ_p_GE(*(ZZ_p_GE*)&other);
  }

  // create new random group element
  virtual GroupElement* newElementRandom() {
    ZZ_p_GE* g = new ZZ_p_GE();
    g->makeRandom();
    return g;
  }

  // number of elements in group
  virtual ZZ size() {
    return ZZ_p::modulus()-1;
  }

};

NTL_CLOSE_NNS;

#endif
