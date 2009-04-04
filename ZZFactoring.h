#ifndef NTL_ZZFactoring__H
#define NTL_ZZFactoring__H

#include <NTL/ZZ.h>
#include "pair_ZZ_long.h"


/* Methods for factoring integers.
 *
 * Written by: Chris Studholme
 * Copyright:  GPL (http://www.fsf.org/copyleft/gpl.html)
 */

NTL_OPEN_NNS;


/* Prove that n is prime.  If n<0, then the absolute value of n is tested.
 * If n is 0 or 1, then 0 is returned.
 *
 * Not implemented; just calls ProbPrime() for now.
 */
long ProvePrime(const ZZ& n);


/* General purpose factoring of integer n using best known methods.  If n<0,
 * the absolute value of n will be factored.  If n==0 or n==1, factors will
 * be set to zero length.
 *
 * When finished, the factors vector will contain the factors and exponents 
 * sorted in increasing order by factor.  
 * 
 * If bnd>0, the method only runs long enough to find prime factors less than
 * bnd with probability ???.  In this case, the method may inadvertantly find
 * larger prime factors.  Furthermore, large composite factors may contain
 * primes <=bnd.
 *
 * If deterministic!=0, all factors will have their primality proved with
 * ProvePrime().  If bnd>0 then only factors <=bnd will have had their 
 * primality proved.
 */
void factor(vec_pair_ZZ_long& factors, const ZZ& n,
            const ZZ& bnd=ZZ::zero(), 
	    long deterministic=0, long verbose=0);


/* Recursively factor n by Pollard Rho method and update factors. 
 *
 * The factors vector is assumed to either be zero length or already contain
 * a properly sorted list of factors.
 *
 * Also, n must satisfy n>1.
 *
 * See factor() for details regarding bnd and deterministic.
 */
void PollardRho(vec_pair_ZZ_long& factors, const ZZ& n,
		const ZZ& bnd=ZZ::zero(), 
		long deterministic=0, long verbose=0);


NTL_CLOSE_NNS;

#endif
