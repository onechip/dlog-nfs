#ifndef NTLX_vec_long__H
#define NTLX_vec_long__H

#include <NTL/vec_long.h>
#include <cstring>

NTL_OPEN_NNS;

inline bool IsZero(long i) {
  return (i==0);
}

inline void clear(long &i) {
  i=0;
}

inline void clear(vec_long& v) {
  memset(v.elts(),0,v.length()*sizeof(long));
}

NTL_CLOSE_NNS;

#endif
