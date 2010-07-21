#ifndef BWTC_LONGSEQUENCES_H_
#define BWTC_LONGSEQUENCES_H_

#include <vector>

#include "../globaldefs.h"

namespace bwtc {

namespace long_sequences {

template <typename T>
int Border(T *source, unsigned length) {
  std::vector<int> border(length + 1);
  border[0] = -1;
  int t = -1;
  for(unsigned j = 1; j <= length; ++j) {
    while(t >= 0 && source[t] != source[j-1]) t = border[t];
    border[j] = ++t;
  }
  return border[length];
}

} //namespace long_sequences

uint64 CompressSequences(byte *from, uint64 length, int memory_constraint);

} 

#endif
