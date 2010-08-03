#ifndef BWTC_LONGSEQUENCES_H_
#define BWTC_LONGSEQUENCES_H_

#include <vector>

#include "../globaldefs.h"

namespace bwtc {

namespace long_sequences {

struct long_seq {
  long_seq() : position(0), length(0), count(0) {}
  long_seq(unsigned c, int l, uint64 p) :
      position(p), length(l), count(c) {}
  
  bool operator<(const long_seq& ll) const {
    return count*length < ll.count*ll.length;
  }

  uint64 position;
  int length;
  unsigned count;
};

/* To avoid costly memory-allocations we allocate the array used in *
 * computations of Borders only once.                               */
template <typename T>
class Border {
 public:
  explicit Border(unsigned length) : length_(length) {
    border_ = new int[length_ + 1];
  }

  ~Border() {
    delete [] border_;
  }

  int operator()(const T *source) {
    border_[0] = -1;
    int t = -1;
    for(unsigned j = 1; j <= length_; ++j) {
      while(t >= 0 && source[t] != source[j-1]) t = border_[t];
      border_[j] = ++t;
    }
    return border_[length_];
  }

 private:
  int *border_;
  unsigned length_;
};

} //namespace long_sequences

uint64 CompressSequences(byte *from, uint64 length, int memory_constraint,
                         unsigned window_size, int threshold);

} 

#endif
