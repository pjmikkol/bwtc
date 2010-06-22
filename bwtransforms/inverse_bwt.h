#ifndef BWTC_INVERSE_BWT_H_
#define BWTC_INVERSE_BWT_H_

#include <vector>

#include "../globaldefs.h"

namespace bwtc {
/**************************************************************************
 * Base class for Inverse Burrows-Wheeler transforms. This class has some *
 * common functionality with BWTransform at the moment (especially memory *
 * allocation mechanism) so merging these two together is an option.      *
 *                                                                        *
 * At the moment biggest differece is that in inverse transform we aren't *
 * prepared to make the transform in pieces, so Connecting to the block   *
 * before transform is not required as it is in forward-transform.        *
 **************************************************************************/

class InverseBWTransform {
 public:
  InverseBWTransform() {}
  virtual ~InverseBWTransform() {}
  virtual uint64 MaxBlockSize(uint64 memory_budget) const = 0;
  virtual std::vector<byte>* DoTransform(const byte* source_bwt,
                                         uint64 bwt_size,
                                         uint64 eob_position) = 0;

 protected:
  /* This is here for the sake of memory management */
  virtual std::vector<byte>* AllocateMemory(uint64 block_size);
};

/**************************************************************************
 * FastInverseBWTransformer: This is the algorithm mergeTL in Sewards     *
 * paper with some additional heuristics to be able to handle very large  *
 * blocks.                                                                *
 **************************************************************************/
class FastInverseBWTransform : public InverseBWTransform {
 public:
  FastInverseBWTransform() {}
  virtual ~FastInverseBWTransform() {}
  virtual uint64 MaxBlockSize(uint64 memory_budget) const;
  virtual std::vector<byte>* DoTransform(const byte* source_bwt,
                                         uint64 bwt_size, uint64 eob_position);
 private:
  static const int64 kMemoryOverhead = 1 << 20;
};

/* For example memory budget would be a good parameter.. */
InverseBWTransform* GiveInverseTransformer();

} //namespace bwtc
#endif
