/***********************************************************************
 * Burrows-Wheeler transform using the Karkkainen-Burkhardt algorithm  *
 * which is based on difference cover sampling.                        *
 *                                                                     *
 * Original implementation of the algorithm can be found at:           *
 * http://code.google.com/p/dcs-bwt-compressor/                        *
 ***********************************************************************/
#ifndef DCBW_TRANSFORM_H_
#define DCBW_TRANSFORM_H_

#include "../block.h"
#include "bw_transform.h"

namespace bwtc {

class DCBWTransform : public BWTransform {
 public:
  DCBWTransform(int period);
  virtual ~DCBWTransform();
  
  virtual std::vector<byte>* DoTransform(uint64* eob_byte);

  virtual uint64 MaxSizeInBytes(uint64 block_size) const;
  virtual uint64 MaxBlockSize(uint64 memory_budget) const;
  virtual uint64 SuggestedBlockSize(uint64 memory_budget) const;

 private:
  static const uint64 kMemoryOverhead = (1 << 20);
  static const int kMinLogPeriod = 3;
  static const int kMaxLogPeriod = 15;
  int log_period_;
  uint32 period_;
};

} //namespace bwtc

#endif
