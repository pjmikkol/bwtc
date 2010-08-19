/*        Part of this code was originally released at        * 
 *        http://code.google.com/p/dcs-bwt-compressor/        */
#include <cassert>

#include <vector>

#include "../globaldefs.h"
#include "bw_transform.h"
#include "dcbwt.h"
#include "sa-is-bwt.h"

namespace bwtc {

/* We assume that the context-block of sentinel char is at the front of *
 * transform.                                                           */
//TODO: Allow take array from preprocessor to parameter since they have
//      already built stats
void BWTransform::BuildStats() {
  std::fill(current_block_->stats_->begin(), current_block_->stats_->end(), 0);
  (*current_block_->stats_)[0] = 1;
  for(uint64 i = 0; i < current_block_->Size(); ++i)
    (*current_block_->stats_)[(*current_block_->block_)[i] + 1]++;
}

std::vector<byte>* BWTransform::AllocateMemory(uint64 size) {
  return new std::vector<byte>(size);
}

BWTransform* GiveTransformer(char transform) {
  /*When there are multiple ways to do transform this is to place to add them*/
  if(transform == 'd')
    return new DCBWTransform(8);
  else
    return new SAISBWTransform();
}

}
