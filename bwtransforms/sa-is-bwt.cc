#include <cassert>

#include <vector>

#include "bw_transform.h"
#include "sa-is-bwt.h"
#include "../block.h"
#include "../globaldefs.h"

namespace bwtc {

SAISBWTransform::SAISBWTransform() {}

std::vector<byte>* SAISBWTransform::DoTransform(uint64 *eob_byte) {
  if(!current_block_) return NULL;

  int64 block_size = current_block_->Size();
  current_block_->Append(0);
  byte *block = current_block_->begin();
  /* Whole transformation is done in single pass. */
  current_block_ = NULL;
  std::vector<byte> *result = AllocateMemory(block_size + 1);
  std::vector<uint32> suffix_array(block_size + 1);
  sa_is::SA_IS_ZeroInclude(block, &suffix_array[0], block_size + 1, 256);

  *eob_byte = BwtFromSuffixArray(block, block_size, &suffix_array[0],
                                 &(*result)[0]);
  return result;
}

} //namespace bwtc
