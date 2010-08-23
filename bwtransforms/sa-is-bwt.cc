#include <cassert>

#include <vector>

#include "bw_transform.h"
#include "sa-is-bwt.h"
#include "sais.hxx"
#include "../block.h"
#include "../globaldefs.h"

namespace bwtc {

SAISBWTransform::SAISBWTransform() {}

std::vector<byte>* SAISBWTransform::DoTransform(uint64 *eob_byte) {
  if(!current_block_) return NULL;

  int block_size = current_block_->Size();
  current_block_->Append(0);
  byte *block = current_block_->begin();
  /* Whole transformation is done in single pass. */
  current_block_ = NULL;

  std::vector<byte> *result = AllocateMemory(block_size + 1);
  std::vector<int> suffix_array(block_size + 1);
  //sa_is::SA_IS_ZeroInclude(block, &suffix_array[0], block_size + 1, 256);
  /*  *eob_byte = sa_is::sais(block, &suffix_array[0], 0, block_size + 1, 256, false);
  for(int i = 0; i < block_size + 1; ++i) {
    (*result)[i] = block[suffix_array[i]];
    }*/
  *eob_byte = saisxx_bwt(block, &(*result)[0], &suffix_array[0], block_size);
  //  *eob_byte = BwtFromSuffixArray(block, block_size, &suffix_array[0],
  //                             &(*result)[0]);
  return result;
}

} //namespace bwtc
