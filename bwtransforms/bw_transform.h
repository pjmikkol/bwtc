/*****************************************************************************
 *                                                                           *
 * Base classes (interfaces) for Burrows-Wheeler transform and its reversal. *
 *                                                                           *
 *****************************************************************************/

#ifndef BWTC_BW_TRANSFORM_H_
#define BWTC_BW_TRANSFORM_H_

#include <cassert>

#include <algorithm> // for reverse
#include <vector>

#include "../block.h"
#include "../globaldefs.h"

namespace bwtc {

/*****************************************************************************
 * For implementing new algorithm for Burrows-Wheeler Transform one needs to *
 * inherit BWTransform.                                                      *
 *                                                                           *
 * Transform will be used in the following way:                              *
 *                                                                           *
 *   BWTransform* tranform = GiveTransform(...);                             *
 *   MainBlock* data; uint64 eob_byte;                                       *
 *   transform->Connect(data);                                               *
 *   transform->BuildStats();                                                *
 *   while( std::vector<byte>* result = transform->DoTransform(&eob_byte) {  *
 *     ...do something with stats and a part of a transform                  *
 *                                                                           *
 *****************************************************************************/
// TODO: If there is need for optimized memor management, then transformer
//       needs to be connected to some manager-object  
class BWTransform {
 public:
  BWTransform() : current_block_(NULL) {}
  virtual ~BWTransform() {}
  
  virtual void Connect(MainBlock* block) {
    current_block_ = block;
    std::reverse(current_block_->begin(), current_block_->end());
  }
  virtual std::vector<byte>* DoTransform(uint64* eob_byte) = 0;
  virtual void BuildStats();

  virtual uint64 MaxSizeInBytes(uint64 block_size) const = 0;
  virtual uint64 MaxBlockSize(uint64 memory_budget) const = 0;
  virtual uint64 SuggestedBlockSize(uint64 memory_budget) const = 0;

 protected:
  /* If some later stage we want to implement external memory manager ...*/
  virtual std::vector<byte>* AllocateMemory(uint64 block_size);
  MainBlock* current_block_;

 private:
  BWTransform(const BWTransform&);
  const BWTransform& operator=(const BWTransform& );
};

/* Block size and memory budget would probably be suitable parameters... */
BWTransform* GiveTransformer(char transform);

/***********************************************************************
 * BwtFromSuffixArray computes the BWT from the suffix array and       *
 * writes it to output. Caller has to be sure that output has at least *
 * length of block_size + 1                                            *
 ***********************************************************************/
template <typename T>
uint64 BwtFromSuffixArray(const byte* block, int64 block_size,
                          const T* suffix_array, byte* output)
{
  byte ch = '\0';
  T eob_position = 0;
  for (T rank = 0; rank < block_size + 1; ++rank) {
    T suffix = suffix_array[rank];
    assert(suffix <= block_size);
    if (suffix != 0) {
      ch = block[suffix - 1];
    } else {
      /* This is the end-of-block (eob) character, which cannot be
       * directly written since it has no code.
       * Instead, write a copy of the previous character here and
       * store this position in order to indicate its location. */
      eob_position = rank;
    }
    *output++ = ch;
  }
  return eob_position;
}

} //namespace bwtc

#endif
