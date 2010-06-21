/*****************************************************************************
 *                                                                           *
 * Base classes (interfaces) for Burrows-Wheeler transform and its reversal. *
 *                                                                           *
 *****************************************************************************/

#ifndef BWTC_BW_TRANSFORM_H_
#define BWTC_BW_TRANSFORM_H_

#include <cassert>


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
 *   transform->SetContextLength(...);                                       *
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
  
  virtual void Connect(MainBlock* block) { current_block_ = block;}
  virtual byte* DoTransform(uint64* eob_byte) = 0;
  virtual void BuildStats() = 0;

  /* If some later stage we want to implement external memory manager ...*/
  virtual byte* AllocateMemory(uint64 block_size);
  virtual uint64 MaxSizeInBytes(uint64 block_size) const = 0;
  virtual uint64 MaxBlockSize(uint64 memory_budget) const = 0;
  virtual uint64 SuggestedBlockSize(uint64 memory_budget) const = 0;

 protected:
  MainBlock* current_block_;

 private:
  BWTransform(const BWTransform&);
  const BWTransform& operator=(const BWTransform& );
};

/* Block size and memory budget would probably be suitable parameters... */
BWTransform* GiveTransformer();

/* Computes BWT to output. Caller has to be sure that output is at least *
 * length of block_size + 1 */
uint32 BwtFromSuffixArray(const byte* block, uint32 block_size,
                          const uint32* suffix_array, byte* output);

} //namespace bwtc

#endif
