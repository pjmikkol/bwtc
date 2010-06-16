/*****************************************************************************
 *                                                                           *
 * Base classes (interfaces) for Burrows-Wheeler transform and its reversal. *
 *                                                                           *
 *****************************************************************************/

#ifndef BWTC_BW_TRANSFORM_H_
#define BWTC_BW_TRANSFORM_H_

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
 *   byte* data; uint64 length_of_data; uint64 eob_byte;                     *
 *   transform->SetContextLength(...);                                       *
 *   transform->Connect(data, length_of_data);                               *
 *   transform->BuildStats(data, length_of_data);                            *
 *   while( std::vector<byte>* result = transform->DoTransform(&eob_byte) {  *
 *     ...do something with stats and a part of a transform                  *
 *                                                                           *
 *****************************************************************************/
class BWTransform {
 public:
  BWTransform() {}
  virtual ~BWTransform() {}
  
  virtual Connect(byte* block, uint64 block_length) = 0;
  uint64 Transform() = 0;
  virtual void SetContextLength(int) = 0;
  virtual void BuildStats();
};

} //namespace bwtc

#endif
