/*****************************************************************************
 *                                                                           *
 * Base classes (interfaces) for Burrows-Wheeler transform and its reversal. *
 *                                                                           *
 *****************************************************************************/

#ifndef BWTC_BW_TRANSFORM_H_
#define BWTC_BW_TRANSFORM_H_

#include <vector>

#include "../block.h"
#include "../globaldefs.h"
//#include "dcbwt.h"

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
  BWTransform() {}
  virtual ~BWTransform() {}
  
  virtual Connect(MainBlock* block) = 0;
  std::vector<byte>* DoTransform(uint64* eob_byte) = 0;
  virtual void SetContextLength(int length) = 0;
  virtual void BuildStats() = 0;
};

/* Block size and memory budget would probably be suitable parameters... */
BWTransform* GiveTransform() {
  /* When there are multiple ways to do transform this is to place to add them */
  return new DCBWTransform();
}

} //namespace bwtc

#endif
