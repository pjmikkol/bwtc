/* 
 *
 */


#ifndef BWTC_BLOCK_H_
#define BWTC_BLOCK_H_

#include <vector>

#include "globaldefs.h"

namespace bwtc {

// TODO: If Block isn't going to be more complex then model it as a struct
class Block {
 public:
  Block(std::vector<char>* block, std::vector<char>::iterator filled);
  ~Block();

  /* Iterators are meant for reading and writing of a block.*/
  std::vector<char>::iterator begin() {
    return block_->begin();
  }
  /* end() returns iterator one past the valid range to a vector.
   * Uses filled_ for deducing the value.*/
  std::vector<char>::iterator end() {
    return block_->end();
  }

 private:
  /* Vector will be initialized to maximum block size given */
  std::vector<char>* block_;
  /* Defines range [0, filled_) in vector, which will hold the relevant data. */
  std::vector<char>::iterator filled_;

  Block& operator=(const Block& b);
  Block(const Block&);
  
};

} // namespace bwtc

#endif
