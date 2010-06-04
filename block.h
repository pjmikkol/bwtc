/* 
 *
 */


#ifndef BWTC_BLOCK_H_
#define BWTC_BLOCK_H_

#include <vector>
#include <boost/cstdint.hpp>

namespace bwtc {

class Block {
 public:
  Block(boost::int64_t max_block_size);
  ~Block();

  /* Iterators are meant for reading and writing of a block.*/
  std::vector<char>::iterator begin();
  /* end() returns iterator one past the valid range to a vector.
   * Uses filled_ for deducing the value.*/
  std::vector<char>::iterator end();

 private:
  /* Vector will be initialized to maximum block size given */
  std::vector<char>* block_;
  /* Defines range [0, filled_) in vector, which will hold the relevant data.
   * Actual size of data that block contains. */
  boost::int64_t filled_;
  
};

} // namespace bwtc

#endif
