#include <vector>

#include <boost/cstdint.hpp>

#include "block.h"

namespace bwtc {

Block::Block(boost::int64_t max_block_size):
    block_(max_block_size), filled_(0) {}

} // namespace bwtc
