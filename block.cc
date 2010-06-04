#include <vector>

#include <boost/cstdint.hpp>

#include "block.h"

namespace bwtc {

Block::Block(boost::int64_t max_block_size) :
    block_(NULL), filled_(0) {}
