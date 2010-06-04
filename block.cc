#include <vector>

#include <boost/cstdint.hpp>

#include "block.h"
#include "globaldefs.h"

namespace bwtc {

Block::Block(int64 max_block_size) :
    block_(NULL), filled_(0) {}
