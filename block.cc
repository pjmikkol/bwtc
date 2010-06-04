#include <vector>

#include <boost/cstdint.hpp>

#include "block.h"
//#include "globaldefs.h"

#include <vector>

namespace bwtc {

Block::Block(std::vector<char>* block) :
    block_(block), filled_(0) {}

} //namespace bwtc
