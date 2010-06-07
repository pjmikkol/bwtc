#include <vector>

#include <boost/cstdint.hpp>

#include "block.h"
//#include "globaldefs.h"

#include <vector>

namespace bwtc {

Block::Block(std::vector<char>* block, std::vector<char>::iterator filled) :
    block_(block), filled_(filled) {}

} //namespace bwtc
