#include <vector>

#include "../globaldefs.h"

//TODO: PostProcessor-class
//      At the moment only algorithms for postprocessing exists

namespace bwtc {
  
uint64 UncompressCommonPairs(std::vector<byte> *from, uint64 length);
uint64 UncompressLongRuns(std::vector<byte> *from, uint64 length);

}
