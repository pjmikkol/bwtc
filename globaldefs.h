#ifndef GLOBAL_DEFS_H_
#define GLOBAL_DEFS_H_

#include <boost/cstdint.hpp>
typedef boost::int64_t int64;
typedef boost::uint32_t int32;
typedef boost::int16_t int16;
typedef unsigned char byte;

// Name of the compressor program
#define COMPRESSOR "compr"
// Name of the decrompressor program
#define DECOMPRESSOR "uncompr"

namespace bwtc {
  extern int verbosity;
}

#endif
