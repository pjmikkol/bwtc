#ifndef GLOBAL_DEFS_H_
#define GLOBAL_DEFS_H_

/* Common typedefs for different sizes of integers */
#include <boost/cstdint.hpp>
typedef boost::int64_t int64;
typedef boost::int32_t int32;
typedef boost::uint32_t uint32;
typedef boost::int16_t int16;
typedef boost::uint16_t uint16;

/* Definitions for probability models and arithmetic coding */
typedef uint16 Probability;
static const int kLogProbabilityScale = 12;
static const Probability kProbabilityScale = (1 << kLogProbabilityScale);

typedef unsigned char byte;

// Name of the compressor program
#define COMPRESSOR "compr"
// Name of the decrompressor program
#define DECOMPRESSOR "uncompr"

namespace bwtc {
  extern int verbosity;
}

#endif
