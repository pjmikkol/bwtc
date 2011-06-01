/**
 * @file globaldefs.h
 * @author Pekka Mikkola <pjmikkol@cs.helsinki.fi>
 *
 * @section LICENSE
 *
 * This file is part of bwtc.
 *
 * bwtc is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * bwtc is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with bwtc.  If not, see <http://www.gnu.org/licenses/>.
 *
 * @section DESCRIPTION
 *
 * Global definitions and typedefs for bwtc.
 */

#ifndef GLOBAL_DEFS_H_
#define GLOBAL_DEFS_H_

#include <limits>
#include <boost/cstdint.hpp>

#ifdef MAIN
#define EXTERN
#else
#define EXTERN extern
#endif

namespace bwtc {

EXTERN int verbosity;

/* Common typedefs for different sizes of integers */
typedef boost::int64_t int64;
typedef boost::uint64_t uint64;
typedef boost::int32_t int32;
typedef boost::uint32_t uint32;
typedef boost::int16_t int16;
typedef boost::uint16_t uint16;
typedef boost::uint8_t uint8;

/* Upper limits for different types */
template <typename Integer>
struct Max {
 public:
  static Integer max; 
};
template <typename Integer>
Integer Max<Integer>::max = std::numeric_limits<Integer>::max();

/* Definitions for probability models and arithmetic coding */
typedef uint16 Probability;
static const int kLogProbabilityScale = 12;
static const Probability kProbabilityScale = (1 << kLogProbabilityScale);

typedef unsigned char byte;

} // namespace bwtc

/**
 * Name of the compressor program
 */
#define COMPRESSOR "compr"

/**
 * Name of the decompressor program
 */
#define DECOMPRESSOR "uncompr"

#endif
