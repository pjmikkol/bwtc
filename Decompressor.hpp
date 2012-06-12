/**
 * @file Decompressor.hpp
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
 * Header for Decompressor-class. The decompressor is an abstraction of
 * an decompression pipeline. It has 3 stages which are illustrated in
 * the following diagram:
 *
 *  input -->  ENTROPY DECODING -> INVERSE BWT -> POSTPROCESSING --> ouput
 *
 * The postprocessing phase decompresses the precompressed data.
 * Memory needed for the decompression is determined by the compressor options.
 *
 * For the description of compressed file format @see Compressor.hpp.
 *
 */
#ifndef BWTC_DECOMPRESSOR_HPP_
#define BWTC_DECOMPRESSOR_HPP_

#include "Compressor.hpp"
#include "EntropyCoders.hpp"
#include "Streams.hpp"

#include <string>

namespace bwtc {

class Decompressor {
 public:
  Decompressor(const std::string& in, const std::string& out);
  Decompressor(RawInStream* in, RawOutStream* out);
  ~Decompressor();

  size_t decompress(size_t threads);
  size_t readGlobalHeader();

 private:
  InStream *m_in;
  OutStream *m_out;
  EntropyDecoder *m_decoder;
  Postprocessor m_postprocessor;
};

} //namespace bwtc

#endif
