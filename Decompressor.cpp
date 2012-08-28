/**
 * @file Decompressor.cpp
 * @author Pekka Mikkola <pmikkol@gmail.com>
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
 * Implementation of Decompressor-class.
 */

#include "Decompressor.hpp"
#include "Streams.hpp"
#include "preprocessors/Postprocessor.hpp"
#include "Profiling.hpp"
#include "bwtransforms/InverseBWT.hpp"

#include <cassert>
#include <string>

namespace bwtc {

Decompressor::Decompressor(const std::string& in, const std::string& out)
    : m_in(new RawInStream(in)), m_out(new RawOutStream(out)),
      m_decoder(0) {}

Decompressor::Decompressor(InStream* in, OutStream* out)
    : m_in(in), m_out(out), m_decoder(0) {}

Decompressor::~Decompressor() {
  delete m_in;
  delete m_out;
  delete m_decoder;
}

size_t Decompressor::readGlobalHeader() {
  char entropyDecoder = static_cast<char>(m_in->readByte());
  delete m_decoder;
  m_decoder = giveEntropyDecoder(entropyDecoder);
  return 1;
}

size_t Decompressor::decompress(size_t threads) {
  PROFILE("Decompressor::decompress");
  if(threads != 1) {
    std::cerr << "Supporting only single thread!" << std::endl;
    return 0;
  }
  InverseBWTransform *ibwt = giveInverseTransformer();

  readGlobalHeader();

  size_t preBlocks = 0, bwtBlocks = 0, decompressedSize = 0;
  while(true) {
    PrecompressorBlock *pb = PrecompressorBlock::readBlockHeader(m_in);
    if(pb->originalSize() == 0) {
      delete pb;
      break;
    }
    ++preBlocks;
    bwtBlocks += pb->slices();

    for(size_t i = 0; i < pb->slices(); ++i) {
      pb->getSlice(i).setBegin(pb->end());
      m_decoder->decodeBlock(pb->getSlice(i), m_in);
      pb->usedAtEnd(pb->getSlice(i).size());
      ibwt->doTransform(pb->getSlice(i));
    }
    // Postprocess pb
    Postprocessor postprocessor(verbosity > 1, pb->grammar());
    size_t postSize = postprocessor.uncompress(pb->begin(), pb->size(), m_out);
    decompressedSize += postSize;
    //m_out->writeBlock(pb->begin(), pb->end());
    assert(postSize == pb->originalSize());
    delete pb;
  }
  delete ibwt;
  return decompressedSize;
}

} //namespace bwtc
