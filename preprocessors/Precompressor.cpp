/**
 * @file Precompressor.hpp
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
 * Implementation of precompressor.
 */

#include "PairReplacer.hpp"
#include "Precompressor.hpp"
#include "../Profiling.hpp"
#include "../Streams.hpp"

namespace bwtc {

Precompressor::Precompressor(const std::string& preprocessing)
    : m_preprocessingOptions(preprocessing) {}

Precompressor::~Precompressor() {}

PrecompressorBlock*
Precompressor::readBlock(size_t blockSize, RawInStream* in) const {
  PrecompressorBlock *result = new PrecompressorBlock(blockSize, in);
  precompress(*result);
  return result;
}

#define PREPROCESS(Type, verb, src, dst) \
  Type r(grammar, (verb));\
  r.analyseData((src), length);\
  r.finishAnalysation();\
  r.decideReplacements();\
  length = r.writeReplacedVersion((src), length, (dst))

void Precompressor::precompress(PrecompressorBlock& block) const {
  PROFILE("Precompressor::precompress");
  //TODO: use less memory!
  byte *dst, *src = block.begin();
  std::vector<byte> tmp;
  Grammar& grammar = block.grammar();
  size_t length = block.size();

  if(m_preprocessingOptions.size() > 0) {
    tmp.resize(length + 1);
    dst = &tmp[0];
    for(size_t i = 0; i < m_preprocessingOptions.size(); ++i) {
      char c = m_preprocessingOptions[i];
      if(c == 'p') {
        PREPROCESS(PairReplacer, verbosity > 1, src, dst);
      }
      std::swap(src, dst);
    }
  }
  if(m_preprocessingOptions.size() & 1) {
    std::copy(src, src + length, dst);
    src = dst;
  }
  block.setSize(length);
  if(verbosity > 0) {
    std::clog << "Size of precompressed block is " << length << std::endl;
  }

}

} //namespace bwtc
