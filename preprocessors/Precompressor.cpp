/**
 * @file Precompressor.hpp
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
 * Implementation of precompressor.
 */

#include "PairReplacer.hpp"
#include "Precompressor.hpp"
#include "../Profiling.hpp"
#include "../Streams.hpp"

namespace bwtc {

Precompressor::Precompressor() {}

Precompressor::Precompressor(const std::string& preprocessing)
    : m_preprocessingOptions(preprocessing) {}

Precompressor::~Precompressor() {}

PrecompressorBlock*
Precompressor::readBlock(size_t blockSize, InStream* in) const {
  PrecompressorBlock *result = new PrecompressorBlock(blockSize, in);
  if(result->originalSize() > 0) precompress(*result);
  return result;
}

#define PREPROCESS(Type, verb, src) \
  Type r(grammar, (verb));\
  r.analyseData((src), length);\
  r.finishAnalysation();\
  r.decideReplacements();\
  length = r.writeReplacedVersion((src), length)

#define PREPROCESSS(Type, verb, src, dst) \
  Type r(grammar, (verb));\
  r.analyseData((src), length);\
  r.finishAnalysation();\
  r.decideReplacements();\
  length = r.writeReplacedVersion((src), length, (dst))

void Precompressor::precompress(PrecompressorBlock& block) const {
  PROFILE("Precompressor::precompress");
  byte *src = block.begin(), *dst;
  std::vector<byte> tmp;

  Grammar& grammar = block.grammar();
  size_t length = block.size();

  size_t roundsWithTempArray = 0;

  if(m_preprocessingOptions.size() > 0) {
    size_t oldLength = length;
    bool useTempArray = false;

    for(size_t i = 0; i < m_preprocessingOptions.size(); ++i) {
      char c = m_preprocessingOptions[i];
      double ratio = (double)length/block.originalSize();
      if(ratio < 0.666 && !useTempArray) {
        tmp.resize(length + 1);
        dst = &tmp[0];
        useTempArray = true;
      }

      if(c == 'p') {
        if(useTempArray) {
          PREPROCESSS(PairReplacer, verbosity > 1, src, dst);
        } else {
          PREPROCESS(PairReplacer, verbosity > 1, src);
        }
      }

      if(length == oldLength) {
        if(verbosity > 1) {
          std::clog << "Aborted precompression because couldn't achieve "
                    << "improvement anymore." << std::endl;
        }
        break;
      }

      if(useTempArray) {
        ++roundsWithTempArray;
        std::swap(src, dst);
      }

    }
  }

  if(roundsWithTempArray % 2 != 0) {
    std::copy(src, src + length, dst);
    src = dst;
  }

  block.setSize(length);
  if(verbosity > 1) {
    std::clog << "Size of precompressed block is "
              << ((double)length/block.originalSize())
              << " of original." << std::endl;
  }

}

} //namespace bwtc
