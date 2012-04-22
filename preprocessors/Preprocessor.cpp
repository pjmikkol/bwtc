/**
 * @file Preprocessor.cpp
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
 * Implementation of 2 preprocessing algorithms. One for replacing the
 * most common pairs and another for replacing long runs of the same byte.
 */

#include "../MainBlock.hpp"
#include "../BlockManager.hpp"
#include "../globaldefs.hpp"
#include "../Streams.hpp"
#include "../Utils.hpp"
#include "Preprocessor.hpp"
#include "FrequencyTable.hpp"
#include "PairReplacer.hpp"
#include "SequenceReplacer.hpp"
#include "../Profiling.hpp"

#include <cassert>
#include <boost/static_assert.hpp>

#include <algorithm>
#include <iostream> /* for std::streamsize*/
#include <map>
#include <string>
#include <utility> /* for pair */
#include <vector>

namespace bwtc {


Preprocessor::Preprocessor(uint64 block_size, const std::string& prepr)
    : m_source(0), m_blockSize(block_size), m_blockManager(0),
      m_preprocessingOptions(prepr) {}

Preprocessor::Preprocessor(uint64 block_size) :
    m_source(0), m_blockSize(block_size), m_blockManager(0)
{}

Preprocessor::~Preprocessor() {
  delete m_source;
}

void Preprocessor::buildStats(std::vector<byte>* data,
                              std::vector<uint64>* stats, uint64 data_size) {
  std::fill(stats->begin(), stats->end(), 0);
  //TODO: at the moment only contexts of length 1 are supported
  for( uint64 i = 0; i < data_size; ++i)
    (*stats)[(*data)[i]]++; 
}

void Preprocessor::connect(const std::string& source_name) {
  m_source = new InStream(source_name);
}

void Preprocessor::addBlockManager(BlockManager* manager) {
  m_blockManager = manager;
}

/* We append sentinel to the block here */
MainBlock* Preprocessor::readBlock() {
  PROFILE("Preprocessor::readBlock");
  assert(m_source);
  assert(m_blockManager);
  std::vector<byte>* to = m_blockManager->getFreeBuffer();
  std::vector<uint64>* stats = m_blockManager->getFreeStats();
  /* TODO:
   * streamsize type has as many bits as long. Since the preprocessor gets
   * blocksize as an uint64 we may end up in problems if user is on 32-bit
   * platform and wants to use too large blocks (amount of bits needed to
   * represent block size is more than 32)?? */

  /* The unused byte is left the block so that SA-IS-transformer can be used.
  *  Also prepare for the worst case with preprocessing algorithms. */
  std::streamsize read = m_source->readBlock(
      &(*to)[0], static_cast<std::streamsize>(
          m_blockSize - 1 - m_preprocessingOptions.size()*5));
  if (!read) return 0;
  
  size_t preprocessedSize = preprocess(*to, read);
  
  return m_blockManager->makeBlock(to, stats, preprocessedSize);
}

#define PREPROCESS(Type, verb, src, dst) \
  Type r(grammar, (verb));      \
  r.analyseData((src), length); \
  r.finishAnalysation(); \
  r.decideReplacements(); \
  size_t comprSize = r.writeReplacedVersion((src), length, (dst)); \
  length = comprSize

size_t Preprocessor::preprocess(std::vector<byte>& original, size_t length) {
  PROFILE("Preprocessor::preprocess");
  byte *dst, *src = &original[0];
  std::vector<byte> tmp;
  /**Stores the replacement rules. */
  Grammar grammar;
  if(m_preprocessingOptions.size() > 0) {
    //TODO: correct limits
    tmp.resize(length + m_preprocessingOptions.size()*5);
    dst = &tmp[0];
    for(size_t i = 0; i < m_preprocessingOptions.size(); ++i) {
      char c = m_preprocessingOptions[i];
      if(c == 'p') {
        PREPROCESS(PairReplacer, verbosity > 1, src, dst);
      } /*else if(c == 's') {
        PREPROCESS(SequenceReplacer, verbosity > 1, src, dst);
        }*/
      std::swap(src, dst);
    }
  }
  if(m_preprocessingOptions.size() & 1) {
    std::copy(src, src + length, dst);
    src = dst;
  }
  //grammar.printRules();
  
  uint32 gSize = grammar.writeGrammar(src+length);
  length += gSize;
  if(verbosity > 0) {
    std::clog << "Size of preprocessed block is " << length <<std::endl;
  }
  return length;
}

#undef PREPROCESS


} //namespace bwtc
