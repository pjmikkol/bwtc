/**
 * @file TestPreprocessor.cpp
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
 * Implementation for TestPreprocessor which enables prototyping of the
 * preprocessing algorithms.
 */

#include <cassert>
#include <vector>

#include "../MainBlock.hpp"
#include "../BlockManager.hpp"
#include "../globaldefs.hpp"
#include "../Streams.hpp"
#include "Preprocessor.hpp"
#include "TestPreprocessor.hpp"
#include "PairReplacer.hpp"
#include "RunReplacer.hpp"
#include "PairAndRunReplacer.hpp"
#include "SequenceReplacer.hpp"

namespace bwtc {

TestPreprocessor::TestPreprocessor(uint64 block_size) :
    Preprocessor(block_size), m_currentBlock(0) {}

TestPreprocessor::~TestPreprocessor() {
  if(m_currentBlock) delete m_currentBlock;
}

/* Testing the case where preprocessing algorithms are interleaved */
size_t TestPreprocessor::pppr() {
  PairReplacer pr1(true,true);
  PairReplacer pr2(true,true);
  PairReplacer pr3(true,true);
  RunReplacer rr1(true, true);

  size_t origSize = m_currentBlock->m_filled;
  byte *src = &(*m_currentBlock->m_block)[0];
  pr1.analyseData(src, origSize);
  pr1.finishAnalysation();
  pr1.decideReplacements();

  byte *tmp = new byte[origSize + 3];
  size_t hSize = pr1.writeAndAnalyseHeader(tmp, pr2);
  size_t compr = pr1.writeAndAnalyseReplacedVersion(src, origSize, tmp+hSize, pr2);
  size_t compressedSize = compr + hSize;

  pr2.finishAnalysation();
  pr2.decideReplacements();
  hSize = pr2.writeAndAnalyseHeader(src, pr3);
  compr = pr2.writeAndAnalyseReplacedVersion(tmp, compressedSize, src+hSize, pr3);
  compressedSize = compr + hSize;

  pr3.finishAnalysation();
  pr3.decideReplacements();
  hSize = pr3.writeAndAnalyseHeader(tmp, rr1);
  compr = pr3.writeAndAnalyseReplacedVersion(src, compressedSize, tmp+hSize, rr1);
  compressedSize = compr + hSize;

  rr1.finishAnalysation();
  rr1.decideReplacements();
  hSize = rr1.writeHeader(src);
  compr = rr1.writeReplacedVersion(tmp, compressedSize, src+hSize);
  compressedSize = compr + hSize;
  
  //std::copy(tmp, tmp + compressedSize, src);
  m_currentBlock->m_filled = compressedSize;
  delete [] tmp;
  return origSize - compressedSize;

}

uint64 TestPreprocessor::compressPairs() {
  PairReplacer pr(true,true);

  size_t origSize = m_currentBlock->m_filled;
  byte *src = &(*m_currentBlock->m_block)[0];
  pr.analyseData(src, origSize);
  pr.finishAnalysation();
  pr.decideReplacements();

  byte *tmp = new byte[origSize + 3];
  size_t hSize = pr.writeHeader(tmp);
  size_t compr = pr.writeReplacedVersion(src, origSize, tmp+hSize);
  size_t compressedSize = compr + hSize;
  std::copy(tmp, tmp + compressedSize, src);
  m_currentBlock->m_filled = compressedSize;
  delete [] tmp;
  return origSize - compressedSize;
  /*
  uint64 filled = compressCommonPairs(&(*m_currentBlock->m_block)[0],
                                      m_currentBlock->m_filled);
  uint64 result = m_currentBlock->m_filled - filled;
  m_currentBlock->m_filled = filled;
  return result;*/
}

uint64 TestPreprocessor::compressRuns() {
  RunReplacer rr(true, true);

  size_t origSize = m_currentBlock->m_filled;
  byte *src = &(*m_currentBlock->m_block)[0];
  rr.analyseData(src, origSize);
  rr.finishAnalysation();
  rr.decideReplacements();

  byte *tmp = new byte[origSize + 3];
  size_t hSize = rr.writeHeader(tmp);
  size_t compr = rr.writeReplacedVersion(src, origSize, tmp+hSize);
  size_t compressedSize = compr + hSize;
  std::copy(tmp, tmp + compressedSize, src);
  m_currentBlock->m_filled = compressedSize;
  delete [] tmp;
  return origSize - compressedSize;
  /*
  uint64 filled = compressLongRuns(&(*m_currentBlock->m_block)[0],
                                   m_currentBlock->m_filled);
  uint64 result = m_currentBlock->m_filled - filled;
  m_currentBlock->m_filled = filled;
  return result;*/
}

size_t TestPreprocessor::compressPairsAndRuns() {
  pairs_and_runs::PairAndRunReplacer pr(true, true);

  size_t origSize = m_currentBlock->m_filled;
  byte *src = &(*m_currentBlock->m_block)[0];
  pr.analyseData(src, origSize);
  pr.finishAnalysation();
  pr.decideReplacements();

  byte *tmp = new byte[origSize + 3];
  size_t hSize = pr.writeHeader(tmp);
  size_t compr = pr.writeReplacedVersion(src, origSize, tmp+hSize);
  size_t compressedSize = compr + hSize;
  std::copy(tmp, tmp + compressedSize, src);
  m_currentBlock->m_filled = compressedSize;
  delete [] tmp;
  return origSize - compressedSize;
}

size_t TestPreprocessor::compressSequences() {
  SequenceReplacer sr(true, true);

  size_t origSize = m_currentBlock->m_filled;
  byte *src = &(*m_currentBlock->m_block)[0];
  sr.analyseData(src, origSize);
  sr.finishAnalysation();
  sr.decideReplacements();

  byte *tmp = new byte[origSize + 2];
  size_t hSize = sr.writeHeader(tmp);
  size_t compr = sr.writeReplacedVersion(src, origSize, tmp+hSize);
  size_t compressedSize = compr + hSize;
  std::copy(tmp, tmp + compressedSize, src);
  m_currentBlock->m_filled = compressedSize;
  delete [] tmp;
  return origSize - compressedSize;
}

void TestPreprocessor::initializeTarget() {
  assert(m_blockManager);
  std::vector<byte>* target = m_blockManager->getFreeBuffer();
  target->resize(m_blockSize);
  std::vector<uint64>* stats = m_blockManager->getFreeStats();
  m_currentBlock = m_blockManager->makeBlock(target, stats, 0UL);
}

uint64 TestPreprocessor::fillBuffer() {
  // TODO: Check the types
  assert(m_source);
  assert(m_currentBlock);
  assert(m_blockSize <= m_currentBlock->m_block->size());
  if (m_currentBlock->m_filled == m_blockSize) return 0;
  std::streamsize read = m_source->readBlock(
      m_currentBlock->begin() + m_currentBlock->size(),
      static_cast<std::streamsize>(m_blockSize - m_currentBlock->size()));
  m_currentBlock->m_filled += read;
  return static_cast<uint64>(read);
}


} //namespace bwtc
