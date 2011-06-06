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

namespace bwtc {

TestPreprocessor::TestPreprocessor(uint64 block_size) :
    Preprocessor(block_size), curr_block_(0) {}

TestPreprocessor::~TestPreprocessor() {
  if(curr_block_) delete curr_block_;
}

uint64 TestPreprocessor::compressPairs() {
  uint64 filled = compressCommonPairs(&(*curr_block_->m_block)[0],
                                      curr_block_->m_filled);
  uint64 result = curr_block_->m_filled - filled;
  curr_block_->m_filled = filled;
  return result;
}

uint64 TestPreprocessor::compressRuns() {
  uint64 filled = compressLongRuns(&(*curr_block_->m_block)[0],
                                   curr_block_->m_filled);
  uint64 result = curr_block_->m_filled - filled;
  curr_block_->m_filled = filled;
  return result;
}

void TestPreprocessor::InitializeTarget() {
  assert(block_manager_);
  std::vector<byte>* target = block_manager_->GetFreeBuffer();
  target->resize(block_size_);
  std::vector<uint64>* stats = block_manager_->GetFreeStats();
  curr_block_ = block_manager_->MakeBlock(target, stats, 0UL);
}

uint64 TestPreprocessor::FillBuffer() {
  // TODO: Check the types
  assert(source_);
  assert(curr_block_);
  assert(block_size_ <= curr_block_->m_block->size());
  if (curr_block_->m_filled == block_size_) return 0;
  std::streamsize read = source_->ReadBlock(
      curr_block_->begin() + curr_block_->size(),
      static_cast<std::streamsize>(block_size_ - curr_block_->size()));
  curr_block_->m_filled += read;
  return static_cast<uint64>(read);
}


} //namespace bwtc
