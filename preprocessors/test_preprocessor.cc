/**************************************************************************
 *  Copyright 2010, Pekka Mikkola, pjmikkol (at) cs.helsinki.fi           *
 *                                                                        *
 *  This file is part of bwtc.                                            *
 *                                                                        *
 *  bwtc is free software: you can redistribute it and/or modify          *
 *  it under the terms of the GNU General Public License as published by  *
 *  the Free Software Foundation, either version 3 of the License, or     *
 *  (at your option) any later version.                                   *
 *                                                                        *
 *  bwtc is distributed in the hope that it will be useful,               *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *  GNU General Public License for more details.                          *
 *                                                                        *
 *  You should have received a copy of the GNU General Public License     *
 *  along with bwtc.  If not, see <http://www.gnu.org/licenses/>.         *
 **************************************************************************/

#include <cassert>

#include <vector>

#include "../block.h"
#include "../block_manager.h"
#include "../globaldefs.h"
#include "../stream.h"
#include "preprocessor.h"
#include "test_preprocessor.h"

namespace bwtc {

TestPreProcessor::TestPreProcessor(uint64 block_size) :
    PreProcessor(block_size), curr_block_(NULL) {}

TestPreProcessor::~TestPreProcessor() {
  if(curr_block_) delete curr_block_;
}

uint64 TestPreProcessor::CompressPairs() {
  uint64 filled = CompressCommonPairs(&(*curr_block_->block_)[0],
                                      curr_block_->filled_);
  uint64 result = curr_block_->filled_ - filled;
  curr_block_->filled_ = filled;
  return result;
}

uint64 TestPreProcessor::CompressRuns() {
  uint64 filled = CompressLongRuns(&(*curr_block_->block_)[0],
                                   curr_block_->filled_);
  uint64 result = curr_block_->filled_ - filled;
  curr_block_->filled_ = filled;
  return result;
}

void TestPreProcessor::InitializeTarget() {
  assert(block_manager_);
  std::vector<byte>* target = block_manager_->GetFreeBuffer();
  target->resize(block_size_);
  std::vector<uint64>* stats = block_manager_->GetFreeStats();
  curr_block_ = block_manager_->MakeBlock(target, stats, 0UL);
}

uint64 TestPreProcessor::FillBuffer() {
  // TODO: Check the types
  assert(source_);
  assert(curr_block_);
  assert(block_size_ <= curr_block_->block_->size());
  if (curr_block_->filled_ == block_size_) return 0;
  std::streamsize read = source_->ReadBlock(
      curr_block_->begin() + curr_block_->Size(),
      static_cast<std::streamsize>(block_size_ - curr_block_->Size()));
  curr_block_->filled_ += read;
  return static_cast<uint64>(read);
}


} //namespace bwtc
