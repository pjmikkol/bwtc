/**
 * @file PrecompressorBlock.cpp
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
 * Implementation of Precompressor-block.
 */

#include "globaldefs.hpp"
#include "PrecompressorBlock.hpp"
#include "Streams.hpp"
#include "Utils.hpp"

#include <cassert>
#include <cstdlib>

namespace bwtc {

PrecompressorBlock::PrecompressorBlock(size_t size)
    : m_data((byte*)malloc(sizeof(byte)*(size+1))), m_used(0),
      m_originalSize(size), m_reserved(size+1) {}

PrecompressorBlock::PrecompressorBlock(size_t maxSize, InStream* in)
    : m_data((byte*)malloc(sizeof(byte)*(maxSize+1))), m_used(0),
      m_originalSize(0), m_reserved(maxSize+1)
{
  uint64 read = in->readBlock(m_data, maxSize);
  m_data = (byte*)realloc(m_data, sizeof(byte)*(read+1));
  m_reserved = read + 1;
  m_originalSize = m_used = read;
}

PrecompressorBlock::~PrecompressorBlock() {
  free(m_data);
}

void PrecompressorBlock::setSize(size_t size) {
  assert(size > 0);
  if(m_used != size) {
    m_used = size;
    m_reserved = size+1;
    m_data = (byte*)realloc(m_data, sizeof(byte)*(m_reserved));
  }
}

size_t PrecompressorBlock::
writeBlockHeader(OutStream* out, MetaData metadata) const {
  size_t bytes = 0;
  size_t integersToPack[2];
  size_t elems;
  if(metadata == PRECOMP_ORIGINAL_SIZE) {
    integersToPack[0] = originalSize();
    integersToPack[1] = slices();
    elems = 2;
  } else if(metadata == PRECOMP_COMPRESSED_SIZE) {
    integersToPack[0] = size();
    elems = 1;
  }

  for(size_t j = 0; j < elems; ++j) {
    int bytesNeeded;
    uint64 packedInteger = utils::packInteger(integersToPack[j], &bytesNeeded);
    for(int i = 0; i < bytesNeeded; ++i) {
      out->writeByte(packedInteger & 0xff);
      packedInteger >>= 8;
    }
    bytes += bytesNeeded;
  }
  
  bytes += m_grammar.writeGrammar(out);
  return bytes;
}

size_t PrecompressorBlock::writeEmptyHeader(OutStream* out) {
  out->writeByte(0);
  return 1;
}

PrecompressorBlock* PrecompressorBlock::
readBlockHeader(InStream* in, MetaData metadata) {
  if(metadata == PRECOMP_ORIGINAL_SIZE) {
    size_t bytes;
    size_t originalSize = utils::readPackedInteger(*in, bytes);
    PrecompressorBlock* pb = new PrecompressorBlock(originalSize);
    if(originalSize != 0) {
      size_t blocks = utils::readPackedInteger(*in, bytes);
      pb->m_bwtBlocks.resize(blocks);
      pb->grammar().readGrammar(in);
    }
    return pb;
  } else if(metadata == PRECOMP_COMPRESSED_SIZE) {
    size_t bytes;
    size_t compressedSize = utils::readPackedInteger(*in, bytes);
    PrecompressorBlock* pb = new PrecompressorBlock(compressedSize);
    if(compressedSize != 0) {
      pb->grammar().readGrammar(in);
      uint64 read = in->readBlock(pb->begin(), compressedSize);
      assert(read == compressedSize);
      pb->setSize(read);
    }
    return pb;
  }
}

void PrecompressorBlock::sliceIntoBlocks(size_t blockSize) {
  //Have at least one additional byte for the sentinel of BWT
  assert(m_used < m_reserved);
  assert(blockSize < (0x80000000 - 2));
  m_bwtBlocks.clear();
  size_t begin = 0;
  while(begin < m_used) {
    uint32 bSize = std::min(blockSize, m_used - begin);
    m_bwtBlocks.push_back(BWTBlock(&m_data[begin], bSize, false));
    begin += bSize;
  }
}


BWTBlock& PrecompressorBlock::getSlice(int i) {
  assert(i >= 0);
  assert((size_t)i < m_bwtBlocks.size());
  return m_bwtBlocks[i];
}

} //namespace bwtc
