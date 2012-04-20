/**
 * @file HuffmanCoders.cpp
 * @author Dominik Kempa <dominik.kempa@cs.helsinki.fi>
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
 * Implementations of Huffman encoder and decoder.
 */

#include <cassert>
#include <cstdio>
#include <ctime>
#include <iterator>
#include <iostream> // For std::streampos
#include <numeric> // for std::accumulate
#include <algorithm> // for sort, reverse, fill
#include <string>
#include <vector>

#include "MainBlock.hpp"
#include "HuffmanCoders.hpp"
#include "globaldefs.hpp"
#include "Utils.hpp"

namespace bwtc {

HuffmanEncoder::HuffmanEncoder(const std::string& destination, char prob_model)
    : m_out(new RawOutStream(destination)),
      m_headerPosition(0), m_compressedBlockLength(0) {
  (void) prob_model;
}

void HuffmanEncoder::writeGlobalHeader(const std::string& preproc, char encoding) {
  // bitpatterns in encoding of preprocessing:
  // 001 -- run, 010 -- pairAndrun, 011 -- pair,  100 -- sequence, 000 -- end 
  uint16 b = 0;
  size_t bits = 0;
  for(size_t i = 0; i < preproc.size(); ++i) {
    if(preproc[i] == 'p') b = (b << 3) | 3;
    else if (preproc[i] == 'r') b = (b << 3) | 1;
    else if (preproc[i] == 'c') b = (b << 3) | 2;
    else if (preproc[i] == 's') b = (b << 3) | 4;
    bits += 3;

    if(bits >= 8) {
      m_out->writeByte(b >> (bits - 8));
      b &= ((1 << (bits - 8)) - 1);
      bits -= 8;
    }
  }
  if(preproc.size() == 0) {
    m_out->writeByte(static_cast<byte>(0));
  } else if(bits == 8) {
    m_out->writeByte(static_cast<byte>(b & 0xff));
    m_out->writeByte(static_cast<byte>(0));
  } else if (bits <= 5){
    m_out->writeByte(static_cast<byte>((b << (8 - bits)) & 0xff));
  } else {
    m_out->writeByte(static_cast<byte>((b << (8 - bits)) & 0xff));
    m_out->writeByte(static_cast<byte>(0));
  }

  m_out->writeByte(static_cast<byte>(encoding));
}

std::string HuffmanDecoder::readGlobalHeader() {
  std::string preproc;
  size_t bitsLeft = 0;
  size_t bitsInCode = 0;
  size_t code = 0;
  byte b = 0;

  while(true) {
    if(bitsInCode == 3) {
      if(code == 0) break;
      else if (code == 2) preproc += 'c';
      else if (code == 1) preproc += 'r';
      else if (code == 3) preproc += 'p';
      else if (code == 4) preproc += 's';
      code = 0;
      bitsInCode = 0;
    }
    if(bitsLeft == 0) {
      b = m_in->readByte();
      bitsLeft = 8;
    }
    code = (code << 1) | ((b >> (bitsLeft - 1)) & 1);
    --bitsLeft;
    ++bitsInCode;
  }
  char probmodel = static_cast<char>(m_in->readByte());
  (void) probmodel;
  return preproc;
}

HuffmanEncoder::~HuffmanEncoder() {
  delete m_out;
}

int HuffmanEncoder::writeTrailer(const std::vector<uint32>& LFpowers) {
  int bytes = 1;
  byte s = (byte)(LFpowers.size()-1);
  m_out->writeByte(s);
  int bitsLeft = 8;
  for(size_t i = 0; i < LFpowers.size(); ++i) {
    for(int j = 30; j >= 0; --j) {
      s = (s << 1) | ((LFpowers[i] >> j) & 0x1);
      --bitsLeft;
      if(bitsLeft == 0) {
        m_out->writeByte(s);
        bitsLeft = 8;
        ++bytes;
      }
    }
  }
  if(bitsLeft < 8) {
    m_out->writeByte(s << bitsLeft);
    ++bytes;
  }
  return bytes;
}

void HuffmanEncoder::serializeShape(uint32 *clen, std::vector<bool> &vec) {
  size_t maxLen = 0;
  byte b = 0;
  std::vector<byte> symbols;
  for(size_t i = 0; i < 256; ++i, ++b) {
    if(clen[i] > 0) {
      symbols.push_back(b);
      if(clen[i] > maxLen) maxLen = clen[i];
    }
  }
  // Largest symbol in alphabet
  utils::pushBits(vec, symbols.back(), 8);
  // Number of distinct symbols
  utils::pushBits(vec, symbols.size(), 8);

  int bytesInLongestCode;
  size_t packedInt = utils::packInteger(maxLen, &bytesInLongestCode);
  utils::pushBits(vec, packedInt, bytesInLongestCode*8);

  utils::binaryInterpolativeCode(symbols, symbols.back(), vec);

  for(size_t i = 0; i < symbols.size(); ++i)
    utils::unaryCode(vec, maxLen - clen[symbols[i]] + 1);
}

size_t HuffmanDecoder::deserializeShape(RawInStream &input, uint32 *clen) {
  size_t maxSym = input.readByte();
  size_t symbols = input.readByte();
  if(symbols == 0) symbols = 256;

  size_t bitsRead = 16;
  size_t maxLen = 0;
  size_t read = 0xff;
  size_t j = 0;
  while(read & 0x80) {
    read = input.readByte();
    maxLen |= ((read & 0x7f) << j);
    j += 7;
    bitsRead += 8;
  }

  std::vector<byte> alphabet;
  bitsRead += utils::binaryInterpolativeDecode(alphabet, input,
                                               maxSym, symbols);

  for(size_t i = 0; i < symbols; ++i) {
    size_t n = utils::unaryDecode(input);
    bitsRead += n;
    size_t len = maxLen - n + 1;
    clen[alphabet[i]] = len;
  }

  input.flushBuffer();
  return (bitsRead + 7) / 8;
}

/******************************************************************************
 *            Encoding and decoding single MainBlock                          *
 *----------------------------------------------------------------------------*
 * Following functions handle encoding and decoding of the main blocks.       *
 * Also block header-format is specified here.                                *
 *                                                                            *
 * Format of the block is following:                                          *
 *  - Block header (no fixed length)                              (1)         *
 *  - Compressed block (no fixed length)                          (2)         *
 *  - Block trailer (coded same way as the context block lengths) (3)         *
 *                                                                            *
 *----------------------------------------------------------------------------*
 *                                                                            *
 * Block header (1) format is following:                                      *
 * a) Length of the (header + compressed block + trailer) in bytes (6 bytes). *
 *    Note that the the length field itself isn't included in total length    *
 * b) List of context block lengths. Lengths are compressed with              *
 *    utils::PackInteger-function.                                            *
 * c) Ending symbol of the block header (2 bytes = 0x8000) which is           *
 *    invalid code for packed integer                                         *
 ******************************************************************************/
void HuffmanEncoder::encodeData(std::vector<byte>* block,
    std::vector<uint64>* stats, uint64 block_size) {
  (void) block_size;
  size_t beg = 0;
  
  // For storing runs data.
  byte *runseq = new byte[block_size];
  uint32 *runlen = new uint32[block_size];
  if (!runseq || !runlen) {
    fprintf(stderr,"Allocation error.\n");
    exit(1);
  }

  byte *block_ptr = &(*block)[0];
  for(size_t i = 0; i < stats->size(); ++i) {
    size_t current_cblock_size = (*stats)[i];
    if(current_cblock_size == 0) continue;

    // Compute lengths of Huffman codes.
    uint32 clen[256];
    std::fill(clen, clen + 256, 0);
    uint64 freqs[256];
    std::fill(freqs, freqs + 256, 0);
    uint64 nRuns = utils::calculateRunFrequenciesAndStoreRuns(freqs,
      runseq, runlen, block_ptr + beg, current_cblock_size);
    std::vector<std::pair<uint64, byte> > codeLengths;
    utils::calculateHuffmanLengths(codeLengths, freqs);
    int32 nCodes = codeLengths.size();
    for (int32 k = 0; k < nCodes; ++k)
      clen[codeLengths[k].second] = codeLengths[k].first;

    // Store the number of runs.
    int bytes = 0;
    uint64 packed_nRuns = utils::packInteger(nRuns, &bytes);
    m_compressedBlockLength += bytes;
    writePackedInteger(packed_nRuns);

    // Store Huffman code lengths.
    std::vector<bool> shape;
    serializeShape(clen, shape);
    for(size_t k = 0; k < shape.size();) {
      byte b = 0; size_t j = 0;
      for(; j < 8 && k < shape.size(); ++k, ++j) {
        b <<= 1;
        b |= (shape[k]) ? 1 : 0;
      }
      if (j < 8) b <<= (8 - j);
      m_out->writeByte(b);
      ++m_compressedBlockLength;
    }

    // Compute Huffman codes.
    uint32 code[256];
    utils::computeHuffmanCodes(clen, code);

    // Encode the data using Huffman code.
    // Assumption: max_code_len <= 47 (roughly).
    uint64 buffer = 0;
    int32 bitsInBuffer = 0;
    for (uint64 k = 0; k < nRuns; ++k) {
      byte c = runseq[k];
      while (bitsInBuffer + clen[c] > 64) {
        bitsInBuffer -= 8;
        m_out->writeByte((buffer >> bitsInBuffer) & 0xff);
        ++m_compressedBlockLength;
      }
      buffer <<= clen[c];
      buffer |= code[c];
      bitsInBuffer += clen[c];
    }

    // Flush the remaining bytes.    
    while (bitsInBuffer >= 8) {
      bitsInBuffer -= 8;
      m_out->writeByte((buffer >> bitsInBuffer) & 0xff);
      ++m_compressedBlockLength;
    }

    // Flush the remaining bits.
    if (bitsInBuffer > 0) {
      buffer <<= (8 - bitsInBuffer);
      m_out->writeByte(buffer & 0xff);
      ++m_compressedBlockLength;
    }

    // Store the lengths of runs.
    buffer = 0;
    bitsInBuffer = 0;
    for (uint64 k = 0; k < nRuns; ++k) {
      int gammaCodeLen = utils::logFloor(runlen[k]) * 2 + 1;
      while (bitsInBuffer + gammaCodeLen > 64) {
        bitsInBuffer -= 8;
        m_out->writeByte((buffer >> bitsInBuffer) & 0xff);
        ++m_compressedBlockLength;
      }
      buffer <<= gammaCodeLen;
      buffer |= runlen[k];
      bitsInBuffer += gammaCodeLen;
    }
    while (bitsInBuffer >= 8) {
      bitsInBuffer -= 8;
      m_out->writeByte((buffer >> bitsInBuffer) & 0xff);
      ++m_compressedBlockLength;
    }
    if (bitsInBuffer > 0) {
      buffer <<= (8 - bitsInBuffer);
      m_out->writeByte(buffer & 0xff);
      ++m_compressedBlockLength;
    }

    beg += current_cblock_size;
  }
  delete[] runseq;
  delete[] runlen;
}

void HuffmanEncoder::finishBlock(const std::vector<uint32>& LFpowers) {
  m_compressedBlockLength +=  writeTrailer(LFpowers);
  m_out->write48bits(m_compressedBlockLength, m_headerPosition);
}

/*********************************************************************
 * The format of header for single main block is the following:      *
 * - 48 bits for the length of the compressed main block, doesn't    *
 *   include 6 bytes used for this                                   *
 * - byte representing the number of separately encoded sections.    *
 *   zero represents 256                                             *
 * - lengths of the sections which are encoded with same wavelet tree*
 *********************************************************************/
void HuffmanEncoder::writeBlockHeader(std::vector<uint64>* stats) {
  uint64 headerLength = 0;
  m_headerPosition = m_out->getPos();
  for (unsigned i = 0; i < 6; ++i) m_out->writeByte(0x00); //fill 48 bits

  /* Deduce sections for separate encoding. At the moment uses not-so-well
   * thought heuristic. */
  std::vector<uint64> temp; std::vector<uint64>& s = *stats;
  size_t sum = 0;
  for(size_t i = 0; i < s.size(); ++i) {
    sum += s[i];
    if(sum >= 10000) {
      temp.push_back(sum);
      sum = 0;
    }
  }
  if (sum != 0) {
    if(temp.size() > 0) temp.back() += sum;
    else temp.push_back(sum);
  }
  s.resize(temp.size());
  std::copy(temp.begin(), temp.end(), s.begin());
  byte len;
  if(temp.size() == 256) len = 0;
  else len = temp.size();
  m_out->writeByte(len);
  headerLength += 1;

  assert(s.size() == temp.size());
  assert(temp.size() <= 256);

  for (size_t i = 0; i < stats->size(); ++i) {
    int bytes;
    uint64 packed_cblock_size = utils::packInteger((*stats)[i], &bytes);
    headerLength += bytes;
    writePackedInteger(packed_cblock_size);
  }
  m_compressedBlockLength = headerLength;
}

/* Integer is written in reversal fashion so that it can be read easier.*/
void HuffmanEncoder::writePackedInteger(uint64 packed_integer) {
  do {
    byte to_written = static_cast<byte>(packed_integer & 0xFF);
    packed_integer >>= 8;
    m_out->writeByte(to_written);
  } while (packed_integer);
}

uint64 HuffmanDecoder::readBlockHeader(std::vector<uint64>* stats) {
  uint64 compressed_length = m_in->read48bits();
  byte sections = m_in->readByte();
  size_t sects = (sections == 0) ? 256 : sections;
  for(size_t i = 0; i < sects; ++i) {
    uint64 value = readPackedInteger();
    stats->push_back(utils::unpackInteger(value));
  }
  return compressed_length;
}

std::vector<byte>* HuffmanDecoder::decodeBlock(std::vector<uint32>& LFpowers) {
  if(m_in->compressedDataEnding()) return 0;

  std::vector<uint64> context_lengths;
  uint64 compr_len = readBlockHeader(&context_lengths);

  if (verbosity > 2) {
    std::clog << "Size of compressed block = " << compr_len << "\n";
  }

  uint64 block_size = std::accumulate(
      context_lengths.begin(), context_lengths.end(), static_cast<uint64>(0));
  
  byte *runseq = new byte[block_size];
  uint32 *runlen = new uint32[block_size];

  std::vector<byte>* data = new std::vector<byte>();
  data->resize(block_size);
  byte *data_ptr = &(*data)[0];
  uint64 beg = 0;

  for(size_t i = 0; i < context_lengths.size(); ++i) {
    if(context_lengths[i] == 0) continue;

    // Get the number of runs withing the current context block.
    uint64 packed_nRuns = readPackedInteger();
    uint64 nRuns = utils::unpackInteger(packed_nRuns);

    // Get Huffman code lengths.
    uint32 clen[256];
    std::fill(clen, clen + 256, 0);
    deserializeShape(*m_in, clen);

    // Compute Huffman codes.
    uint32 code[256];
    utils::computeHuffmanCodes(clen, code);

    uint64 buffer = 0;
    uint32 bitsInBuffer = 0;
    for (uint64 k = 0; k < nRuns; ++k) {
      buffer = 0;
      bitsInBuffer = 0;
      byte next_byte;
      bool found = false;
      while (!found) {
        int32 bit = m_in->readBit();
        buffer = (buffer << 1) | bit;
        ++bitsInBuffer;
        for (int32 j = 0; j < 256; ++j) {
          if (clen[j] == bitsInBuffer && (uint64)code[j] == buffer) {
            found = true;
            next_byte = static_cast<byte>(j);
            break;
          }
        }
      }
      runseq[k] = next_byte;
    }
    m_in->flushBuffer();

    // Now read gamma codes that store lenghts of runs.
    for (uint64 k = 0; k < nRuns; ++k) {
      int zeros = 0;
      while (!m_in->readBit())
        ++zeros;
      uint64 value = 0;        
      for (int32 t = 0; t < zeros; ++t) {
        int32 bit = m_in->readBit();
        value = (value << 1) | bit;
      }
      value |= (1 << zeros);
      runlen[k] = value;
    }
    m_in->flushBuffer();
    
    // Fill the block with runs data.
    for (uint64 k = 0; k < nRuns; ++k) {
      for (uint32 t = 0; t < runlen[k]; ++t)
        data_ptr[beg++] = runseq[k];
    }
  }

  uint32 LFpows = m_in->readByte()+1;
  LFpowers.resize(LFpows);
  for(uint32 i = 0; i < LFpows; ++i) {
    uint32 pos = 0;
    for(uint32 j = 0; j < 31; ++j)
      pos = (pos << 1) | (m_in->readBit() ? 1 : 0);
    LFpowers[i] = pos;
  }
  m_in->flushBuffer();
  delete[] runseq;
  delete[] runlen;
  return data;
}

uint64 HuffmanDecoder::readPackedInteger() {
  static const uint64 kEndSymbol = static_cast<uint64>(1) << 63;
  static const uint64 kEndMask = static_cast<uint64>(1) << 7;

  uint64 packed_integer = 0;
  bool bits_left = true;
  int i;
  for(i = 0; bits_left; ++i) {
    uint64 read = static_cast<uint64>(m_in->readByte());
    bits_left = (read & kEndMask) != 0;
    packed_integer |= (read << i*8);
  }
  if (packed_integer == 0x80) return kEndSymbol;
  return packed_integer;
}
/*********** Encoding and decoding single MainBlock-section ends ********/

HuffmanDecoder::HuffmanDecoder(const std::string& source) :
    m_in(new RawInStream(source)) {}

HuffmanDecoder::~HuffmanDecoder() {
  delete m_in;
}


} // namespace bwtc

