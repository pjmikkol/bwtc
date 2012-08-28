/**
 * @file HuffmanCoders.cpp
 * @author Dominik Kempa <dominik.kempa@cs.helsinki.fi>
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
#include <map> // for entropy profiling

#include "HuffmanCoders.hpp"
#include "globaldefs.hpp"
#include "Utils.hpp"
#include "Profiling.hpp"

namespace bwtc {

HuffmanEncoder::HuffmanEncoder()
    : m_headerPosition(0), m_compressedBlockLength(0) {}

HuffmanEncoder::~HuffmanEncoder() {}

size_t HuffmanEncoder::
transformAndEncode(BWTBlock& block, BWTManager& bwtm, OutStream* out) {
  std::vector<uint32> characterFrequencies(256, 0);
  //TODO: Also gather information about the runs  during BWT
  bwtm.doTransform(block, &characterFrequencies[0]); 

  writeBlockHeader(block, characterFrequencies, out);
  encodeData(block.begin(), characterFrequencies, block.size(), out);
  finishBlock(out);
  return m_compressedBlockLength + 6;
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

size_t HuffmanDecoder::deserializeShape(InStream &input, uint32 *clen) {
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

void HuffmanEncoder::
encodeData(const byte* block, const std::vector<uint32>& stats,
           uint32 blockSize, OutStream* out) {
  PROFILE("HuffmanEncoder::encodeData");
  size_t beg = 0;
  
  // For storing runs data.
  byte *runseq = new byte[blockSize];
  uint32 *runlen = new uint32[blockSize];
  if (!runseq || !runlen) {
    fprintf(stderr,"Allocation error.\n");
    exit(1);
  }

  const byte *block_ptr = block;
  for(size_t i = 0; i < stats.size(); ++i) {
    size_t current_cblock_size = stats[i];
    if(current_cblock_size == 0) continue;

    // Compute lengths of Huffman codes.
    uint32 clen[256];
    std::fill(clen, clen + 256, 0);
    uint64 freqs[256];
    std::fill(freqs, freqs + 256, 0);
    uint64 nRuns = utils::calculateRunFrequenciesAndStoreRuns(freqs,
      runseq, runlen, block_ptr + beg, current_cblock_size);

#ifdef ENTROPY_PROFILER
    {
      std::map<uint32, uint32> runDistribution, charDistribution;
      for(size_t j = 0; j < nRuns; ++j) {
        ++runDistribution[runlen[j]];
        ++charDistribution[runseq[j]];
      }
      
      for(std::map<uint32, uint32>::const_iterator it = runDistribution.begin();
          it != runDistribution.end(); ++it)
        std::cout << it->first << ":" <<  it->second << std::endl;

      std::cout << "----" << std::endl;

      for(std::map<uint32, uint32>::const_iterator it = charDistribution.begin();
          it != charDistribution.end(); ++it)
        std::cout << it->first << ":" <<  it->second << std::endl;

      std::cout << "####" << std::endl;
    }
#endif

    std::vector<std::pair<uint64, uint32> > codeLengths;
    utils::calculateHuffmanLengths(codeLengths, freqs);
    int32 nCodes = codeLengths.size();
    for (int32 k = 0; k < nCodes; ++k)
      clen[codeLengths[k].second] = codeLengths[k].first;

    // Store the number of runs.
    int bytes = 0;
    uint64 packed_nRuns = utils::packInteger(nRuns, &bytes);
    m_compressedBlockLength += bytes;
    writePackedInteger(packed_nRuns, out);

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
      out->writeByte(b);
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
        out->writeByte((buffer >> bitsInBuffer) & 0xff);
        ++m_compressedBlockLength;
      }
      buffer <<= clen[c];
      buffer |= code[c];
      bitsInBuffer += clen[c];
    }

    // Flush the remaining bytes.    
    while (bitsInBuffer >= 8) {
      bitsInBuffer -= 8;
      out->writeByte((buffer >> bitsInBuffer) & 0xff);
      ++m_compressedBlockLength;
    }

    // Flush the remaining bits.
    if (bitsInBuffer > 0) {
      buffer <<= (8 - bitsInBuffer);
      out->writeByte(buffer & 0xff);
      ++m_compressedBlockLength;
    }

    // Store the lengths of runs.
    buffer = 0;
    bitsInBuffer = 0;
    for (uint64 k = 0; k < nRuns; ++k) {
      int gammaCodeLen = utils::logFloor(runlen[k]) * 2 + 1;
      while (bitsInBuffer + gammaCodeLen > 64) {
        bitsInBuffer -= 8;
        out->writeByte((buffer >> bitsInBuffer) & 0xff);
        ++m_compressedBlockLength;
      }
      buffer <<= gammaCodeLen;
      buffer |= runlen[k];
      bitsInBuffer += gammaCodeLen;
    }
    while (bitsInBuffer >= 8) {
      bitsInBuffer -= 8;
      out->writeByte((buffer >> bitsInBuffer) & 0xff);
      ++m_compressedBlockLength;
    }
    if (bitsInBuffer > 0) {
      buffer <<= (8 - bitsInBuffer);
      out->writeByte(buffer & 0xff);
      ++m_compressedBlockLength;
    }

    beg += current_cblock_size;
  }
  delete[] runseq;
  delete[] runlen;
}

void HuffmanEncoder::finishBlock(OutStream* out) {
  out->write48bits(m_compressedBlockLength, m_headerPosition);
}

/*********************************************************************
 * The format of header for single main block is the following:      *
 * - 48 bits for the length of the compressed main block, doesn't    *
 *   include 6 bytes used for this                                   *
 * - byte representing the number of separately encoded sections.    *
 *   zero represents 256                                             *
 * - lengths of the sections which are encoded with same wavelet tree*
 *********************************************************************/
void HuffmanEncoder::
writeBlockHeader(const BWTBlock& block, std::vector<uint32>& stats,
                 OutStream* out) {
  uint64 headerLength = 0;
  m_headerPosition = out->getPos();
  for (unsigned i = 0; i < 6; ++i) out->writeByte(0x00); //fill 48 bits

  headerLength += block.writeHeader(out);
  
  /* Deduce sections for separate encoding. At the moment uses not-so-well
   * thought heuristic. */
  std::vector<uint32> temp; std::vector<uint32>& s = stats;
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
  out->writeByte(len);
  headerLength += 1;

  assert(s.size() == temp.size());
  assert(temp.size() <= 256);

  for (size_t i = 0; i < stats.size(); ++i) {
    int bytes;
    uint64 packed_cblock_size = utils::packInteger(stats[i], &bytes);
    headerLength += bytes;
    writePackedInteger(packed_cblock_size, out);
  }
  m_compressedBlockLength = headerLength;
}

/* Integer is written in reversal fashion so that it can be read easier.*/
void HuffmanEncoder::writePackedInteger(uint64 packed_integer, OutStream* out) {
  do {
    byte to_written = static_cast<byte>(packed_integer & 0xFF);
    packed_integer >>= 8;
    out->writeByte(to_written);
  } while (packed_integer);
}

uint64 HuffmanDecoder::
readBlockHeader(BWTBlock& block, std::vector<uint64>* stats, InStream* in) {
  uint64 compressed_length = in->read48bits();
  block.readHeader(in);
  
  byte sections = in->readByte();
  size_t sects = (sections == 0) ? 256 : sections;
  for(size_t i = 0; i < sects; ++i) {
    uint64 value = readPackedInteger(in);
    stats->push_back(utils::unpackInteger(value));
  }
  return compressed_length;
}

void HuffmanDecoder::decodeBlock(BWTBlock& block, InStream* in) {
  PROFILE("HuffmanDecoder::decodeBlock");
  if(in->compressedDataEnding()) return;

  std::vector<uint64> context_lengths;
  uint64 compr_len = readBlockHeader(block, &context_lengths, in);

  if (verbosity > 2) {
    std::clog << "Size of compressed block = " << compr_len << "\n";
  }

  uint64 block_size = std::accumulate(
      context_lengths.begin(), context_lengths.end(), static_cast<uint64>(0));
  
  byte *runseq = new byte[block_size];
  uint32 *runlen = new uint32[block_size];

  byte *data_ptr = block.begin();

  for(size_t i = 0; i < context_lengths.size(); ++i) {
    if(context_lengths[i] == 0) continue;

    // Get the number of runs withing the current context block.
    uint64 packed_nRuns = readPackedInteger(in);
    uint64 nRuns = utils::unpackInteger(packed_nRuns);

    // Get Huffman code lengths.
    uint32 clen[256];
    std::fill(clen, clen + 256, 0);
    deserializeShape(*in, clen);

    // Compute Huffman codes.
    uint32 code[256];
    utils::computeHuffmanCodes(clen, code);

    // Comput lookup arrays used in Huffman decoding.
    //   1. For each byte k, leadingZeros[p][k] tells what is the number of
    //      leading zeros on p least significant bits of k (the other bits
    //      doesn't matter).
    byte leadingZeros[9][256];
    
    //   2. For each length t and byte k lookupWhich[t][k] tells which
    //      Huffman code consists of t zeros followed by some prefix of k.
    //      In case there is no such code, the value 256 is stored.
    uint16 lookupWhich[50][256];
    
    //   3. For each length t and byte k such that lookupWhich[t][k] < 256,
    //      lookupLength[t][k] tells the length of corresponding code prefix.
    byte lookupLength[50][256];
    
    // Compute lookup arrays.
    for (int32 p = 0; p <= 8; ++p) {
      for (int32 k = 0; k < 256; ++k) {
        leadingZeros[p][k] = p;
        for (int32 t = 0; t < p; ++t)
          if (k & (1 << t))
            leadingZeros[p][k] = p - t - 1;
      }
    }
    for (int32 k = 0; k < 50; ++k)
      std::fill(&lookupWhich[k][0], &lookupWhich[k][256], 256);
    for (int k = 0; k < 256; ++k) if (clen[k] && code[k]) {
      int32 max_non0_bit = 0;
      for (int32 j = 0; j < 16; ++j)
        if (code[k] & (1 << j))
          max_non0_bit = j;
      ++max_non0_bit;
      int32 lead0 = clen[k] - max_non0_bit;
      int32 rest = 8 - max_non0_bit;
      for (int32 any_set = 0; any_set < (1 << rest); ++any_set) {
        lookupWhich[lead0][(code[k] << rest) + any_set] = k;
        lookupLength[lead0][(code[k] << rest) + any_set] = max_non0_bit;
      }
    }
    
    // Store the symbol (and its length) that have code with zeros only.
    int32 len0 = 0;
    byte b = 0, code0 = 0;
    for (int32 k = 0; k < 256; ++k, ++b) {
      if (clen[k] > 0 && code[k] == 0) {
        len0 = clen[k];
        code0 = b;
      }
    }

    // Decode Huffman codes.
    byte buffer = 0;
    int32 bitsInBuffer = 0, zeroCount = 0;
    uint32 haveDecoded = 0;
    while (haveDecoded < nRuns) {
      byte b = in->readByte();
      if (zeroCount > 0 && !bitsInBuffer) {
        buffer = b;
        if (!buffer) {
          zeroCount += 8;
          
          // Extract codes consisting of zeros (code0s).
          while (zeroCount >= len0 && haveDecoded < nRuns) {
            runseq[haveDecoded++] = code0;
            zeroCount -= len0;
          }
        } else {
        
          // Extract code0s.
          zeroCount += leadingZeros[8][buffer];
          bitsInBuffer = 8 - leadingZeros[8][buffer];
          while (zeroCount >= len0 && haveDecoded < nRuns) {
            runseq[haveDecoded++] = code0;
            zeroCount -= len0;
          }
          
          // Extract code on the boundary of previous byte (buffer) and b.
          if (zeroCount + bitsInBuffer > 8 && haveDecoded < nRuns) {
            byte sbuffer = (buffer << (8 - bitsInBuffer));
            if (lookupWhich[zeroCount][sbuffer] < 256 &&
                bitsInBuffer >= lookupLength[zeroCount][sbuffer]) {
              runseq[haveDecoded++] = lookupWhich[zeroCount][sbuffer];
              bitsInBuffer -= lookupLength[zeroCount][sbuffer];
              zeroCount = leadingZeros[bitsInBuffer][buffer];
              bitsInBuffer -= zeroCount;
              buffer &= (1 << bitsInBuffer) - 1;
            }
            while (zeroCount >= len0 && haveDecoded < nRuns) {
              runseq[haveDecoded++] = code0;
              zeroCount -= len0;
            }
          }

          // Extract codes from b.
          if (zeroCount + bitsInBuffer <= 8 && haveDecoded < nRuns) {
            while (haveDecoded < nRuns) { bool decoded = false;
              byte sbuffer = (buffer << (8 - bitsInBuffer));
              if (lookupWhich[zeroCount][sbuffer] < 256
                  && bitsInBuffer >= lookupLength[zeroCount][sbuffer]) {
                runseq[haveDecoded++] = lookupWhich[zeroCount][sbuffer];
                bitsInBuffer -= lookupLength[zeroCount][sbuffer];
                zeroCount = leadingZeros[bitsInBuffer][buffer];
                bitsInBuffer -= zeroCount; decoded = true;
                buffer &= (1 << bitsInBuffer) - 1;
              }
              while (zeroCount >= len0 && haveDecoded < nRuns) {
                runseq[haveDecoded++] = code0;
                zeroCount -= len0; decoded = true;
              }
              if (!decoded) break;
            }
            continue;
          }

          // Tricky case, the code starts in byte before b (buffer) and ends in
          // byte after b (nextb), but since b is nonzero, we are guaranteed
          // that this code will end in the very next byte after b (nextb).
          if (zeroCount + bitsInBuffer > 8 && haveDecoded < nRuns) {
            
            // Read nextb.
            byte nextb = in->readByte();
            byte cbuffer = (buffer << (8 - bitsInBuffer)) |
              (nextb >> bitsInBuffer);
            assert(lookupWhich[zeroCount][cbuffer] < 256);
            runseq[haveDecoded++] = lookupWhich[zeroCount][cbuffer];
            bitsInBuffer = bitsInBuffer + (8 - lookupLength[zeroCount][cbuffer]);
            zeroCount = leadingZeros[bitsInBuffer][nextb];
            bitsInBuffer -= zeroCount;
            buffer = nextb;
            nextb &= (1 << bitsInBuffer) - 1;

            // Extract codes from nextb.
            while (haveDecoded < nRuns) { bool decoded = false;
              byte sbuffer = (buffer << (8 - bitsInBuffer));
              if (lookupWhich[zeroCount][sbuffer] < 256  &&
                  bitsInBuffer >= lookupLength[zeroCount][sbuffer]) {
                runseq[haveDecoded++] = lookupWhich[zeroCount][sbuffer];
                bitsInBuffer -= lookupLength[zeroCount][sbuffer];
                zeroCount = leadingZeros[bitsInBuffer][buffer];
                bitsInBuffer -= zeroCount; decoded = true;
                buffer &= (1 << bitsInBuffer) - 1;
              }
              while (zeroCount >= len0 && haveDecoded < nRuns) {
                runseq[haveDecoded++] = code0;
                zeroCount -= len0; decoded = true;
              }
              if (!decoded) break;
            }
            continue;
          }
        }
      } else if (!zeroCount && !bitsInBuffer) {
      
        // Extract codes from b.
        buffer = b;
        bitsInBuffer = 8;
        zeroCount = leadingZeros[bitsInBuffer][buffer];
        bitsInBuffer -= zeroCount;
        while (haveDecoded < nRuns) { bool decoded = false;
          byte sbuffer = (buffer << (8 - bitsInBuffer));
          if (lookupWhich[zeroCount][sbuffer] < 256
              && bitsInBuffer >= lookupLength[zeroCount][sbuffer]) {
            runseq[haveDecoded++] = lookupWhich[zeroCount][sbuffer];
            bitsInBuffer -= lookupLength[zeroCount][sbuffer];
            zeroCount = leadingZeros[bitsInBuffer][buffer];
            bitsInBuffer -= zeroCount;
            decoded = true;
            buffer &= (1 << bitsInBuffer) - 1;
          }
          while (zeroCount >= len0 && haveDecoded < nRuns) {
            runseq[haveDecoded++] = code0;
            zeroCount -= len0;
            decoded = true;
          }
          if (!decoded) break;
        }
        continue;
      } else { // bitsInBuffer > 0

        // Previous byte (buffer) contains nonzero non-consumed bits, hence
        // this code will end in b.
        byte cbuffer = (buffer << (8 - bitsInBuffer)) | (b >> bitsInBuffer);
        assert(lookupWhich[zeroCount][cbuffer] < 256);
        runseq[haveDecoded++] = lookupWhich[zeroCount][cbuffer];
        bitsInBuffer = bitsInBuffer + (8 - lookupLength[zeroCount][cbuffer]);
        zeroCount = leadingZeros[bitsInBuffer][b];
        bitsInBuffer -= zeroCount;
        buffer = b;
        buffer &= (1 << bitsInBuffer) - 1;

        // Extract codes from b.
        while (haveDecoded < nRuns) {
          bool decoded = false;
          byte sbuffer = (buffer << (8 - bitsInBuffer));
          if (lookupWhich[zeroCount][sbuffer] < 256 &&
              bitsInBuffer >= lookupLength[zeroCount][sbuffer]) {
            runseq[haveDecoded++] = lookupWhich[zeroCount][sbuffer];
            bitsInBuffer -= lookupLength[zeroCount][sbuffer];
            zeroCount = leadingZeros[bitsInBuffer][buffer];
            bitsInBuffer -= zeroCount;
            decoded = true;
            buffer &= (1 << bitsInBuffer) - 1;
          }
          while (zeroCount >= len0 && haveDecoded < nRuns) {
            runseq[haveDecoded++] = code0;
            zeroCount -= len0;
            decoded = true;
          }
          if (!decoded) break;
        }
        continue;
      }
    }
    in->flushBuffer();

    // Now read gamma codes that store lenghts of runs.
    for (uint64 k = 0; k < nRuns; ++k) {
      int zeros = 0;
      while (!in->readBit())
        ++zeros;
      uint64 value = 0;        
      for (int32 t = 0; t < zeros; ++t) {
        int32 bit = in->readBit();
        value = (value << 1) | bit;
      }
      value |= (1 << zeros);
      runlen[k] = value;
    }
    in->flushBuffer();
    
    // Fill the block with runs data.
    byte *runseq_ptr = &runseq[0];
    for (uint64 k = 0; k < nRuns; ++k) {
      for (uint32 t = 0; t < runlen[k]; ++t)
        *data_ptr++ = *runseq_ptr;
      ++runseq_ptr;
    }
  }

  block.setSize(block_size);
  
  delete[] runseq;
  delete[] runlen;
}

uint64 HuffmanDecoder::readPackedInteger(InStream* in) {
  static const uint64 kEndSymbol = static_cast<uint64>(1) << 63;
  static const uint64 kEndMask = static_cast<uint64>(1) << 7;

  uint64 packed_integer = 0;
  bool bits_left = true;
  int i;
  for(i = 0; bits_left; ++i) {
    uint64 read = static_cast<uint64>(in->readByte());
    bits_left = (read & kEndMask) != 0;
    packed_integer |= (read << i*8);
  }
  if (packed_integer == 0x80) return kEndSymbol;
  return packed_integer;
}
/*********** Encoding and decoding single BWTBlock-section ends ********/

HuffmanDecoder::HuffmanDecoder() {}

HuffmanDecoder::~HuffmanDecoder() {}


} // namespace bwtc

