/**
 * @file Preprocessor.hpp
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
 * Header for preprocessor.
 */


#ifndef BWTC_PREPROCESSOR_HPP_
#define BWTC_PREPROCESSOR_HPP_

#include <iostream> /* for std::streamsize*/
#include <string>

#include "../MainBlock.hpp"
#include "../BlockManager.hpp"
#include "../globaldefs.hpp"
#include "../Streams.hpp"

namespace bwtc {

//TODO: Make this class to abstract class (interface for real implementations)
class Preprocessor {
 public:
  Preprocessor(uint64 block_size);
  virtual ~Preprocessor();
  virtual void connect(std::string source_name);
  virtual void addBlockManager(BlockManager* bm);
  /* Reads and preprocesses data to byte array provided by block_manager_*/
  virtual MainBlock* readBlock();

 protected:
  InStream* m_source;
  uint64 m_blockSize;
  BlockManager* m_blockManager;

 private:
  /* This should be done during preprocessing*/
  void buildStats(std::vector<byte>* data, std::vector<uint64>* stats,
                  uint64 data_size);
  Preprocessor& operator=(const Preprocessor& p);
  Preprocessor(const Preprocessor&);
};

/* This function returns chosen preprocessor */ 
Preprocessor* givePreprocessor(
    char choice, uint64 block_size, const std::string& input);

uint64 compressCommonPairs(byte *from, uint64 length);
uint64 compressLongRuns(byte *from, uint64 length);

} // namespace bwtc


#endif
