/**
 * @file bwt_and_preproctest.cc
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
 * File for testing the speed of preprocessing algorithms and their
 * effect to the speed of Burrows-Wheeler transform.
 */
#include <cassert>
#include <cstdlib>
#include <ctime>

#include <iostream>
#include <string>
#include <vector>

#include "../preprocessors/TestPreprocessor.hpp"
#include "../preprocessors/Postprocessor.hpp"
#include "../globaldefs.hpp"
#include "../bwtransforms/BWTransform.hpp"

namespace bwtc {
int verbosity = 0;
}

using bwtc::uint64;
using bwtc::byte;

namespace tests {

const int kTimes = 1;

void PreprocBWTSpeed(char *preprocs, int threshold, const std::string& input,
                     uint64 block_size, unsigned mem_constr)
{
  (void) threshold;
  (void) mem_constr;
  bwtc::BlockManager bm(block_size, 1);
  bwtc::TestPreprocessor pp(block_size);
  bwtc::BWTransform *transformer = bwtc::giveTransformer('s');
  pp.addBlockManager(&bm);
  pp.connect(input);
  pp.initializeTarget();
  uint64 orig_size = pp.fillBuffer();
  uint64 compressed_size = orig_size;

  /* Preprocessing */
  int str_index = 0;
  clock_t prepr_start = clock();

  do {
    switch (preprocs[str_index]) {
      /*
      case 'l': // Replace long sequences 
        compressed_size = bwtc::compressSequences(pp.m_currentBlock->begin()
                                                  ,pp.m_currentBlock->filled_
                                                  ,mem_constr, 16 //Window size
                                                  ,threshold);
        pp.m_currentBlock->filled_ = compressed_size;
        break;*/
      case 'c': /* Replcae pairs and runs.*/
        compressed_size -= pp.compressPairsAndRuns();
        break;
      case 'r': /* Replace runs of same byte */
        compressed_size -= pp.compressRuns();
        break;
      case 'p': /* Replace common pairs */
        compressed_size -= pp.compressPairs();
        break;
      default:
        break;
    }
  } while(preprocs[++str_index]);
  clock_t prepr_end = clock();

  /* Burrows-Wheeler Transform */
  uint64 eob;
  transformer->connect(pp.m_currentBlock);
  clock_t bwt_start = clock();
  std::vector<byte> *result = transformer->doTransform(&eob);
  clock_t bwt_end = clock();
  delete result;
  delete transformer;

  std::cout << "##################################################\n"
            << "Time spent on Burrows-Wheeler Transform: "
            << (bwt_end - bwt_start)/static_cast<double>(CLOCKS_PER_SEC)
            << "\n";
  std::cout << "Time spent on preprocessing: "
            << (prepr_end - prepr_start)/static_cast<double>(CLOCKS_PER_SEC)
            << "\n";
  std::cout << "############################\n";
  std::cout << "Size of data was " << orig_size << "B\n";
  std::cout << "Result of preprocessing was " << compressed_size
            << "B which is "
            << (compressed_size/static_cast<double>(orig_size))*100.0
            << "% of the original data\n";
}

void ValidatePreproc(char *preprocs, int threshold, const std::string& input,
                     uint64 block_size, unsigned mem_constr)
{
  (void) threshold;
  (void) mem_constr;
  bwtc::verbosity = 3;
  bwtc::BlockManager bm(block_size, 1);
  bwtc::TestPreprocessor pp(block_size);
  pp.addBlockManager(&bm);
  pp.connect(input);
  pp.initializeTarget();
  uint64 orig_size = pp.fillBuffer();
  uint64 compressed_size = orig_size;

  std::vector<byte> original(orig_size);
  std::copy(pp.m_currentBlock->begin(), pp.m_currentBlock->end(), original.begin());
  
  /* Preprocessing */
  std::cout << "############################\n";
  std::cout << "Preprocessing\n";
  int str_index = 0;
  do {
    switch (preprocs[str_index]) {
      /*
      case 'l': // Replace long sequences 
        compressed_size = bwtc::compressSequences(pp.m_currentBlock->begin()
                                                  ,pp.m_currentBlock->filled_
                                                  ,mem_constr, 16 //Window size
                                                  ,threshold);
        pp.m_currentBlock->filled_ = compressed_size;
        break;*/
      case 'c': /* Replcae pairs and runs.*/
        compressed_size -= pp.compressPairsAndRuns();
        break;
      case 'r': /* Replace runs of same byte */
        compressed_size -= pp.compressRuns();
        break;
      default: /* Replace common pairs */
        compressed_size -= pp.compressPairs();
        break;
    }
  } while(preprocs[++str_index]);

  --str_index;
  uint64 uncompr_size = 0;
  std::cout << "############################\n";
  std::cout << "Postprocessing\n";
  do {
    switch(preprocs[str_index]) {
      /*
      case 'l':
        uncompr_size = bwtc::uncompressSequences(pp.m_currentBlock->block_,
                                                 pp.m_currentBlock->filled_);
        pp.m_currentBlock->filled_ = uncompr_size;
        break;*/
      case 'c': /* Replcae pairs and runs.*/
        uncompr_size = bwtc::postprocessor::
            uncompressPairsAndRuns(pp.m_currentBlock->m_block,
                                   pp.m_currentBlock->m_filled);
        pp.m_currentBlock->m_filled = uncompr_size;
        break;
      case 'r':
        uncompr_size = bwtc::postprocessor::
            uncompressLongRuns(pp.m_currentBlock->m_block,
                               pp.m_currentBlock->m_filled);
        pp.m_currentBlock->m_filled = uncompr_size;
        break;
      default:
        uncompr_size = bwtc::postprocessor::
            uncompressCommonPairs(pp.m_currentBlock->m_block,
                                  pp.m_currentBlock->m_filled);
        pp.m_currentBlock->m_filled = uncompr_size; 
    }
      
  } while(--str_index >= 0);

  std::vector<byte>& uncompressed = *pp.m_currentBlock->m_block;
  std::cout << uncompr_size << " " <<  original.size() << std::endl;
  for(uint64 j = 0; j < original.size(); ++j) {
    assert(uncompressed[j] == original[j]);
  }
  std::cout << "############################\n";
  std::cout << "Size of data was " << orig_size << "B\n";
  std::cout << "Result of preprocessing was " << compressed_size
            << "B which is "
            << (compressed_size/static_cast<double>(orig_size))*100.0
            << "% of the original data\n############################\n";

}

} // namespace tests

/****************************************************************************
 * The syntax for using the program is the following:                       *
 * ./preprocspeedtest <preproc algorithms string> <t> <input file> <c> <b>  *
 * where t is threshold (integer value) for replacing long strings.         *
 * Preprocessor algorithms string consists of three different characters,   *
 * namely 'l','p' and 'r', where each of them denotes single preprocessor   *
 * algorithm (long sequences, pairs and long runs). For example if given    *
 * string "lrpp" then we first replace long sequences, then long runs and   *
 * finally we replace common pairs twice. Parameters c and b are optional.  *
 * c is memory constraint (c bytes per single input byte for use) and b is  *
 * buffer size.                                                             *
 ****************************************************************************/
int main(int argc, char **argv) {
  uint64 block_size = 259715200;
  unsigned mem_constr = 2;
  if (argc < 4) {
    std::cout << "Usage: " << argv[0] << " [clpr]+ <n> input_file "
              << "[mem_constr] [buffer]" << std::endl;
    return 0;
  }
  if (argc > 4) mem_constr = atoi(argv[4]);
  if (argc > 5) block_size = atoi(argv[5]);
  tests::ValidatePreproc(argv[1], atoi(argv[2]), std::string(argv[3]),
                         block_size, mem_constr);
  tests::PreprocBWTSpeed(argv[1], atoi(argv[2]), std::string(argv[3]),
                         block_size, mem_constr);
}
