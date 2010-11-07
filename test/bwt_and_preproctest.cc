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

#include "../preprocessors/longsequences.h"
#include "../preprocessors/test_preprocessor.h"
#include "../preprocessors/postprocessor.h"
#include "../globaldefs.h"
#include "../bwtransforms/dcbwt.h"
#include "../bwtransforms/bw_transform.h"

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
  bwtc::BlockManager bm(block_size, 1);
  bwtc::TestPreProcessor pp(block_size);
  bwtc::BWTransform *transformer = bwtc::GiveTransformer('s');
  pp.AddBlockManager(&bm);
  pp.Connect(input);
  pp.InitializeTarget();
  uint64 orig_size = pp.FillBuffer();
  uint64 compressed_size = orig_size;

  /* Preprocessing */
  int str_index = 0;
  clock_t prepr_start = clock();
  do {
    switch (preprocs[str_index]) {
      case 'l': /* Replace long sequences */
        compressed_size = bwtc::CompressSequences(pp.curr_block_->begin()
                                                  ,pp.curr_block_->filled_
                                                  ,mem_constr, 16 /*Window size*/
                                                  ,threshold);
        pp.curr_block_->filled_ = compressed_size;
        break;
      case 'r': /* Replace runs of same byte */
        compressed_size -= pp.CompressRuns();
        break;
      case 'p': /* Replace common pairs */
        compressed_size -= pp.CompressPairs();
        break;
      default:
        break;
    }
  } while(preprocs[++str_index]);
  clock_t prepr_end = clock();

  /* Burrows-Wheeler Transform */
  uint64 eob;
  transformer->Connect(pp.curr_block_);
  clock_t bwt_start = clock();
  std::vector<byte> *result = transformer->DoTransform(&eob);
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
  bwtc::verbosity = 3;
  bwtc::BlockManager bm(block_size, 1);
  bwtc::TestPreProcessor pp(block_size);
  pp.AddBlockManager(&bm);
  pp.Connect(input);
  pp.InitializeTarget();
  uint64 orig_size = pp.FillBuffer();
  uint64 compressed_size = orig_size;

  std::vector<byte> original(orig_size);
  std::copy(pp.curr_block_->begin(), pp.curr_block_->end(), original.begin());
  
  /* Preprocessing */
  std::cout << "############################\n";
  std::cout << "Preprocessing\n";
  int str_index = 0;
  do {
    switch (preprocs[str_index]) {
      case 'l': /* Replace long sequences */
        compressed_size = bwtc::CompressSequences(pp.curr_block_->begin()
                                                  ,pp.curr_block_->filled_
                                                  ,mem_constr, 16 /*Window size*/
                                                  ,threshold);
        pp.curr_block_->filled_ = compressed_size;
        break;
      case 'r': /* Replace runs of same byte */
        compressed_size -= pp.CompressRuns();
        break;
      default: /* Replace common pairs */
        compressed_size -= pp.CompressPairs();
        break;
    }
  } while(preprocs[++str_index]);

  --str_index;
  uint64 uncompr_size = 0;
  std::cout << "############################\n";
  std::cout << "Postprocessing\n";
  do {
    switch(preprocs[str_index]) {
      case 'l':
        uncompr_size = bwtc::UncompressSequences(pp.curr_block_->block_,
                                                 pp.curr_block_->filled_);
        pp.curr_block_->filled_ = uncompr_size;
        break;
      case 'r':
        uncompr_size = bwtc::UncompressLongRuns(pp.curr_block_->block_,
                                                pp.curr_block_->filled_);
        pp.curr_block_->filled_ = uncompr_size;
        break;
      default:
        uncompr_size = bwtc::UncompressCommonPairs(pp.curr_block_->block_,
                                                   pp.curr_block_->filled_);
        pp.curr_block_->filled_ = uncompr_size; 
    }
      
  } while(--str_index >= 0);

  std::vector<byte>& uncompressed = *pp.curr_block_->block_;
  assert(uncompr_size == original.size());
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
  if (argc < 4) return 0;
  if (argc > 4) mem_constr = atoi(argv[4]);
  if (argc > 5) block_size = atoi(argv[5]);
  tests::PreprocBWTSpeed(argv[1], atoi(argv[2]), std::string(argv[3]),
                         block_size, mem_constr);
  tests::ValidatePreproc(argv[1], atoi(argv[2]), std::string(argv[3]),
                         block_size, mem_constr);

}
