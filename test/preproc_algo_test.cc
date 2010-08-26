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
#include <cstdlib>

#include <iostream>
#include <string>
#include <vector>

#include "../preprocessors/test_preprocessor.h"
#include "../preprocessors/postprocessor.h"
#include "../globaldefs.h"
#include "../bwtransforms/dcbwt.h"
#include "../bwtransforms/bw_transform.h"

namespace bwtc {
int verbosity = 2;
}

namespace tests {

const int kTimes = 1;

void TestRunUncompression(std::string source, int times, uint64 block_size)
{
  uint64 total_data, total_reduction;

  bwtc::BlockManager bm(block_size, 1);
  for(int i = 0; i < kTimes; ++i) {
    bwtc::TestPreProcessor pp(block_size);
    pp.AddBlockManager(&bm);
    pp.Connect(source);
    pp.InitializeTarget();
    total_data = pp.FillBuffer();
    std::vector<byte> original(pp.curr_block_->filled_);
    std::copy(pp.curr_block_->begin(), pp.curr_block_->end(), original.begin());
    assert(total_data == pp.curr_block_->filled_);
    assert(total_data == original.size());
    uint64 reduction = 0;
    for(int j = 0; j < times; ++j) {
      reduction += pp.CompressRuns();
    }
    
    total_reduction = reduction;
    assert(total_reduction == total_data - pp.curr_block_->filled_);
    uint64 uncompressed_size = 0;
    /* Make sure that we can also uncompress the thing */
    for(int j = 0; j < times; ++j) {
      uncompressed_size = bwtc::UncompressLongRuns(pp.curr_block_->block_,
                                                   pp.curr_block_->filled_);
      /* This one should be done in same kind of wrapper than
       * Preprocessor::CompressPairs*/
      pp.curr_block_->filled_ = uncompressed_size; 
    }
    std::vector<byte>& uncompressed = *pp.curr_block_->block_;

    assert(uncompressed_size == original.size());
    for(uint64 j = 0; j < uncompressed_size; ++j) {
      assert(uncompressed[j] == original[j]);
    }
    if (bwtc::verbosity < 2) {
      std::cout << ".";
      std::cout.flush();
    }
  }
  if (bwtc::verbosity < 2) std::cout << "\n";
  std::cout << "Run compression:\n";
  std::cout << "Size of data was " << total_data << "B\n";
  std::cout << "Result of preprocessing was " << total_data - total_reduction
            << "B which is "
            << (1.0 - (total_reduction/static_cast<double>(total_data)))*100.0
            << "% of the original data\n############################\n";
}

void TestPairUncompression(std::string source, int times, uint64 block_size)
{
  bwtc::BlockManager bm(block_size, 1);
  for(int i = 0; i < kTimes; ++i) {
    bwtc::TestPreProcessor pp(block_size);
    pp.AddBlockManager(&bm);
    pp.Connect(source);
    pp.InitializeTarget();
    uint64 data_size = 0;
    data_size = pp.FillBuffer();
    std::vector<byte> original(pp.curr_block_->filled_);
    std::copy(pp.curr_block_->begin(), pp.curr_block_->end(), original.begin());
    assert(data_size == pp.curr_block_->filled_);
    assert(data_size == original.size());
    uint64 compressed = 0;
    for(int j = 0; j < times; ++j) {
      compressed += pp.CompressPairs();
    }
    assert(compressed == data_size - pp.curr_block_->filled_);
    uint64 uncompressed_size = 0;
    /* Make sure that we can also uncompress the thing */
    for(int j = 0; j < times; ++j) {
      uncompressed_size = bwtc::UncompressCommonPairs(pp.curr_block_->block_,
                                                      pp.curr_block_->filled_);
      /* This one should be done in same kind of wrapper than
       * Preprocessor::CompressPairs*/
      pp.curr_block_->filled_ = uncompressed_size; 
    }
    std::vector<byte>& uncompressed = *pp.curr_block_->block_;

    assert(uncompressed_size == original.size());
    for(uint64 j = 0; j < data_size; ++j) {
      assert(uncompressed[j] == original[j]);
    }
    if (bwtc::verbosity < 2) {
      std::cout << ".";
      std::cout.flush();
    }
  }
  if (bwtc::verbosity < 2) std::cout << "\n";
}

void TestPairCompression(std::string source_name, int times, uint64 block_size)
{
  uint64 total_data, total_reduction;

  bwtc::BlockManager bm(block_size, 1);
  //bwtc::BWTransform *transformer = bwtc::GiveTransformer(); 
  for(int i = 0; i < kTimes; ++i) {
    bwtc::TestPreProcessor pp(block_size);
    pp.AddBlockManager(&bm);
    pp.Connect(source_name);
    pp.InitializeTarget();
    uint64 data_size = 0;
    uint64 data_reduction = 0;
    for(int j = 0; j < times; ++j) {
      data_size += pp.FillBuffer();
      data_reduction += pp.CompressPairs();
    }
    total_data = data_size;
    total_reduction = data_reduction;
    /*    uint64 eob;
    transformer->Connect(pp.curr_block_);
    transformer->BuildStats();
    std::vector<byte>* result = transformer->DoTransform(&eob);
    delete result;*/
    /* Make sure that we can also uncompress the thing */
    if (bwtc::verbosity < 2) {
      std::cout << ".";
      std::cout.flush();
    }
  }
  if (bwtc::verbosity < 2) std::cout << "\n";
  std::cout << "Pair compression:\n";
  std::cout << "Size of data was " << total_data << "B\n";
  std::cout << "Result of preprocessing was " << total_data - total_reduction
            << "B which is "
            << (1.0 - (total_reduction/static_cast<double>(total_data)))*100.0
            << "% of the original data\n############################\n";
  /*  std::cout << "Average CPU-cycles spent on Burrows-Wheerer Transform: "
            << total_cycles_bwt << "\n";
  std::cout << "Ratio of (bwt cycles)/(preproc cycles) is "
            << static_cast<double>(total_cycles_bwt)/total_cycles_preproc
            << "\n";
            delete transformer;*/
}

void TestComboCompression(std::string source_name, int times, uint64 block_size)
{
  uint64 total_data, total_reduction;

  bwtc::BlockManager bm(block_size, 1);
  for(int i = 0; i < kTimes; ++i) {
    bwtc::TestPreProcessor pp(block_size);
    pp.AddBlockManager(&bm);
    pp.Connect(source_name);
    pp.InitializeTarget();
    uint64 data_size = pp.FillBuffer();
    std::vector<byte> original(data_size);
    std::copy(pp.curr_block_->begin(), pp.curr_block_->end(), original.begin());
    uint64 data_reduction = 0;
    for(int j = 0; j < times; ++j) {
      data_reduction += pp.CompressPairs();
      data_reduction += pp.CompressRuns();
    }
    total_data = data_size;
    total_reduction = data_reduction;

    uint64 uncompressed_size = 0;
    for(int j = 0; j < times; ++j) {
      uncompressed_size = bwtc::UncompressLongRuns(pp.curr_block_->block_,
                                                   pp.curr_block_->filled_);
      pp.curr_block_->filled_ = uncompressed_size; 
      uncompressed_size = bwtc::UncompressCommonPairs(pp.curr_block_->block_,
                                                      pp.curr_block_->filled_);
      pp.curr_block_->filled_ = uncompressed_size; 
    }
    std::vector<byte>& uncompressed = *pp.curr_block_->block_;

    std::cout << uncompressed_size << " " << data_size << "\n";
    assert(uncompressed_size == original.size());
    for(uint64 j = 0; j < uncompressed_size; ++j) {
      assert(uncompressed[j] == original[j]);
    }
    
    if (bwtc::verbosity < 2) {
      std::cout << ".";
      std::cout.flush();
    }
  }
  if (bwtc::verbosity < 2) std::cout << "\n";
  std::cout << "Combo compression:\n";
  std::cout << "Size of data was " << total_data << "B\n";
  std::cout << "Result of preprocessing was " << total_data - total_reduction
            << "B which is "
            << (1.0 - (total_reduction/static_cast<double>(total_data)))*100.0
            << "% of the original data\n############################\n";
}

} // namespace tests

int main(int argc, char **argv) {
  int times = 1;
  uint64 block_size = 259715200;
  if(argc > 2) times =  atoi(argv[2]);
  if (argc > 3) block_size = atoi(argv[3]);
  if(argc == 1) return 0;
  tests::TestRunUncompression(std::string(argv[1]), times, block_size);
  tests::TestComboCompression(std::string(argv[1]), times, block_size);
  tests::TestPairCompression(std::string(argv[1]), times, block_size);
  tests::TestPairUncompression(std::string(argv[1]), times, block_size);
}
