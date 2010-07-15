/* Measurement of spent clock cycles is originally from *
 * http://www.agner.org/optimize/ */
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

namespace {
/* Code in this namespace is from http://www.agner.org/optimize/ */
#if ! defined(_MSC_VER) 
// define type __int64 if not defined
typedef int64 __int64;
#endif

#ifdef USE_ALIB    // use alibXXX.lib
// for compilers with insufficient inline assembly support, use external
// library for function ReadTSC()
#include "alib.h"

#else
// for compilers that understand the inline assembly code you can define
// the function ReadTSC() here:

// Read time stamp counter
// The return value is the internal clock count
__int64 ReadTSC() {
   int res[2];                              // store 64 bit result here
   
   #if defined(__GNUC__) && !defined(__INTEL_COMPILER)
   // Inline assembly in AT&T syntax

   #if defined (_LP64)                      // 64 bit mode
      __asm__ __volatile__  (               // serialize (save rbx)
      "xorl %%eax,%%eax \n push %%rbx \n cpuid \n"
       ::: "%rax", "%rcx", "%rdx");
      __asm__ __volatile__  (               // read TSC, store edx:eax in res
      "rdtsc\n"
       : "=a" (res[0]), "=d" (res[1]) );
      __asm__ __volatile__  (               // serialize again
      "xorl %%eax,%%eax \n cpuid \n pop %%rbx \n"
       ::: "%rax", "%rcx", "%rdx");
   #else                                    // 32 bit mode
      __asm__ __volatile__  (               // serialize (save ebx)
      "xorl %%eax,%%eax \n pushl %%ebx \n cpuid \n"
       ::: "%eax", "%ecx", "%edx");
      __asm__ __volatile__  (               // read TSC, store edx:eax in res
      "rdtsc\n"
       : "=a" (res[0]), "=d" (res[1]) );
      __asm__ __volatile__  (               // serialize again
      "xorl %%eax,%%eax \n cpuid \n popl %%ebx \n"
       ::: "%eax", "%ecx", "%edx");
   #endif
   #else
   // Inline assembly in MASM syntax
      __asm {
         xor eax, eax
         cpuid                              // serialize
         rdtsc                              // read TSC
         mov dword ptr res, eax             // store low dword in res[0]
         mov dword ptr res+4, edx           // store high dword in res[1]
         xor eax, eax
         cpuid                              // serialize again
      };
   #endif   // __GNUC__
   
   return *(__int64*)res;                   // return result
}

#endif   // USE_ALIB
}

void TestRunUncompression(std::string source, int times, uint64 block_size)
{
  uint64 total_cycles_postproc = 0;
  uint64 total_cycles_preproc = 0;
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
    __int64 beginPre = ReadTSC();
    for(int j = 0; j < times; ++j) {
      reduction += pp.CompressRuns();
    }
    
    __int64 endPre = ReadTSC();
    total_reduction = reduction;
    total_cycles_preproc += (endPre - beginPre);
    assert(total_reduction == total_data - pp.curr_block_->filled_);
    uint64 uncompressed_size;
    /* Make sure that we can also uncompress the thing */
    __int64 beginPost = ReadTSC();
    for(int j = 0; j < times; ++j) {
      uncompressed_size = bwtc::UncompressLongRuns(pp.curr_block_->block_,
                                                   pp.curr_block_->filled_);
      /* This one should be done in same kind of wrapper than
       * Preprocessor::CompressPairs*/
      pp.curr_block_->filled_ = uncompressed_size; 
    }
    __int64 endPost = ReadTSC();
    total_cycles_postproc += (endPost - beginPost);
    std::vector<byte>& uncompressed = *pp.curr_block_->block_;

    std::cout << uncompressed_size << " " << original.size() << "\n";
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
  total_cycles_postproc /= kTimes;
  total_cycles_preproc /= kTimes;
  std::cout << "Average CPU-cycles spent on preprocessing: "
            << total_cycles_preproc << "\n";
  std::cout << "Size of data was " << total_data << "B\n";
  std::cout << "Result of preprocessing was " << total_data - total_reduction
            << "B which is "
            << (1.0 - (total_reduction/static_cast<double>(total_data)))*100.0
            << "% of the original data\n";
  std::cout << "Average CPU-cycles spent on postprocessing: "
            << total_cycles_postproc << "\n";
}

void TestPairUncompression(std::string source, int times, uint64 block_size)
{
  uint64 total_cycles_postproc = 0;
  
  bwtc::BlockManager bm(block_size, 1);
  //bwtc::BWTransform *transformer = bwtc::GiveTransformer(); 
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
    uint64 uncompressed_size;
    /* Make sure that we can also uncompress the thing */
    __int64 beginPost = ReadTSC();
    for(int j = 0; j < times; ++j) {
      uncompressed_size = bwtc::UncompressCommonPairs(pp.curr_block_->block_,
                                                      pp.curr_block_->filled_);
      /* This one should be done in same kind of wrapper than
       * Preprocessor::CompressPairs*/
      pp.curr_block_->filled_ = uncompressed_size; 
    }
    __int64 endPost = ReadTSC();
    total_cycles_postproc += (endPost - beginPost);
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
  total_cycles_postproc /= kTimes;
  std::cout << "Average CPU-cycles spent on postprocessing: "
            << total_cycles_postproc << "\n";
}

void TestPairCompression(std::string source_name, int times, uint64 block_size)
{
  uint64 total_cycles_preproc = 0;
  uint64 total_cycles_bwt = 0;
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
    __int64 beginPre = ReadTSC();
    for(int j = 0; j < times; ++j) {
      data_size += pp.FillBuffer();
      data_reduction += pp.CompressPairs();
    }
    __int64 endPre = ReadTSC();
    total_cycles_preproc += (endPre - beginPre);
    total_data = data_size;
    total_reduction = data_reduction;
    /*    uint64 eob;
    __int64 beginBWT = ReadTSC();
    transformer->Connect(pp.curr_block_);
    transformer->BuildStats();
    std::vector<byte>* result = transformer->DoTransform(&eob);
    __int64 endBWT = ReadTSC();
    total_cycles_bwt += (endBWT - beginBWT);
    delete result;*/
    /* Make sure that we can also uncompress the thing */
    if (bwtc::verbosity < 2) {
      std::cout << ".";
      std::cout.flush();
    }
  }
  if (bwtc::verbosity < 2) std::cout << "\n";
  total_cycles_preproc /= kTimes;
  total_cycles_bwt /= kTimes;
  std::cout << "Average CPU-cycles spent on preprocessing: "
            << total_cycles_preproc << "\n";
  std::cout << "Size of data was " << total_data << "B\n";
  std::cout << "Result of preprocessing was " << total_data - total_reduction
            << "B which is "
            << (1.0 - (total_reduction/static_cast<double>(total_data)))*100.0
            << "% of the original data\n";
  /*  std::cout << "Average CPU-cycles spent on Burrows-Wheerer Transform: "
            << total_cycles_bwt << "\n";
  std::cout << "Ratio of (bwt cycles)/(preproc cycles) is "
            << static_cast<double>(total_cycles_bwt)/total_cycles_preproc
            << "\n";
            delete transformer;*/
}

void TestComboCompression(std::string source_name, int times, uint64 block_size)
{
  uint64 total_cycles_preproc = 0;
  uint64 total_cycles_postproc = 0;
  uint64 total_data, total_reduction;

  bwtc::BlockManager bm(block_size, 1);
  for(int i = 0; i < kTimes; ++i) {
    bwtc::TestPreProcessor pp(block_size);
    pp.AddBlockManager(&bm);
    pp.Connect(source_name);
    pp.InitializeTarget();
    uint64 data_size = pp.FillBuffer();
    std::vector<byte> original(pp.curr_block_->filled_);
    std::copy(pp.curr_block_->begin(), pp.curr_block_->end(), original.begin());

    uint64 data_reduction = 0;
    __int64 beginPre = ReadTSC();
    for(int j = 0; j < times; ++j) {
      data_reduction += pp.CompressPairs();
      data_reduction += pp.CompressRuns();
    }
    __int64 endPre = ReadTSC();
    total_cycles_preproc += (endPre - beginPre);
    total_data = data_size;
    total_reduction = data_reduction;

    uint64 uncompressed_size;
    __int64 beginPost = ReadTSC();
    for(int j = 0; j < times; ++j) {
      uncompressed_size = bwtc::UncompressLongRuns(pp.curr_block_->block_,
                                                   pp.curr_block_->filled_);
      pp.curr_block_->filled_ = uncompressed_size; 
      uncompressed_size = bwtc::UncompressCommonPairs(pp.curr_block_->block_,
                                                      pp.curr_block_->filled_);
      pp.curr_block_->filled_ = uncompressed_size; 
    }
    __int64 endPost = ReadTSC();
    total_cycles_postproc += (endPost - beginPost);
    std::vector<byte>& uncompressed = *pp.curr_block_->block_;
    std::cout << uncompressed_size << " " << original.size() << "\n";
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
  total_cycles_preproc /= kTimes;
  total_cycles_postproc /= kTimes;
  std::cout << "Combo compression:\n";
  std::cout << "Average CPU-cycles spent on preprocessing: "
            << total_cycles_preproc << "\n";
  std::cout << "Size of data was " << total_data << "B\n";
  std::cout << "Result of preprocessing was " << total_data - total_reduction
            << "B which is "
            << (1.0 - (total_reduction/static_cast<double>(total_data)))*100.0
            << "% of the original data\n";
  std::cout << "Average CPU-cycles spent on postprocessing: "
            << total_cycles_postproc << "\n";
}

} // namespace tests

int main(int argc, char **argv) {
  int times = 1;
  uint64 block_size = 209715200;
  if(argc > 2) times =  atoi(argv[2]);
  if (argc > 3) block_size = atoi(argv[3]);
  if(argc == 1) return 0;
  tests::TestRunUncompression(std::string(argv[1]), times, block_size);
  tests::TestComboCompression(std::string(argv[1]), times, block_size);
  tests::TestPairCompression(std::string(argv[1]), times, block_size);
  tests::TestPairUncompression(std::string(argv[1]), times, block_size);
}
