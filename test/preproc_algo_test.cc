/* Measurement of spent clock cycles is originally from *
 * http://www.agner.org/optimize/ */
#include <cassert>
#include <cstdlib>

#include <iostream>
#include <string>

#include "../preprocessors/test_preprocessor.h"
#include "../globaldefs.h"
#include "../bwtransforms/dcbwt.h"
#include "../bwtransforms/bw_transform.h"

namespace bwtc {
int verbosity = 0;
}

namespace tests {

const int kTimes = 20;

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

void TestPairCompression(std::string source_name, int times, uint64 block_size)
{
  uint64 total_cycles_preproc = 0;
  uint64 total_cycles_bwt = 0;
  uint64 total_data, total_reduction;

  bwtc::BlockManager bm(block_size, 1);
  bwtc::BWTransform *transformer = bwtc::GiveTransformer(); 
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
    uint64 eob;
    __int64 beginBWT = ReadTSC();
    transformer->Connect(pp.curr_block_);
    transformer->BuildStats();
    std::vector<byte>* result = transformer->DoTransform(&eob);
    __int64 endBWT = ReadTSC();
    total_cycles_bwt += (endBWT - beginBWT);
    delete result;
    std::cout << ".";
    std::cout.flush();
  }
  std::cout << "\n";
  total_cycles_preproc /= kTimes;
  total_cycles_bwt /= kTimes;
  std::cout << "Average CPU-cycles spent on preprocessing: "
            << total_cycles_preproc << "\n";
  std::cout << "Size of data was " << total_data << "B\n";
  std::cout << "Result of preprocessing was " << total_data - total_reduction
            << "B which is "
            << (1.0 - (total_reduction/static_cast<double>(total_data)))*100.0
            << "% of the original data\n";
  std::cout << "Average CPU-cycles spent on Burrows-Wheerer Transform: "
            << total_cycles_bwt << "\n";
  std::cout << "Ratio of (bwt cycles)/(preproc cycles) is "
            << static_cast<double>(total_cycles_bwt)/total_cycles_preproc
            << "\n";
  delete transformer;
}

} // namespace tests

int main(int argc, char **argv) {
  if(argc > 2)
    tests::TestPairCompression(std::string(argv[1]), atoi(argv[2]),
                               1024*100000);
}
