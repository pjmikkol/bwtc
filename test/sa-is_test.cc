#include <cstdlib>
#include <cstring>
#include <ctime>

#include <iostream>

#include "../globaldefs.h"
#include "../bwtransforms/sa-is-bwt.h"
#include "../bwtransforms/sais.hxx"

namespace tests {

int SufCmp(byte *str, uint32 s1, uint32 s2, unsigned n) {
  unsigned limit = n - ((s1 > s2)? s1 : s2);
  for(unsigned i = 0; i < limit; ++i) {
    if(str[s1 + i] < str[s2 + i]) return -1;
    else if (str[s1 + i] > str[s2 + i]) return 1;
  }
  if(s1 > s2) return -1;
  else return 1;
}

void testSA_IS() {
  unsigned size = (rand() & 0x000FFFFF) + 1;
  byte *str = new byte[size+1];
  for(unsigned i = 0; i < size + 1; ++i) {
    str[i] = rand() & 0xFF;
  }
  if(size > 255) {
    for(unsigned i = 0; i < 256; ++i)
      str[i] = i;
  }
  str[size] = 0;
  uint32 *SA = new uint32[size+1];
  sa_is::SA_IS_ZeroInclude(str, SA, size + 1, 256);
  for(unsigned i = 1; i < size; ++i) {
    assert(SufCmp(str, SA[i], SA[i+1], size) < 0);
    assert(SA[i] < size);
  }
  std::cout << "." << std::flush;
  assert(SA[0] == size);
}

void test_sais() {
  int size = (rand() & 0x000FFFFF) + 1;
  byte *str = new byte[size+1];
  for(unsigned i = 0; i < size; ++i) {
    str[i] = rand() & 0xFF;
  }
  str[size] = 0;
  int *SA = new int[size + 1];
  saisxx(str, SA, size + 1, 256);
  for(unsigned i = 1; i < size; ++i) {
    assert(SufCmp(str, SA[i], SA[i+1], size) < 0);
    assert(SA[i] < size);
  }
  std::cout << "." << std::flush;
  assert(SA[0] == size);
}

void SimpleTest(char *arg) {
  uint32 len = strlen(arg);
  byte *str = new byte[len+1];
  strcpy((char*)str, arg);
  uint32 *SA = new uint32[len+1];
  sa_is::SA_IS(str, SA, len + 1, 256);
  for(uint32 i = 0; i < len + 1; ++i) {
    if(SA[i] > 0) std::cout << str[SA[i] - 1];
    else std::cout << '$';
  }
  std::cout << "\n";
  delete [] SA;
  delete [] str;
}

void BwtSaisTest(char *arg) {
  int len = strlen(arg);
  byte *str = new byte[len+1];
  strcpy((char*)str, arg);
  int *SA = new int[len+1];
  byte *res = new byte[len+1];
  int val = saisxx_bwt(str, res, SA, len + 1, 256);
  for(int i = 0; i < len + 1; ++i) {
    if(i != val) std::cout << res[i];
    else std::cout << '$';
  }
  std::cout << "\nReturn value: " << val << "\n";
  delete [] SA;
  delete [] str;
}

}


int main(int argc, char** argv) {
  using namespace tests;
  if(argc > 1) {
    SimpleTest(argv[1]);
    BwtSaisTest(argv[1]);
  }
  srand( time(NULL));
  for(int i = 0; i < 10; ++i) {
    testSA_IS();
  }
  for(int i = 0; i < 10; ++i) {
    test_sais();
  }
  std::cout << "\nSA_IS passed all tests.\n";
}
