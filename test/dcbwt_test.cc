#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include "../block.h"
#include "../bwtransforms/bw_transform.h"
#include "../bwtransforms/dcbwt.h"

namespace bwtc {
int verbosity = 7;
}

namespace tests {

void SimpleBWTtest(const char* input, uint64 length) {
  bwtc::BWTransform* bwt = bwtc::GiveTransformer();
  std::vector<byte> data(&input[0], &input[length]);
  std::vector<uint64> dummy;
  bwtc::MainBlock block(&data, &dummy, length);
  bwt->Connect(&block);
  std::reverse(data.begin(), data.end());
  uint64 eob;
  std::vector<byte>* result = bwt->DoTransform(&eob);
  for(unsigned i = 0; i < length + 1; ++i) {
    if(i == eob) std::cout << '$';
    else std::cout << (*result)[i];
  }
  std::cout << "\n" << "eob-byte position : " << eob << "\n";
  delete bwt;
  delete result;
}

void BWTfromFile(const char *fname) {
  bwtc::BWTransform* bwt = bwtc::GiveTransformer();
  std::ifstream file(fname);
  std::vector<byte> data((std::istream_iterator<byte>(file)),
                         std::istream_iterator<byte>());
  std::vector<uint64> dummy;
  uint64 length = data.size();
  bwtc::MainBlock block(&data, &dummy, length);
  bwt->Connect(&block);
  std::reverse(data.begin(), data.end());
  uint64 eob;
  std::vector<byte>* result = bwt->DoTransform(&eob);
  for(unsigned i = 0; i < length + 1; ++i) {
    if(i == eob) std::cout << '$';
    else std::cout << (*result)[i];
  }
  delete bwt;
  delete result;
}


} //namespace tests

int main(int argc, char **argv) {
  std::string tst("mississippi");
  tests::SimpleBWTtest(tst.c_str(), tst.size());
  if (argc > 1) tests::BWTfromFile(argv[1]);
}
