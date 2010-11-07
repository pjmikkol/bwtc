/**
 * @file long_sequences_test.cc
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
 * Testing for detecting and replacing long repetitive sequences.
 */


#include "../preprocessors/long_sequences.h"
#include "../preprocessors/sequence_detector.h"
#include "../preprocessors/sequence_detector-inl.h"
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

namespace bwtc {
int verbosity = 3;

namespace tests {

using namespace long_sequences;
typedef hash_functions::MaskHasher Hasher;

void ValidateBucketSort(const std::vector<uint32>& bucket_limits,
                        const std::vector<chunk>& positions,
                        const std::vector<bucket_struct>& buckets)
{
  uint32 lo = 0, hi = long_sequences::NextValidPos(bucket_limits, 0);
  while(lo < positions.size()){
    uint32 index_in_bl = positions[GetPosition(&buckets[lo])].hash_value;
    for(int i = lo; i < bucket_limits[hi]; ++i) {
      assert(index_in_bl == positions[GetPosition(&buckets[i])].hash_value);
    }
    lo = bucket_limits[hi];
    hi = long_sequences::NextValidPos(bucket_limits, hi+1);
  }
}

void TestValidity(byte *from, unsigned length, unsigned win_size) {
  uint32 freqs[256];
  std::vector<chunk> pos_ordered;
  std::vector<uint32> hash_counts;
  uint32 size_recommendation = length/(2*win_size);
  SequenceDetector<Hasher> seq_det(from, size_recommendation, freqs,
                                   &pos_ordered, &hash_counts, win_size);
  seq_det.ScanAndStore(length);
  std::vector<bucket_struct> buckets;
  SortIntoBuckets(&hash_counts, &pos_ordered, &buckets);
  assert(buckets.size() == pos_ordered.size());
  assert(buckets.size() < (Max<uint32>::max & 0x7FFFFFFF));
  ValidateBucketSort(hash_counts, pos_ordered, buckets);
  /* TODO: Validation of sorting */
  SortBuckets(from, win_size, hash_counts, &buckets, pos_ordered);

  //HandleBuckets(from, win_size, &buckets, &pos_ordered);
  
}


}//namespace tests
}//namespace bwtc

using bwtc::byte;

void ReadFile(std::vector<byte>& to, const char *from) {
  std::ifstream infile(from);
  infile.seekg(0, std::ios::end);
  unsigned length = infile.tellg();
  infile.seekg(0, std::ios::beg);
  to.resize(length);
  infile.read((char*)&to[0], length);  
}

int main(int argc, char **argv) {
  unsigned win_size = 16;
  if(argc < 2) return 0;
  if(argc > 2) win_size = atoi(argv[2]);
  std::vector<byte> input;
  ReadFile(input, argv[1]);
  bwtc::tests::TestValidity(&input[0], input.size(), win_size);
}
