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
#include <algorithm>
#include <numeric>
#include <cassert>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <utility>
#include <vector>

namespace bwtc {
int verbosity = 3;

namespace tests {

using namespace long_sequences;
typedef hash_functions::MaskHasher Hasher;
/** The first vector has number of elements in bucket. The second has number
 *  of distinct strings in bucket. */
typedef std::pair<std::vector<uint32>*, std::vector<uint32>* > Stats;

/** Check that each entry in single bucket had same hash value before
 *  sorting the buckets with stringsort. */
void ValidateBucketSort(const std::vector<uint32>& bucket_limits,
                        const std::vector<chunk>& positions,
                        const std::vector<bucket_struct>& buckets)
{
  uint32 lo = 0, hi = NextValidPos(bucket_limits, 0);
  while(lo < positions.size()){
    uint32 index_in_bl = positions[GetPosition(&buckets[lo])].hash_value;
    for(uint32 i = lo; i < bucket_limits[hi]; ++i) {
      assert(index_in_bl == positions[GetPosition(&buckets[i])].hash_value);
    }
    lo = bucket_limits[hi];
    hi = NextValidPos(bucket_limits, hi+1);
  }
}

bool StringEq(uint32 f, uint32 s, byte *src, uint32 length) {
  for(uint32 i = 0; i < length; ++i)
    if(src[f+i] != src[s+i]) return false;
  return true;
}

/** Validate that strings are sorted correctly */
uint32 ValidateStrSort(const std::vector<bucket_struct>& buckets,
                       const std::vector<chunk>& pos_ordered, byte *from,
                       uint32 length)
{
  uint32 lo = 0, hi = 1;
  uint32 count = 0;;
  assert(StartOfBucket(&buckets[lo]));
  while(lo < buckets.size()) {
    ++count;
    assert(GetPosition(&buckets[lo]) < pos_ordered.size());
    while(hi < buckets.size() && !StartOfBucket(&buckets[hi])) ++hi;
    for(uint32 i = lo+1; i < hi; ++i) {
      assert(GetPosition(&buckets[i]) < pos_ordered.size());
      assert(StringEq(pos_ordered[GetPosition(&buckets[lo])].position,
                      pos_ordered[GetPosition(&buckets[i])].position,
                      from, length));
    }
    if(hi < buckets.size())
      assert(!StringEq(pos_ordered[GetPosition(&buckets[lo])].position,
                       pos_ordered[GetPosition(&buckets[hi])].position,
                       from, length));
    lo = hi++;
  }
  return count;
}

Stats BucketStatistics(const std::vector<bucket_struct>& buckets,
                       const std::vector<uint32>& hash_values)
{
  std::vector<uint32> *total = new std::vector<uint32>();
  std::vector<uint32> *distincts = new std::vector<uint32>();
  std::vector<uint32> &t = *total, &d = *distincts;
  t.push_back(0); d.push_back(0);
  uint32 curr_buck = 0, h_lo = 0, h_hi = 0;
  h_hi = NextValidPos(hash_values, 0);
  for(uint32 i = 0; i < hash_values[h_hi]; ++i) {
    if(StartOfBucket(&buckets[i])) ++d[curr_buck];
    ++t[curr_buck];
  }
  h_lo = h_hi++;
  h_hi = NextValidPos(hash_values, h_hi);
  while(h_hi < hash_values.size()) {
    t.push_back(0); d.push_back(0); ++curr_buck;
    for(uint32 i = hash_values[h_lo]; i < hash_values[h_hi]; ++i) {
      if(StartOfBucket(&buckets[i])) ++d[curr_buck];
      ++t[curr_buck];
    }
    h_lo = h_hi++;
    h_hi = NextValidPos(hash_values, h_hi);
  }
  return std::make_pair(total, distincts);
}

void TestValidity(byte *from, unsigned length, unsigned win_size) {
  uint32 freqs[256];
  std::vector<chunk> pos_ordered;
  std::vector<uint32> hash_counts;
  uint32 size_recommendation = length/(win_size*2);
  SequenceDetector<Hasher> seq_det(from, size_recommendation, freqs,
                                   &pos_ordered, &hash_counts, win_size);
  std::cout << "Size of hash table: " << hash_counts.size() << std::endl;
  seq_det.ScanAndStore(length);
  std::vector<bucket_struct> buckets;
  SortIntoBuckets(&hash_counts, &pos_ordered, &buckets);
  assert(buckets.size() == pos_ordered.size());
  assert(buckets.size() < (Max<uint32>::max & 0x7FFFFFFF));
  ValidateBucketSort(hash_counts, pos_ordered, buckets);

  SortBuckets(from, win_size, hash_counts, &buckets, pos_ordered);
  Stats st = BucketStatistics(buckets, hash_counts);
  assert(st.first->size() == st.second->size());
  uint32 real_buckets = ValidateStrSort(buckets, pos_ordered, from, win_size);
  std::cout << "Buckets found in scan-phase: " << st.second->size() << "\n";
  std::cout << "Real number of buckets: " << real_buckets << "\n";
  std::vector<uint32>::iterator max_el =
      std::max_element(st.second->begin(), st.second->end());
  std::cout << "Maximum number of collisions: " << *max_el << "\n";
  std::cout << "On average there were " <<
      ((1.0*std::accumulate(st.second->begin(), st.second->end(),0))/
       st.second->size()) << " distinct elements in single bucket\n";
  std::vector<uint32>::iterator max_c =  
      std::max_element(st.first->begin(), st.first->end());
  std::cout << "Largest bucket held " << *max_c << " elements, where were "
            << st.second->at(max_c - st.first->begin()) << " separate\n";
  delete st.first; delete st.second;
  //HandleBuckets(from, win_size, &buckets, &pos_ordered);
  
}

void TestScanTime(byte *from, unsigned length, unsigned win_size) {
  clock_t start = clock();
  uint32 freqs[256];
  std::vector<chunk> pos_ordered;
  std::vector<uint32> hash_counts;
  uint32 size_recommendation = length/(win_size*2);
  SequenceDetector<Hasher> seq_det(from, size_recommendation, freqs,
                                   &pos_ordered, &hash_counts, win_size);
  seq_det.ScanAndStore(length);
  clock_t end = clock();
  std::cout << "Time spent on scanning: "
            << (end-start)/static_cast<double>(CLOCKS_PER_SEC) << "\n";
  std::vector<bucket_struct> buckets;
  SortIntoBuckets(&hash_counts, &pos_ordered, &buckets);
  clock_t t = clock();
  std::cout << "Time spent on sorting into buckets: "
            << (t-end)/static_cast<double>(CLOCKS_PER_SEC) << "\n";
  end = t;
  SortBuckets(from, win_size, hash_counts, &buckets, pos_ordered);
  t = clock();
  std::cout << "Time spent on sorting the buckets: "
            << (t-end)/static_cast<double>(CLOCKS_PER_SEC) << "\n";
  end = t;
  std::cout << "Time spent on scanning and sorting buckets: "
            << (end-start)/static_cast<double>(CLOCKS_PER_SEC) << "\n";
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
  bwtc::tests::TestScanTime(&input[0], input.size(), win_size);
}
