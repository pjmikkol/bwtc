#include <cstring>

#include <iostream>
#include <vector>

#include "../block.h"
#include "dcbwt.h"
#include "difference_cover.h"
#include "difference_cover-inl.h"
#include "prefix_doubling.h"
#include "prefix_doubling-inl.h"
#include "stringsort.h"
#include "stringsort-inl.h"

namespace bwtc {

using dcsbwt::DifferenceCover;
using bwtc::verbosity;

DCBWTransform::DCBWTransform() {}

DCBWTransform::~DCBWTransform() {}
  
void DCBWTransform::Connect(MainBlock* block) {}

void DCBWTransform::BuildStats() {}


uint64 DCBWTransform::MaxSizeInBytes(uint64 block_size) const {
  uint64 period = static_cast<uint64>(1) << log_period_;
  uint64 size_per_period = sizeof(byte) * period
      + sizeof(uint32) * (period
                          + DifferenceCover::CoverSizeForPeriod(period));
  return kMemoryOverhead
      + ((block_size * size_per_period) / period)
      + DifferenceCover::SizeInBytesForPeriod(period);
}

uint64 DCBWTransform::MaxBlockSize(uint64 memory_budget) const {
  uint64 period = static_cast<uint64>(1) << log_period_;
  uint64 size_per_period = sizeof(byte) * period
      + sizeof(uint32) * (period
                          + DifferenceCover::CoverSizeForPeriod(period));
  memory_budget -= (kMemoryOverhead
                    + DifferenceCover::SizeInBytesForPeriod(period));
  return (period * memory_budget) / size_per_period;
}

uint64 DCBWTransform::SuggestedBlockSize(uint64 memory_budget) const {
  return MaxBlockSize(memory_budget);
}

/* BWT with difference cover sampling */
namespace { 


template <typename Iterator>
class IteratorGreater : public std::binary_function<Iterator, Iterator, bool> {
 public:
  IteratorGreater() {}
  bool operator() (Iterator a, Iterator b) const {
    return *a > *b;
  }
};



//////////////////////////////////////////////////////////////////
// Optimize the order of the alphabet.
//////////////////////////////////////////////////////////////////

class EdgeGreater : std::binary_function<int, int, bool> {
 public:
  EdgeGreater(const uint32* weights) : weights_(weights) {}
  bool operator() (int a, int b) const {
    return weights_[a] > weights_[b];
  }
 private:
  const uint32* const weights_;
};

uint32 OptimizeAlphabetOrder(const uint32* bucket_size, uint8* char_of_rank) {

  class Graph {
   public:
    Graph(const uint32* weights) : kRemovedEdgeWeight(0xF0FFFFFF) {
      memcpy(weights_, weights, 0x10000 * sizeof(uint32));
    }
    int Edge(int i, int j) const {
      assert(i >= 0);
      assert(i < 256);
      assert(j >= 0);
      assert(i < 256);
      return (i << 8) + j;
    }
    int Source(int edge) const {
      assert(edge >= 0);
      assert(edge < 0x10000);
      return edge >> 8;
    }
    int Target(int edge) const {
      assert(edge >= 0);
      assert(edge < 0x10000);
      return edge & 255;
    }
    uint32 Weight(int edge) const { return weights_[edge]; }
    void ReduceWeight(int edge, uint32 weight) {
      assert(weight >= 0);
      assert(weight <= weights_[edge]);
      weights_[edge] -= weight;
    }
    void Remove(int edge) { weights_[edge] = kRemovedEdgeWeight; }
    bool Exists(int edge) const {
      return weights_[edge] != kRemovedEdgeWeight;
    }

    EdgeGreater GetEdgeGreater() const { return EdgeGreater(weights_); }
   private:
    uint32 weights_[0x10000];
    const uint32 kRemovedEdgeWeight;
  };

  Graph graph(bucket_size);

  uint32 total_weight = 0;
  for (int i = 0; i < 256; ++i) {
    for (int j = 0; j < 256; ++j) {
      int edge = graph.Edge(i, j);
      total_weight += graph.Weight(edge);
    }
  }

  uint32 commited_weight = 0;
  uint32 avoided_weight = 0;

  // Remove self-loops.
  uint32 one_cycle_weight = 0;
  for (int i = 0; i < 256; ++i) {
    int edge = graph.Edge(i, i);
    one_cycle_weight += graph.Weight(edge);
    graph.Remove(edge);
  }
  avoided_weight += one_cycle_weight;

  // '\0' is always the first character.
  // Remove edges adjacent to it.
  for (int i = 1; i < 256; ++i) {
    int edge = graph.Edge(0, i);
    commited_weight += graph.Weight(edge);
    graph.Remove(edge);
    edge = graph.Edge(i, 0);
    avoided_weight += graph.Weight(edge);
    graph.Remove(edge);
  }

  if (verbosity > 4) {
    std::clog << " total_weight=" << total_weight
              << " commited_weight=" << commited_weight
              << " avoided_weight=" << avoided_weight
              << " one_cycle_weight=" << one_cycle_weight
              << std::endl;
  }

  uint32 two_cycle_weight = 0;
  for (int i = 1; i < 255; ++i) {
    for (int j = i + 1; j < 256; ++j) {
      int edge = graph.Edge(i, j);
      int reverse_edge = graph.Edge(j, i);
      uint32 min_weight = std::min(graph.Weight(edge),
                                   graph.Weight(reverse_edge));
      two_cycle_weight += min_weight;
      graph.ReduceWeight(edge, min_weight);
      graph.ReduceWeight(reverse_edge, min_weight);
    }
  }
  commited_weight += two_cycle_weight;
  avoided_weight += two_cycle_weight;

  if (verbosity > 4) {
    std::clog << " commited_weight=" << commited_weight
              << " avoided_weight=" << avoided_weight
              << " two_cycle_weight=" << two_cycle_weight
              << std::endl;
  }

  std::vector<int> positive_edges;
  for (int i = 1; i < 256; ++i) {
    for (int j = 1; j < 256; ++j) {
      int edge = graph.Edge(i, j);
      if (graph.Exists(edge) && graph.Weight(edge) > 0)
        positive_edges.push_back(edge);
    }
  }
  std::sort(positive_edges.begin(), positive_edges.end(),
            graph.GetEdgeGreater());

  uint32 selected_weight = 0;
  uint32 discarded_weight = 0;

  // Dynamic transitive closure
  std::vector<uint8> predecessors[256];
  std::vector<uint8> successors[256];
  for (std::vector<int>::iterator it = positive_edges.begin();
       it != positive_edges.end(); ++it) {
    int discarded_edge = *it;
    if (!graph.Exists(discarded_edge)) continue;
    int source = graph.Target(discarded_edge);
    int target = graph.Source(discarded_edge);
    std::vector<uint8> new_target_predecessors;
    new_target_predecessors.push_back(source);
    std::vector<uint8> new_source_successors;
    new_source_successors.push_back(target);
    std::vector<uint8>::iterator it1, it2;
    for (it1 = predecessors[source].begin();
         it1 != predecessors[source].end(); ++it1) {
      int node = *it1;
      int edge = graph.Edge(node, target);
      if (graph.Exists(edge)) {
        new_target_predecessors.push_back(node);
      }
    }
    for (it2 = successors[target].begin();
         it2 != successors[target].end(); ++it2) {
      int node = *it2;
      int edge = graph.Edge(source, node);
      if (graph.Exists(edge)) {
        new_source_successors.push_back(node);
      }
    }
    for (it1 = new_target_predecessors.begin();
         it1 != new_target_predecessors.end(); ++it1) {
      int node1 = *it1;
      for (it2 = new_source_successors.begin();
           it2 != new_source_successors.end(); ++it2) {
        int node2 = *it2;
        int selected_edge = graph.Edge(node1, node2);
        if (graph.Exists(selected_edge)) {
          int discarded_edge = graph.Edge(node2, node1);
          selected_weight += graph.Weight(selected_edge);
          graph.Remove(selected_edge);
          discarded_weight += graph.Weight(discarded_edge);
          graph.Remove(discarded_edge);
          predecessors[node2].push_back(node1);
          successors[node1].push_back(node2);
        }
      }
    }
  }

  if (verbosity > 4) {
    std::clog << " selected_weight=" << commited_weight
              << " discarded_weight=" << avoided_weight
              << std::endl;
  }

  // Recover the node order.
  uint8 node_weight[256];
  std::vector<uint8*> node_ptr;
  for (int i = 1; i < 256; ++i) {
    node_weight[i] = successors[i].size();
    node_ptr.push_back(&node_weight[i]);
  }
  std::sort(node_ptr.begin(), node_ptr.end(),
            IteratorGreater<uint8*>());
  char_of_rank[0] = 0;
  for (int i = 1; i < 256; ++i) {
    char_of_rank[i] = node_ptr[i-1] - node_weight;
  }
  if (discarded_weight < selected_weight) {
    std::reverse(char_of_rank + 1, char_of_rank + 256);
    std::swap(discarded_weight, selected_weight);
  }
  commited_weight += selected_weight;
  avoided_weight += discarded_weight;

  if (verbosity > 4) {
    std::clog << " total_weight=" << total_weight
              << " commited_weight=" << commited_weight
              << " avoided_weight=" << avoided_weight;
  }
  assert(total_weight == commited_weight + avoided_weight);

  return commited_weight;
}











/*****************************************************************
 * Compute the suffix cover sample ranks and provide methods     *
 * for using it.                                                 *
 *****************************************************************/
class DifferenceCoverSuffixComparer {
 public:
  DifferenceCoverSuffixComparer(byte* block, uint32 block_size,
                                uint32 period,
                                uint32* suffix_area, uint32* suffix_area_end)
      : dcsampler_(period, block_size),
        dcranks_(dcsampler_.size() + 1) {
    uint32 num_samples = dcsampler_.size();
    assert(suffix_area_end - suffix_area >= 2 * num_samples + 1);
    uint32* sample_suffixes_in = suffix_area_end - num_samples;
    uint32* sample_suffixes_end = suffix_area + 1 + num_samples;

    uint32* end = dcsampler_.Fill(suffix_area);
    assert(end == sample_suffixes_end - 1);

    // Compute the sizes of the buckets.
    uint32 bucket_size[0x10000];
    memset(bucket_size, 0, 0x10000 * sizeof(uint32));
    // The suffix of length one needs special treatment if it is in the sample.
    // If it is not, this will remain 0x10000.
    int bucket_of_length_one_suffix = 0x10000;
    for (uint32 i = 0; i < num_samples; ++i) {
      uint32 suffix = suffix_area[i];
      assert(suffix < block_size);
      int bucket = static_cast<byte>(block[suffix]) << 8;
      if (suffix + 1 < block_size)
        bucket += static_cast<byte>(block[suffix + 1]);
      else
        // The suffix of length one is treated as if there was
        // an extra '\0' at the end of the block.
        bucket_of_length_one_suffix = bucket;
      ++bucket_size[bucket];
    }

    // Compute the bucket positions.
    uint32 bucket_begin[0x10000];
    uint32 sum = 0;
    for (int bucket = 0; bucket < 0x10000; ++bucket) {
      bucket_begin[bucket] = sum;
      sum += bucket_size[bucket];
    }
    assert(sum == num_samples);

    // Distribute the suffixes to the buckets.
    if (bucket_of_length_one_suffix != 0x10000) {
      // The suffix of length one is always the first in its bucket.
      sample_suffixes_in[bucket_begin[bucket_of_length_one_suffix]++]
          = block_size - 1;
    }
    for (uint32 i = 0; i < num_samples; ++i) {
      uint32 suffix = suffix_area[i];
      assert(suffix < block_size);
      if (suffix < block_size - 1) {
        int bucket = static_cast<byte>(block[suffix]) << 8;
        bucket += static_cast<byte>(block[suffix + 1]);
        sample_suffixes_in[bucket_begin[bucket]++] = suffix;
      }
    }

    uint32* dcranks = &dcranks_[0];
    uint32* dcranks_end = dcranks + num_samples + 1;
    // This is the empty suffix needed by prefix doubling.
    dcranks[num_samples] = 0;
    suffix_area[0] = num_samples;

    dcsbwt::PrefixDoubler<uint32> doubler(num_samples);
    RankRecorder rank_recorder(&dcsampler_, dcranks, suffix_area,
                               doubler.FinishedSuffixMarker());

    // Sort the buckets.
    end = suffix_area + 1;  // Skip the empty suffix.
    uint32 shift = sample_suffixes_in - end;
    for (int bucket = 0; bucket < 0x10000; ++bucket) {
      uint32* begin = end;
      end += bucket_size[bucket];
      if (bucket == bucket_of_length_one_suffix) {
        // The suffix of length one is already in its place
        // at the beginning of the bucket. Process it separately and
        // then sort the rest of the bucket. We cannot sort it with
        // others, because StringsortSuffixes assumes that all suffixes
        // have at least 2 characters.
        assert(*(begin + shift) == block_size - 1);
        *begin = block_size - 1;
        rank_recorder(begin, begin + 1);
        ++begin;
      }
      uint32* bucket_end =
          dcsbwt::StringsortSuffixes(
              block, block + block_size,
              begin, begin + shift, end + shift,
              2, period,
              rank_recorder);
      assert(end == bucket_end);
    }

    // Final sorting by prefix doubling.
    uint32 period_interval = dcsampler_.PeriodInterval();
    doubler.SortSuffixes(dcranks, dcranks_end,
                         suffix_area, sample_suffixes_end, period_interval);
  }

  // This functor can be used for comparing two suffixes in constant time
  // if they have a common prefix of length period - 1.
  class Less : public std::binary_function<uint32, uint32, bool> {
   public:
    Less(const dcsbwt::DifferenceCoverSample<uint32>* dcsample,
         const uint32* dcranks)
        : dcsample_(dcsample), dcranks_(dcranks) {}
    bool operator() (uint32 a, uint32 b) const {
      uint32 shift = dcsample_->Shift(a, b);
      uint32 apos = dcsample_->Rank(a + shift);
      uint32 bpos = dcsample_->Rank(b + shift);
      return dcranks_[apos] < dcranks_[bpos];
    }
   private:
    const dcsbwt::DifferenceCoverSample<uint32>* const dcsample_;
    const uint32* const dcranks_;
  };
  Less GetLess() const { return Less(&dcsampler_, &dcranks_[0]); }

  uint32 Size() const { return dcranks_.size(); }

 private:
  dcsbwt::DifferenceCoverSample<uint32> dcsampler_;
  std::vector<uint32> dcranks_;

  class RankRecorder {
   public:
    RankRecorder(const dcsbwt::DifferenceCoverSample<uint32>* dcsample,
                 uint32* rank_array, uint32* suffix_array,
                 uint32 finished_suffix_marker)
        : dcsample_(dcsample),
          rank_array_(rank_array), suffix_array_(suffix_array),
          finished_suffix_marker_(finished_suffix_marker) {}
    void operator() (uint32* begin, uint32* end) const {
      uint32 rank = end - suffix_array_ - 1;
      if (end - begin == 1) {
        rank_array_[dcsample_->Rank(*begin)] = rank;
        *begin = finished_suffix_marker_;
      } else {
        for (uint32* it = begin; it != end; ++it) {
          uint32 pos = dcsample_->Rank(*it);
          rank_array_[pos] = rank;
          *it = pos;
        }
      }
    }
   private:
    const dcsbwt::DifferenceCoverSample<uint32>* const dcsample_;
    uint32* const rank_array_;
    const uint32* const suffix_array_;
    uint32 finished_suffix_marker_;
  };
};











class DCSorter {
 public:
  DCSorter(const DifferenceCoverSuffixComparer::Less& dcless)
      : dcless_(dcless) {}
  void operator() (uint32* begin, uint32* end) const {
    std::sort(begin, end, dcless_);
  }
 private:
  DifferenceCoverSuffixComparer::Less dcless_;
};

} //empty namespace








byte* DCBWTransform::DoTransform(uint64* eob_byte) {
  if (current_block_) return NULL;

  uint64 block_size = current_block_->Size();
  byte* block = current_block_->begin();
  
  assert(block_size > 0);
  assert(block_size < 0xFFFFFFFF - period_);  // TODO: Check this
  uint32 num_suffixes = block_size + 1;

  byte* result = AllocateMemory(block_size);

  // Use basic string sorting for small blocks.
  if (block_size < 1024) {
    if (verbosity > 2) {
      std::clog << "DifferenceCoverBWTransformer::DoTransform"
                << " block_size=" << block_size
                << " (using plain string sorting)"
                << std::endl;
    }
    std::vector<uint32> suffix_array(2 * (num_suffixes));
    uint32* begin = &suffix_array[0];
    uint32* end = begin + suffix_array.size();
    uint32* suffixes_in = end - num_suffixes;
    for (int i = 0; i < num_suffixes; ++i) suffixes_in[i] = i;
    dcsbwt::NullFinishedGroupReporter<uint32> null_reporter;
    uint32* suffixes_end =
        dcsbwt::StringsortSuffixes(
            block, block + block_size,
            begin, suffixes_in, end,
            0, block_size,
            null_reporter);
    assert(suffixes_end == begin + num_suffixes);

    *eob_byte =  BwtFromSuffixArray(block, block_size, begin, result);
    return result;
  }
  // For short blocks and long periods, we may have to reduce the period.
  // The period is restored at the end.
  uint32 saved_period = period_;
  while (block_size < period_ - 1) period_ >>= 1;

  if (verbosity > 2) {
    std::clog << "DifferenceCoverBWTransformer::DoTransform"
              << " block_size=" << block_size
              << " period=" << period_
              << std::endl;
  }
  assert(block_size >= period_ - 1);

  // The extra +4 might be needed by DifferenceCoverSuffixComparer
  // if period_==8.
  std::vector<uint32> suffix_array(num_suffixes + 4);
  uint32* suffix_area = &suffix_array[0];
  uint32* suffix_area_end = suffix_area + suffix_array.size();

  // Construct the difference cover sample.
  DifferenceCoverSuffixComparer dcscomp(block, block_size, period_,
                                        suffix_area, suffix_area_end);
  if (verbosity > 3) {
    std::clog << "Using difference cover sample of size " << dcscomp.Size()
              << std::endl;
  }

  // Compute the sizes of the buckets.
  if (verbosity > 4) {
    std::clog << "Computing bucket sizes" << std::endl;
  }
  uint32 bucket_size[0x10000];
  memset(bucket_size, 0, 0x10000*sizeof(uint32));
  int bucket = static_cast<byte>(block[0]) << 8;
  for (uint32 i = 1; i < block_size; ++i) {
    bucket += static_cast<byte>(block[i]);
    ++bucket_size[bucket];
    bucket = (bucket << 8) & 0xFF00;
  }
  ++bucket_size[bucket];  // suffix of length one
  ++bucket_size[0];       // empty suffix

  if (verbosity > 4) {
    std::clog << "Optimizing alphabet order" << std::endl;
  }
  byte char_of_rank[256];
  uint32 num_suffixes_to_sort
      = OptimizeAlphabetOrder(bucket_size, char_of_rank);
  uint8 rank_of_char[256];
  for (int i = 0; i < 256; ++i) {
    rank_of_char[char_of_rank[i]] = i;
  }

  // Compute the positions of the selected buckets.
  if (verbosity > 4) {
    std::clog << "Computing selected bucket positions" << std::endl;
  }
  uint32 selected_bucket_begin[0x10000];
  uint32 num_selected_suffixes = 0;
  for (int ch1 = 0; ch1 < 256; ++ch1) {
    for (int ch2 = 0; ch2 < 256; ++ch2) {
      if (rank_of_char[ch1] < rank_of_char[ch2]) {
        int bucket = (ch1 << 8) + ch2;
        selected_bucket_begin[bucket] = num_selected_suffixes;
        num_selected_suffixes += bucket_size[bucket];
      }
    }
  }
  assert(num_selected_suffixes == num_suffixes_to_sort);

  uint32* suffixes_end = suffix_area + num_suffixes;
  uint32* selected_suffixes = suffix_area_end - num_selected_suffixes;

  // Distribute the suffixes to the selected buckets.
  if (verbosity > 4) {
    std::clog << "Distributing into selected buckets" << std::endl;
  }
  int ch1 = static_cast<byte>(block[0]);
  for (uint32 i = 1; i < block_size; ++i) {
    int ch2 = static_cast<byte>(block[i]);
    if (rank_of_char[ch1] < rank_of_char[ch2]) {
      bucket = (ch1 << 8) + ch2;
      selected_suffixes[selected_bucket_begin[bucket]++] = i - 1;
    }
    ch1 = ch2;
  }

  // Compute the positions of all buckets.
  uint32 bucket_begin[0x10001];
  bucket_begin[0] = 0;
  uint32 sum = bucket_size[0];
  for (int bucket = 1; bucket < 0x10000; ++bucket) {
    bucket_begin[bucket] = sum;
    sum += bucket_size[bucket];
  }
  assert(sum == num_suffixes);
  bucket_size[0x10000] = sum;

  // Sort the selected buckets.
  if (verbosity > 3) {
    std::clog << "Stringsorting " << num_suffixes_to_sort
              << " out of " << num_suffixes << " suffixes"
              << std::endl;
  }
  DCSorter dcsorter(dcscomp.GetLess());
  uint32* source_end = selected_suffixes;
  for (int ch1 = 0; ch1 < 256; ++ch1) {
    for (int ch2 = 0; ch2 < 256; ++ch2) {
      if (rank_of_char[ch1] < rank_of_char[ch2]) {
        int bucket = (ch1 << 8) + ch2;
        uint32* source_begin = source_end;
        source_end += bucket_size[bucket];
        uint32* target_begin = suffix_area + bucket_begin[bucket];
        assert(target_begin <= source_begin);
        uint32* target_end =
            dcsbwt::StringsortSuffixes(
                block, block + block_size,
                target_begin, source_begin, source_end,
                2, period_ - 1,
                dcsorter);
        assert(target_end == target_begin + bucket_size[bucket]);
      }
    }
  }

  // Move other suffixes into their final positions.
  if (verbosity > 4) {
    std::clog << "Setting other suffixes" << std::endl;
  }

  suffix_area[bucket_begin[0]] = block_size;
  for (int i = 0; i < 256; ++i) {
    int chi = char_of_rank[i];
    uint32 sub_bucket_begin[256];
    uint32 sub_bucket_end[256];
    for (int ch = 0; ch < 256; ++ch) {
      int bucket = (ch << 8) + chi;
      sub_bucket_begin[ch] = bucket_begin[bucket];
      sub_bucket_end[ch] = bucket_begin[bucket] + bucket_size[bucket];
    }
    if (chi == 0) ++sub_bucket_begin[0];
    uint32 rank = bucket_begin[chi << 8];
    while (rank < sub_bucket_begin[chi]) {
      uint32 suffix = suffix_area[rank++];
      if (suffix > 0) {
        int ch = static_cast<byte>(block[suffix - 1]);
        if (rank_of_char[ch] >= rank_of_char[chi])
          suffix_area[sub_bucket_begin[ch]++] = suffix - 1;
      }
    }
    rank = bucket_begin[(chi + 1) << 8];
    while (rank > sub_bucket_end[chi]) {
      uint32 suffix = suffix_area[--rank];
      if (suffix > 0) {
        int ch = static_cast<byte>(block[suffix - 1]);
        if (rank_of_char[ch] >= rank_of_char[chi])
          suffix_area[--sub_bucket_end[ch]] = suffix - 1;
      }
    }
  }

  period_ = saved_period;
  *eob_byte =  BwtFromSuffixArray(block, block_size, suffix_area, result);
  return result;
}







#if 0

class DifferenceCoverBWTransformer : public BWTransformer {


 public:
  class Options : public BWTransformer::AlgorithmSpecificOptions {
   public:
    Options() { Set(std::string("6")); }
    virtual ~Options() {}

    virtual bool Set(const std::string& options) {
      if (0 == options.length()) return true;
      const char* begin = options.c_str();
      const char* end = begin + options.length();
      char* next;
      int64 log_period = strtoll(begin, &next, 0);
      if (next != end) return false;
      if (log_period < kMinLogPeriod || log_period > kMaxLogPeriod)
        return false;
      if (verbosity > 4) {
        std::clog << "Setting DifferenceCoverBWTransform::Options"
                  << " period=" << (1 << log_period) << std::endl;
      }
      log_period_ = log_period;
      return true;
    }
    virtual std::string Get() const {
      assert(log_period_ >= 0);
      assert(log_period_ <= kMaxLogPeriod);
      char buffer[4];
      sprintf(buffer, "%d", log_period_);
      std::string options(buffer);
      return options;
    }

    virtual int64 MaxSizeInBytes(int64 block_size) const {
      uint64 period = 1LL << log_period_;
      uint64 size_per_period = sizeof(char) * period
          + sizeof(uint32) * (period
                              + DifferenceCover::CoverSizeForPeriod(period));
      return kMemoryOverhead
          + ((block_size * size_per_period) / period)
          + DifferenceCover::SizeInBytesForPeriod(period);
    }
    virtual int64 MaxBlockSize(int64 memory_budget) const {
      uint64 period = 1LL << log_period_;
      uint64 size_per_period = sizeof(char) * period
          + sizeof(uint32) * (period
                              + DifferenceCover::CoverSizeForPeriod(period));
      memory_budget -= (kMemoryOverhead
                        + DifferenceCover::SizeInBytesForPeriod(period));
      return (period * memory_budget) / size_per_period;
    }
    virtual int64 SuggestedBlockSize(int64 memory_budget) const {
      return MaxBlockSize(memory_budget);
    }

    virtual BWTransformer* GetTransformer(int64 memory_budget) {
      return new DifferenceCoverBWTransformer(1 << log_period_);
    }
   private:
    static const int64 kMemoryOverhead = (1 << 20);
    static const int kMinLogPeriod = 3;
    static const int kMaxLogPeriod = 15;
    int log_period_;
  };

  DifferenceCoverBWTransformer(uint32 period) : period_(period) {}
  virtual ~DifferenceCoverBWTransformer() {}

 private:
  uint32 period_;

  virtual int64 DoTransform(char* block, size_t block_size,
                            OutStream* output);

  class DCSorter {
   public:
    DCSorter(const DifferenceCoverSuffixComparer::Less& dcless)
        : dcless_(dcless) {}
    void operator() (uint32* begin, uint32* end) const {
      std::sort(begin, end, dcless_);
    }
   private:
    DifferenceCoverSuffixComparer::Less dcless_;
  };

  DifferenceCoverBWTransformer(const DifferenceCoverBWTransformer&);
  DifferenceCoverBWTransformer& operator=(const DifferenceCoverBWTransformer&);
};
#endif

} //namespace bwtc
