#include <cassert>

#include <algorithm>
#include <iostream> /* for std::streamsize*/
#include <map>
#include <string>
#include <utility> /* for pair */

#include "../block.h"
#include "../block_manager.h"
#include "../globaldefs.h"
#include "../stream.h"
#include "preprocessor.h"

namespace bwtc {

PreProcessor* GivePreProcessor(char choice, uint64 block_size,
                               const std::string& input)
{
  PreProcessor* pp;
  /* Expand this to conform for different PreProcessing algorithms */
  switch (choice) {
    case 'n':
    default:
      pp = new PreProcessor(block_size);
  }
  pp->Connect(input);
  return pp;
}

PreProcessor::PreProcessor(uint64 block_size) :
    source_(NULL), block_size_(block_size), block_manager_(NULL) { }

PreProcessor::~PreProcessor() {
  delete source_;
}

void PreProcessor::BuildStats(std::vector<byte>* data,
                              std::vector<uint64>* stats, uint64 data_size) {
  std::fill(stats->begin(), stats->end(), 0);
  //TODO: at the moment only contexts of length 1 are supported
  for( uint64 i = 0; i < data_size; ++i)
    (*stats)[(*data)[i]]++; 
}

void PreProcessor::Connect(std::string source_name) {
  source_ = new InStream(source_name);
}

void PreProcessor::AddBlockManager(BlockManager* manager) {
  block_manager_ = manager;
}

MainBlock* PreProcessor::ReadBlock() {
  assert(source_);
  assert(block_manager_);
  std::vector<byte>* to = block_manager_->GetFreeBuffer();
  std::vector<uint64>* stats = block_manager_->GetFreeStats();
  /* TODO:
   * streamsize type has as many bits as long. Since the preprocessor gets
   * blocksize as an uint64 we may end up in problems if user is on 32-bit
   * platform and wants to use too large blocks (amount of bits needed to
   * represent block size is more than 32)?? */
  /*** Stub implementation ***/
  std::streamsize read = source_->ReadBlock(
      &(*to)[0], static_cast<std::streamsize>(block_size_) );
  if (!read) return NULL;
  return block_manager_->MakeBlock(to, stats, static_cast<uint64>(read));
}


/*#################### Preprocessing algorithms #############################*/

/* Common utility-functions for preprocessing algorithms */

/* In C++0x these two could be implemented more naturally with the use of
 * lambda-functions. */
template <typename F, typename S>
bool ComparePairSecondAsc(std::pair<F,S> p1, std::pair<F,S> p2) {
  return (p1.second < p2.second);
}

template <typename F, typename S>
bool ComparePairSecondDesc(std::pair<F,S> p1, std::pair<F,S> p2)
{
  return (p1.second > p2.second);
}

template<typename Key, typename Value>
void InitPairsWithValue(std::pair<Key, Value> *pairs, Value v, uint64 length)
{
  for(uint64 i = 0; i < length; ++i)
    pairs[i] = std::make_pair(static_cast<Key>(i), v);
}

#define FREQUPDATE 1

class FreqTable {
 public:
  FreqTable() {
    for(int i = 0; i < 256; ++i) {
      freq_[i] = std::make_pair(static_cast<byte>(i), 0U);
    }
    InitLocations();
    assert(Test());
  }

  FreqTable(uint64* frequencies) {
    /* Assumes that frequencies has length of 256 */
    for(int i = 0; i < 256; ++i) {
      freq_[i] = std::make_pair(static_cast<byte>(i), frequencies[i]);
    }
    std::sort(freq_, freq_ + 256, ComparePairSecondAsc<byte, uint64>);
    InitLocations();
  }

  const uint64& operator[](unsigned i) {
    assert(i <= 255);
    return freq_[i].second;
  }

  byte Key(unsigned i) {
    assert(i < 256);
    return freq_[i].first;
  }

  //  void Decrease(unsigned freq_index, uint64 value) {
  bool Decrease(unsigned key, uint64 value) {
    //assert(freq_index <= 255);
    assert(key <= 255);
    unsigned freq_index = location_[key];
    //assert(freq_[freq_index].second >= value);
    if(freq_[freq_index].second < value) return false;
    uint64 new_value = freq_[freq_index].second - value;
    std::pair<byte, uint64> new_pair =
        std::make_pair(freq_[freq_index].first, new_value);

    while (freq_index > 0 && new_value < freq_[freq_index - 1].second)
    {
      ++location_[freq_[freq_index - 1].first];
      freq_[freq_index] = freq_[freq_index - 1];
      --freq_index;
    }
    freq_[freq_index] = new_pair;
    location_[new_pair.first]= freq_index;
    assert(freq_[location_[new_pair.first]].first == new_pair.first);
    return true;
  }

  void Increase(unsigned key, uint64 value) {
    //    assert(freq_index <= 255);
    assert(key <= 255);
    unsigned freq_index = location_[key];

    uint64 new_value = freq_[freq_index].second + value;
    std::pair<byte, uint64> new_pair = 
        std::make_pair(freq_[freq_index].first, new_value);

    while (freq_index < 255 && new_value > freq_[freq_index + 1].second)
    {
      --location_[freq_[freq_index + 1].first];
      freq_[freq_index] = freq_[freq_index + 1];
      ++freq_index;
    }
    freq_[freq_index] = new_pair;
    location_[new_pair.first]= freq_index;
    assert(freq_[location_[new_pair.first]].first == new_pair.first);
 }

 private:
  std::pair<byte, uint64> freq_[256];
  byte location_[256];

  void InitLocations() {
    for(unsigned i = 0; i < 256; ++i) {
      location_[freq_[i].first] = i;
    }
  }

  bool Test() {
    for(int i = 0; i < 256; ++i) {
      assert(freq_[location_[i]].first == i );
    }
    return true;
  }
};

/*##################### Replacing the most common pairs ######################*/
  /**************************************************************************
   *    Reasoning behind choosing the replaceable pairs                     *
   *------------------------------------------------------------------------*
   * For each pair we replace with some symbol we have to write that symbol *
   * and the correspongind pair to the header. We also make some symbols    *
   * 'free' from the original data, which means that we escape these symbols*
   * with special escape byte. Then we can use these freed symbols in pair  *
   * replacements.                                                          *
   *                                                                        *
   * Denote the frequency of symbol x_i with f(x_i) and frequency of pair   *
   * P_i with f(P_i). Writing the replacement info for single replacement   *
   * takes 3 bytes. So if P_i is going to be replaced with x_i we require   *
   * that                                                                   *
   *             f(x_i) + 3 < f(P_i)                              (p1)      *
   * On the left side is amount of additional bytes after replacement and on*
   * the right side is amount of bytes saved when writing pairs.            *
   *                                                                        *
   * If we free symbols the effect of making one byte to escape byte has to *
   * be also notified. Let x_i,...,x_j be the freed symbols, P_i,...,P_j    *
   * the pairs replaced with x_k's and x be the escape byte.                *
   * Total improvement from replacements is                                 *
   *      sum from k=i to j: f(P_k) - f(x_k) - 3                            *
   * The penalty from escape byte is f(x) so we require that                *
   *   sum from k=i to j: f(P_k) - f(x_k) - 3 > f(x)              (p2)      *
   **************************************************************************/
namespace commonpairs {

void ComputePairFrequencies(byte *data, uint64* freqs,
                            std::pair<uint16, uint32>* pair_freqs, uint64 len)
{
  uint16 index = data[0];
  ++freqs[index];
  for (uint64 i = 1; i < len; ++i) {
    ++freqs[data[i]];
    index <<= 8;
    index |= data[i];
    ++pair_freqs[index].second;
  }
}

/*****************************************************************************
 * FindReplaceablePairs - Finds the candidates for pairs which are going to  *
 *                        be replaced by single symbols. It is required that *
 *                        if pair = <p1,p2> is in result, then there is no   *
 *                        pair where p2 would be first symbol of pair or     *
 *                        where p1 would be second symbol. Optimal solution  *
 *                        for the problem is NP-hard (max-cut).              *
 *                                                                           *
 * Takes the following arguments:                                            *
 * replaceable_pairs  - empty vector where pairs are stored                  *
 * pair_freqs         - unsorted array of <pair, frequency>-pairs, size of   *
 *                      65536                                                *
 * frequencies        - array of <byte, frequency>-pairs sorted by frequency *
 *                      in descending order                                  *
 *****************************************************************************/
void FindReplaceablePairs(std::vector<std::pair<uint16,uint32> >*
                          replaceable_pairs,
                          std::pair<uint16, uint32>* pair_freqs,
                          FreqTable *freqs)
{
  const unsigned kStep = 256;
  unsigned current_pair = 0, current_symbol = 0;
  unsigned limit = 0;

  while(replaceable_pairs->size() < 254) {
    if(current_pair == limit) {
      limit += kStep;
      assert(limit < 65536);
      std::partial_sort(pair_freqs + current_pair, pair_freqs + limit,
                        pair_freqs  + 65536,
                        ComparePairSecondDesc<uint16, uint32>);
    }
    byte fst = static_cast<byte>((pair_freqs[current_pair].first & 0xFF00) >> 8);
    byte snd = static_cast<byte>(pair_freqs[current_pair].first & 0x00FF);
#ifdef FREQUPDATE
    if (!freqs->Decrease(fst, pair_freqs[current_pair].second)) {
      ++current_pair;
      continue;
    }
    if(fst != snd) freqs->Decrease(snd, pair_freqs[current_pair].second);
#endif
    // TODO: make the control flow more reasonable 
    
    /* Condition (p1) */
    if((*freqs)[current_symbol] + 3 >= pair_freqs[current_pair].second) {
#ifdef FREQUPDATE
      freqs->Increase(fst, pair_freqs[current_pair].second);
      freqs->Increase(snd, pair_freqs[current_pair].second);
#endif      
      break; /* We won't benefit from any changes any more*/
    }
    /* Reject pairs which have conflicting symbols
     * This one NEEDS better heuristics. Now it uses greedy heuristic */
    if (fst != snd) { /* Do not approve pairs of the same char */
      bool valid = true;
      for(std::vector<std::pair<uint16, uint32> >::iterator it =
              replaceable_pairs->begin();
          it != replaceable_pairs->end(); ++it)
      {
        uint16 current_fst = static_cast<byte>((it->first & 0xFF00) >> 8);
        uint16 current_sec = static_cast<byte>((it->first & 0x00FF));
        if (current_fst == snd || current_sec == fst) {
          valid = false;
#ifdef FREQUPDATE
          freqs->Increase(fst, pair_freqs[current_pair].second);
          freqs->Increase(snd, pair_freqs[current_pair].second);
#endif
          break;
        }
      }
      if(valid) {
        replaceable_pairs->push_back(pair_freqs[current_pair]);
#ifdef LAZYFREQ
        freqs->Decrease(fst, pair_freqs[current_pair].second);
        freqs->Decrease(snd, pair_freqs[current_pair].second);
#endif
        ++current_symbol;
      }
    }
    else {
#ifdef FREQUPDATE
      freqs->Increase(fst, pair_freqs[current_pair].second);
#endif
    }
    ++current_pair;
  }
  assert(current_symbol == replaceable_pairs->size());
}

/* Returns the index for the 'escape_char' in freqs-array or value of
 * free_symbols if freeing of the chars is not profitable */
unsigned EscapeCharIndex(FreqTable* freqs,
                         const std::vector<std::pair<uint16, uint32> >&
                         suitable_pairs,
                         unsigned free_symbols)
{
  if (suitable_pairs.size() <= free_symbols) return free_symbols;
  int64 utility = 0; 
  unsigned i; 
  for(i = free_symbols; i < suitable_pairs.size(); ++i) {
    utility += (suitable_pairs[i].second - (*freqs)[i] - 3);
  }
  /* Condition (p2) */
  while(utility <= static_cast<int64>((*freqs)[i]) &&
        i > free_symbols)
  {
    --i;
    byte fst = static_cast<byte>((suitable_pairs[i].first & 0xFF00) >> 8);
    byte snd = static_cast<byte>(suitable_pairs[i].first & 0x00FF);
#ifdef FREQUPDATE
    freqs->Increase(fst, suitable_pairs[i].second);
    freqs->Increase(snd, suitable_pairs[i].second);
#endif
    utility -= (suitable_pairs[i].second - (*freqs)[i] - 3);
  }
  return i;
}

/* Writes uint16 to given address. Used for writing the pairs to result-array*/
void WriteBytes(uint16 value, byte *address) {
  *address = (value >> 8);
  *(address + 1) = value & 0xFF;
}

/*******************************************************************************
 * WriteReplacements - Writes replacements of pairs to result array.           *
 *                                                                             *
 * replacements - array of size 65536 where index of element is interpreted    *
 *                as a pair of two bytes where their bit-representations are   *
 *                concatenated. Let p be the value of pair. We require the     *
 *                following conditions:                                        *
 *                if replacements[p] == common_byte there is no replacement    *
 *                                      for p                                  *
 *                if replacements[p] == escape_byte the first element of p     *
 *                                      needs escaping because it is made free *
 *                in other cases replacement for p is replacements[p]          *
 * to           - is an array where we write the replacements                  *
 * from         - source of an original data                                   *
 * length       - length of from                                               *
 * common_byte and escape_byte are described above                             *
 *******************************************************************************/
uint64 WriteReplacements(byte *replacements, byte *to, byte *from, uint64 length,
                         byte common_byte, byte escape_byte)
{
  unsigned escapes = 0, pairs = 0;

  uint64 result_index = 0;
  uint16 pair = from[0];
  uint64 i = 1;
  while(1) {
    pair <<= 8;
    pair |= from[i];
    if(replacements[pair] == common_byte) {
      to[result_index++] = from[i-1];
    }
    else if (replacements[pair] == escape_byte) {
      to[result_index++] = escape_byte;
      to[result_index++] = from[i-1];
      escapes++;
    }
    else { // pair will be replaced
      to[result_index++] = replacements[pair];
      pairs++;
      if( i == length - 1) break;
      pair = static_cast<byte>(from[++i]);
    }

    if ( i >= length - 1) {
      assert(i == length - 1);
      if(from[i] == escape_byte && escape_byte != common_byte) {
        to[result_index++] = escape_byte; escapes++;}
      to[result_index++] = from[i];
      break;
    }
    ++i;
  }
  return result_index;
}

} //namespace commonpairs

/************************************************************************
 * CompressCommonPairs - Replaces common pairs of bytes with single     *
 *                       byte values. Writes the result to source array *
 *                       (from). Requires that from has at least length *
 *                       of length+3, where the actual data is in range *
 *                       [0, length).                                   *
 ************************************************************************/
 
uint64 CompressCommonPairs(byte *from, uint64 length)
{
  using namespace commonpairs;

  assert(length > 0);
  assert(from);
  uint64 freq[256] = {0U};
  //std::pair<byte, uint64> freq[256];
  std::pair<uint16, uint32> pair_freq[65536];
  //  InitPairsWithValue<byte, uint64>(freq, 0, 256);
  InitPairsWithValue<uint16, uint32>(pair_freq, 0, 65536);
  ComputePairFrequencies(from, freq, pair_freq, length);

  //std::sort(freq, freq + 256);
  //std::sort(&freq[0], &freq[0] + 256, ComparePairSecondAsc<byte, uint64>);

  FreqTable freqs(freq);
  unsigned free_symbols = 0;
  while(freqs[free_symbols] == 0) ++free_symbols;
  //while(freq[free_symbols].second == 0) ++free_symbols;

  std::vector<std::pair<uint16, uint32> > replaceable_pairs;
  FindReplaceablePairs(&replaceable_pairs, pair_freq, &freqs);

  unsigned escape_index = free_symbols;
  if(replaceable_pairs.size() > free_symbols) {
    escape_index = EscapeCharIndex(&freqs, replaceable_pairs, free_symbols);
  }

  byte common_byte = (replaceable_pairs.size() > 0) ?
      static_cast<byte>(replaceable_pairs[0].first >> 8) : 0;
  byte escape_byte;
  if(escape_index > free_symbols) escape_byte = freqs.Key(escape_index);
  else escape_byte = common_byte;

  byte replacements[65536];
  std::fill(replacements, replacements + 65536, common_byte);
  byte *temp = new byte[length + 3];

  /* Initialize replacement table and write header */
  int symbols_in_use = 0;
  uint64 position = 0;
  unsigned candidates = static_cast<unsigned>(replaceable_pairs.size());
  unsigned k;
  for(k = 0; k < std::min(free_symbols, candidates); ++k) {
    replacements[replaceable_pairs[k].first] = freqs.Key(k);
    temp[position++] = freqs.Key(k);
    WriteBytes(replaceable_pairs[k].first, temp + position);
    position += 2;
  }
  if (k > 0) symbols_in_use = k;
  else symbols_in_use = 0;

  if ( free_symbols < escape_index) {
    for(unsigned i = free_symbols; i < escape_index; ++i) {
      replacements[replaceable_pairs[i].first] = freqs.Key(i);
      temp[position++] = freqs.Key(i);
      WriteBytes(replaceable_pairs[i].first, temp + position);
      position += 2;
      uint16 pair_value = freqs.Key(i) << 8;
      for(unsigned j = 0; j < 256; ++j, ++pair_value) {
        if (replacements[pair_value] == common_byte)
          replacements[pair_value] = escape_byte;
      }
    }
    symbols_in_use += escape_index - free_symbols + 1;
    uint16 pair_value = escape_byte << 8;
    for(unsigned j = 0; j < 256; ++j, ++pair_value) {
        if (replacements[pair_value] == common_byte)
          replacements[pair_value] = escape_byte;
    }
  }
  unsigned new_symbols = (free_symbols == escape_index)?
      0 : symbols_in_use - free_symbols;
  /* Find the dummy byte for header and write the end of header */
  byte dummy = escape_byte + 1;
  if (new_symbols > 0) dummy = freqs.Key(escape_index - 1);
  else if (symbols_in_use > 0) dummy = freqs.Key(symbols_in_use - 1);
  temp[position++] = dummy;
  if (free_symbols < escape_index) temp[position++] = escape_byte;
  else {
    temp[position++] = dummy;
    if( symbols_in_use == 0) temp[position++] = dummy;
  }

  if (verbosity > 1) {
    std::clog << "Replacing " << ((symbols_in_use)?(symbols_in_use - 1):0)
              << " pairs. ";
    if (new_symbols > 0)
      std::clog << "Made " << new_symbols << " symbols free.\n";
    else
      std::clog << "No symbols made free.\n";
  }

  uint64 total_size = position;
  total_size += WriteReplacements(replacements, temp + position, from, length,
                                  common_byte, escape_byte);
  assert(total_size <= length + 3);
  std::copy(temp, temp + total_size, from);
  delete [] temp;
  return total_size;
}




/*##################### Replacing the runs of same byte ######################*/
namespace longruns {

const unsigned kMaxLenOfSeq = 1 << 15;

/* Floor of logarithm of base two */
unsigned LogFloor(unsigned n) {
  assert(n > 0);
  unsigned log = 0;
  while(n > 1) {
    n >>= 1;
    ++log;
  }
  return log;
}

unsigned MostSignificantBit(unsigned n) {
  assert(n < (1 << 16));
  n |= (n >> 1);
  n |= (n >> 2);
  n |= (n >> 4);
  n |= (n >> 8);
  return n & (~(n >> 1));
}

void UpdateFreqs(std::map<unsigned, uint32> *run_freq, byte symbol, unsigned length)
{
  assert(length <= kMaxLenOfSeq);
  assert(length > 1);
  length -= (length % 2);
  unsigned original = length;
  while(length) {
    unsigned longest = MostSignificantBit(length);
    if (run_freq[symbol].count(longest))
      run_freq[symbol][longest] += original/longest;
    else run_freq[symbol][longest] = original/longest;
    length -= longest;
  }
}

void ComputeRunFrequencies(byte *from, std::pair<byte, uint64> *freq,
                           std::map<unsigned, uint32> *run_freq, uint64 length)
{
  byte prev = from[0];
  unsigned run_length = 1;
  ++freq[prev].second;
  for(uint64 i = 1; i < length; ++i) {
    if (from[i] == prev && run_length < kMaxLenOfSeq)
      ++run_length;
    else {
      if (run_length > 1)
        UpdateFreqs(run_freq, prev, run_length);
      prev = from[i];
      run_length = 1;
    }
    ++freq[prev].second;
  }
}

struct triple {
  triple(byte sym, byte len, uint32 freq) :
      symbol(sym), length(len), frequency(freq) {}
  byte symbol;
  unsigned length;
  uint32 frequency;
};

/* Answers the question: which one of the sequences is more profitable *
 * to replace with single symbol */
bool CompareTripleDesc(triple t1, triple t2) {
  return (t1.length - 1)*t1.frequency > (t2.length - 1)*t2.frequency;
}

triple ProfitableSequence(const std::map<unsigned, uint32>& seq_freqs,
                          unsigned bval)
{
  unsigned max_length = 1;
  uint32 max_freq = 0;
  for(std::map<unsigned, uint32>::const_iterator it = seq_freqs.begin();
      it != seq_freqs.end(); ++it)
  {
    if ((max_length - 1)*max_freq < (it->first - 1)*it->second) {
      max_length = it->first;
      max_freq = it->second;
    }
  }
  return triple(static_cast<byte>(bval), max_length, max_freq);
}

void FindReplaceableRuns(std::vector<triple> *longest_runs,
                         const std::map<unsigned,uint32> *run_freq,
                         const std::pair<byte, uint64> *freq)
{
  /* We choose only one run for each different character */
  std::vector<triple> &runs = *longest_runs;
  assert(runs.size() == 0);

  for(unsigned i = 0; i < 256; ++i) {
    if(run_freq[i].size() > 0)
      runs.push_back(ProfitableSequence(run_freq[i], i));
  }
  std::sort(runs.begin(), runs.end(), CompareTripleDesc);
  /* Remove the */
  unsigned i = runs.size() - 1;
  while(freq[i].second >= (runs[i].length-1)*runs[i].frequency - 3) {
    runs.pop_back();
    --i;
  }

  /*
  while(runs.size() < 254) {
    if()
      break;
    runs.push_back(sequences[i]);
    ++i;
    }*/
  assert(i == runs.size());
}

} // namespace longruns

uint64 CompressLongRuns(byte *from, uint64 length)
{
  using namespace longruns;

  assert(length > 0);
  assert(from);
  std::pair<byte, uint64> freq[256];
  std::map<unsigned, uint32> run_freq[256];
  InitPairsWithValue<byte, uint64>(freq, 0, 256);
  ComputeRunFrequencies(from, freq, run_freq, length);

  std::sort(&freq[0], &freq[0] + 256, ComparePairSecondAsc<byte, uint64>);

  unsigned free_symbols = 0;
  while(freq[free_symbols].second == 0) ++free_symbols;

  std::vector<triple> longest_runs;
  FindReplaceableRuns(&longest_runs, run_freq, freq);
  
}

} //namespace bwtc
