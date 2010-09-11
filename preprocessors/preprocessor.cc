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

/* We append sentinel to the block here */
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
  /* We leave on unused byte to the block so that we can use SA-IS-transformer */
  std::streamsize read = source_->ReadBlock(
      &(*to)[0], static_cast<std::streamsize>(block_size_ - 1));
  if (!read) return NULL;
  return block_manager_->MakeBlock(to, stats, static_cast<uint64>(read));
}

/* Floor of logarithm of base two */
byte LogFloor(unsigned n) {
  assert(n > 0);
  byte log = 0;
  while(n > 1) {
    n >>= 1;
    ++log;
  }
  return log;
}

unsigned MostSignificantBit16(unsigned n) {
  assert(n < (1 << 16));
  n |= (n >> 1);
  n |= (n >> 2);
  n |= (n >> 4);
  n |= (n >> 8);
  return n & (~(n >> 1));
}

unsigned MostSignificantBit(unsigned n) {
  assert(sizeof(unsigned) == 4);
  n |= (n >> 1);
  n |= (n >> 2);
  n |= (n >> 4);
  n |= (n >> 8);
  n |= (n >> 16);
  return n & (~(n >> 1));
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

FreqTable::FreqTable() {
  for(int i = 0; i < 256; ++i) {
    freq_[i] = std::make_pair(static_cast<byte>(i), 0U);
  }
  InitLocations();
  assert(Test());
}

FreqTable::FreqTable(uint64* frequencies) {
  /* Assumes that frequencies has length of 256 */
  for(int i = 0; i < 256; ++i) {
    freq_[i] = std::make_pair(static_cast<byte>(i), frequencies[i]);
  }
  std::sort(freq_, freq_ + 256, ComparePairSecondAsc<byte, uint64>);
  InitLocations();
}

const uint64& FreqTable::operator[](unsigned i) const {
  assert(i <= 255);
  return freq_[i].second;
}

byte FreqTable::Key(unsigned i) const {
  assert(i < 256);
  return freq_[i].first;
}

bool FreqTable::Decrease(unsigned key, uint64 value) {
  assert(key <= 255);
  unsigned freq_index = location_[key];
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

void FreqTable::Increase(unsigned key, uint64 value) {
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

void FreqTable::InitLocations() {
  for(unsigned i = 0; i < 256; ++i) {
    location_[freq_[i].first] = i;
  }
}

bool FreqTable::Test() {
  for(int i = 0; i < 256; ++i) {
    assert(freq_[location_[i]].first == i );
  }
  return true;
}

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
/**
 * Implementation  for replacing the most common pairs.
 */
namespace commonpairs {

void ComputePairFrequencies(byte *data, uint64* freqs,
                            std::pair<uint16, unsigned>* pair_freqs, uint64 len)
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

/**
 * Finds the candidates for pairs which are going to be replaced by
 * single symbols.
 *
 * If pair = <p1,p2> is in result, then isn't any pair where p2 is the first
 * symbol or p1 is the second symbol. Greedy heuristic is used. Optimal
 * solution for choosing the pairs is NP-hard (max-cut).
 *
 * @param replaceable_pairs empty vector where pairs are stored
 * @param pair_freqs unsorted array of <pair, frequency>-pairs with
 *                   65536 entries
 * @param freqs pointer to FreqTable-object
 * @return result is stored into replaceable_pairs
 */
void FindReplaceablePairs(std::vector<std::pair<uint16,unsigned> >*
                          replaceable_pairs,
                          std::pair<uint16, unsigned>* pair_freqs,
                          FreqTable *freqs)
{
  const unsigned kStep = 256;
  unsigned current_pair = 0, current_symbol = 0;
  unsigned limit = 0;

  while(replaceable_pairs->size() < 254) {
    if(current_pair + 1 >= limit) {
      limit += kStep;
      assert(limit < 65536);
      std::partial_sort(pair_freqs + current_pair, pair_freqs + limit,
                        pair_freqs  + 65536,
                        ComparePairSecondDesc<uint16, unsigned>);
    }
    byte fst = static_cast<byte>((pair_freqs[current_pair].first & 0xFF00) >> 8);
    byte snd = static_cast<byte>(pair_freqs[current_pair].first & 0x00FF);
    if(fst == snd) {
      ++current_pair;
      continue;
    }
    if (!freqs->Decrease(fst, pair_freqs[current_pair].second)) {
      ++current_pair;
      continue;
    }
    if(!freqs->Decrease(snd, pair_freqs[current_pair].second)) {
      freqs->Increase(fst, pair_freqs[current_pair].second);      
      ++current_pair;
      continue;
    }
    // TODO: make the control flow more readable 
    
    /* Condition (p1) */
    if((*freqs)[current_symbol] + 3 >= pair_freqs[current_pair].second) {
      freqs->Increase(fst, pair_freqs[current_pair].second);
      freqs->Increase(snd, pair_freqs[current_pair].second);
      break; /* We won't benefit from any changes any more*/
    }
    /* Reject pairs which have conflicting symbols. */
    bool valid = true;
    for(std::vector<std::pair<uint16, unsigned> >::iterator it =
            replaceable_pairs->begin();
        it != replaceable_pairs->end(); ++it)
    {
      uint16 current_fst = static_cast<byte>((it->first & 0xFF00) >> 8);
      uint16 current_sec = static_cast<byte>((it->first & 0x00FF));
      if (current_fst == snd || current_sec == fst) {
        valid = false;
        freqs->Increase(fst, pair_freqs[current_pair].second);
        freqs->Increase(snd, pair_freqs[current_pair].second);
        break;
      }
    }
    if(valid) {
      replaceable_pairs->push_back(pair_freqs[current_pair]);
      ++current_symbol;
    }
    ++current_pair;
  }
  assert((size_t)current_symbol == replaceable_pairs->size());
}

/* Returns the index for the 'escape_char' in freqs-array or value of
 * free_symbols if freeing of the chars is not profitable */
unsigned EscapeCharIndex(FreqTable* freqs,
                         const std::vector<std::pair<uint16, unsigned> >&
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
    freqs->Increase(fst, suitable_pairs[i].second);
    freqs->Increase(snd, suitable_pairs[i].second);
    utility -= (suitable_pairs[i].second - (*freqs)[i] - 3);
  }
  return i;
}

/**
 * Writes pair to given address. Used for writing the pair-replacement info
 * to result-array.
 *
 * @param pair Two bytes represented with uint16. First 8 bits represent the
 *             first byte of pair. The latter 8 bits are the second byte of pair.
 * @param address pointer to where the pair is written
 */
void WritePair(uint16 pair, byte *address) {
  *address = (pair >> 8);
  *(address + 1) = pair & 0xFF;
}

/**
 * Writes replacements of pairs to result array.
 *
 * @param replacements Array of size 65536 where index of element is
 *                     interpreted as a pair of two bytes where their
 *                     bit-representations are concatenated.
 *                     Let p be the value of pair. We require the following
 *                     conditions: if replacements[p] == common_byte there is
 *                     no replacement for p if replacements[p] == escape_byte
 *                     the first element of p needs escaping because it is made
 *                     free in other cases replacement for p is value at
 *                     replacements[p].
 * @param to an array where the replacements are written
 * @param from source of an original data
 * @param length length of from (source)
 * @param common_byte explained above
 * @param escape_byte explained above
 * @return Bytes used at result array when replaced the pairs.
 */
uint64 WriteReplacements(byte *replacements, byte *to, byte *from, uint64 length,
                         byte common_byte, byte escape_byte)
{
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
    }
    else { // pair will be replaced
      to[result_index++] = replacements[pair];
      if( i == length - 1) break;
      pair = from[++i];
    }

    if ( i == length - 1) {
      pair = (from[i] << 8) | 0;
      if(replacements[pair] == escape_byte && escape_byte != common_byte)
        to[result_index++] = escape_byte;
      to[result_index++] = from[i];
      break;
    }
    ++i;
  }
  return result_index;
}

/**
 * Constructs replacement table.
 *
 * @param symbols Number of symbols used in writing the replacements,
 * @param replacements Array of size 65536 indexed by concatenated pairs.
 *                     The replacements for pairs are writen to this. 
 * @param freqs {@link FreqTable} containing frequencies of characters.
 * @param escaping true if escaping is in use.
 * @param free_symbols Number of symbols which aren't present in the source.
 * @return Result is written to replacements.
 */
void ConstructReplacementTable(unsigned symbols, byte *replacements,
                               const std::vector<std::pair<uint16, unsigned> >&
                               replaceable_pairs, const FreqTable& freqs,
                               bool escaping, unsigned free_symbols)
{
  assert(symbols > 0);
  byte common_byte = replacements[0];
  unsigned limit = (escaping)?symbols-1:symbols;
  for(unsigned i = 0; i < limit; ++i) {
    replacements[replaceable_pairs[i].first] = freqs.Key(i);
  }
  byte escape_byte = freqs.Key(symbols-1);
  /* Prepare escaped characters for writing of the replacements */
  for(unsigned i = free_symbols; i < symbols; ++i) {
    uint16 pair_value = freqs.Key(i) << 8;
    for(unsigned j = 0; j < 256; ++j, ++pair_value)
      if(replacements[pair_value] == common_byte)
        replacements[pair_value] = escape_byte;
  }
}

/**
 * Writes header for pair replacements and gathers replacement table for pairs.
 */
unsigned WriteReplaceablePairs(byte *to, byte *replacements, byte escape_byte,
                               byte common_byte)
{
  unsigned position = 0;
  byte last_value = escape_byte;
  for(unsigned i = 0; i < 65536; ++i) {
    if(replacements[i] != escape_byte && replacements[i] != common_byte) {
      last_value = replacements[i];
      *to++ = last_value;
      WritePair(i&0xFFFF, to);
      to += 2;
      position += 3;
    }
  }
  *to++ = last_value;
  *to = (escape_byte != common_byte)?escape_byte : last_value;
  return position + 2;
}

} //namespace commonpairs

/**
 * Replaces common pairs of bytes with single byte values. Writes the
 * result to source array.
 *
 *  Requires that from has at least 3 additional
 * bytes reserved of length+3, where the actual data is in range [0, length).
 *
 * @param from Source array. It is required that from points to memory area
 *             where at least length+3 bytes is allocated.
 * @param length Input data is in the range [from[0], from[length])
 * @return Result is written to input array (from).
 */ 
uint64 CompressCommonPairs(byte *from, uint64 length)
{
  using namespace commonpairs;

  assert(length > 0);
  assert(from);
  uint64 freq[256] = {0U};
  std::pair<uint16, unsigned> pair_freq[65536];
  InitPairsWithValue<uint16, unsigned>(pair_freq, 0, 65536);
  ComputePairFrequencies(from, freq, pair_freq, length);

  FreqTable freqs(freq);
  unsigned free_symbols = 0;
  while(freqs[free_symbols] == 0) ++free_symbols;

  std::vector<std::pair<uint16, unsigned> > replaceable_pairs;
  FindReplaceablePairs(&replaceable_pairs, pair_freq, &freqs);

  unsigned escape_index = free_symbols;
  if(replaceable_pairs.size() > free_symbols) {
    escape_index = EscapeCharIndex(&freqs, replaceable_pairs, free_symbols);
  }
  byte common_byte = freqs.Key(255);
  byte escape_byte;
  if(escape_index > free_symbols) escape_byte = freqs.Key(escape_index);
  else escape_byte = common_byte;

  /* Initialize replacement table and write header */
  unsigned candidates = static_cast<unsigned>(replaceable_pairs.size());
  unsigned symbols = (candidates < free_symbols)?candidates:
      ((free_symbols < escape_index)?escape_index+1:free_symbols); 

  if (verbosity > 1) {
    std::clog << "Replacing ";
    if (symbols > free_symbols)
      std::clog << (symbols - 1) << " pairs. Made " << (symbols - free_symbols)
                << " symbols free.\n";
    else
      std::clog << symbols << " pairs. No symbols made free.\n";
  }

  byte replacements[65536];
  std::fill(replacements, replacements + 65536, common_byte);
  byte *temp = new byte[length + 3];
  uint64 total_size = 0;
  if(symbols > 0) {
    ConstructReplacementTable(symbols, replacements, replaceable_pairs,
                              freqs, escape_byte != common_byte,
                              free_symbols);
    total_size = WriteReplaceablePairs(temp, replacements, escape_byte,
                                       common_byte);
  } else {
    for(int i = 0; i < 3; ++i)
      temp[total_size++] = 0;
  }
  total_size += WriteReplacements(replacements, temp + total_size, from,
                                  length, common_byte, escape_byte);
  assert(total_size <= length + 3);
  std::copy(temp, temp + total_size, from);
  assert(temp[total_size-1] != escape_byte || escape_byte == common_byte);
  delete [] temp;
  return total_size;
}




/*##################### Replacing the runs of same byte ######################*/

namespace longruns {

const unsigned kMaxLenOfSeq = 1 << 15;

struct triple {
  triple(byte sym, unsigned len, uint32 freq) :
      symbol(sym), length(len), frequency(freq) {}

  triple(const triple& t) {
    symbol = t.symbol;
    length = t.length;
    frequency = t.frequency;
  }

  byte symbol;
  unsigned length;
  uint32 frequency;
};

/* Answers the question: which one of the sequences is more profitable *
 * to replace with single symbol */
bool CompareTripleDesc(triple t1, triple t2) {
  return (t1.length - 1)*t1.frequency > (t2.length - 1)*t2.frequency;
}

bool operator>(triple t1, triple t2) {
  return CompareTripleDesc(t1, t2);
}

bool operator<(triple t1, triple t2) {
  return CompareTripleDesc(t2, t1);
}

class SequenceHeap {
  /* Class for storing the counts of runs, picking the most profitable run and *
   * updating the frequencies of runs. Essentially this is heap with additional*
   * keys, which provide fast search. This class doesn't do any memory since   *
   * the idea is to use this only during the CompressLongRuns. */

#define parent(x) (((x)-1)/2)
#define left(x) (2*(x) + 1)
#define right(x) (2*(x) + 2)

 public:
  SequenceHeap(std::vector<triple> &sequences) : seqs_(sequences)
  {
    InitLocations();
    last_ = seqs_.size() - 1;
    BuildMaxHeap();
  }

  triple DeleteMax() {
    triple max = seqs_[0];
    /* Update or delete the sequences of same byte */
    for(std::map<unsigned, unsigned>::iterator it =
            locations_[max.symbol].begin();
        it != locations_[max.symbol].end();
        ++it)
    {
      if(it->first >= max.length) {
        Delete(it->second);

      }
      else {
        triple target = seqs_[it->second];
        Decrease(it->second, (max.length/target.length)*max.frequency);
      }
    }
    return max;
  }

 private:
  std::vector<triple> &seqs_;
  /* Pairs in maps are <length of sequence, location in seqs_ >. They
   * are indexed by their byte-value */
  std::map<unsigned, unsigned> locations_[256];
  int last_;

  void Decrease(unsigned index, unsigned value) {
    if (static_cast<int>(index) <= last_) return;
    assert(seqs_[index].frequency >= value);
    seqs_[index].frequency -= value;
    Heapify(index);
  }

  void Delete(unsigned index) {
    if (static_cast<int>(index) > last_) return;
    locations_[seqs_[index].symbol][seqs_[index].length] = last_;
    locations_[seqs_[last_].symbol][seqs_[last_].length] = index;    
    std::swap(seqs_[index], seqs_[last_]);
    --last_;
    Heapify(index);
  }

  void InitLocations() {
    for(unsigned i = 0; i < seqs_.size(); ++i) {
      locations_[seqs_[i].symbol].insert(
          std::pair<unsigned, unsigned>(seqs_[i].length, i));
      assert(seqs_[i].length > 0);
    }
  }

  void Heapify(int i) {
    int l = left(i), r = right(i);
    while(r <= last_) {
      int largest = (seqs_[l] > seqs_[r]) ? l : r;
      if(seqs_[i] < seqs_[largest]) {
        assert(locations_[seqs_[largest].symbol].count(seqs_[largest].length) > 0);
        assert(locations_[seqs_[i].symbol].count(seqs_[i].length) > 0);
        std::swap(locations_[seqs_[i].symbol][seqs_[i].length],
                  locations_[seqs_[largest].symbol][seqs_[largest].length]);
        std::swap(seqs_[i], seqs_[largest]);
        i = largest;
        l = left(i), r = right(i);
      } else return;
    }
    if( l == last_ && seqs_[i] < seqs_[l]) {
        std::swap(locations_[seqs_[i].symbol][seqs_[i].length],
                  locations_[seqs_[l].symbol][seqs_[l].length]);
        std::swap(seqs_[i], seqs_[l]);
    }
  } 


  void BuildMaxHeap() {
    for(int i = parent(last_); i >= 0; --i) {
      Heapify(i);
    }
  }

#undef parent
#undef left
#undef right
};

void UpdateFreqs(std::map<unsigned, uint32> *run_freq, byte symbol,
                 unsigned length)
{
  assert(length <= kMaxLenOfSeq);
  assert(length > 1);
  length -= (length % 2);
  unsigned original = length;
  while(length) { /* Compute the number of sequences of length 2^k for some k */
    unsigned longest = MostSignificantBit16(length);
    if (run_freq[symbol].count(longest))
      run_freq[symbol][longest] += original/longest;
    else run_freq[symbol][longest] = original/longest;
    length -= longest;
  }
}

void ComputeRunFrequencies(byte *from, uint64 *freq,
                           std::map<unsigned, uint32> *run_freq, uint64 length)
{
  byte prev = from[0];
  unsigned run_length = 1;
  ++freq[prev];
  for(uint64 i = 1; i < length; ++i) {
    if (from[i] == prev && run_length < kMaxLenOfSeq)
      ++run_length;
    else {
      if (run_length > 1)
        UpdateFreqs(run_freq, prev, run_length);
      prev = from[i];
      run_length = 1;
    }
    ++freq[prev];
  }
}

void FindReplaceableRuns(std::vector<triple> *runs,
                         std::vector<triple> *longest_runs,
                         FreqTable *freqs)
{
  SequenceHeap seq_heap(*runs);
  assert(longest_runs->size() == 0);

  unsigned current_symbol = 0;
  while(longest_runs->size() < std::min(static_cast<size_t>(254),runs->size()))
  {
    triple best = seq_heap.DeleteMax();
    freqs->Decrease(best.symbol, best.length*best.frequency);
    if((*freqs)[current_symbol] + 3 >= (best.length - 1)*(best.frequency)) {
      freqs->Increase(best.symbol, best.length*best.frequency);
      break;
    }
    longest_runs->push_back(best);
    ++current_symbol;
  }
  assert(current_symbol == longest_runs->size());
}

unsigned EscapeCharIndex(FreqTable *freqs, const std::vector<triple>& runs,
                         unsigned free_symbols)
{
  if(runs.size() <= free_symbols) return free_symbols;
  int64 utility = 0;
  unsigned i;
  for(i = free_symbols; i < runs.size(); ++i) {
    utility += ((runs[i].length - 1)*runs[i].frequency - (*freqs)[i] - 3);
  }
  while(utility <= static_cast<int64>((*freqs)[i]) &&
        i > free_symbols)
  {
    --i;
    freqs->Increase(runs[i].symbol, (runs[i].length - 1)*runs[i].frequency);
    utility -= ((runs[i].length - 1)*runs[i].frequency - (*freqs)[i]- 3);
  }
  return i;
}

struct runlist_elem {
  runlist_elem(int l, byte s, int n) : length(l), symbol(s), next_elem(n) {}
  int length;
  byte symbol;
  int next_elem;
};

class ReplacementTable {

 public:
  ReplacementTable(byte escape_byte) {
    std::fill(table_, table_ + 256, -1);
    rpls_.push_back(runlist_elem(1, escape_byte, -1));
  }

  bool Empty(byte key) {
    return table_[key] == -1;
  }

  void PushBack(byte run_symbol, unsigned length, byte replacement) {
    int list_index = table_[run_symbol];
    if(list_index == 0) {
      table_[run_symbol] = rpls_.size();
      rpls_.push_back(runlist_elem(length, replacement, 0));
      return;
    }
    runlist_elem *e = 0;
    while(list_index != -1 && rpls_[list_index].length > (int)length) {
      e = &rpls_[list_index];
      list_index = e->next_elem;
    }
    int next = -1;
    if(length > 1) {
      if(!e)
        table_[run_symbol] = rpls_.size();
      else {
        next = e->next_elem;
        e->next_elem = rpls_.size();
      }
      rpls_.push_back(runlist_elem(length, replacement, next));
    } else {
      if(!e)
        table_[run_symbol] = 0;
      else {
        e->next_elem = 0;
      }
    }
  }

  int ListBegin(byte symbol) {
    return table_[symbol];
  }

  runlist_elem ListElement(int index) {
    assert(index >= 0);
    assert(index < static_cast<int>(rpls_.size()));
    return rpls_[index];
  }

 private:
  int table_[256];
  std::vector<runlist_elem> rpls_;
};


unsigned WriteRunReplacement(ReplacementTable *repl, unsigned run_length,
                             byte escape, byte symbol, byte *to)
{
  int tbl_index = repl->ListBegin(symbol);
  unsigned j = 0;
  do {
    if (tbl_index == -1) {
      std::fill(to + j, to + j + run_length, symbol);
      return run_length + j;
    }
    runlist_elem el = repl->ListElement(tbl_index);
    unsigned times = run_length/el.length;
    if(el.length == 1) {
      for(unsigned k = 0; k < times; ++k) {
        to[j++] = escape; to[j++] = symbol;
      }
    } else if (times > 0) {
      std::fill(to + j, to + j + times, el.symbol);
      j += times;
    }
    run_length -= times*el.length;
    tbl_index = el.next_elem;
  } while(run_length > 0);
  return j;
}

uint64 WriteReplacements(ReplacementTable *replacements, byte *to, byte *from,
                         uint64 length, byte escape)
{
  uint64 j = 0; /* Index of target */
  byte prev = from[0];
  unsigned run_length = 1;
  for(uint64 i = 1; i < length; ++i) {
    if(prev == from[i] && run_length < kMaxLenOfSeq)
      ++run_length;
    else {
      j += WriteRunReplacement(replacements, run_length, escape,
                               prev, to + j);
      prev = from[i];
      run_length = 1;
    }
  }
  j += WriteRunReplacement(replacements, run_length, escape,
                           prev, to + j);
  return j;
}

} // namespace longruns


/* Use of map for replacments may be too slow. */
uint64 CompressLongRuns(byte *from, uint64 length)
{
  using namespace longruns;

  assert(length > 0);
  assert(from);
  uint64 freq[256] = { 0 };
  std::map<unsigned, uint32> run_freq[256];
  ComputeRunFrequencies(from, freq, run_freq, length);

  FreqTable freqs(freq);

  unsigned free_symbols = 0;
  while(freqs[free_symbols] == 0) ++free_symbols;

  std::vector<triple> runs, longest_runs;
  runs.reserve(256);
  for(unsigned i = 0; i < 256; ++i) {
    for(std::map<unsigned, uint32>::const_iterator it = run_freq[i].begin();
        it != run_freq[i].end(); ++it)
    {
      assert(it->second > 0);
      assert(it->first > 1);
      runs.push_back(triple(static_cast<byte>(i), it->first, it->second));
    }
  }
  FindReplaceableRuns(&runs,&longest_runs, &freqs);

  unsigned escape_index = free_symbols;
  if(longest_runs.size() > free_symbols) {
    escape_index = EscapeCharIndex(&freqs, longest_runs, free_symbols);
  }
  unsigned new_symbols =  (escape_index == free_symbols)?
      0 : escape_index - free_symbols + 1;
  unsigned symbols_in_use = (new_symbols > 0) ?
      (escape_index + 1) : std::min((size_t)free_symbols, longest_runs.size());

  unsigned run_replacements = (new_symbols > 0)?
      (symbols_in_use - 1) : symbols_in_use;

  if (verbosity > 1) {
    std::clog << "Replacing " << run_replacements << " runs. ";
    if (new_symbols > 0)
      std::clog << "Made " << new_symbols << " symbols free.\n";
    else
      std::clog << "No symbols made free.\n";
  }

  byte *temp = new byte[length + 2];
  uint64 position = 0;
  byte escape_byte = freqs.Key(escape_index);

  if(symbols_in_use > 0) {
    //TODO: Store replacement in somewhere so they can be merged with other
    unsigned limit = run_replacements - (run_replacements % 2);

    for(unsigned i = 0; i < limit; i += 2) {
      temp[position++] = freqs.Key(i);
      assert(longest_runs[i].length <= kMaxLenOfSeq);
      assert(longest_runs[i+1].length <= kMaxLenOfSeq);
      byte lengths = (LogFloor(longest_runs[i].length) << 4) |
          LogFloor(longest_runs[i+1].length);
      temp[position++] = lengths;
      temp[position++] = longest_runs[i].symbol;
      temp[position++] = freqs.Key(i+1);
      temp[position++] = longest_runs[i+1].symbol;
    }
    byte sentinel = (escape_index != free_symbols) ?
        escape_byte : freqs.Key(symbols_in_use - 1);
    
    if( run_replacements != limit )
    {
      assert(limit == symbols_in_use - 1 || limit == symbols_in_use - 2 );
      temp[position++] = freqs.Key(limit);
      temp[position++] = (LogFloor(longest_runs[limit].length) << 4);
      temp[position++] = longest_runs[limit].symbol;
      temp[position++] = sentinel;
    } else {
      temp[position++] = sentinel;
      temp[position++] = 0;
    }
  } else {
    /* No replacements being made */
    temp[position++] = 0;
    temp[position++] = 0;
  }
  /* Replacement table */
  ReplacementTable replacements(escape_byte);
  /* Escaped characters*/
  if( new_symbols > 0) {
    for(unsigned i = free_symbols; i <= escape_index; ++i) {
      replacements.PushBack(freqs.Key(i), 1, escape_byte);
    }
  }
  //  for(unsigned i = 0; i < 256; ++i) {
  //  if(replacements.Empty(i))
  //   replacements.PushBack(i, 1, );
  //}
  for (unsigned i = 0; i < run_replacements; ++i) {
    assert( i < longest_runs.size() );
    replacements.PushBack(longest_runs[i].symbol, longest_runs[i].length,
                          freqs.Key(i));
  }
  uint64 total_size = position;
  //  for(unsigned i = 0; i < 256; ++i) {
  //std::sort(replacements[i].rbegin(), replacements[i].rend());
    //std::reverse(replacements[i].begin(), replacements[i].end());
  //}
  total_size += WriteReplacements(&replacements, temp + position, from, length,
                                  escape_byte);
  std::copy(temp, temp + total_size, from);
  delete [] temp;
  assert(total_size <= length + 2);
  return total_size;
}

} //namespace bwtc
