/**
 * @file Preprocessor.cpp
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
 * Implementation of 2 preprocessing algorithms. One for replacing the
 * most common pairs and another for replacing long runs of the same byte.
 */

#include "../MainBlock.hpp"
#include "../BlockManager.hpp"
#include "../globaldefs.hpp"
#include "../Streams.hpp"
#include "../Utils.hpp"
#include "Preprocessor.hpp"

#include <cassert>
#include <boost/static_assert.hpp>

#include <algorithm>
#include <iostream> /* for std::streamsize*/
#include <map>
#include <string>
#include <utility> /* for pair */


namespace bwtc {

Preprocessor* givePreprocessor(char choice, uint64 block_size,
                               const std::string& input)
{
  Preprocessor* pp;
  /* Expand this to conform for different PreProcessing algorithms */
  switch (choice) {
    case 'n':
    default:
      pp = new Preprocessor(block_size);
  }
  pp->connect(input);
  return pp;
}

Preprocessor::Preprocessor(uint64 block_size) :
    m_source(0), m_blockSize(block_size), m_blockManager(0) { }

Preprocessor::~Preprocessor() {
  delete m_source;
}

void Preprocessor::buildStats(std::vector<byte>* data,
                              std::vector<uint64>* stats, uint64 data_size) {
  std::fill(stats->begin(), stats->end(), 0);
  //TODO: at the moment only contexts of length 1 are supported
  for( uint64 i = 0; i < data_size; ++i)
    (*stats)[(*data)[i]]++; 
}

void Preprocessor::connect(std::string source_name) {
  m_source = new InStream(source_name);
}

void Preprocessor::addBlockManager(BlockManager* manager) {
  m_blockManager = manager;
}

/* We append sentinel to the block here */
MainBlock* Preprocessor::readBlock() {
  assert(m_source);
  assert(m_blockManager);
  std::vector<byte>* to = m_blockManager->getFreeBuffer();
  std::vector<uint64>* stats = m_blockManager->getFreeStats();
  /* TODO:
   * streamsize type has as many bits as long. Since the preprocessor gets
   * blocksize as an uint64 we may end up in problems if user is on 32-bit
   * platform and wants to use too large blocks (amount of bits needed to
   * represent block size is more than 32)?? */
  /*** Stub implementation ***/
  /* We leave on unused byte to the block so that we can use SA-IS-transformer */
  std::streamsize read = m_source->readBlock(
      &(*to)[0], static_cast<std::streamsize>(m_blockSize - 1));
  if (!read) return 0;
  return m_blockManager->makeBlock(to, stats, static_cast<uint64>(read));
}

/*#################### Preprocessing algorithms #############################*/

/* Common utility-functions for preprocessing algorithms */

/* In C++0x these two could be implemented more naturally with the use of
 * lambda-functions. */
template <typename F, typename S>
bool comparePairSecondAsc(std::pair<F,S> p1, std::pair<F,S> p2) {
  return (p1.second < p2.second);
}

template <typename F, typename S>
bool comparePairSecondDesc(std::pair<F,S> p1, std::pair<F,S> p2)
{
  return (p1.second > p2.second);
}

template<typename Key, typename Value>
void initPairsWithValue(std::pair<Key, Value> *pairs, Value v, uint64 length)
{
  for(uint64 i = 0; i < length; ++i)
    pairs[i] = std::make_pair(static_cast<Key>(i), v);
}

FreqTable::FreqTable() {
  for(int i = 0; i < 256; ++i) {
    m_frequencies[i] = std::make_pair(static_cast<byte>(i), 0U);
  }
  initLocations();
  assert(test());
}

FreqTable::FreqTable(uint64* frequencies) {
  /* Assumes that frequencies has length of 256 */
  for(int i = 0; i < 256; ++i) {
    m_frequencies[i] = std::make_pair(static_cast<byte>(i), frequencies[i]);
  }
  std::sort(m_frequencies, m_frequencies + 256, comparePairSecondAsc<byte, uint64>);
  initLocations();
}

const uint64& FreqTable::operator[](uint32 i) const {
  assert(i <= 255);
  return m_frequencies[i].second;
}

byte FreqTable::key(uint32 i) const {
  assert(i < 256);
  return m_frequencies[i].first;
}

bool FreqTable::decrease(uint32 key, uint64 value) {
  assert(key <= 255);
  uint32 freq_index = m_location[key];
  if(m_frequencies[freq_index].second < value) return false;
  uint64 new_value = m_frequencies[freq_index].second - value;
  std::pair<byte, uint64> new_pair =
      std::make_pair(m_frequencies[freq_index].first, new_value);
  
  while (freq_index > 0 && new_value < m_frequencies[freq_index - 1].second)
  {
    ++m_location[m_frequencies[freq_index - 1].first];
    m_frequencies[freq_index] = m_frequencies[freq_index - 1];
    --freq_index;
  }
  m_frequencies[freq_index] = new_pair;
  m_location[new_pair.first]= freq_index;
  assert(m_frequencies[m_location[new_pair.first]].first == new_pair.first);
  return true;
}

void FreqTable::increase(uint32 key, uint64 value) {
  assert(key <= 255);
  uint32 freq_index = m_location[key];
  
  uint64 new_value = m_frequencies[freq_index].second + value;
  std::pair<byte, uint64> new_pair = 
      std::make_pair(m_frequencies[freq_index].first, new_value);
  
  while (freq_index < 255 && new_value > m_frequencies[freq_index + 1].second)
  {
    --m_location[m_frequencies[freq_index + 1].first];
    m_frequencies[freq_index] = m_frequencies[freq_index + 1];
    ++freq_index;
  }
  m_frequencies[freq_index] = new_pair;
  m_location[new_pair.first]= freq_index;
  assert(m_frequencies[m_location[new_pair.first]].first == new_pair.first);
}

void FreqTable::initLocations() {
  for(uint32 i = 0; i < 256; ++i) {
    m_location[m_frequencies[i].first] = i;
  }
}

bool FreqTable::test() {
  for(int i = 0; i < 256; ++i) {
    assert(m_frequencies[m_location[i]].first == i );
  }
  return true;
}

/*##################### Replacing the most common pairs ######################*/

/**
 * Implementation  for replacing the most common pairs.
 *
 * Reasoning behind choosing the replaceable pairs:
 *
 * For each pair we replace with some symbol we have to write that symbol
 * and the correspongind pair to the header. We also make some symbols
 * 'free' from the original data, which means that we escape these symbols
 * with special escape byte. Then we can use these freed symbols in pair
 * replacements.
 *
 * Denote the frequency of symbol x_i with f(x_i) and frequency of pair
 * P_i with f(P_i). Writing the replacement info for single replacement
 * takes 3 bytes. So if P_i is going to be replaced with x_i we require
 * that
 *             f(x_i) + 3 < f(P_i)                              (p1)
 * On the left side is amount of additional bytes after replacement and on
 * the right side is amount of bytes saved when writing pairs.
 *
 * If we free symbols the effect of making one byte to escape byte has to
 * be also notified. Let x_i,...,x_j be the freed symbols, P_i,...,P_j
 * the pairs replaced with x_k's and x be the escape byte.
 * Total improvement from replacements is
 *      sum from k=i to j: f(P_k) - f(x_k) - 3
 * The penalty from escape byte is f(x) so we require that
 *   sum from k=i to j: f(P_k) - f(x_k) - 3 > f(x)              (p2)
 */
namespace commonpairs {

void computePairFrequencies(byte *data, uint64* freqs,
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
void findReplaceablePairs(std::vector<std::pair<uint16,uint32> >*
                          replaceable_pairs,
                          std::pair<uint16, uint32>* pair_freqs,
                          FreqTable *freqs)
{
  const uint32 kStep = 256;
  uint32 current_pair = 0, current_symbol = 0;
  uint32 limit = 0;

  while(replaceable_pairs->size() < 254) {
    if(current_pair + 1 >= limit) {
      limit += kStep;
      assert(limit < 65536);
      std::partial_sort(pair_freqs + current_pair, pair_freqs + limit,
                        pair_freqs  + 65536,
                        comparePairSecondDesc<uint16, uint32>);
    }
    byte fst = static_cast<byte>((pair_freqs[current_pair].first & 0xFF00) >> 8);
    byte snd = static_cast<byte>(pair_freqs[current_pair].first & 0x00FF);
    if(fst == snd) {
      ++current_pair;
      continue;
    }
    if (!freqs->decrease(fst, pair_freqs[current_pair].second)) {
      ++current_pair;
      continue;
    }
    if(!freqs->decrease(snd, pair_freqs[current_pair].second)) {
      freqs->increase(fst, pair_freqs[current_pair].second);      
      ++current_pair;
      continue;
    }
    // TODO: make the control flow more readable 
    
    /* Condition (p1) */
    if((*freqs)[current_symbol] + 3 >= pair_freqs[current_pair].second) {
      freqs->increase(fst, pair_freqs[current_pair].second);
      freqs->increase(snd, pair_freqs[current_pair].second);
      break; /* We won't benefit from any changes any more*/
    }

    /* Reject pairs which have conflicting symbols. */
    bool valid = true;
    for(std::vector<std::pair<uint16, uint32> >::iterator it =
            replaceable_pairs->begin();
        it != replaceable_pairs->end(); ++it)
    {
      uint16 current_fst = static_cast<byte>((it->first & 0xFF00) >> 8);
      uint16 current_sec = static_cast<byte>((it->first & 0x00FF));
      if (current_fst == snd || current_sec == fst) {
        valid = false;
        freqs->increase(fst, pair_freqs[current_pair].second);
        freqs->increase(snd, pair_freqs[current_pair].second);
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
// TODO: uint32 -> possibly 64-bit
uint32 escapeCharIndex(FreqTable* freqs,
                       const std::vector<std::pair<uint16, uint32> >&
                       suitable_pairs,
                       uint32 free_symbols)
{
  if (suitable_pairs.size() <= free_symbols) return free_symbols;
  int64 utility = 0; 
  uint32 i; 
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
    freqs->increase(fst, suitable_pairs[i].second);
    freqs->increase(snd, suitable_pairs[i].second);
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
void writePair(uint16 pair, byte *address) {
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
uint64 writeReplacements(byte *replacements, byte *to, byte *from, uint64 length,
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
void constructReplacementTable(uint32 symbols, byte *replacements,
                               const std::vector<std::pair<uint16, uint32> >&
                               replaceable_pairs, const FreqTable& freqs,
                               bool escaping, uint32 free_symbols)
{
  assert(symbols > 0);
  byte common_byte = replacements[0];
  uint32 limit = (escaping)?symbols-1:symbols;
  for(uint32 i = 0; i < limit; ++i) {
    replacements[replaceable_pairs[i].first] = freqs.key(i);
  }
  byte escape_byte = freqs.key(symbols-1);
  /* Prepare escaped characters for writing of the replacements */
  for(uint32 i = free_symbols; i < symbols; ++i) {
    uint16 pair_value = freqs.key(i) << 8;
    for(uint32 j = 0; j < 256; ++j, ++pair_value)
      if(replacements[pair_value] == common_byte)
        replacements[pair_value] = escape_byte;
  }
}

/**
 * Writes header for pair replacements and gathers replacement table for pairs.
 */
uint32 writeReplaceablePairs(byte *to, byte *replacements, byte escape_byte,
                               byte common_byte)
{
  uint32 position = 0;
  byte last_value = escape_byte;
  for(uint32 i = 0; i < 65536; ++i) {
    if(replacements[i] != escape_byte && replacements[i] != common_byte) {
      last_value = replacements[i];
      *to++ = last_value;
      writePair(i&0xFFFF, to);
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
uint64 compressCommonPairs(byte *from, uint64 length)
{
  using namespace commonpairs;

  assert(length > 0);
  assert(from);
  uint64 freq[256] = {0U};
  std::pair<uint16, uint32> pair_freq[65536];
  initPairsWithValue<uint16, uint32>(pair_freq, 0, 65536);
  computePairFrequencies(from, freq, pair_freq, length);

  FreqTable freqs(freq);
  uint32 free_symbols = 0;
  while(freqs[free_symbols] == 0) ++free_symbols;

  std::vector<std::pair<uint16, uint32> > replaceable_pairs;
  findReplaceablePairs(&replaceable_pairs, pair_freq, &freqs);

  uint32 escape_index = free_symbols;
  if(replaceable_pairs.size() > free_symbols) {
    escape_index = escapeCharIndex(&freqs, replaceable_pairs, free_symbols);
  }
  byte common_byte = freqs.key(255);
  byte escape_byte;
  if(escape_index > free_symbols) escape_byte = freqs.key(escape_index);
  else escape_byte = common_byte;

  /* Initialize replacement table and write header */
  uint32 candidates = static_cast<uint32>(replaceable_pairs.size());
  uint32 symbols = (candidates < free_symbols)?candidates:
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
    constructReplacementTable(symbols, replacements, replaceable_pairs,
                              freqs, escape_byte != common_byte,
                              free_symbols);
    total_size = writeReplaceablePairs(temp, replacements, escape_byte,
                                       common_byte);
  } else {
    for(int i = 0; i < 3; ++i)
      temp[total_size++] = 0;
  }
  total_size += writeReplacements(replacements, temp + total_size, from,
                                  length, common_byte, escape_byte);
  assert(total_size <= length + 3);
  std::copy(temp, temp + total_size, from);
  assert(temp[total_size-1] != escape_byte || escape_byte == common_byte);
  delete [] temp;
  return total_size;
}




/*##################### Replacing the runs of same byte ######################*/

namespace longruns {

// this exists because of faster computation of log
const size_t kMaxLenOfSeq = 1 << 15;

struct triple {
  triple(byte sym, uint32 len, uint32 freq) :
      symbol(sym), length(len), frequency(freq) {}

  triple(const triple& t) {
    symbol = t.symbol;
    length = t.length;
    frequency = t.frequency;
  }

  byte symbol;
  uint32 length;
  uint32 frequency;
};

/* Answers the question: which one of the sequences is more profitable *
 * to replace with single symbol */
bool compareTripleDesc(triple t1, triple t2) {
  return (t1.length - 1)*t1.frequency > (t2.length - 1)*t2.frequency;
}

bool operator>(triple t1, triple t2) {
  return compareTripleDesc(t1, t2);
}

bool operator<(triple t1, triple t2) {
  return compareTripleDesc(t2, t1);
}

class SequenceHeap {
  /* Class for storing the counts of runs, picking the most profitable run and *
   * updating the frequencies of runs. Essentially this is heap with additional*
   * keys, which provide fast search. This class doesn't do any memory since   *
   * the idea is to use this only during the compressLongRuns. */

#define parent(x) (((x)-1)/2)
#define left(x) (2*(x) + 1)
#define right(x) (2*(x) + 2)

 public:
  SequenceHeap(std::vector<triple> &sequences) : seqs_(sequences)
  {
    initLocations();
    last_ = seqs_.size() - 1;
    buildMaxHeap();
  }

  triple removeMax() {
    triple max = seqs_[0];
    /* Update or delete the sequences of same byte */
    for(std::map<uint32, uint32>::iterator it =
            locations_[max.symbol].begin();
        it != locations_[max.symbol].end();
        ++it)
    {
      if(it->first >= max.length) {
        remove(it->second);

      }
      else {
        triple target = seqs_[it->second];
        decrease(it->second, (max.length/target.length)*max.frequency);
      }
    }
    return max;
  }

 private:
  std::vector<triple> &seqs_;
  /* Pairs in maps are <length of sequence, location in seqs_ >. They
   * are indexed by their byte-value */
  std::map<uint32, uint32> locations_[256];
  int last_;

  void decrease(uint32 index, uint32 value) {
    if (static_cast<int>(index) <= last_) return;
    assert(seqs_[index].frequency >= value);
    seqs_[index].frequency -= value;
    heapify(index);
  }

  void remove(uint32 index) {
    if (static_cast<int>(index) > last_) return;
    locations_[seqs_[index].symbol][seqs_[index].length] = last_;
    locations_[seqs_[last_].symbol][seqs_[last_].length] = index;    
    std::swap(seqs_[index], seqs_[last_]);
    --last_;
    heapify(index);
  }

  void initLocations() {
    for(uint32 i = 0; i < seqs_.size(); ++i) {
      locations_[seqs_[i].symbol].insert(
          std::pair<uint32, uint32>(seqs_[i].length, i));
      assert(seqs_[i].length > 0);
    }
  }

  void heapify(int i) {
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


  void buildMaxHeap() {
    for(int i = parent(last_); i >= 0; --i) {
      heapify(i);
    }
  }

#undef parent
#undef left
#undef right
};

void updateFreqs(std::vector<size_t> *run_freq, byte symbol, uint32 length)
{
  assert(length <= kMaxLenOfSeq);
  assert(length > 1);
  length &= 0xFFFFFFFE;
  while(length) { /* Compute the number of sequences of length 2^k for some k */
    uint32 loglongest = utils::logFloor(length);
    ++run_freq[symbol][loglongest-1];
    length ^= (1 << loglongest); 
  }
}

void computeRunFrequencies(byte *from, uint64 *freq,
                           std::vector<size_t> *run_freq, uint64 length)
{
  byte prev = from[0];
  size_t run_length = 1;
  ++freq[prev];
  for(size_t i = 1; i < length; ++i) {
    if (from[i] == prev && run_length < kMaxLenOfSeq)
      ++run_length;
    else {
      if (run_length > 1)
        updateFreqs(run_freq, prev, run_length);
      prev = from[i];
      run_length = 1;
    }
    ++freq[prev];
  }
}

void findReplaceableRuns(std::vector<triple> *runs,
                         std::vector<triple> *longest_runs,
                         FreqTable *freqs)
{
  SequenceHeap seq_heap(*runs);
  assert(longest_runs->size() == 0);

  uint32 current_symbol = 0;
  while(longest_runs->size() < std::min(static_cast<size_t>(254),runs->size()))
  {
    triple best = seq_heap.removeMax();
    freqs->decrease(best.symbol, best.length*best.frequency);
    if((*freqs)[current_symbol] + 3 >= (best.length - 1)*(best.frequency)) {
      freqs->increase(best.symbol, best.length*best.frequency);
      break;
    }
    longest_runs->push_back(best);
    ++current_symbol;
  }
  assert(current_symbol == longest_runs->size());
}

uint32 escapeCharIndex(FreqTable *freqs, const std::vector<triple>& runs,
                         uint32 free_symbols)
{
  if(runs.size() <= free_symbols) return free_symbols;
  int64 utility = 0;
  uint32 i;
  for(i = free_symbols; i < runs.size(); ++i) {
    utility += ((runs[i].length - 1)*runs[i].frequency - (*freqs)[i] - 3);
  }
  while(utility <= static_cast<int64>((*freqs)[i]) &&
        i > free_symbols)
  {
    --i;
    freqs->increase(runs[i].symbol, (runs[i].length - 1)*runs[i].frequency);
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

  bool empty(byte key) {
    return table_[key] == -1;
  }

  void pushBack(byte run_symbol, uint32 length, byte replacement) {
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

  int listBegin(byte symbol) const {
    return table_[symbol];
  }

  const runlist_elem& listElement(int index) const {
    assert(index >= 0);
    assert(index < static_cast<int>(rpls_.size()));
    return rpls_[index];
  }

 private:
  int table_[256];
  std::vector<runlist_elem> rpls_;
};

class RepleacableRun {
 public:
  RepleacableRun(uint32 length, byte symbol, byte replacement)
      : m_length(length), m_symbol(symbol), m_replacement(replacement) {}
  byte symbol() const { return m_symbol;} 
  byte replacement() const { return m_replacement;} 
  uint32 length() const { return m_length;}

  bool operator<(const RepleacableRun& other) const {
    return symbol() < other.symbol() || length() < other.length();
  }

 private:
  uint32 m_length;
  byte m_symbol;
  byte m_replacement;
};


class ReplacementTable2 {
 public:
 
  class RTIterator {
   public:
    RTIterator(const ReplacementTable2& rt, byte symbol, int index)
        : m_rt(rt), m_index(index), m_symbol(symbol) {}
    
    void operator++() { ++m_index; }
    void operator--() { --m_index; }

    bool atEnd() const {
      return m_index >= m_rt.m_buckets.size() || m_rt.m_buckets[m_index].symbol() != m_symbol;
    }

    bool onCorrectBucket() const {
      return m_index >= 0 && m_rt.m_buckets[m_index].symbol() == m_symbol;
    }

    const RepleacableRun& operator*() const {
      return m_rt.m_buckets[m_index];
    }
    
   private:
    const ReplacementTable2& m_rt;
    int m_index;
    const byte m_symbol;
  };

  typedef RTIterator iterator;
  
  ReplacementTable2() {
    std::fill(m_bucketBegins, m_bucketBegins + 256, -1);
  }

  void addReplacement(byte symbol, uint32 length, byte replacement) {
    m_buckets.push_back(RepleacableRun(length, symbol, replacement));
  }

  void addByteToBeEscaped(byte symbol) {
    m_bucketBegins[symbol] = -2;
  }

  bool replaceable(byte symbol) const {
    return m_bucketBegins[symbol] >= 0;
  }
  
  bool toBeEscaped(byte symbol) const {
    return m_bucketBegins[symbol] == -2;
  }

  RTIterator listBegin(byte symbol) const {
    return RTIterator(*this, symbol, m_bucketBegins[symbol]);
  }

  void prepareForQueries() {
    std::sort(m_buckets.begin(), m_buckets.end());
    for(size_t i = 0; i < m_buckets.size(); ++i) {
      byte symbol = m_buckets[i].symbol();
      assert(m_bucketBegins[symbol] != -2);
      if(m_bucketBegins[symbol] == -1) m_bucketBegins[symbol] = i;
    }
  }
  
 private:
  /* There exists 256 (one for each byte) buckets (some of them may be empty)
   * packed in a single vector. Single bucket is ordered in the ascending order
   * of the length of runs. The starting indices of buckets are stored in
   * m_bucketBegins. If there isn't bucket for character the value is -1.
   * If the character is to be escaped the value is -2. */
  int m_bucketBegins[256];
  std::vector<RepleacableRun> m_buckets;
};

uint32 writeRunReplacement(const ReplacementTable& repl, uint32 run_length,
                           byte escape, byte symbol, byte *to)
{
  int tbl_index = repl.listBegin(symbol);
  uint32 j = 0;
  do {
    if (tbl_index == -1) {
      std::fill(to + j, to + j + run_length, symbol);
      return run_length + j;
    }
    const runlist_elem& el = repl.listElement(tbl_index);
    uint32 times = run_length/el.length; //el.length == 2^k for some k
    if(el.length == 1) {
      for(uint32 k = 0; k < times; ++k) {
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

uint64 writeReplacements(const ReplacementTable& replacements, byte *to,
                         byte *from, uint64 length, byte escape)
{
  uint64 j = 0; /* Index of target */
  byte prev = from[0];
  uint32 run_length = 1;
  for(uint64 i = 1; i < length; ++i) {
    if(prev == from[i] && run_length < kMaxLenOfSeq)
      ++run_length;
    else {
      j += writeRunReplacement(replacements, run_length, escape,
                               prev, to + j);
      prev = from[i];
      run_length = 1;
    }
  }
  j += writeRunReplacement(replacements, run_length, escape,
                           prev, to + j);
  return j;
}

size_t writeReplacements(const ReplacementTable2& rt, byte *to, byte *from,
                         size_t length, byte escape)
{

  size_t j = 0, i = 1;
  byte prev = from[0];
  uint32 runLength = 1;
  for(; i < length; ++i) {
    if(!rt.replaceable(prev)) {
      if(rt.toBeEscaped(prev)) to[j++] = escape;
      to[j++] = prev;
      prev = from[i];
      runLength = 1;
    } else {
      ReplacementTable2::iterator it = rt.listBegin(prev);
      while(i < length) {
        if(prev == from[i]) {
          ++runLength;
          if(runLength == (*it).length())
            ++it;
          if(it.atEnd()) {
            --it;
            break;
          }
          ++i;
        } else {
          --i;
          break;
        }
      }
      while(runLength > 0) {
        bool corrList = it.onCorrectBucket();
        while(corrList && runLength < (*it).length()) {
          --it;
          corrList = it.onCorrectBucket();
        }
        if(corrList) {
          size_t times = runLength/(*it).length();
          std::fill(to + j, to + j + times, (*it).replacement());
          j += times;
          runLength -= times * (*it).length();
        } else {
          std::fill(to + j, to + j + runLength, prev);
          j += runLength;
          break;
        }
      }
      if(i < length+1) prev = from[++i];
      runLength = 1;
    }
  }
  if(i == length) to[j++] = prev;
  return j;
}


} // namespace longruns


uint64 compressLongRuns(byte *from, uint64 length)
{
  using namespace longruns;

  assert(length > 0);
  assert(from);
  uint64 freq[256] = { 0 };
  /* Table for runs. Vector is indexed by log(runlength) - 1 */
  std::vector<size_t> run_freq[256];
  size_t vsize = utils::logFloor(std::min(length, kMaxLenOfSeq));
  for(size_t i = 0; i < 256; ++i) {
    run_freq[i].resize(vsize-1);
    std::fill(run_freq[i].begin(), run_freq[i].end(), 0);
  }

  computeRunFrequencies(from, freq, run_freq, length);

  FreqTable freqs(freq);

  uint32 free_symbols = 0;
  while(freqs[free_symbols] == 0) ++free_symbols;

  std::vector<triple> runs, longest_runs;
  runs.reserve(256);
  for(uint32 i = 0; i < 256; ++i) {
    for(size_t j = 0; j < run_freq[i].size(); ++j) 
    {
      if(run_freq[i][j] > 0) {
        runs.push_back(triple(static_cast<byte>(i), 1 << (j+1), run_freq[i][j]));
      }
    }
  }
  findReplaceableRuns(&runs,&longest_runs, &freqs);

  uint32 escape_index = free_symbols;
  if(longest_runs.size() > free_symbols) {
    escape_index = escapeCharIndex(&freqs, longest_runs, free_symbols);
  }
  uint32 new_symbols =  (escape_index == free_symbols)?
      0 : escape_index - free_symbols + 1;
  uint32 symbols_in_use = (new_symbols > 0) ?
      (escape_index + 1) : std::min((size_t)free_symbols, longest_runs.size());

  uint32 run_replacements = (new_symbols > 0)?
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
  byte escape_byte = freqs.key(escape_index);

  if(symbols_in_use > 0) {
    //TODO: Store replacement in somewhere so they can be merged with other
    uint32 limit = run_replacements & 0xFFFFFFFE;

    for(uint32 i = 0; i < limit; i += 2) {
      temp[position++] = freqs.key(i);
      assert(longest_runs[i].length <= kMaxLenOfSeq);
      assert(longest_runs[i+1].length <= kMaxLenOfSeq);
      byte lengths = (utils::logFloor(longest_runs[i].length) << 4) |
          utils::logFloor(longest_runs[i+1].length);
      temp[position++] = lengths;
      temp[position++] = longest_runs[i].symbol;
      temp[position++] = freqs.key(i+1);
      temp[position++] = longest_runs[i+1].symbol;
    }
    byte sentinel = (escape_index != free_symbols) ?
        escape_byte : freqs.key(symbols_in_use - 1);
    
    if( run_replacements != limit )
    {
      assert(limit == symbols_in_use - 1 || limit == symbols_in_use - 2 );
      temp[position++] = freqs.key(limit);
      temp[position++] = (utils::logFloor(longest_runs[limit].length) << 4);
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
  ReplacementTable replacements(escape_byte);
  //ReplacementTable2 replacements;
  /* Escaped characters*/
  if( new_symbols > 0) {
    for(uint32 i = free_symbols; i <= escape_index; ++i) {
      replacements.pushBack(freqs.key(i), 1, escape_byte);
      //replacements.addByteToBeEscaped(freqs.key(i));
    }
  }
  for (uint32 i = 0; i < run_replacements; ++i) {
    assert( i < longest_runs.size() );
    replacements.pushBack(longest_runs[i].symbol, longest_runs[i].length,
                          freqs.key(i));
    //replacements.addReplacement(longest_runs[i].symbol, longest_runs[i].length,
    //                            freqs.key(i));
  }
  //replacements.prepareForQueries();
  uint64 total_size = position;
  total_size += writeReplacements(replacements, temp + position, from, length,
                                  escape_byte);
  std::copy(temp, temp + total_size, from);
  delete [] temp;
  assert(total_size <= length + 2);
  return total_size;
}

} //namespace bwtc
