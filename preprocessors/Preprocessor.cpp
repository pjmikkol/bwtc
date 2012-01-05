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
#include "FrequencyTable.hpp"
#include "PairReplacer.hpp"
#include "PairAndRunReplacer.hpp"
#include "RunReplacer.hpp"
#include "../Profiling.hpp"

#include <cassert>
#include <boost/static_assert.hpp>

#include <algorithm>
#include <iostream> /* for std::streamsize*/
#include <map>
#include <string>
#include <utility> /* for pair */
#include <vector>

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

Preprocessor::Preprocessor(uint64 block_size, const std::string& prepr) :
    m_source(0), m_blockSize(block_size), m_blockManager(0),
    m_preprocessingOptions(prepr) {}

Preprocessor::Preprocessor(uint64 block_size) :
    m_source(0), m_blockSize(block_size), m_blockManager(0) {}

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

void Preprocessor::connect(const std::string& source_name) {
  m_source = new InStream(source_name);
}

void Preprocessor::addBlockManager(BlockManager* manager) {
  m_blockManager = manager;
}

/* We append sentinel to the block here */
MainBlock* Preprocessor::readBlock() {
  PROFILE("Preprocessor::readBlock");
  assert(m_source);
  assert(m_blockManager);
  std::vector<byte>* to = m_blockManager->getFreeBuffer();
  std::vector<uint64>* stats = m_blockManager->getFreeStats();
  /* TODO:
   * streamsize type has as many bits as long. Since the preprocessor gets
   * blocksize as an uint64 we may end up in problems if user is on 32-bit
   * platform and wants to use too large blocks (amount of bits needed to
   * represent block size is more than 32)?? */

  /* The unused byte is left the block so that SA-IS-transformer can be used.
  *  Also prepare for the worst case with preprocessing algorithms. */
  std::streamsize read = m_source->readBlock(
      &(*to)[0], static_cast<std::streamsize>(
          m_blockSize - 1 - m_preprocessingOptions.size()*5));
  if (!read) return 0;
  
  size_t preprocessedSize = preprocess(&(*to)[0], read);
  
  return m_blockManager->makeBlock(to, stats, preprocessedSize);
}

#define PREPROCESS(Type, verb, src, dst) \
  Type r((verb)); \
  r.analyseData((src), length); \
  r.finishAnalysation(); \
  r.decideReplacements(); \
  size_t hSize = r.writeHeader((dst)); \
  size_t comprSize = r.writeReplacedVersion((src), length, (dst)+hSize); \
  length = comprSize + hSize
  

size_t Preprocessor::preprocess(byte *src, size_t length) {
  PROFILE("Preprocessor::preprocess");
  std::vector<byte> tmp;
  tmp.resize(length + 5);
  byte *dst = &tmp[0];
  for(size_t i = 0; i < m_preprocessingOptions.size(); ++i) {
    char c = m_preprocessingOptions[i];
    if(c == 'p') {
      PREPROCESS(PairReplacer, verbosity > 1, src, dst);
    } else if(c == 'r') {
      PREPROCESS(RunReplacer, verbosity > 1, src, dst);
    } else if(c == 'c') {
      PREPROCESS(pairs_and_runs::PairAndRunReplacer, verbosity > 1, src, dst);
    }
    if(length + 5 > tmp.size()) tmp.resize(length + 5);
    std::swap(src, dst);
  }
  if(m_preprocessingOptions.size() & 1) {
    std::copy(src, src + length, dst);
  }
  return length;
}

#undef PREPROCESS

/*#################### Preprocessing algorithms #############################*/

/* Common utility-functions for preprocessing algorithms */

/* In C++0x these two could be implemented more naturally with the use of
 * lambda-functions. */
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
 * @param freqs pointer to FrequencyTable-object
 * @return result is stored into replaceable_pairs
 */
void findReplaceablePairs(std::vector<std::pair<uint16,uint32> >*
                          replaceable_pairs,
                          std::pair<uint16, uint32>* pair_freqs,
                          FrequencyTable *freqs)
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
    if(freqs->getFrequency(current_symbol) + 3 >= pair_freqs[current_pair].second) {
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
uint32 escapeCharIndex(FrequencyTable* freqs,
                       const std::vector<std::pair<uint16, uint32> >&
                       suitable_pairs,
                       uint32 free_symbols)
{
  if (suitable_pairs.size() <= free_symbols) return free_symbols;
  int64 utility = 0; 
  uint32 i; 
  for(i = free_symbols; i < suitable_pairs.size(); ++i) {
    utility += (suitable_pairs[i].second - freqs->getFrequency(i) - 3);
  }
  /* Condition (p2) */
  while(utility <= static_cast<int64>(freqs->getFrequency(i)) &&
        i > free_symbols)
  {
    --i;
    byte fst = static_cast<byte>((suitable_pairs[i].first & 0xFF00) >> 8);
    byte snd = static_cast<byte>(suitable_pairs[i].first & 0x00FF);
    freqs->increase(fst, suitable_pairs[i].second);
    freqs->increase(snd, suitable_pairs[i].second);
    utility -= (suitable_pairs[i].second - freqs->getFrequency(i) - 3);
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
 * @param freqs {@link FrequencyTable} containing frequencies of characters.
 * @param escaping true if escaping is in use.
 * @param free_symbols Number of symbols which aren't present in the source.
 * @return Result is written to replacements.
 */
void constructReplacementTable(uint32 symbols, byte *replacements,
                               const std::vector<std::pair<uint16, uint32> >&
                               replaceable_pairs, const FrequencyTable& freqs,
                               bool escaping, uint32 free_symbols)
{
  assert(symbols > 0);
  byte common_byte = replacements[0];
  uint32 limit = (escaping)?symbols-1:symbols;
  for(uint32 i = 0; i < limit; ++i) {
    replacements[replaceable_pairs[i].first] = freqs.getKey(i);
  }
  byte escape_byte = freqs.getKey(symbols-1);
  /* Prepare escaped characters for writing of the replacements */
  for(uint32 i = free_symbols; i < symbols; ++i) {
    uint16 pair_value = freqs.getKey(i) << 8;
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

  FrequencyTable freqs(freq);
  uint32 free_symbols = 0;
  while(freqs.getFrequency(free_symbols) == 0) ++free_symbols;

  std::vector<std::pair<uint16, uint32> > replaceable_pairs;
  findReplaceablePairs(&replaceable_pairs, pair_freq, &freqs);

  uint32 escape_index = free_symbols;
  if(replaceable_pairs.size() > free_symbols) {
    escape_index = escapeCharIndex(&freqs, replaceable_pairs, free_symbols);
  }
  byte common_byte = freqs.getKey(255);
  byte escape_byte;
  if(escape_index > free_symbols) escape_byte = freqs.getKey(escape_index);
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

  size_t maxUtility() const {
    triple max = seqs_[0];
    return (max.length-1)*max.frequency;
  }

  bool empty() const { return last_ < 0; }

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
                         FrequencyTable *freqs)
{
  SequenceHeap seq_heap(*runs);
  assert(longest_runs->size() == 0);

  uint32 current_symbol = 0;
  while(longest_runs->size() < std::min(static_cast<size_t>(254),runs->size()))
  {
    triple best = seq_heap.removeMax();
    freqs->decrease(best.symbol, best.length*best.frequency);
    if(freqs->getFrequency(current_symbol) + 3 >= (best.length - 1)*(best.frequency)) {
      freqs->increase(best.symbol, best.length*best.frequency);
      break;
    }
    longest_runs->push_back(best);
    ++current_symbol;
  }
  assert(current_symbol == longest_runs->size());
}

uint32 escapeCharIndex(FrequencyTable *freqs, const std::vector<triple>& runs,
                         uint32 free_symbols)
{
  if(runs.size() <= free_symbols) return free_symbols;
  int64 utility = 0;
  uint32 i;
  for(i = free_symbols; i < runs.size(); ++i) {
    utility += ((runs[i].length - 1)*runs[i].frequency - freqs->getFrequency(i) - 3);
  }
  while(utility <= static_cast<int64>(freqs->getFrequency(i)) &&
        i > free_symbols)
  {
    --i;
    freqs->increase(runs[i].symbol, (runs[i].length - 1)*runs[i].frequency);
    utility -= ((runs[i].length - 1)*runs[i].frequency - freqs->getFrequency(i)- 3);
  }
  return i;
}

struct runlist_elem {
  runlist_elem(int l, byte s, int n) : length(l), symbol(s), next_elem(n) {}

  bool operator<(const runlist_elem& other) const {
    return next_elem < other.next_elem || (next_elem == other.next_elem && length < other.length);
  }

  int length;
  byte symbol;
  int next_elem;
};

class ReplacementTable {

 public:
  ReplacementTable(byte escape_byte) {
    rpls_.push_back(runlist_elem(1, escape_byte, -1));
    std::fill(table_, table_ + 256, -1);
  }

  bool empty(byte key) {
    return table_[key] == -1;
  }

  void pushBack(byte run_symbol, uint32 length, byte replacement) {
    rpls_.push_back(runlist_elem(length, replacement, run_symbol));
  }

  void addEscaping(byte symbol) {
    table_[symbol] = -2;
  }

  bool isEscaped(byte symbol) const {
    return table_[symbol] == -2;
  }

  void prepare() {
    std::sort(rpls_.begin(), rpls_.end());
    for(size_t i = 0; i < rpls_.size(); ++i) {
      byte symbol = (byte)rpls_[i].next_elem;
      assert(table_[symbol] != -2);
      int old = table_[symbol];
      table_[symbol] = i;
      rpls_[i].next_elem = old;
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

size_t writeRunReplacement(const ReplacementTable& repl, uint32 run_length,
                           byte escape, byte symbol, byte *to)
{
  size_t j = 0;
  if(repl.isEscaped(symbol)) {
    for(size_t i = 0; i < run_length; ++i) {
      to[j++] = escape; to[j++] = symbol;
    }
    return j;
  }
  int tbl_index = repl.listBegin(symbol);
  do {
    if (tbl_index == -1) {
      std::fill(to + j, to + j + run_length, symbol);
      return run_length + j;
    }
    const runlist_elem& el = repl.listElement(tbl_index);
    uint32 times = (el.length == 0)?0:run_length/el.length;
    if (times > 0) {
      std::fill(to + j, to + j + times, el.symbol);
      j += times;
      run_length -= times*el.length;
    }
    tbl_index = el.next_elem;
  } while(run_length > 0);
  return j;
}

uint64 writeReplacements(const ReplacementTable& replacements, byte *to,
                         byte *from, size_t length, byte escape)
{
  size_t j = 0; /* Index of target */
  byte prev = from[0];
  for(size_t i = 1; i < length; ++i) {
    if(prev != from[i]) {
      if(replacements.isEscaped(prev)) to[j++] = escape;
      to[j++] = prev;
    } else {
      uint32 run_length = 2;
      ++i;
      while(i < length && prev == from[i]) { ++i; ++run_length;}
      j += writeRunReplacement(replacements, run_length, escape,
                               prev, to + j);
    }
    prev = from[i];
  }
  if(prev != from[length - 2]) {
    if(replacements.isEscaped(prev)) to[j++] = escape;
    to[j++] = prev;
  }
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

  FrequencyTable freqs(freq);

  uint32 free_symbols = 0;
  while(freqs.getFrequency(free_symbols) == 0) ++free_symbols;

  std::vector<triple> runs, longest_runs;
  runs.reserve(256);
  for(size_t i = 0; i < 256; ++i) {
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
  byte escape_byte = freqs.getKey(escape_index);

  if(symbols_in_use > 0) {
    //TODO: Store replacement in somewhere so they can be merged with other
    uint32 limit = run_replacements & 0xFFFFFFFE;

    for(uint32 i = 0; i < limit; i += 2) {
      temp[position++] = freqs.getKey(i);
      assert(longest_runs[i].length <= kMaxLenOfSeq);
      assert(longest_runs[i+1].length <= kMaxLenOfSeq);
      byte lengths = (utils::logFloor(longest_runs[i].length) << 4) |
          utils::logFloor(longest_runs[i+1].length);
      temp[position++] = lengths;
      temp[position++] = longest_runs[i].symbol;
      temp[position++] = freqs.getKey(i+1);
      temp[position++] = longest_runs[i+1].symbol;
    }
    byte sentinel = (escape_index != free_symbols) ?
        escape_byte : freqs.getKey(symbols_in_use - 1);
    
    if( run_replacements != limit )
    {
      assert(limit == symbols_in_use - 1 || limit == symbols_in_use - 2 );
      temp[position++] = freqs.getKey(limit);
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
  /* Escaped characters*/
  if( new_symbols > 0) {
    for(uint32 i = free_symbols; i <= escape_index; ++i) {
      replacements.addEscaping(freqs.getKey(i));
      //replacements.pushBack(freqs.key(i), 1, escape_byte);
    }
  }
  for (uint32 i = 0; i < run_replacements; ++i) {
    assert( i < longest_runs.size() );
    replacements.pushBack(longest_runs[i].symbol, longest_runs[i].length,
                          freqs.getKey(i));
  }
  replacements.prepare();
  uint64 total_size = position;
  total_size += writeReplacements(replacements, temp + position, from, length,
                                  escape_byte);
  std::copy(temp, temp + total_size, from);
  delete [] temp;
  assert(total_size <= length + 2);
  return total_size;
}

namespace pr {

void gatherFrequencies(byte *from, size_t length, size_t *freq,
                       std::pair<uint16, size_t> *pairFreqs,
                       std::vector<size_t> *runFreqs)
{
  byte prev = from[0];
  ++freq[prev];
  for(size_t i = 1; i < length; ++i) {
    if(prev != from[i]) {
      ++pairFreqs[(prev << 8) | from[i]].second;
      prev = from[i];
    } else {
      size_t runLength = 2;
      do {
        ++freq[prev];
        ++i;
      } while(i < length && prev == from[i] && runLength <= longruns::kMaxLenOfSeq);
      longruns::updateFreqs(runFreqs, prev, runLength);
      if(i < length) prev = from[i];
    }
    ++freq[prev];
  }
}

class Replacements {
 public:
  Replacements(byte escape) : runReplacements(escape), escapeSymbol(escape) {
    for(size_t i = 0; i < (1 << 16); ++i) pairReplacements[i] = escapeSymbol;
    for(size_t i = 0; i < 256; ++i) escapedCharacters[i] = false;
  }

  void addReplaceablePair(uint16 pair, byte symbol) {
    pairReplacements[pair] = symbol;
  }

  void addEscape(byte symbol) { escapedCharacters[symbol] = true; }

  bool isEscaped(byte symbol) const { return escapedCharacters[symbol]; }

  void addReplaceableRun(byte runSymbol, size_t length, byte replacement) {
    runReplacements.pushBack(runSymbol, length, replacement);
  }
  
  void prepare() {
    runReplacements.prepare();
  }
  
 private:
  longruns::ReplacementTable runReplacements;
  /**If there isn't replacement for the pair, the value is escapeSymbol.
   * Otherwise the value is the replacement symbol. */
  byte pairReplacements[1 << 16];
  bool escapedCharacters[256];
  byte escapeSymbol;
};

void findReplaceablePairsAndRuns(FrequencyTable& freqs,
                                 std::vector<longruns::triple>& allRuns,
                                 std::vector<longruns::triple>& replaceableRuns,
                                 std::pair<uint16, size_t>* pairFreqs,
                                 std::vector<std::pair<uint16, size_t> >& replaceablePairs)
{
  // 1 if character used in run replacements, 2 if in pair repl as the first
  // character and 4 if in pair repl as the secind character
  byte alphabetPartition[256] = {0};

  // Should we sort this only partially?
  std::sort(pairFreqs, pairFreqs + (1 << 16), comparePairSecondDesc<uint16, size_t>);
  int pairPtr = 0;

  longruns::SequenceHeap seqHeap(allRuns);
  int currentSymbol = 0;
  bool searchPair = true, searchRuns = true;
  while(currentSymbol < 254 && (searchPair || searchRuns)) {
    // TODO: fix the formula to take written metadata in to account
    if(seqHeap.maxUtility() > pairFreqs[pairPtr].second && searchRuns) {
      longruns::triple best = seqHeap.removeMax();
      if(alphabetPartition[best.symbol] != 1 && alphabetPartition[best.symbol] != 0) continue;

      freqs.decrease(best.symbol, best.length*best.frequency);

      if(freqs.getFrequency(currentSymbol) + 3 >= (best.frequency-1)*best.length) {
        freqs.increase(best.symbol, best.length*best.frequency);
        searchRuns = false;
        continue;
      }

      alphabetPartition[best.symbol] = 1;
      replaceableRuns.push_back(best);
    } else if (searchPair) {
      byte fst = static_cast<byte>((pairFreqs[pairPtr].first & 0xFF00) >> 8);
      byte snd = static_cast<byte>(pairFreqs[pairPtr].first & 0x00FF);
      byte &fstStatus = alphabetPartition[fst], &sndStatus = alphabetPartition[snd];

      if(fst == snd || (fstStatus != 2 && fstStatus != 0) ||
         (sndStatus != 3 && sndStatus != 0) )
      {
        ++pairPtr;
        continue;
      }

      if(!freqs.decrease(fst, pairFreqs[pairPtr].second)) {
        ++pairPtr; continue;
      }
      if(!freqs.decrease(snd, pairFreqs[pairPtr].second)) {
        freqs.increase(fst, pairFreqs[pairPtr].second);
        ++pairPtr; continue;
      }
      // Check if it is wise to make replacement
      if(freqs.getFrequency(currentSymbol) + 3 >= pairFreqs[pairPtr].second) {
        freqs.increase(fst, pairFreqs[pairPtr].second);
        freqs.increase(snd, pairFreqs[pairPtr].second);
        searchPair = false;
      }
      
      fstStatus = 2; sndStatus = 3;
      replaceablePairs.push_back(pairFreqs[pairPtr]);
      ++pairPtr;
    }
    ++currentSymbol;
  }
}

size_t findEscapeIndex(FrequencyTable& freqs,
                       std::vector<std::pair<uint16, size_t> >& pairs,
                       std::vector<longruns::triple>& runs)
{
  size_t escapeIndex = pairs.size() + runs.size();
  while(true) {
    if(escapeIndex == 0) break;
    const longruns::triple& tr = runs.back();
    int runUtil = (tr.length - 1)*tr.frequency;
    if(runUtil < pairs.back().second) {
      if(runUtil >= freqs.getFrequency(escapeIndex) + 3) break;
      else {
        runs.pop_back();
        freqs.increase(tr.symbol, runUtil + tr.frequency);
      }
    } else {
      if(pairs.back().second >= freqs.getFrequency(escapeIndex) + 3) break;
      else {
        freqs.increase(pairs.back().first & 0xff, pairs.back().second);        
        freqs.increase((pairs.back().first >> 8)& 0xff, pairs.back().second);
        pairs.pop_back();
      }
    }
    --escapeIndex;
  }
  return escapeIndex;
}

} //namespace pr

size_t compressPairsAndRuns(byte *from, size_t length) {
  using namespace pr;
  assert(length > 0);
  assert(from);

  // Initialize vector for runs
  std::vector<size_t> runFreqs[256];
  size_t vsize = utils::logFloor(std::min(length, longruns::kMaxLenOfSeq));
  for(size_t i = 0; i < 256; ++i) {
    runFreqs[i].resize(vsize-1);
    std::fill(runFreqs[i].begin(), runFreqs[i].end(), 0);
  }

  //Initialize array for pairs
  std::pair<uint16, size_t> pairFreqs[1 << 16];
  for(size_t i = 0; i < 256; ++i) {
    for(size_t j = 0; j < 256; ++j) {
      uint16 p = (i << 8) | j;
      pairFreqs[p] = std::make_pair(p, 0);
    }
  }

  size_t freq[256] = {0};
  gatherFrequencies(from, length, freq, pairFreqs, runFreqs);

  FrequencyTable freqs(freq);
  size_t freeSymbols = 0;
  while(freq[freeSymbols++] == 0) ;

  std::vector<longruns::triple> allRuns, replaceableRuns;
  allRuns.reserve(256);
  for(size_t i = 0; i < 256; ++i) {
    for(size_t j = 0; j < runFreqs[i].size(); ++j) {
      if(runFreqs[i][j] > 0) {
        allRuns.push_back(longruns::triple(static_cast<byte>(i), 1 << (j+1), runFreqs[i][j]));
      }
    }
  }

  std::vector<std::pair<uint16, size_t> > replaceablePairs;
  findReplaceablePairsAndRuns(freqs, allRuns, replaceableRuns,
                              pairFreqs, replaceablePairs);

  size_t candidates = replaceablePairs.size() + replaceableRuns.size();
  size_t escapeIndex = freeSymbols;
  if(candidates > freeSymbols) {
    escapeIndex = findEscapeIndex(freqs, replaceablePairs, replaceableRuns);
    candidates = replaceablePairs.size() + replaceableRuns.size();
  } else {
    escapeIndex = 255;
  }

  if(verbosity > 1) {
    std::clog << "Replacing " << replaceablePairs.size() << " pairs and "
              << replaceableRuns.size() << " runs. ";
    if(escapeIndex > freeSymbols) {
      std::clog << "Made " << (escapeIndex - freeSymbols + 1) << " symbols free.\n";
    } else {
      std::clog << "No symbol made free\n.";
    }
  }
  
  byte escape = freqs.getKey(escapeIndex);

  //TODO: write metadata for replacements
  byte *temp = new byte[length + 2];
  size_t position = 0;

  Replacements replacements(escape);
  int repIndex = 0;
  for(size_t i = 0; i < replaceablePairs.size(); ++i, ++repIndex) {
    replacements.addReplaceablePair(replaceablePairs[i].first, freqs.getKey(repIndex));
  }
  for(size_t i = 0; i < replaceableRuns.size(); ++i, ++repIndex) {
    longruns::triple& t = replaceableRuns[i];
    replacements.addReplaceableRun(t.symbol, t.length, freqs.getKey(repIndex));
  }
  assert(escapeIndex == 255 || escapeIndex == repIndex);
  if(escapeIndex != 255) {
    for(size_t i = freeSymbols; i <= escapeIndex; ++i)
      replacements.addEscape(freqs.getKey(i));
  }  
}


} //namespace bwtc
