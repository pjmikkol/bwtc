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

namespace {

/* In C++0x these two could be implemented more naturally with the use of
 * lambda-functions. */
template <typename F, typename S>
bool ComparePairSecondAsc(std::pair<F,S> p1, std::pair<F,S> p2) {
  return (p1.second < p2.second);
}

template <typename F, typename S>
bool ComparePairSecondDesc(std::pair<F,S> p1, std::pair<F,S> p2) {
  return (p1.second > p2.second);
}

}

template<typename Key, typename Value>
void InitPairsWithValue(std::pair<Key, Value> *pairs, Value v, uint64 length)
{
  for(uint64 i = 0; i < length; ++i)
    pairs[i] = std::make_pair(static_cast<Key>(i), v);
}

void ComputePairFrequencies(byte *data, std::pair<byte, uint64>* freqs,
                            std::pair<uint16, uint32>* pair_freqs, uint64 len)
{
  uint16 index = data[0];
  ++freqs[index].second;
  for (uint64 i = 1; i < len; ++i) {
    ++freqs[data[i]].second;
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
 * free_symbols       - count of pairs in above array where frequency == 0   *
 *****************************************************************************/
void FindReplaceablePairs(std::vector<std::pair<uint16,uint32> >*
                          replaceable_pairs,
                          std::pair<uint16, uint32>* pair_freqs,
                          const std::pair<byte, uint64>* frequencies)
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
    } /* Conditions (p1) and (p2) */
    if(frequencies[current_symbol].second+3 >= pair_freqs[current_pair].second)
      break; /* We won't benefit from any changes any more*/

    /* Reject pairs which have conflicting symbols
     * This one NEEDS better heuristics. Now it uses greedy heuristic */
    byte fst = static_cast<byte>((pair_freqs[current_pair].first & 0xFF00) >> 8);
    byte snd = static_cast<byte>(pair_freqs[current_pair].first & 0x00FF);
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
          break;
        }
      }
      if(valid) {
        replaceable_pairs->push_back(pair_freqs[current_pair]);
        ++current_symbol;
      }
    }
    ++current_pair;
  }
  assert(current_symbol == replaceable_pairs->size());
}

/* Returns the index for the 'escape_char' in freqs-array or value of
 * free_symbols if freeing of the chars is not profitable */
unsigned EscapeCharIndex(const std::pair<byte, uint64>* freqs,
                         const std::vector<std::pair<uint16, uint32> >&
                         suitable_pairs,
                         unsigned free_symbols)
{
  if (suitable_pairs.size() <= free_symbols) return free_symbols;
  int64 utility = 0; 
  unsigned i; 
  for(i = free_symbols; i < suitable_pairs.size(); ++i) {
    utility += (suitable_pairs[i].second - freqs[i].second - 3);
  }

  while(utility <= static_cast<int64>(freqs[i].second) &&
        i > free_symbols)
  {
    --i;
    utility -= (suitable_pairs[i].second - freqs[i].second - 3);
  }
  return i;
}

void WriteBytes(uint16 value, byte *address) {
  *address = (value >> 8);
  *(address + 1) = value & 0xFF;
}

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

uint64 CompressCommonPairs(byte *from, uint64 length)
{
  assert(length > 0);
  assert(from);

  std::pair<byte, uint64> freq[256];
  std::pair<uint16, uint32> pair_freq[65536];
  InitPairsWithValue<byte, uint64>(freq, 0, 256);
  InitPairsWithValue<uint16, uint32>(pair_freq, 0, 65536);
  ComputePairFrequencies(from, freq, pair_freq, length);

  std::sort(&freq[0], &freq[0] + 256, ComparePairSecondAsc<byte, uint64>);

  unsigned free_symbols = 0;
  while(freq[free_symbols].second == 0) ++free_symbols;

  std::vector<std::pair<uint16, uint32> > replaceable_pairs;
  FindReplaceablePairs(&replaceable_pairs, pair_freq, freq);

  unsigned escape_index = free_symbols;
  if(replaceable_pairs.size() > free_symbols) {
    escape_index = EscapeCharIndex(freq, replaceable_pairs, free_symbols);
  }

  byte common_byte = (replaceable_pairs.size() > 0) ?
      static_cast<byte>(replaceable_pairs[0].first >> 8) : 0;
  byte escape_byte;
  if(escape_index > free_symbols) escape_byte = freq[escape_index].first;
  else escape_byte = common_byte;

  byte replacements[65536];
  std::fill(replacements, replacements + 65536, common_byte);
  byte *temp = new byte[length + 2];
  /* Initialize replacement table and write header */
  unsigned symbols_in_use = 0;
  uint64 position = 0;
  unsigned candidates = static_cast<unsigned>(replaceable_pairs.size());
  unsigned k;
  for(k = 0; k < std::min(free_symbols, candidates); ++k) {
    replacements[replaceable_pairs[k].first] = freq[k].first;
    temp[position++] = freq[k].first;
    WriteBytes(replaceable_pairs[k].first, temp + position);
    position += 2;
  }
  if (k > 0) symbols_in_use = k - 1;
  else symbols_in_use = 0;

  if ( free_symbols < escape_index) {
    for(unsigned i = free_symbols; i < escape_index; ++i) {
      replacements[replaceable_pairs[i].first] = freq[i].first;
      temp[position++] = freq[i].first;
      WriteBytes(replaceable_pairs[i].first, temp + position);
      position += 2;
      uint16 pair_value = freq[i].first << 8;
      for(unsigned j = 0; j < 256; ++j, ++pair_value)
        replacements[pair_value] = escape_byte;
    }
    symbols_in_use += escape_index - free_symbols + 1;
    uint16 pair_value = escape_byte << 8;
    for(unsigned j = 0; j < 256; ++j, ++pair_value)
      replacements[pair_value] = escape_byte;
  }
  byte dummy = escape_byte + 0;
  if (escape_index > free_symbols) dummy = freq[escape_index - 1].first;
  else if (symbols_in_use) dummy = freq[symbols_in_use - 1].first;

  temp[position++] = dummy;
  if (free_symbols < escape_index) temp[position++] = escape_byte;
  else temp[position++] = dummy;

  unsigned new_symbols = std::max(0U, symbols_in_use - free_symbols);

  if (verbosity > 1) {
    std::clog << "Replacing " << symbols_in_use - 1 << " pairs. ";
    if (new_symbols > 0)
      std::clog << "Made " << new_symbols << " symbols free.\n";
    else
      std::clog << "No symbols made free.\n";
  }


  uint64 total_size = position;
  total_size += WriteReplacements(replacements, temp + position, from, length,
                                  common_byte, escape_byte);
  assert(total_size <= length + 2);
  std::copy(temp, temp + total_size, from);
  delete [] temp;
  return total_size;
}

/*#################### Preprocessing algorithms #############################*/

  /**************************************************************************
   *    Reasoning behind finding the replaceable pairs   .                  *
   *------------------------------------------------------------------------*
   *(p1)                                                                    *
   * For symbols which are not present in original data, it is required     *
   * that the pair it will replaced appears at least 4 times in original    *
   * data, because otherwise we won't benefit from the replacement.         *
   *                                                                        *
   *(p2)                                                                    *
   * We replace some bytes which are present in original data with two      *
   * bytes. Idea is to use these replaced bytes as coding symbols.          *
   * Denote the frequency of symbol 'x' with f_x and denote the  frequency  *
   * of pair 'p' with P_p. If we are going to replace 'x' with two symbols  *
   * and use it to represent 'p' we require                                 *
   *        2*f_x + P_p < f_x + 2*P_p    <=>     f_x < P_p                  *
   * since otherwise we wouldn't  benefit from the replacement.             *
   *                                                                        *
   *(p3)                                                                    *
   * In addition to (2) we need to satisfy the following condition:         *
   * Let x_j, x_j+1, ..., x_k be the frequencies of freed symbols 'x_i'.    *
   * Let P_j, P_j+1, ..., P_k be the frequencies of pairs 'P_i' which       *
   * are going to be replaced with 'x_i'. Let 'x' be the special symbol     *
   * which is going to be used in each replacements for 'x_i'. We need      *
   * also use 2 symbols for representing 'x' so the following inequality    *
   * must hold for replacements. Denote the frequency of 'x' with x.        *
   *  sum from i = j to k : [2*x_i + P_i] + 2*x <                           *
   *  sum from i = j to k : [x_i + 2*P_i] +   x                             *
   *                 <=>                                                    *
   *  sum from i = j to k : [P_i - x_i] > x                                 *
   **************************************************************************/
#if 0
/* Empty namespace for implementing the common functions needed in different
 * preprocessing algorithms. */
namespace {

/* In C++0x these two could be implemented more naturally with the use of
 * lambda-functions. */
template <typename F, typename S>
bool ComparePairSecondAsc(std::pair<F,S> p1, std::pair<F,S> p2) {
  return (p1.second < p2.second);
}

template <typename F, typename S>
bool ComparePairSecondDesc(std::pair<F,S> p1, std::pair<F,S> p2) {
  return (p1.second > p2.second);
}

/*****************************************************************************
 * FindReplaceablePairs - Finds the candidates for pairs which are going to  *
 *                        be replaced by single symbols. It is required that *
 *                        if pair = <p1,p2> is in result, then there is no   *
 *                        pair where p2 would be first symbol of pair or     *
 *                        where p1 would be second symbol.                   *
 *                                                                           *
 * Takes the following arguments:                                            *
 * replaceable_pairs  - empty vector where pairs are stored                  *
 * pair_freqs         - unsorted array of <pair, frequency>-pairs, size of   *
 *                      65536                                                *
 * frequencies        - array of <byte, frequency>-pairs sorted by frequency *
 *                      in descending order                                  *
 * free_symbols       - count of pairs in above array where frequency == 0   *
 *****************************************************************************/
void FindReplaceablePairs(std::vector<uint16>* replaceable_pairs,
                          std::pair<uint16, uint32>* pair_freqs,
                          const std::pair<byte, uint64>* frequencies)
{
  const unsigned kStep = 256;
  unsigned current_pair = 0, current_symbol = 0;
  unsigned limit = 0;

  while(replaceable_pairs->size() < 254) {
    if(current_pair == limit) {
      limit += kStep;
      assert(limit < 65536);
      std::partial_sort(pair_freqs + current_pair, pair_freqs + limit,
                        pair_freqs + 65536,
                        ComparePairSecondDesc<uint16, uint32>);
    } /* Conditions (p1) and (p2) */
    if(frequencies[current_symbol].second + 3 >= pair_freqs[current_pair].second)
      break; /* We won't benefit from any changes any more*/

    /* Reject pairs which have conflicting symbols
     * This one NEEDS better heuristics. Now it uses greedy heuristic */
    byte fst = static_cast<byte>((pair_freqs[current_pair].first & 0xFF00) >> 8);
    byte snd = static_cast<byte>(pair_freqs[current_pair].first & 0x00FF);
    if (fst != snd) { /* Do not approve pairs of the same char */
      bool valid = true;
      for(std::vector<uint16>::iterator it = replaceable_pairs->begin();
          it != replaceable_pairs->end(); ++it)
      {
        uint16 current_fst = static_cast<byte>(((*it) & 0xFF00) >> 8);
        uint16 current_sec = static_cast<byte>(((*it) & 0x00FF));
        if (current_fst == snd || current_sec == fst) {
          valid = false;
          break;
        }
      }
      if(valid) {
        replaceable_pairs->push_back(pair_freqs[current_pair].first);
        ++current_symbol;
      }
    }
    ++current_pair;
  }
  assert(current_symbol == replaceable_pairs->size());
}

} //empty namespace

/* Assumes that from-array has in reality length of length+2, since in the *
 * worst case we increase the size of the input by two                     */
uint64 CompressCommonPairs(byte *from, uint64 length)
{
  assert(length > 0);
  assert(from);

  /* First compute frequencies of characters and pairs */
  std::pair<byte, uint64> frequencies[256];
  for(unsigned i = 0; i < 256; ++i) frequencies[i] = std::make_pair(i, 0);
  std::pair<uint16, uint32> pair_freqs[65536];
  for(unsigned i = 0; i < 65536; ++i) pair_freqs[i] = std::make_pair(i, 0);

  uint16 index = from[0];
  ++frequencies[from[0]].second;
  for (uint64 i = 1; i < length; ++i) {
    ++frequencies[from[i]].second;
    index <<= 8;
    index |= from[i];
    ++pair_freqs[index].second;
  }

  std::sort(&frequencies[0], &frequencies[0] + 256,
            ComparePairSecondAsc<byte, uint64>);

  unsigned free_symbols = 0;
  while(frequencies[free_symbols].second == 0) ++free_symbols;

  std::vector<uint16> replaceable_pairs;
  FindReplaceablePairs(&replaceable_pairs, pair_freqs, frequencies);

  byte escape_char;
  unsigned new_symbols;
  if (replaceable_pairs.size() > free_symbols) {
    /* We have found at least one symbol worth of freeing. *
     * Find the suitable symbol for 'x' in (p3). */
    // TODO: Find out if this could be optimised by calculating the sum in
    //       FindReplaceablePairs

  int64 utility = 0; 
  unsigned i;
  for(i = free_symbols; i < replaceable_pairs.size(); ++i)
    utility += (pair_freqs[i].second - frequencies[i].second - 3);
  while(utility <= static_cast<int64>(frequencies[i].second) &&
        i > free_symbols)
  {
    --i;
    utility -= (pair_freqs[i].second - frequencies[i].second - 3);
  }

  /*
  int64 check_sum = 0;
    unsigned i;
    for (i = free_symbols; i < replaceable_pairs.size(); ++i) 
      check_sum += (pair_freqs[i].second - frequencies[i].second + 3);
    while(check_sum >= static_cast<int64>(frequencies[i].second) &&
          i > free_symbols) {
      --i;
      check_sum -= (pair_freqs[i].second - frequencies[i].second + 3);
    }
  */  //    if (i >= free_symbols) {
      escape_char = frequencies[i].first;
      new_symbols = i - free_symbols;
      // }
      //else {
      //new_symbols = 0;
      //}
  }
  else {
    new_symbols = 0;
  }
  unsigned symbols_in_use = std::min(
      free_symbols + new_symbols, static_cast<unsigned>(
          replaceable_pairs.size()));
  if (verbosity > 1) {
    std::clog << "Replacing " << symbols_in_use << " pairs.\n";
    if (new_symbols > 0)
      std::clog << "Made " << new_symbols + 1 << " symbols free.\n";
    else
      std::clog << "No symbols made free.\n";
  }
  /* Store the replacements in map<byte, uint16> where the latter is the
   * pair to be replaced */
  std::map<byte, uint16> pairs;
  
  /* failbyte can't be any of the symbols used in replacements so it
   * is a good choice for signaling "no replacement for this pair" */
  byte failbyte = 0;
  if(replaceable_pairs.size() > 0)
    failbyte = static_cast<byte>(replaceable_pairs[0] & 0x00FF);

  byte replacements[65536];
  std::fill(replacements, replacements + 65536, failbyte);
  if(new_symbols == 0) escape_char = failbyte; /* No need for escaping */
  /*else {
    uint16 pair_value = (escape_char << 8);
    for(unsigned j = 0; j < 256; ++j, ++pair_value) {
      replacements[pair_value] = escape_char;
    }
  }*/
  for(unsigned i = 0; i < std::min(free_symbols,symbols_in_use); ++i) {
    /* Here replaceable_pairs[i] is the pair which will be replaced with
     * frequencies[i] */
    replacements[replaceable_pairs[i]] = frequencies[i].first;
    pairs[frequencies[i].first] = replaceable_pairs[i];
  }
  for(unsigned i = free_symbols; i < symbols_in_use; ++i) {
    /* Mark every pair where the first byte is one of the new symbols with
     * escape_char*/
    replacements[replaceable_pairs[i]] = frequencies[i].first;
    pairs[frequencies[i].first] = replaceable_pairs[i];
    uint16 pair_value = (frequencies[i].first << 8);
    for(unsigned j = 0; j < 256; ++j, ++pair_value) {
      replacements[pair_value] = escape_char;
    }
  }
  /* Allocate array for result. Prepare for the worst */
  byte *temp = new byte[length+2];

  /* Write result to temporary-array */
  uint64 result_index = 0;
  /* First write replacements */
  std::pair<byte, uint16> curr_pair;
  for(std::map<byte, uint16>::const_iterator it = pairs.begin();
      it != pairs.end(); ++it) {
    curr_pair = *it;
    temp[result_index++] = curr_pair.first;
    temp[result_index++] = static_cast<byte>((curr_pair.second & 0xFF00) >> 8);
    temp[result_index++] = static_cast<byte>(curr_pair.second & 0x00FF);
  }
  /* End of pair replacements- list */
  temp[result_index++] = curr_pair.first;
  /* Write the "escape char" from character replacements or end of *
   * sequence. */
  if (new_symbols > 0) temp[result_index++] = escape_char;
  else temp[result_index++] = curr_pair.first;

  uint16 pair = static_cast<uint16>(from[0]);
  uint64 i = 1;
   /*   while(1) {
    pair <<= 8;
    pair |= from[i];
    if(replacements[pair] == failbyte) {
      temp[result_index++] = from[i-1];
    }
    else if (replacements[pair] == escape_char) {
      temp[result_index++] = escape_char;
      temp[result_index++] = from[i-1];
    }
    else { // pair will be replaced
      temp[result_index++] = replacements[pair];
      if( i == length - 1) break;
      pair = static_cast<byte>(from[++i]);
    }

    if ( i >= length - 1) {
      assert(i == length - 1);
      if(from[i] != escape_char && escape_char != failbyte)
      temp[result_index++] = from[i];
      break;
    }
    ++i;
  }*/
    
  for(i = 1; i < length; ++i, ++result_index) {
    pair <<= 8;
    pair |= from[i];
    if(replacements[pair] == failbyte) {
      temp[result_index] = from[i-1];
    }
    else if (replacements[pair] == escape_char) {
      temp[result_index++] = escape_char;
      temp[result_index] = from[i-1];
    }
    else { // pair will be replaced
      temp[result_index] = replacements[pair];
      pair = static_cast<byte>(from[++i]);
    }
  }
  if (replacements[pair] == failbyte || replacements[pair] == escape_char) {
    pair <<= 8;
    if (replacements[pair] == escape_char && escape_char != failbyte) {
      temp[result_index++] = escape_char;
    }
    temp[result_index++] = from[length - 1];  
  }
  assert(result_index <= length + 2);

  std::copy(temp, temp + result_index, from);
  delete [] temp;
  return result_index;
}
#endif

} //namespace bwtc
