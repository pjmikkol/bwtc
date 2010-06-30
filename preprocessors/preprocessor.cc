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
/* In C++0x these could be implemented more naturally with the use of
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

/* Assumes that frequencies is already sorted and the size of frequencies is
 * 256 
void FreeSymbols(std::pair<byte, uint64> *frequencies, byte *data,
                 uint64 length)
{
  
}*/

// TODO: What if we need more than length bytes for result
uint64 CompressCommonPairs(byte *from, uint64 length) {
  assert(length > 0);
  assert(from);

  /* First compute frequencies of characters */
  std::pair<byte, uint64> frequencies[256];
  for(unsigned i = 0; i < 256; ++i) frequencies[i] = std::make_pair(i, 0);
  for (uint64 i = 0; i < length; ++i) ++frequencies[from[i]].second;

  std::sort(&frequencies[0], &frequencies[0] + 256,
            ComparePairSecondAsc<byte, uint64>);

  unsigned free_symbols = 0;
  while(frequencies[free_symbols].second == 0) {
    ++free_symbols;
  }

  /* Calculate the frequencies of pairs. */
  std::pair<uint16, uint32> pair_freqs[65536];
  for(unsigned i = 0; i < 65536; ++i) pair_freqs[i] = std::make_pair(i, 0);

  uint16 index = from[0];
  for(uint64 i = 1; i < length; ++i ) {
    index <<= 8;
    index |= from[i];
    ++pair_freqs[index].second;
  }

  
  /**************************************************************************
   *    Deduce the replaceable pairs from frequency table.                  *
   *------------------------------------------------------------------------*
   * Reasoning behind the choices made in loop:                             *
   *                                                                        *
   *(1)                                                                     *
   * For symbols which are not present in original data, it is required     *
   * that the pair it will replaced appears at least 4 times in original    *
   * data, because otherwise we won't benefit from the replacement.         *
   *                                                                        *
   *(2)                                                                     *
   * We replace some bytes which are present in original data with two      *
   * bytes. Idea is to use these replaced bytes as coding symbols.          *
   * Denote the frequency of symbol 'x' with f_x and denote the  frequency  *
   * of pair 'p' with P_p. If we are going to replace 'x' with two symbols  *
   * and use it to represent 'p' we require                                 *
   *        2*f_x + P_p < f_x + 2*P_p    <=>     f_x < P_p                  *
   * since otherwise we wouldn't  benefit from the replacement.             *
   *                                                                        *
   *(3)                                                                     *
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
   *  sum from i = j to k : [x_i - P_i] < x                                 *
   **************************************************************************/

  /* Step size for sorting the pairs */
  const unsigned kStep = 256; 

  unsigned current_pair = 0, current_symbol = 0;
  unsigned limit = 0;
  std::vector<uint16> replaceable_pairs;
  /* Here some more clever heuristic would be a good idea */ 
  while(replaceable_pairs.size() < 254) {
    if(current_pair == limit) {
      limit += kStep;
      assert(limit < 65536);
      std::partial_sort(&pair_freqs[current_pair], &pair_freqs[limit],
                        &pair_freqs[65536],
                        ComparePairSecondDesc<uint16, uint32>);
    } /* Conditions (1) and (2) */
    if((current_symbol < free_symbols && pair_freqs[current_pair].second < 4) ||
       (current_symbol >= free_symbols && frequencies[current_symbol].second >=
        pair_freqs[current_pair].second)) {
      break; /* We won't benefit from any changes any more*/
    }

    byte fst = static_cast<byte>((pair_freqs[current_pair].first & 0xFF00)>> 8);
    byte snd = static_cast<byte>(pair_freqs[current_pair].first & 0x00FF);
    if (fst != snd) { /* Do not mind pairs of the of the same char */
      bool valid = true;
      for(std::vector<uint16>::iterator it = replaceable_pairs.begin();
          it != replaceable_pairs.end(); ++it)
      {
        uint16 current_fst = static_cast<byte>(((*it) & 0xFF00) >> 8);
        uint16 current_sec = static_cast<byte>(((*it) & 0x00FF));
        if (current_fst == snd || current_sec == fst) {
          valid = false;
          break;
        }
      }
      if(valid) {
        replaceable_pairs.push_back(pair_freqs[current_pair].first);
        ++current_symbol;
      }
    }
    ++current_pair;
  }
  assert(current_symbol == replaceable_pairs.size());
  byte escape_char;
  unsigned new_symbols;
  if (replaceable_pairs.size() > free_symbols) {
    /* We have found at least one symbol worth of freeing. *
     * Find the suitable symbol for 'x' in (3). */
    // TODO: Find out if this could be optimised by calculating the sum in
    //       'if valid'-step in previous loop
    int64 check_sum = 0;
    unsigned i;
    for (i = free_symbols; i < replaceable_pairs.size(); ++i) 
      check_sum += (frequencies[i].second - pair_freqs[i].second);
    while(check_sum >= static_cast<int64>(frequencies[i].second) &&
          i > free_symbols) {
      --i;
      check_sum -= (frequencies[i].second - pair_freqs[i].second);
    }
    if (i > free_symbols) {
      escape_char = frequencies[i].first;
      new_symbols = i - free_symbols;
    }
    else {
      new_symbols = 0;
    }
  }
  else {
    new_symbols = 0;
  }
  unsigned symbols_in_use = std::min(
      free_symbols + new_symbols,static_cast<unsigned>(
          replaceable_pairs.size()));
  if (verbosity > 1) {
    std::clog << "Replacing " << symbols_in_use << " pairs.\n"
              << "Using " << new_symbols << " freed symbols.\n";
  }
  /* Store the replacements in map<byte, uint16> where the latter is the
   * pair to be replaced */
  std::map<byte, uint16> pairs;
  
  /* failbyte can't be any of the symbols used in replacements so it
   * is a good choice for signaling "no replacement for this pair" */
  byte failbyte = 0;
  if(replaceable_pairs.size() > 0)
    failbyte = static_cast<byte>(replaceable_pairs[0] & 0x00FF);

  assert(replaceable_pairs.size() >= symbols_in_use);
  byte replacements[65536];
  std::fill(replacements, replacements + 65536, failbyte);
  if(new_symbols == 0) escape_char = failbyte; /* No need for escaping */
  else replacements[escape_char] = escape_char;
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
  /* Allocate array for result */
  byte *temp = new byte[length];

  /* WHAT IF: there won't be any replacements */
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
  /* Write the "escape char" from character replacements*/
  if (new_symbols > 0) temp[result_index++] = escape_char;
  /* Write the "end of sequence" */
  temp[result_index++] = curr_pair.first;

  uint16 pair = static_cast<uint16>(from[0]);
  for(uint64 i = 1; i < length; ++i, ++result_index) {
    pair <<= 8;
    pair |= from[i];
    if(replacements[pair] == failbyte) {
      temp[result_index] = from[i-1];
    }
    else if (replacements[pair] != escape_char) {
      temp[result_index] = replacements[pair];
      ++i;
      pair = static_cast<byte>(from[i]);
    }
    else { // replacements[pair] == escape_char
      temp[result_index++] = escape_char;
      temp[result_index] = from[i-1];
    }
  }
  if (replacements[pair] == failbyte) {
    pair <<= 8;
    if (replacements[pair] == escape_char) {
      temp[result_index++] = escape_char;
    }
    temp[result_index] = from[length - 1];
  }
  assert(result_index < length + 3);
  std::copy(temp, temp + result_index, from);
  delete [] temp;
  return result_index;
}

} //namespace bwtc
