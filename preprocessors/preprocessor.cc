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

  unsigned symbols = 0;
  while(frequencies[symbols].second == 0) {
    ++symbols;
  }

  /* Lower bound for free symbols */
  const unsigned kLimit = 0; 
  if( symbols <=  kLimit) {
    /* We need to free some symbols */
    //TODO: Free symbols and replace some symbols, with pairs
  }
  assert(symbols > 0);

  /* Calculate the frequencies of pairs. It is done here because freeing of
   * symbols modifies the input */
  std::pair<uint16, uint32> pair_freqs[65536];
  for(unsigned i = 0; i < 65536; ++i) pair_freqs[i] = std::make_pair(i, 0);

  uint16 index = from[0];
  for(uint64 i = 1; i < length; ++i ) {
    index <<= 8;
    index |= from[i];
    ++pair_freqs[index].second;
  }
  /* Deduce the replaceable pairs from frequency table */
  unsigned current = 0;
  unsigned limit = 0;
  std::vector<uint16> replaceable_pairs;
  /* Here some more clever heuristic would be a good idea */ 
  while(replaceable_pairs.size() < symbols) {
    if(current == limit) {
      limit += 2*symbols;
      assert(limit < 65536);
      std::partial_sort(&pair_freqs[current], &pair_freqs[limit],
                        &pair_freqs[65536],
                        ComparePairSecondDesc<uint16, uint32>);
    }
    if(pair_freqs[current].second <= 3) {
      break;
    }

    byte fst = static_cast<byte>((pair_freqs[current].first & 0xFF00) >> 8);
    byte snd = static_cast<byte>(pair_freqs[current].first & 0x00FF);
    if (fst == snd) goto failure;

    for(std::vector<uint16>::iterator it = replaceable_pairs.begin();
        it != replaceable_pairs.end(); ++it)
    {
      uint16 current_fst = static_cast<byte>(((*it) & 0xFF00) >> 8);
      uint16 current_sec = static_cast<byte>(((*it) & 0x00FF));
      if (current_fst == snd || current_sec == fst) {
        pair_freqs[current].second = 0; /* For the later phases.. */
        goto failure;
      }
    }
    replaceable_pairs.push_back(pair_freqs[current].first);
 failure:
    ++current;
  }
  unsigned symbols_in_use = replaceable_pairs.size();
  /* Store the replacements in map<byte, uint16> where the latter is the
   * pair to be replaced */
  std::map<byte, uint16> pairs;
  
  /* failbyte can't be any of the symbols used in replacements so it
   * is a good choice for signaling "no replacement for this pair" */
  byte failbyte = static_cast<byte>(replaceable_pairs[0] & 0x00FF);
  byte replacements[65536];
  std::fill(replacements, replacements + 65536, failbyte);
  for(unsigned i = 0; i < symbols_in_use; ++i) {
    /* Here replaceable_pairs[i] is the pair which will be replaced with
     * frequencies[i] */
    replacements[replaceable_pairs[i]] = frequencies[i].first;
    pairs[frequencies[i].first] = replacements[replaceable_pairs[i]];
  }
  assert(pairs.size() == symbols_in_use);
  /* Allocate array for result */
  byte *temp = new byte[length];

  /* Write result to temporary-array */
  uint64 result_index = 0;
  std::pair<byte, uint16> curr_pair;
  /* Write replacements first */

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
  // TODO: character replacements
  /* Write the "end of sequence" */
  temp[result_index++] = curr_pair.first;

  uint16 pair = static_cast<uint16>(from[0]);
  for(uint64 i = 1; i < length; ++i, ++result_index) {
    pair <<= 8;
    pair |= from[i];
    if(replacements[pair] != failbyte) {
      temp[result_index] = replacements[pair];
      ++i;
      pair = static_cast<byte>(from[i]);
    }
    else {
      temp[result_index] = from[i-1];
    }
  }
  if (replacements[pair] == failbyte) temp[result_index] = from[length - 1];

  assert(result_index < length);
  std::copy(temp, temp + result_index, from);
  delete [] temp;
  return length - result_index;
}

} //namespace bwtc
