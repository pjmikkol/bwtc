#ifndef BWTC_CODERS_H_
#define BWTC_CODERS_H_

#include <iostream>
#include <string>
#include <vector>

#include "block.h"
#include "globaldefs.h"
#include "rl_compress.h"
#include "probmodels/base_prob_model.h"

namespace bwtc {

/*************************************************************************
 * PackInteger and UnpackInteger                                         *
 * Packs integer to bytefields, so that the first bit of the byte        *
 * tells if the number has succeeding bytes. Few examples:               *
 * 0xF0   -> 0x01F0                                                      *
 *   -- the last byte is F0, because of the continuation-bit             *
 * 0x2    -> 0x2                                                         *
     -- no overhead here                                                 *
 * 0x142A -> 0x28AA                                                      *
 *   -- no overhead here because the most significant byte had leading   *
 *      zeroes                                                           *
 *************************************************************************/
uint64 PackInteger(uint64 integer, int* bytes_needed);

uint64 UnpackInteger(uint64 packed_integer);

/**********************************************************************
 * Encoder and decoder are pretty similar in structure.               *
 *                                                                    *
 * Both have a field for ProbabilityModel-object which ultimately     *
 * decides how the encoding or decoding is done.                      *
 *                                                                    *
 * Both have also destination or source field which is for arithmetic *
 * encoder/decoder (objects of a type BitEncoder/BitDecoder).         *
 **********************************************************************/

class Encoder {
 public:
  Encoder(const std::string& destination, char prob_model);
  ~Encoder();
  void WriteGlobalHeader(char preproc, char encoding);
  void EncodeByte(byte b);
  void EncodeData(std::vector<byte>* data, std::vector<uint64>* stats,
                  uint64 data_size);
  //void EncodeMainBlock(MainBlock* block, uint64 trailer);
  void EncodeRange(const byte* begin, const byte* end);
  void WriteBlockHeader(std::vector<uint64>* stats);
  void WritePackedInteger(uint64 packed_integer, int bytes);
  int FinishBlockHeader();
  void EndContextBlock();
  int WriteTrailer(uint64 trailer_value);
  void FinishBlock(uint64 eob_byte); //TODO: this calls write trailer

 private:
  OutStream* out_;
  dcsbwt::BitEncoder* destination_;
  ProbabilityModel* pm_;
  std::streampos header_position_;
  uint64 compressed_block_length_;
  /* We may have to encode result of the transformation in pieces so we
   * have track down the progress of handling single MainBlock. */
  uint64 current_stat_handled_;
  unsigned current_stat_index_;


  Encoder(const Encoder&);
  Encoder& operator=(const Encoder&);
};

class Decoder {
 public:
  Decoder(const std::string& source, char prob_model);
  Decoder(const std::string& source);
  ~Decoder();
  /* ReadGlobalHeader returns char denoting the preprocessing algorithm.
   * It changes the used probability model automatically. */
  char ReadGlobalHeader();
  byte DecodeByte();
  void Start() { source_->Start(); }
  /* If end symbol is encountered, then the most significant bit is activated */
  uint64 ReadPackedInteger();
  /* Allocates memory for block, reads and decodes it. */
  std::vector<byte>* DecodeBlock(uint64* eof_byte_in_bwt);
  /* Returns length of the compressed sequence and stores lengths of the context
   * blocks into stats-array.*/
  uint64 ReadBlockHeader(std::vector<uint64>* stats);
  void EndContextBlock();

 private:
  InStream* in_;
  dcsbwt::BitDecoder* source_;
  ProbabilityModel* pm_;

  Decoder(const Decoder&);
  Decoder& operator=(const Decoder&);
};

ProbabilityModel* GiveProbabilityModel(char choice);
  
} // namespace bwtc

#endif
