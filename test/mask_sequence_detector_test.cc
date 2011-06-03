/**
 * @file mask_seq_det_test.cc
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
 * Unit tests for bwtc::long_sequences::SequenceDetector-class. Detection of
 * repeating sequencies.
 */

#define BOOST_TEST_MODULE 
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <cassert>
#include <iostream>
#include "../globaldefs.hpp"

#include "../preprocessors/sequence_detector.h"
#include "../preprocessors/sequence_detector-inl.h"

namespace bwtc {
int verbosity = 0;

namespace tests {

const uint32 kWindowSize = 9;

using namespace long_sequences;

/**
 * Arrays and vectors needed for SequenceDetector.
 */
struct Vectors {
 public:
  Vectors() { std::fill(freqs, freqs + 256, 0); }
  std::vector<uint32> uints;
  std::vector<chunk> chunks;
  uint32 freqs[256];
};

BOOST_AUTO_TEST_SUITE(FirstPhase)

BOOST_AUTO_TEST_CASE(MaskPeriodRecognition1) {
  byte *t_string = (byte*)
      "abcdefghi" "abcdefghiabcdefghiabcdefghi";
  Vectors v;
  long_sequences::SequenceDetector<hash_functions::MaskHasher>
      seq_det(t_string, 1030, v.freqs, &v.chunks, &v.uints, kWindowSize);
  seq_det.ScanAndStore(36);
  BOOST_CHECK_EQUAL(4, seq_det.Count((byte*)"abcdefghi", 9));
}

BOOST_AUTO_TEST_CASE(MaskPeriodRecognition2) {
  byte *t_string = (byte*)
      "abcdefgab" "cdefghabcdefghabc";
  Vectors v;
  long_sequences::SequenceDetector<hash_functions::MaskHasher>
      seq_det(t_string, 1030, v.freqs, &v.chunks, &v.uints, kWindowSize);
  seq_det.ScanAndStore(26);
  BOOST_CHECK_EQUAL(2, seq_det.Count((byte*)"cdefghabc", 9));
}

BOOST_AUTO_TEST_CASE(MaskPeriodRecognition3) {
  byte *t_string = (byte*)
      "abcdeabcd" "eabcdeabcdeabcdeab";
  Vectors v;
  long_sequences::SequenceDetector<hash_functions::MaskHasher>
      seq_det(t_string, 1030, v.freqs, &v.chunks, &v.uints, kWindowSize);
  seq_det.ScanAndStore(27);
  BOOST_CHECK_EQUAL(2, seq_det.Count((byte*)"eabcdeabc", 9));
}

BOOST_AUTO_TEST_CASE(MaskChoosingRightValues) {
  byte *t_string = (byte*)
      "aaaAffffg" "sfsaaseaaaAffffgwqcaffaaaAffffgwqvnaaaAffffg";
  Vectors v;
  long_sequences::SequenceDetector<hash_functions::MaskHasher>
      seq_det(t_string, 1030, v.freqs, &v.chunks, &v.uints, kWindowSize);
  seq_det.ScanAndStore(53);
  BOOST_CHECK_EQUAL(4, seq_det.Count((byte*)"aaaAffffg", 9));
}

BOOST_AUTO_TEST_CASE(MaskDifferentValues) {
  byte *t_string = (byte*)
      "abcdefghhjabcdllkksnwlsdfiwpasdklnvnca2qlgfp43y2sdrc"
      "amnwa-dfsasd2lk4rfdshbdfdbsfdhbmna2epldmapspbvnsdagt";
  Vectors v;
  long_sequences::SequenceDetector<hash_functions::MaskHasher>
      seq_det(t_string, 1030, v.freqs, &v.chunks, &v.uints, kWindowSize);
  seq_det.ScanAndStore(104);
  BOOST_CHECK_EQUAL(7, seq_det.ChunksCount());
  BOOST_CHECK_EQUAL(1, seq_det.Count((byte*)"abcdefghh", 9));
  BOOST_CHECK_EQUAL(1, seq_det.Count((byte*)"ksnwlsdfi", 9));
  BOOST_CHECK_EQUAL(1, seq_det.Count((byte*)"vnca2qlgf", 9));
}

BOOST_AUTO_TEST_CASE(CorrectFreqCounts) {
  byte *t_string = (byte*)
      "a4cdefghhjabcdllkksnwlsdfiwpwsdklrvnca2qlgfp4302sdf"
      "cxmnwaedfsasd2lk4rfdshbdfdbsfdhbmna2epldmapspbvnsda";
  Vectors v;
  long_sequences::SequenceDetector<hash_functions::MaskHasher>
      seq_det(t_string, 1030, v.freqs, &v.chunks, &v.uints, kWindowSize);
  seq_det.ScanAndStore(102);
  BOOST_CHECK_EQUAL(v.freqs[(byte)'a'], 8);
  BOOST_CHECK_EQUAL(v.freqs[(byte)'c'], 4);
  BOOST_CHECK_EQUAL(v.freqs[(byte)'d'], 13);
  BOOST_CHECK_EQUAL(v.freqs[(byte)'m'], 3);
}

BOOST_AUTO_TEST_SUITE_END()

} //namespace tests
} //namespace long_sequences

