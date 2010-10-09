/**
 * @file seq_det_test.cc
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
#include "../globaldefs.h"

#define BWTC_SEQ_DET_SCAN_CONSTS_
namespace bwtc {
namespace long_sequences {
const uint32 kMinPeriod = 4;
const uint32 kWindowSize = 8;
}
}

#include "../preprocessors/sequence_detector.h"
#include "../preprocessors/sequence_detector-inl.h"

namespace bwtc {
int verbosity = 0;

namespace tests {

using namespace long_sequences;

BOOST_AUTO_TEST_SUITE(FirstPhase)

void DebPrint(const char* str, int len,
              SequenceDetector<hash_functions::PrimeHasher>& seq_det)
{
  std:: cout << str << " " << seq_det.Count((byte*)str,len) << std::endl;
}


BOOST_AUTO_TEST_CASE(PrimePeriodRecognition1) {
  /* Testing recognition of periodical string */
  byte *t_string = (byte*)
      "abcdabcd" "abcdabcdabcd" "abcdabcdab";
  uint32 freqs[256];
  long_sequences::SequenceDetector<hash_functions::PrimeHasher>
      seq_det(t_string, 330, freqs);
  seq_det.ScanAndStore(30);
  BOOST_CHECK_EQUAL(4, seq_det.Count((byte*)"abcdabcd", 8));
}

BOOST_AUTO_TEST_CASE(PrimePeriodRecognition2) {
  byte *t_string = (byte*)
      "abcdefga" "bcdefghabcdefgha";
  uint32 freqs[256];
  long_sequences::SequenceDetector<hash_functions::PrimeHasher>
      seq_det(t_string, 330, freqs);
  seq_det.ScanAndStore(24);
  BOOST_CHECK_EQUAL(2, seq_det.Count((byte*)"bcdefgha", 8));
}

BOOST_AUTO_TEST_CASE(PrimePeriodRecognition3) {
  byte *t_string = (byte*)
      "abcdeabc" "deabcdeabcdeabcde";
  uint32 freqs[256];
  long_sequences::SequenceDetector<hash_functions::PrimeHasher>
      seq_det(t_string, 330, freqs);
  seq_det.ScanAndStore(25);
  BOOST_CHECK_EQUAL(3, seq_det.Count((byte*)"abcdeabc", 8));
}

BOOST_AUTO_TEST_CASE(PrimeChoosingRightValues) {
  byte *t_string = (byte*)
      "aaaaffff" "sfsaaseaaaaffffwqcaffaaaaffffwqvnaaaaffff";
  uint32 freqs[256];
  long_sequences::SequenceDetector<hash_functions::PrimeHasher>
      seq_det(t_string, 330, freqs);
  seq_det.ScanAndStore(49);
  BOOST_CHECK_EQUAL(4, seq_det.Count((byte*)"aaaaffff", 8));
}

/* Note that the parameters hera are chosen so that there won't be
 * collisions*/
BOOST_AUTO_TEST_CASE(PrimeDifferentValues) {
  byte *t_string = (byte*)
      "a4cdefghhjabcdllkksnwlsdfiwpwsdklrvnca2qlgfp4302sdf"
      "cxmnwaedfsasd2lk4rfdshbdfdbsfdhbmna2epldmapspbvnsda";
  uint32 freqs[256];
  long_sequences::SequenceDetector<hash_functions::PrimeHasher>
      seq_det(t_string, 630, freqs);
  seq_det.ScanAndStore(102);
  BOOST_CHECK_EQUAL(7, seq_det.ChunksCount());
  BOOST_CHECK_EQUAL(1, seq_det.Count((byte*)"a4cdefgh", 8));
  BOOST_CHECK_EQUAL(1, seq_det.Count((byte*)"lkksnwls", 8));
  BOOST_CHECK_EQUAL(1, seq_det.Count((byte*)"dklrvnca", 8));
  BOOST_CHECK_EQUAL(1, seq_det.Count((byte*)"sasd2lk4", 8));
}

BOOST_AUTO_TEST_SUITE_END()

} //namespace tests
} //namespace long_sequences
