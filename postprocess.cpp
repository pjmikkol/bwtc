/**
 * @file postprocess.cpp
 * @author Pekka Mikkola <pjmikkol@cs.helsinki.fi>
 * @author Dominik Kempa <dominik.kempa@cs.helsinki.fi>
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
 * Main program for postprocessing.
 */

/**Needed for allocating space for the bwtc::verbosity. */
#define MAIN
#define POSTPROCESSOR "postprocess"

/* bwtc-postprocessor main program */
#include <string>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <algorithm>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "preprocessors/Postprocessor.hpp"
#include "Streams.hpp"
#include "Utils.hpp"
#include "globaldefs.hpp"
#include "Profiling.hpp"
#include "PrecompressorBlock.hpp"

using bwtc::verbosity;

uint64 readPackedInteger(bwtc::InStream *in) {
  static const uint64 kEndSymbol = static_cast<uint64>(1) << 63;
  static const uint64 kEndMask = static_cast<uint64>(1) << 7;

  uint64 packed_integer = 0;
  bool bits_left = true;
  int i;
  for(i = 0; bits_left; ++i) {
    uint64 read = static_cast<uint64>(in->readByte());
    bits_left = (read & kEndMask) != 0;
    packed_integer |= (read << i*8);
  }
  if (packed_integer == 0x80) return kEndSymbol;
  return packed_integer;
}

void postprocess(const std::string& input_name, const std::string& output_name,
                int verbosity) {
  PROFILE("TOTAL_postprocessing_time");
  if (verbosity > 0) {
    if (input_name != "") std::clog << "Input: " << input_name << std::endl;
    else std::clog << "Input: stdin" << std::endl;
    if (output_name != "") std::clog << "Output: " << output_name << std::endl;
    else std::clog << "Output: stdout" << std::endl;
  }

  bwtc::RawInStream in(input_name);

  bwtc::RawOutStream out(output_name);

  unsigned blocks = 0;
  while(true) {
    bwtc::PrecompressorBlock *pb =
        bwtc::PrecompressorBlock::readBlockHeader(&in);
    if(pb->originalSize() == 0) {
      delete pb;
      break;
    }
    ++blocks;

    // Postprocess pb
    bwtc::Postprocessor postprocessor(verbosity > 1, pb->grammar());
    size_t postSize = postprocessor.uncompress(pb->begin(), pb->size(), &out);
    assert(postSize == pb->originalSize());
    delete pb;
  }

  if (verbosity > 0) {
    std::clog << "Read " << blocks << " block" << ((blocks < 2)?"":"s") << "\n";
  }
}

int main(int argc, char** argv) {
  std::string input_name, output_name;
  bool stdout, stdin;

  try {
    po::options_description description(
        "usage: "POSTPROCESSOR" [options] inputfile outputfile\n\nOptions");
    description.add_options()
        ("help,h", "print help message")
        ("stdin,i", "input from standard in")
        ("stdout,c", "output to standard out")
        ("verb,v", po::value<int>(&verbosity)->default_value(0),
         "verbosity level")
        ("input-file", po::value<std::string>(&input_name),
         "file to decompress")
        ("output-file", po::value<std::string>(&output_name),
         "target file")
        ;

    /* Allow input files given in user friendly form (without "--*put-file") */
    po::positional_options_description pos;
    pos.add("input-file", 1);
    pos.add("output-file",1);

    po::variables_map varmap;
    po::store(po::command_line_parser(argc, argv).options(description).
	      positional(pos).run(), varmap);
    po::notify(varmap);

    if (varmap.count("help")) {
      std::cout << description << std::endl;
      return 1;
    }

    stdout = varmap.count("stdout") != 0;
    stdin  = varmap.count("stdin") != 0;
  } /* try-block */

  catch(std::exception& e) {
    std::cerr << "error: " << e.what() << std::endl;
    return 1;
  }
  catch(...) {
    std::cerr << "Exception of unknown type!" << std::endl;
    return 1;    
  }

  if (stdout) output_name = "";
  if (stdin)  input_name = "";

  postprocess(input_name, output_name, verbosity);
  
  PRINT_PROFILE_DATA
  return 0;
}
