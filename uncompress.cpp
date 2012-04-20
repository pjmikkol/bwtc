/**
 * @file uncompress.cpp
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
 * Main program for decompression.
 */

/**Needed for allocating space for the bwtc::verbosity. */
#define MAIN

/* bwtc-decompressor main program */
#include <iostream>
#include <iterator>
#include <algorithm>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "preprocessors/Postprocessor.hpp"
#include "EntropyCoders.hpp"
#include "Streams.hpp"
#include "globaldefs.hpp"
#include "bwtransforms/InverseBWT.hpp"
#include "Profiling.hpp"

using bwtc::verbosity;

void decompress(const std::string& input_name, const std::string& output_name,
                int verbosity) {
  PROFILE("TOTAL_decompression_time");
  if (verbosity > 0) {
    if (input_name != "") std::clog << "Input: " << input_name << std::endl;
    else std::clog << "Input: stdin" << std::endl;
    if (output_name != "") std::clog << "Output: " << output_name << std::endl;
    else std::clog << "Output: stdout" << std::endl;
  }
  bwtc::EntropyDecoder *decoder = bwtc::giveEntropyDecoder(input_name);
  std::string postproc = decoder->readGlobalHeader();
  std::reverse(postproc.begin(), postproc.end());
  bwtc::PostProcessor postProcessor(postproc);

  if(verbosity > 1) {
    std::clog << "Postprocessor initiated with parameter: " << postproc
              << std::endl;
  }

  
  bwtc::RawOutStream out(output_name);

  bwtc::InverseBWTransform* transformer = bwtc::giveInverseTransformer();

  unsigned blocks = 0;

  std::vector<uint32> LFpowers;
  while (std::vector<byte>* bwt_block = decoder->decodeBlock(LFpowers)) {

    if(verbosity > 1) {
      std::clog << "Read " << LFpowers.size() << " starting points for "
                << "inverse transform." << std::endl;
    }
    ++blocks;
    std::vector<byte>* unbwt_block = transformer->doTransform(&(*bwt_block)[0],
                                                              bwt_block->size(),
                                                              LFpowers);
    delete bwt_block;

    postProcessor.postProcess(unbwt_block);
    //out.writeBlock(unbwt_block->begin(), unbwt_block->end());
    byte *unbwt_block_ptr = &(*unbwt_block)[0];
    out.writeBlock(unbwt_block_ptr, unbwt_block_ptr + ((size_t)unbwt_block->size()));

    out.flush();
    delete unbwt_block;
  }
  if (verbosity > 0) {
    std::clog << "Read " << blocks << " block" << ((blocks < 2)?"":"s") << "\n";
  }

  delete decoder;
  delete transformer;
}

int main(int argc, char** argv) {
  std::string input_name, output_name;
  bool stdout, stdin;

  try {
    po::options_description description(
        "usage: "DECOMPRESSOR" [options] inputfile outputfile\n\nOptions");
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

  decompress(input_name, output_name, verbosity);
  
  PRINT_PROFILE_DATA
  return 0;
}
