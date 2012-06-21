/**
 * @file preprocess.cpp
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
 * Main program for preprocessor.
 */

/**Needed for allocating space for the bwtc::verbosity. */
#define MAIN
#define PREPROCESSOR "preprocess"

/* bwtc-preprocessor main program */
#include <iostream>
#include <string>
#include <iterator>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "preprocessors/Precompressor.hpp"
#include "Streams.hpp"
#include "Utils.hpp"
#include "globaldefs.hpp"

#include "Profiling.hpp"

using bwtc::verbosity;

void writePackedInteger(bwtc::uint64 packed_integer, bwtc::OutStream *out) {
  do {
    bwtc::byte to_written = static_cast<bwtc::byte>(packed_integer & 0xFF);
    packed_integer >>= 8;
    out->writeByte(to_written);
  } while (packed_integer);
}

void preprocess(const std::string& input_name, const std::string& output_name,
              bwtc::uint64 block_size, const std::string& preproc, bool escaping)
{
  PROFILE("TOTAL_preprocessing_time");

  if (verbosity > 1) {
    if (input_name != "") std::clog << "Input: " << input_name << std::endl;
    else std::clog << "Input: stdin" << std::endl;
    if (output_name != "") std::clog << "Output: " << output_name << std::endl;
    else std::clog << "Output: stdout" << std::endl;
  }
  bwtc::RawInStream in(input_name);
  
  bwtc::Precompressor preprocessor(preproc);

  bwtc::RawOutStream out(output_name);

  unsigned blocks = 0;
  bwtc::uint64 last_s = 0;
  while(bwtc::PrecompressorBlock* block = preprocessor.readBlock(block_size, &in)) {
    if(block->originalSize() == 0) {
      delete block;
      break;
    }

    ++blocks;
    // TODO:
    // Write pb->size() to be able to know the data needed to read after
    // header. Read that in postprocess.cpp and modify PrecompressorBlock a bit

    block->writeBlockHeader(&out); //write Grammar
    out.writeBlock(block->begin(), block->end());

    last_s = block->originalSize();
    delete block;
  }

  if (verbosity > 0) {
    std::clog << "Read " << blocks << " block" << ((blocks < 2)?"":"s") << "\n";
    std::clog << "Total size: " << (blocks-1)*block_size + last_s << "B\n";
  }
}

/* Notifier function for preprocessing option choice */
void validatePreprocOption(const std::string& p) {
  class PreprocException : public std::exception {
    virtual const char* what() const throw() {
      return "Invalid choice for preprocessing.";
    }
  } exc;

  for(size_t i = 0; i < p.size(); ++i) {
    char c = p[i];
    if(c != 'c' && c != 'p' && c != 'r' && c != 's') throw exc;
  }
}

int main(int argc, char** argv) {
  uint64 block_size;
  std::string input_name, output_name, preprocessing;
  bool stdout, stdin, escaping;

  try {
    po::options_description description(
        "usage: "PREPROCESSOR" [options] inputfile outputfile\n\nOptions");
    description.add_options()
        ("help,h", "print help message")
        ("stdin,i", "input from standard in")
        ("stdout,c", "output to standard out")
        ("block,b", po::value<uint64>(&block_size)->default_value(100000),
         "Block size for compression (in kB)")
        ("verb,v", po::value<int>(&verbosity)->default_value(0),
         "verbosity level")
        ("escape", po::value<bool>(&escaping)->default_value(true),
         "are preprocessing algorithms using escaping (0 to disable)")
        ("input-file", po::value<std::string>(&input_name),
         "file to compress, defaults to stdin")
        ("output-file", po::value<std::string>(&output_name),"target file")
        ("prepr", po::value<std::string>(&preprocessing)->default_value("")->
         notifier(&validatePreprocOption),
         "preprocessor options:\n"
         "  p -- pair replacer\n"
         "  r -- run replacer\n"
         "  c -- pair and run replacer\n"
         "  s -- long recurring sequences replacer\n"
         "For example \"ppr\" would run pair replacer twice and run replacer once")
        ;

    /* Allow input and output files given in user friendly form,
     * depending on their positions */
    po::positional_options_description pos;
    pos.add("input-file", 1);
    pos.add("output-file", 1);


    po::variables_map varmap;
    po::store(po::command_line_parser(argc, argv).options(description).
	      positional(pos).run(), varmap);
    po::notify(varmap);

    if (varmap.count("help")) {
      std::cout << description << std::endl;
      return 0;
    }

    stdout = varmap.count("stdout") != 0;
    stdin  = varmap.count("stdin") != 0;
    // TODO: Check that the block-size is OK
  } /* try-block */

  catch(std::exception& e) {
    std::cerr << "error: " << e.what() << std::endl;
    return 1;
  }
  catch(...) {
    std::cerr << "Exception of unknown type!" << std::endl;
    return 1;
  }

  if (verbosity > 0) {
    std::clog << "Block size = " << block_size <<  "kB" << std::endl;
  }
  if (block_size <= 0) block_size = 1;

  if (stdout) output_name = "";
  if (stdin)  input_name = "";

  preprocess(input_name, output_name, block_size*1024, preprocessing, escaping);

  PRINT_PROFILE_DATA
  return 0;
}
