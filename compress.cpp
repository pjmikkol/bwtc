/**
 * @file compress.cpp
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
 * Main program for compression.
 */

/**Needed for allocating space for the bwtc::verbosity. */
#define MAIN

/* bwtc-compressor main program */
#include <iostream>
#include <string>
#include <iterator>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "MainBlock.hpp"
#include "BlockManager.hpp"
#include "EntropyCoders.hpp"
#include "preprocessors/Preprocessor.hpp"
#include "Streams.hpp"
#include "globaldefs.hpp"
#include "bwtransforms/BWTransform.hpp"

#include "Profiling.hpp"

using bwtc::verbosity;

void compress(const std::string& input_name, const std::string& output_name,
              uint64 block_size, const std::string& preproc, char encoding,
              uint32 startingPoints, char bwtAlgo)
{
  PROFILE("TOTAL_compression_time");

  if (verbosity > 1) {
    if (input_name != "") std::clog << "Input: " << input_name << std::endl;
    else std::clog << "Input: stdin" << std::endl;
    if (output_name != "") std::clog << "Output: " << output_name << std::endl;
    else std::clog << "Output: stdout" << std::endl;
  }
  bwtc::Preprocessor preprocessor(block_size, preproc);
  preprocessor.connect(input_name);

  bwtc::BlockManager block_manager(block_size + 5*preproc.size() + 1, 1);
  preprocessor.addBlockManager(&block_manager);

  bwtc::BWTransform* transformer = bwtc::giveTransformer(bwtAlgo);

  bwtc::EntropyEncoder *encoder = bwtc::giveEntropyEncoder(output_name, encoding);
  encoder->writeGlobalHeader(encoding);

  unsigned blocks = 0;
  uint64 last_s = 0;
  while( bwtc::MainBlock* block = preprocessor.readBlock() ) {
    ++blocks;
    //Transformer could have some memory manager..
    transformer->connect(block);
    transformer->buildStats();

    /* The following way enables the calculation of transformation in 
     * several phases */
    std::vector<uint32> LFpowers;
    LFpowers.resize(startingPoints);
    if(verbosity > 1) {
      std::clog << "Writing " << startingPoints
                << " starting points for inverse transform." << std::endl;
    }
    transformer->doTransform(LFpowers);

    encoder->writeBlockHeader(block->m_stats);
    encoder->encodeData(block->begin(), block->m_stats, block->size());


    encoder->finishBlock(LFpowers);
    last_s = block->m_filled;
    delete block;
  }

#ifdef ENTROPY_PROFILER  
  if(encoding != 'H') {
    std::cout << "kBytes for characters " << (encoder->m_bytesForCharacters/1024.0)
              << std::endl << "kBytes for runs " << (encoder->m_bytesForRuns/1024.0)
              << std::endl;
  }
#endif
  
  if (verbosity > 0) {
    std::clog << "Read " << blocks << " block" << ((blocks < 2)?"":"s") << "\n";
    std::clog << "Total size: " << (blocks-1)*block_size + last_s << "B\n";
  }
  delete encoder;
  delete transformer;
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
    //if(c != 'p' && c != 's') throw exc;
    if(c != 'p' ) throw exc;
  }
}

/* Notifier function for starting points option choice */
void validateStartingPoints(uint32 sp) {
  if (sp > 0 && sp <= 256) return;

  class StartPExc : public std::exception {
    virtual const char* what() const throw() {
      return "Invalid choice for starting points of decompression.";
    }
  } exc;

  throw exc;
}

/* Notifier function for encoding option choice */
void validateEncodingOption(char c) {
  if (c == 'H' || c == 'm' || c == 'M' || c == 'u' || c == 'b' || c == 'B'
      /* || c == <other option> */) return;

  class EncodingExc : public std::exception {
    virtual const char* what() const throw() {
      return "Invalid choice for entropy encoding.";
    }
  } exc;

  throw exc;
}

/* Notifier function for encoding option choice */
void validateBWTchoice(char c) {
  if (c == 'd' || c == 's') return;

  class EncodingExc : public std::exception {
    virtual const char* what() const throw() {
      return "Invalid choice for BWT-algorithm.";
    }
  } exc;

  throw exc;
}


int main(int argc, char** argv) {
  uint64 block_size;
  char encoding, bwtAlgo;
  std::string input_name, output_name, preprocessing;
  bool stdout, stdin;
  uint32 startingPoints;

  try {
    po::options_description description(
        "usage: "COMPRESSOR" [options] inputfile outputfile\n\nOptions");
    description.add_options()
        ("help,h", "print help message")
        ("stdin,i", "input from standard in")
        ("stdout,c", "output to standard out")
        ("block,b", po::value<uint64>(&block_size)->default_value(100000),
         "Block size for compression (in kB)")
        ("starts,s", po::value<uint32>(&startingPoints)->default_value(8)->
         notifier(&validateStartingPoints),
         "Starting points for decompression (more means faster decompression).")
        ("verb,v", po::value<int>(&verbosity)->default_value(0),
         "verbosity level")
        ("input-file", po::value<std::string>(&input_name),
         "file to compress, defaults to stdin")
        ("output-file", po::value<std::string>(&output_name),"target file")
        ("bwt", po::value<char>(&bwtAlgo)->default_value('d')->notifier(&validateBWTchoice),
         "BWT-algorithm to use:\n"
         "  d -- Yuta Mori's libdivsufsort\n"
         "  s -- Yuta Mori's sais")
        ("prepr", po::value<std::string>(&preprocessing)->default_value("")->
         notifier(&validatePreprocOption),
         "preprocessor options:\n"
         "  p -- pair replacer\n"
         "  s -- long recurring sequences replacer\n"
         "For example \"pps\" would run pair replacer twice and sequence replacer once")
        ("enc,e", po::value<char>(&encoding)->default_value('B')->
         notifier(&validateEncodingOption),
         "entropy encoding scheme, options:\n"
         "  H -- Huffman coding with run-length encoding\n"
         "  M -- Remembering 16 previous bits (Wavelet tree)\n"
         "  m -- Remembering 8 previous bits (Wavelet tree)\n"
         "  b -- Finite State Machine with unbiased and equal predictors in each state (Wavelet tree)\n"
         "  B -- Slightly optimised version of above (Wavelet tree)\n"
         "  u -- Simple predictor with 4 states. These are used in FSM's states (Wavelet tree)")
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

  compress(input_name, output_name, block_size*1024, preprocessing, encoding, startingPoints, bwtAlgo);

  PRINT_PROFILE_DATA
  return 0;
}
