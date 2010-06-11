#include <iostream>
#include <string>
#include <iterator>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "block.h"
#include "block_manager.h"
#include "coders.h"
#include "preprocessor.h"
#include "stream.h"
#include "globaldefs.h"

int bwtc::verbosity;
using bwtc::verbosity;

void compress(const std::string& input_name, const std::string& output_name,
              int64 block_size, char preproc, char encoding)
{
  if (verbosity > 1) {
    if (input_name != "") std::clog << "Input: " << input_name << std::endl;
    else std::clog << "Input: stdin" << std::endl;
    if (output_name != "") std::clog << "Output: " << output_name << std::endl;
    else std::clog << "Output: stdout" << std::endl;
  }
  bwtc::PreProcessor* preprocessor = bwtc::GivePreProcessor(
      preproc,block_size, input_name);
  bwtc::BlockManager block_manager(block_size);
  preprocessor->AddBlockManager(&block_manager);

  bwtc::Encoder encoder(output_name, encoding);

  unsigned blocks = 0;
  int64 last_s = 0;

  while( bwtc::MainBlock* block = preprocessor->ReadBlock() ) {
    blocks++;
    /* while (Block* b = DoTransform()) 
       encodeBlock(b)
     */
    last_s = block->Size();
    delete block;
  }

  if (verbosity > 0) {
    std::clog << "Read " << blocks << " block" << ((blocks < 2)?"":"s") << "\n";
    std::clog << "Total size: " << (blocks-1)*block_size + last_s << "B\n";
  }

  delete preprocessor;
}

/* Notifier function for preprocessing option choice */
void ValidatePreprocOption(char c) {
  if (c == 'n' /* || c == <other option> */) return;

  class PreprocException : public std::exception {
    virtual const char* what() const throw() {
      return "Invalid choice for preprocessing.";
    }
  } exc;

  throw exc;
}


/* Notifier function for encoding option choice */
void ValidateEncodingOption(char c) {
  if (c == 'n' /* || c == <other option> */) return;

  class EncodingExc : public std::exception {
    virtual const char* what() const throw() {
      return "Invalid choice for entropy encoding.";
    }
  } exc;

  throw exc;
}


int main(int argc, char** argv) {
  int64 block_size;
  char preproc, encoding;
  std::string input_name, output_name;
  bool stdout, stdin;

  try {
    po::options_description description(
        "usage: "COMPRESSOR" [options] inputfile outputfile\n\nOptions");
    description.add_options()
        ("help,h", "print help message")
        ("stdin,i", "input from standard in")
        ("stdout,c", "output to standard out")
        ("block,b", po::value<int64>(&block_size)->default_value(1000),
         "Block size for compression (in kB)")
        ("verb,v", po::value<int>(&verbosity)->default_value(0),
         "verbosity level")
        ("input-file", po::value<std::string>(&input_name),
         "file to compress, defaults to stdin")
        ("output-file", po::value<std::string>(&output_name),"target file")
        ("pre,p", po::value<char>(&preproc)->default_value('n')->
         notifier(&ValidatePreprocOption),
         "pre-processing algorithm, options:\n"
         "  n -- does nothing")
        ("enc,e", po::value<char>(&encoding)->default_value('n')->
         notifier(&ValidateEncodingOption),
         "entropy encoding scheme, options:\n"
         "  n -- does nothing")
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

  compress(input_name, output_name, block_size*1024, preproc, encoding);

  return 0;
}
