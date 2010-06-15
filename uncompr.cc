#include <iostream>
#include <iterator>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "coders.h"
#include "stream.h"
#include "globaldefs.h"

int bwtc::verbosity;
using bwtc::verbosity;

void decompress(const std::string& input_name, const std::string& output_name,
                int verbosity) {
  if (verbosity > 0) {
    if (input_name != "") std::clog << "Input: " << input_name << std::endl;
    else std::clog << "Input: stdin" << std::endl;
    if (output_name != "") std::clog << "Output: " << output_name << std::endl;
    else std::clog << "Output: stdout" << std::endl;
  }
  bwtc::Decoder decoder(input_name);
  char preproc = decoder.ReadGlobalHeader();
  bwtc::OutStream out(output_name);

  unsigned blocks = 0;
  uint64 eof_byte; // eof_byte position in BWT
  while (std::vector<byte>* block = decoder.DecodeBlock(&eof_byte)) {
    // UnBWT(block, eof_byte);
    // PostProcess
    out.WriteBlock(block->begin(), block->end());
    out.Flush();
    delete block;
  }
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
  
  return 0;
}
