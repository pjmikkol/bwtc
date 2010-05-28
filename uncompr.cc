#include<iostream>
#include<iterator>

#include<boost/program_options.hpp>
namespace po = boost::program_options;

#include "globaldefs.h"


int main(int argc, char** argv) {
  int verbosity;
  std::string input_name, output_name;
  bool output_stdout;

  try {
    po::options_description description(
        "usage: "DECOMPRESSOR" [options] inputfile outputfile\n\nOptions");
    description.add_options()
        ("help,h", "print help message")
        ("stdout,c", "output to standard out")
        ("verb,v", po::value<int>(&verbosity)->default_value(0),
         "verbosity level")
        ("input-file", po::value<std::string>(&input_name),
         "file to decompress")
        ("output-file", po::value<std::string>(&output_name),
         "target file")
        ;

    // Allow input files given in user friendly form (without "--input-file")
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

    output_stdout = static_cast<bool>(varmap.count("stdout"));
  } // try-block

  catch(std::exception& e) {
    std::cerr << "error: " << e.what() << std::endl;
    return 1;
  }
  catch(...) {
    std::cerr << "Exception of unknown type!" << std::endl;
    return 1;    
  }

  if (input_name != "") {
    // TODO: Create inputfile object
    if (verbosity > 0) 
      std::clog << "Input: " << input_name << std::endl;
  } else {
    // TODO: Create outputstream for cin
    if (verbosity > 0) 
      std::clog << "Input: stdin" << std::endl;
  }
      
  if (output_name != "" && !output_stdout) {
    // TODO: Create outputfile object
    if (verbosity > 0)
      std::clog << "Output: " << output_name << std::endl;
  } else {
    // TODO: Create file object for cin
    if (verbosity > 0)
      std::clog << "Output: stdout" << std::endl;
  }
  
  return 0;
}
