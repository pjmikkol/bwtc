#include<iostream>
#include<iterator>

#include<boost/program_options.hpp>
namespace po = boost::program_options;

#include "globaldefs.h"


int main(int argc, char** argv) {
  int verbosity;
  const std::vector<std::string> *input_files;

  try {
    po::options_description description(
        "usage: "DECOMPRESSOR" [options] input-files\n\nOptions");
    description.add_options()
        ("help,h", "print help message")
        ("verb,v", po::value<int>(&verbosity)->default_value(0), "verbosity level")
        ("input-file", po::value< std::vector<std::string> >(), "files to decompress")
        ;

    // Allow input files given in user friendly form (without "--input-file")
    po::positional_options_description pos;
    pos.add("input-file", -1);

    po::variables_map varmap;
    po::store(po::command_line_parser(argc, argv).options(description).
	      positional(pos).run(), varmap);
    po::notify(varmap);

    if (varmap.count("help")) {
      std::cout << description << std::endl;
      return 1;
    }

    if (varmap.count("input-file")) {
      input_files = &varmap["input-file"].as< std::vector<std::string> >();

      /* Create file objects*/

      if (verbosity) {
	std::cout << "Input files: ";
	for(std::vector<std::string>::const_iterator it = input_files->begin();
	    it != input_files->end(); ++it)
	  std::cout << *it << " ";
	std::cout << std::endl;
      }
    } else {
      // Create file object for cin

    }
    
  } // try-block

  catch(std::exception& e) {
    std::cerr << "error: " << e.what() << std::endl;
    return 1;
  }
  catch(...) {
    std::cerr << "Exception of unknown type!" << std::endl;
    return 1;    
  }

  return 0;
}
