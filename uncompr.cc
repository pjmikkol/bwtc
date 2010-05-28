#include<iostream>
#include<iterator>

#include<boost/program_options.hpp>
namespace po = boost::program_options;

#include "globaldefs.h"

using std::endl;
using std::cout;
using std::exception;
using std::cerr;
using std::vector;
using std::string;


int main(int argc, char** argv) {
  int verbosity;
  const vector<string> *input_files;

  try {
    po::options_description description(
        "usage: "DECOMPRESSOR" [options] input-files\n\nOptions");
    description.add_options()
        ("help,h", "print help message")
        ("verb,v", po::value<int>(&verbosity)->default_value(0), "verbosity level")
        ("input-file", po::value< vector<string> >(), "files to decompress")
        ;

    // Allow input files given in user friendly form (without "--input-file")
    po::positional_options_description pos;
    pos.add("input-file", -1);

    po::variables_map varmap;
    po::store(po::command_line_parser(argc, argv).options(description).
	      positional(pos).run(), varmap);
    po::notify(varmap);

    if (varmap.count("help")) {
      cout << description << endl;
      return 1;
    }

    if (varmap.count("input-file")) {
      input_files = &varmap["input-file"].as< vector<string> >();

      /* Create file objects*/

      if (verbosity) {
	cout << "Input files: ";
	for(vector<string>::const_iterator it = input_files->begin();
	    it != input_files->end(); ++it)
	  cout << *it << " ";
	cout << endl;
      }
    } else {
      // Create file object for cin

    }
    
  } // try-block

  catch(exception& e) {
    cerr << "error: " << e.what() << endl;
    return 1;
  }
  catch(...) {
    cerr << "Exception of unknown type!" << endl;
    return 1;    
  }

  return 0;
}
