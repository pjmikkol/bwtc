#include<iostream>
#include<vector>
#include<string>
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


// Validates preprocessing options
void validate_preproc_opt(const char c) {
  if (c == 'n' /* || c == <other option> */) return;

  class preproc_exception : public exception {
    virtual const char* what() const throw() {
      return "Invalid choice for preprocessing.";
    }
  } exc;

  throw exc;
}

// Validate entropy encoding option
void validate_entropy_end(const char c) {
  if (c == 'n' /* || c == <other option> */) return;

  class entropy_exc : public exception {
    virtual const char* what() const throw() {
      return "Invalid choice for entropy encoding.";
    }
  } exc;

  throw exc;
}


int main(int argc, char** argv) {

  int block_size, verbosity;
  char preproc, encoding;
  const vector<string> *input_files;

  try {
    po::options_description description("usage: "COMPRESSOR" [options] input-files\n\nOptions");
    description.add_options()
      ("help,h", "print help message")
      ("stdout,c", "output to standard out")
      ("block,b", po::value<int>(&block_size)->default_value(1000),
       "Block size for compression (in kB)")
      ("verb,v", po::value<int>(&verbosity)->default_value(0), "verbosity level")
      ("input-file", po::value< vector<string> >(), "files to compress, defaults to stdin")
      ("pre,p", po::value<char>(&preproc)->default_value('n')->notifier(&validate_preproc_opt), 
       "pre-processing algorithm, options:\n"
       "  n -- \tdoes nothing")
      ("enc,e", po::value<char>(&encoding)->default_value('n'),
       "entropy encoding scheme, options:\n"
       "  n -- \tdoes nothing")
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
      return 0;
    }

    // TODO: Check that the block-size is reasonable
    if (verbosity < 0) verbosity = 0;
    if (verbosity > MAXVERB) verbosity = MAXVERB; // artificial upper bound for verbosity level

    cout << preproc << endl;

    //if ()


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

    if (varmap.count("stdout") && verbosity) {
      cout << "Outputting to standard out." << endl;
    }

    if (verbosity)
      cout << "Block size = " << block_size <<  "kB" << endl;

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
