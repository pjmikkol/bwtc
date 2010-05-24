#include<iostream>
#include<iterator>
using namespace std;

#include<boost/program_options.hpp>
namespace po = boost::program_options;

int main(int argc, char** argv) {
  try {

    po::options_description description("Options");
    description.add_options()
      ("help,h", "print help message")
      ("keep,k", "keep (don't remove) input files")
      ("stdout,c", "output to standard out")
      ;
    
    po::variables_map varmap;
    po::store(po::parse_command_line(argc, argv, description), varmap);
    po::notify(varmap);

    if (varmap.count("help")) {
      cout << description << endl;
      return 1;
    }

    if (varmap.count("keep")) {
      cout << "Keeping input files." << endl; 
    }

    if (varmap.count("stdout")) {
      cout << "Outputting to standard out." << endl;
    }
  }

  catch(exception& e) {
    cerr << "error: " << e.what() << endl;
    return 1;
  }
  catch(...) {
    cerr << "Exception of unknown type!" << endl;
  }

  return 0;

}
