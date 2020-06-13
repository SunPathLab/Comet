/* (c) 2020 - Sun Ruping
   ruping@umn.edu
   Allocate Passenger mutations to each cell after tumopp-passenger run
   clang++ -o passenger passenger.cpp -lz -lboost_iostreams -std=c++14 */

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <map>
#include <string>
#include <cstring>
#include <sstream>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include "passenger.h"
#include <random>

using namespace std;

inline void splitstring(const string &str, vector<string> &elements, const string &delimiter);

int main ( int argc, char *argv[] ) {

  struct parameters *param = 0;
  param = interface(param, argc, argv);

  float urate = std::stof(param->urate);
  //std::poisson_distribution<int> poisson_distribution (urate);
  //std::uniform_int_distribution<int> uniform_distribution(1000001, 50000000); // define the range of mutational space, here a broad exome
  
  //uint_fast32_t seed = std::random_device{}();
  //uint_fast32_t seed2 = std::random_device{}();

  //using urbg_t = std::mt19937_64;
  //std::unique_ptr<urbg_t> engine_ = std::make_unique<urbg_t>(seed);
  //std::unique_ptr<urbg_t> engine2_ = std::make_unique<urbg_t>(seed2);

  std::cout << "id\tcoor\n";    // print out header
  
  std::ifstream passenger_f(param->passenger_f, std::ios_base::in | std::ios_base::binary);
  boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf_p;
  inbuf_p.push(boost::iostreams::gzip_decompressor());
  inbuf_p.push(passenger_f);
  //Convert streambuf to istream
  std::istream instream_p(&inbuf_p);
  //Iterate lines
  std::string line;
  bool firstline = true;
  
  while(std::getline(instream_p, line)) {  // each line of $passenger file

    if (firstline) {
      firstline = false;
      continue;
    }
    
    vector <string> line_content;
    splitstring(line, line_content, "\t");
    vector <string>::iterator iter = line_content.begin();
    unsigned int i;

    unsigned int cellid = 0;
    string passengers;
 
    for(i = 1; iter != line_content.end(); iter++, i++) {
      switch (i) {
      case 1:
        cellid = atoi((*iter).c_str());
        continue;
      default:
        break;
      }
    }
 
    // number of passenger mutations accordin to urate
    //unsigned int n_passengers = poisson_distribution(*engine_);
    
    //for(int n=0; n < n_passengers; ++n)
      // sample a coordinate
      //passengers = passengers + std::to_string(uniform_distribution(*engine_)) + ",";
    
    // print passenger mutations
    std::cout << cellid << "\t" << passengers << std::endl;
  }
  //Cleanup
  passenger_f.close();

  return 0;

} //main


inline void splitstring(const string &str, vector<string> &elements, const string &delimiter) {
  string::size_type lastPos = str.find_first_not_of(delimiter, 0);
  string::size_type pos     = str.find_first_of(delimiter, lastPos);

  while (string::npos != pos || string::npos != lastPos) {
    elements.push_back(str.substr(lastPos, pos - lastPos));
    lastPos = str.find_first_not_of(delimiter, pos);
    pos = str.find_first_of(delimiter, lastPos);
  }
}
