#include "BTAnalyzer.h"
#include <iostream>
#include <string>

int main(int argc, char* argv[]) {
  if(argc < 3) {
    std::cout << "usage: ./btanalysis <inFilename> <sensorType> <outFilename>\n";
    std::exit(EXIT_FAILURE);  
  }

  std::string iFile = argv[1];
  std::cout << " File name " << iFile << std::endl;
  //unsigned int nmodules = std::stoi(argv[2]);
  //std::string sensorType = argv[2];
  std::string oFile = argv[2];


  //BTAnalyzer *analyzer = new  BTAnalyzer(iFile.c_str(), nmodules, sensorType.c_str(), oFile.c_str());
  //BTAnalyzer *analyzer = new  BTAnalyzer(iFile.c_str(), sensorType.c_str(), oFile.c_str());
  BTAnalyzer *analyzer = new  BTAnalyzer(iFile.c_str(), oFile.c_str());
  std::cout << "Executing event loop!" << std::endl;
  analyzer->Loop();
  std::cout << "Saving Histos!" << std::endl;
  analyzer->SaveHistos();

  if(analyzer)   delete analyzer;

  return 0;
}
