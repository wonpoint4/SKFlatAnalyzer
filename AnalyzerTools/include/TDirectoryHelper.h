#ifndef TDirectoryHelperg_h
#define TDirectoryHelperg_h

#include <iostream>
#include <fstream>
#include <sstream>

#include "TROOT.h"
#include "TDirectory.h"

namespace TDirectoryHelper{

inline TDirectory* GetTempDirectory(std::string dirName){

  gROOT->cd();
  TDirectory* histDir = NULL;
  int counter = 0;
  while (!histDir) {
    //==== First, let's find a directory name that doesn't exist yet
    std::stringstream dirname;
    dirname << dirName << counter;
    if (gROOT->GetDirectory((dirname.str()).c_str())) {
      ++counter;
      continue;
    }
    //==== Let's try to make this directory
    histDir = gROOT->mkdir((dirname.str()).c_str());
  }
  std::cout << "[TDirectoryHelper::GetTempDirectory()] histDir name = " << histDir->GetName() << std::endl;
  return histDir;

}

} // END namespace TDirectoryHelper

#endif
