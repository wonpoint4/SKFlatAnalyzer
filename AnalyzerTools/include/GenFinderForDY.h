#ifndef GenFinderForDY_h
#define GenFinderForDY_h

#include "TString.h"
#include <iostream>
#include "Gen.h"

using namespace std;

class GenFinderForDY{

public:

  GenFinderForDY();
  ~GenFinderForDY();

  bool FoundGenZ;
  Particle GenZ;
  void Find(const vector<Gen>& gens);
  int MethodUsed;
  int LeptonPID;
  bool IsLepton(int pid);

};

#endif
