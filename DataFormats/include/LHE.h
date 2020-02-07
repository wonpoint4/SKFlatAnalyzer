#ifndef LHE_h
#define LHE_h

#include "Particle.h"

class LHE: public Particle {
public:

  LHE();
  ~LHE();

  void SetIndexIDStatus(int i, int p, int t);
  inline int Index() const { return j_index; }
  inline int ID() const { return j_ID; }
  inline int Status() const { return j_status; }

  virtual void Print();

private:

  int j_index, j_ID, j_status;

  ClassDef(LHE,1)
};

#endif
