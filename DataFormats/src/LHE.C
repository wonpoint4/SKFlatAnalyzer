#include "LHE.h"

ClassImp(LHE)

LHE::LHE() : Particle() {
  j_index = -1;
  j_ID = 0;
  j_status = 0;
}

LHE::~LHE(){
}

void LHE::SetIndexIDStatus(int i, int p, int t){
  j_index = i;
  j_ID = p;
  j_status = t;
}

void LHE::Print(){

  cout << "(Index, ID, Status, Pt, Eta, Phi, M) = " << Index() << "\t" << ID() << "\t" << Status() << "\t" << Pt() << "\t" << Eta() << "\t" << Phi() << "\t" << M() << endl;

}
