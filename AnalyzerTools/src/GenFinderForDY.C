#include "GenFinderForDY.h"

GenFinderForDY::GenFinderForDY(){

  FoundGenZ = false;
  GenZ = Particle();

}

GenFinderForDY::~GenFinderForDY(){
}

void GenFinderForDY::Find(const vector<Gen>& gens){

  GenZ = Particle();

  //==== Method 1) Following EXO-19-016, AN2019_019_v3

  vector<Gen> ZLeptons;
  MethodUsed = -1;
  for(unsigned int i=2; i<gens.size(); i++){

    if(ZLeptons.size()==2) break;

    Gen gen = gens.at(i);
    int pid = abs(gen.PID());
    if( pid==11 || pid==13 ){
      if( gen.Status() == 1 && gen.fromHardProcessFinalState() ){
        ZLeptons.push_back( gen );
      }
    }
    else if( pid==15 ){
      if( gen.Status() == 2 ){
        ZLeptons.push_back( gen );
      }
    }

  }
  if(ZLeptons.size()!=2){

    cout << "===========================================================" << endl;
    cout << "ZLeptons.size() = " << ZLeptons.size() << endl;
    cout << "index\tPID\tStatus\tMIdx\tMPID\tStart\tPt\tEta\tPhi\tM" << endl;
    for(unsigned int i=2; i<gens.size(); i++){
      Gen gen = gens.at(i);
      cout << i << "\t" << gen.PID() << "\t" << gen.Status() << "\t" << gen.MotherIndex() << "\t" << gens.at(gen.MotherIndex()).PID()<< "\t" << "-" << "\t";
      printf("%.2f\t%.2f\t%.2f\t%.2f\n",gen.Pt(), gen.Eta(), gen.Phi(), gen.M());
    }

  }
  else{

    //==== Validate

    bool IsSFOS = ( ZLeptons.at(0).PID()+ZLeptons.at(1).PID() ) == 0;
    if(IsSFOS){

      FoundGenZ = true;
      LeptonPID = abs(ZLeptons.at(0).PID());
      GenZ = ZLeptons.at(0)+ZLeptons.at(1);

    }

  }

}

bool GenFinderForDY::IsLepton(int pid){

  return ( abs(pid)==11 || abs(pid)==13 || abs(pid)==15 );

}


