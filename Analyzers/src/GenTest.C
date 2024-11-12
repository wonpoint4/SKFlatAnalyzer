#include "GenTest.h"

void GenTest::executeEvent(){
  genfid_b0=NULL;
  for(auto& lhe:lhes){
    if(lhe.Status()!=1) continue;
    if(abs(lhe.ID())!=5) continue;
    if(lhe.Pt()<40) continue;
    if(fabs(lhe.Eta())>2.4) continue;
    if(genfid_b0&&genfid_b0->Pt()>lhe.Pt()) continue;
    genfid_b0=&lhe;
  }  
  executeEventGen();
}
void GenTest::GetAFBGenParticles(const vector<Gen>& gens,Gen& parton0,Gen& parton1,Gen& l0,Gen& l1,int mode){
  //mode 0:bare 1:dressed01 2:dressed04 3:beforeFSR
  if(!IsDYSample&&!MCSample.Contains("GamGamToLL")&&!MCSample.Contains("TTLL")){
    cout <<"[SMPAnalyzerCore::GetAFBGenParticles] this is only for dilepton event"<<endl;
    exit(EXIT_FAILURE);
  }
  parton0=Gen();
  parton1=Gen();
  l0=Gen();
  l1=Gen();
  vector<const Gen*> leptons;
  vector<const Gen*> photons;
  int ngen=gens.size();
  for(int i=0;i<ngen;i++){
    if(!gens.at(i).isPrompt()) continue;
    int genpid=gens.at(i).PID();
    if(gens.at(i).isHardProcess()){
      if(abs(genpid)<7||genpid==21||genpid==22){
        if(parton0.IsEmpty()) parton0=gens[i];
        else if(parton1.IsEmpty()) parton1=gens[i];
      }
    }
    if(gens.at(i).Status()==1){
      if(abs(genpid)==11||abs(genpid)==13) leptons.push_back(&gens[i]);
      else if(gens.at(i).PID()==22) photons.push_back(&gens[i]);
    }
  }
  int nlepton=leptons.size();
  const double maxdr=0.4;
  double lhe_dimass=(lhe_l0+lhe_l1).M();
  if(MCSample.Contains("MiNNLO")&&lhe_dimass>1){
    double diff=999999;
    for(int i=0;i<nlepton;i++){
      for(int j=0;j<nlepton;j++){
	if(i==j) continue;
	if(leptons[i]->PID()!=lhe_l0.ID()) continue;
	if(leptons[j]->PID()!=lhe_l1.ID()) continue;
	double dimass=((*leptons[i])+(*leptons[j])).M();
	if(fabs(dimass-lhe_dimass)<diff){
	  l0=*leptons[i];
	  l1=*leptons[j];
	  diff=fabs(dimass-lhe_dimass);
	}
      }
    }
    //cout<<"diff:"<<diff<<endl;
  }else{
    for(int i=0;i<nlepton;i++){
      if(leptons[i]->PID()!=lhe_l0.ID()) continue;
      if(leptons[i]->DeltaR(lhe_l0)>maxdr) continue;
      if( fabs(leptons[i]->E()-lhe_l0.E()) < fabs(l0.E()-lhe_l0.E()) ){
	l0=*leptons[i];
      }
    }
    if(l0.PID()==0){
      for(int i=0;i<nlepton;i++){
	if(leptons[i]->PID()!=lhe_l0.ID()) continue;
	if(l0.PID()==0 || leptons[i]->DeltaR(lhe_l0)<l0.DeltaR(lhe_l0)){
	  l0=*leptons[i];
	}
      }
    }
    for(int i=0;i<nlepton;i++){
      if(leptons[i]->PID()!=lhe_l1.ID()) continue;
      if(leptons[i]->DeltaR(lhe_l1)>maxdr) continue;
      if( fabs(leptons[i]->E()-lhe_l1.E()) < fabs(l1.E()-lhe_l1.E()) ){
	l1=*leptons[i];
      }
    }
    if(l1.PID()==0){
      for(int i=0;i<nlepton;i++){
	if(leptons[i]->PID()!=lhe_l1.ID()) continue;
	if(l1.PID()==0 || leptons[i]->DeltaR(lhe_l1)<l1.DeltaR(lhe_l1)){
	  l1=*leptons[i];
	}
      }
    }
  }
  if(l0.Pt()<l1.Pt()){
    Gen tmp=l0;
    l0=l1;
    l1=tmp;
  }
  if(mode>=3){
    if(nlepton>=4){
      for(int i=0;i<nlepton;i++){
        if(leptons[i]->Index()==l0.Index()||leptons[i]->Index()==l1.Index()) continue;
        for(int j=i+1;j<nlepton;j++){
          if(leptons[j]->Index()==l0.Index()||leptons[j]->Index()==l1.Index()) continue;
          if(!(leptons[i]->PID()+leptons[j]->PID()==0)) continue;
          vector<int> history_i=TrackGenSelfHistory(*leptons[i],gens);
          vector<int> history_j=TrackGenSelfHistory(*leptons[j],gens);
          if(history_i.at(1)==history_j.at(1)){
            photons.push_back(leptons[i]);
            photons.push_back(leptons[j]);
          }
        }
      }
    }
    for(const auto& photon:photons){
      vector<int> history=TrackGenSelfHistory(*photon,gens);
      if(gens[history.at(1)].PID()==l0.PID()) l0+=*photon;
      else if(gens[history.at(1)].PID()==l1.PID()) l1+=*photon;
      else if(gens[history.at(1)].PID()==23){ // for minnlo+photos
        if(photon->DeltaR(l0)<photon->DeltaR(l1)) l0+=*photon;
        else l1+=*photon;
      }
    }
  }else if(mode>=1){
    double delr=mode==1?0.1:0.4;
    vector<const Gen*> toadd0,toadd1;
    for(const auto& photon:photons){
      if(l0.DeltaR(*photon)<delr) toadd0.push_back(photon);
      if(l1.DeltaR(*photon)<delr) toadd1.push_back(photon);
    }
    for(const auto& photon:toadd0) l0+=*photon;
    for(const auto& photon:toadd1) l1+=*photon;
  }
}

GenTest::GenTest(){}
GenTest::~GenTest(){}
