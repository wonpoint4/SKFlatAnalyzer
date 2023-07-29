#include "L1PrefireAnalyzer.h"

L1PrefireAnalyzer::L1PrefireAnalyzer(){
}
L1PrefireAnalyzer::~L1PrefireAnalyzer(){
}
void L1PrefireAnalyzer::initializeAnalyzer(){
  SMPAnalyzerCore::initializeAnalyzer(); //setup zpt roc z0
  fChain->SetBranchStatus("pfMET_*",false);
  fChain->SetBranchStatus("jet_*",false);
  fChain->SetBranchStatus("fatjet_*",false);
  fChain->SetBranchStatus("muon_*",false);
  fChain->SetBranchStatus("tau_*",false);
}
void L1PrefireAnalyzer::executeEvent(){
  _event.PassTrigger({"HLT_Ele27_WPTight_Gsf_v","HLT_Ele32_WPTight_Gsf_v"});
  vector<Photon> photons=GetAllPhotons();
  FillHist("photons",photons.size(),1.,10,0,10);

  vector<Electron> electrons=GetAllElectrons();
  FillHist("electrons",electrons.size(),1.,10,0,10);
  for(auto tag:electrons){
    if(!(tag.PassPath("HLT_Ele27_WPTight_Gsf_v")||tag.PassPath("HLT_Ele32_WPTight_Gsf_v"))) continue;
    for(auto probe:photons){
      if(probe.DeltaR(tag)<0.3) continue;
      FillHist("etapt",probe.Eta(),probe.Pt(),1.,netabin,etabins,nptbin,ptbins);
      if(IsUnprefirableEvent()){
	FillHist("etapt_unprefirable",probe.Eta(),probe.Pt(),1.,netabin,etabins,nptbin,ptbins);
      }
    }
  }
}
bool L1PrefireAnalyzer::IsUnprefirableEvent(){
  if(!fUnprefirableEvents.size()){
    TFile f("/gv0/Users/hsseo/L1Prefiring/UnprefirableEventList_SingleElectron_Run2017BtoF.root");
    TTree* tree=(TTree*)f.Get("tree");
    Long64_t this_run,this_event;
    tree->SetBranchAddress("run",&this_run);
    tree->SetBranchAddress("event",&this_event);
    for(int i=0,n=tree->GetEntries();i<n;i++){
      tree->GetEntry(i);
      fUnprefirableEvents.insert(make_pair(this_run,this_event));
    }
  }
  if(fUnprefirableEvents.find(make_pair(run,event))!=fUnprefirableEvents.end()) return true;
  return false;
}
    
