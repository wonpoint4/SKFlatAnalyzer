#include "SkimTree_LJL.h"

void SkimTree_LJL::initializeAnalyzer(){

  outfile->cd();
  cout << "[SkimTree_LJL::initializeAnalyzer()] gDirectory = " << gDirectory->GetName() << endl;
  newtree = fChain->CloneTree(0);

  triggers.clear();
  if(DataYear==2016){
    triggers = {
      "HLT_IsoMu24_v",
      "HLT_IsoTkMu24_v",
      "HLT_Ele27_WPTight_Gsf_v"
    };
  }else if(DataYear==2017){
    triggers = {
      "HLT_IsoMu27_v",
      "HLT_Ele32_WPTight_Gsf_L1DoubleEG_v"
    };
  }else if(DataYear==2018){
    triggers = {
      "HLT_IsoMu24_v",
      "HLT_Ele32_WPTight_Gsf_v"
    };
  }else{
    cout<<"[SkimTree_LJL::initializeAnalyzer] DataYear is wrong : " << DataYear << endl;
  }

  cout << "[SkimTree_LJL::initializeAnalyzer] triggers to skim = " << endl;
  for(unsigned int i=0; i<triggers.size(); i++){
    cout << "[SkimTree_LJL::initializeAnalyzer]   " << triggers.at(i) << endl;
  }  
}

void SkimTree_LJL::executeEvent(){

  Event ev;
  ev.SetTrigger(*HLT_TriggerName);

  if( ev.PassTrigger(triggers) ){
    vector<Muon> isomuons=SMPGetMuons("POGTightWithTightIso",25,2.4);
    vector<Muon> nonisomuons=SMPGetMuons("POGTightWithAntiMediumIso",15,2.4);
    if(isomuons.size()&&nonisomuons.size()){
      newtree->Fill();
    }
  }
}

void SkimTree_LJL::executeEventFromParameter(AnalyzerParameter param){

}

SkimTree_LJL::SkimTree_LJL(){
  newtree=NULL;
}

SkimTree_LJL::~SkimTree_LJL(){

}

void SkimTree_LJL::WriteHist(){

  outfile->mkdir("recoTree");
  outfile->cd("recoTree");
  newtree->Write();
  outfile->cd();

}


