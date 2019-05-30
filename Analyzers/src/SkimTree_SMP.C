#include "SkimTree_SMP.h"

void SkimTree_SMP::initializeAnalyzer(){

  outfile->cd();
  cout << "[SkimTree_SMP::initializeAnalyzer()] gDirectory = " << gDirectory->GetName() << endl;
  newtree = fChain->CloneTree(0);

  triggers.clear();
  if(DataYear==2016){
    triggers = {
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v",
      "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
      "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"
    };
  }else if(DataYear==2017){
    triggers = {
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v",
      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"
    };
  }else if(DataYear==2018){
    triggers = {
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v",
      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"
    };
  }else{
    cout<<"[SkimTree_SMP::initializeAnalyzer] DataYear is wrong : " << DataYear << endl;
  }

  cout << "[SkimTree_SMP::initializeAnalyzer] triggers to skim = " << endl;
  for(unsigned int i=0; i<triggers.size(); i++){
    cout << "[SkimTree_SMP::initializeAnalyzer]   " << triggers.at(i) << endl;
  }

}

void SkimTree_SMP::executeEvent(){

  Event ev;
  ev.SetTrigger(*HLT_TriggerName);

  if( ev.PassTrigger(triggers) ){
    newtree->Fill();
  }

}

void SkimTree_SMP::executeEventFromParameter(AnalyzerParameter param){

}

SkimTree_SMP::SkimTree_SMP(){
  newtree=NULL;
}

SkimTree_SMP::~SkimTree_SMP(){

}

void SkimTree_SMP::WriteHist(){

  outfile->mkdir("recoTree");
  outfile->cd("recoTree");
  newtree->Write();
  outfile->cd();

}


