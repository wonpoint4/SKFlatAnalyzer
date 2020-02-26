#include "SkimTree_HNFake.h"

void SkimTree_HNFake::initializeAnalyzer(){

  outfile->cd();
  cout << "[SkimTree_HNFake::initializeAnalyzer()] gDirectory = " << gDirectory->GetName() << endl;
  newtree = fChain->CloneTree(0);

  triggers.clear();
  if(DataYear==2016){
    triggers = {
      "HLT_Mu3_PFJet40_v",                             // DoubleMuon
      "HLT_Mu8_TrkIsoVVL_v",                           // DoubleMuon
      "HLT_Mu17_TrkIsoVVL_v",                          // DoubleMuon
      "HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v",     // DoubleEG
      "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v",    // DoubleEG
      "HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v",          // DoubleEG
      "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v"     // DoubleEG
    };
  }else if(DataYear==2017){
    triggers = {
      "HLT_Mu3_PFJet40_v",                             // SingleMuon
      "HLT_Mu8_TrkIsoVVL_v",                           // DoubleMuon
      "HLT_Mu17_TrkIsoVVL_v",                          // DoubleMuon
      "HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v",     // SingleElectron
      "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v",    // SingleElectron
      "HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v",          // SingleElectron
      "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v"     // SingleElectron
    };
  }else if(DataYear==2018){
    triggers = {
      "HLT_Mu3_PFJet40_v",                             // SingleMuon
      "HLT_Mu8_TrkIsoVVL_v",                           // DoubleMuon
      "HLT_Mu17_TrkIsoVVL_v",                          // DoubleMuon
      "HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v",     // EGamma
      "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v",    // EGamma
      "HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v",          // EGamma
      "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v"     // EGamma
    };
  }else{
    cout<<"[SkimTree_HNFake::initializeAnalyzer] DataYear is wrong : " << DataYear << endl;
  }

  cout << "[SkimTree_HNFake::initializeAnalyzer] triggers to skim = " << endl;
  for(unsigned int i=0; i<triggers.size(); i++){
    cout << "[SkimTree_HNFake::initializeAnalyzer]   " << triggers.at(i) << endl;
  }

}

void SkimTree_HNFake::executeEvent(){

  Event ev;
  ev.SetTrigger(*HLT_TriggerName);

  if( ev.PassTrigger(triggers) ){
    newtree->Fill();
  }

}

void SkimTree_HNFake::executeEventFromParameter(AnalyzerParameter param){

}

SkimTree_HNFake::SkimTree_HNFake(){
  newtree=NULL;
}

SkimTree_HNFake::~SkimTree_HNFake(){

}

void SkimTree_HNFake::WriteHist(){

  outfile->mkdir("recoTree");
  outfile->cd("recoTree");
  newtree->Write();
  outfile->cd();

}
