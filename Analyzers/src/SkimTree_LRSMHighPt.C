#include "SkimTree_LRSMHighPt.h"

void SkimTree_LRSMHighPt::initializeAnalyzer(){

  outfile->cd();
  cout << "[SkimTree_LRSMHighPt::initializeAnalyzer()] gDirectory = " << gDirectory->GetName() << endl;
  newtree = fChain->CloneTree(0);

  triggers.clear();
  if(DataYear==2016){
    triggers = {
      "HLT_Mu50_v",
      "HLT_TkMu50_v",
      "HLT_Ele27_WPTight_Gsf_v",
      "HLT_Photon175_v",
      "HLT_Ele115_CaloIdVT_GsfTrkIdT_v",
      "HLT_DoublePhoton60_v",
    };
  }
  else if(DataYear==2017){
    triggers = {
      "HLT_Mu50_v",
      "HLT_OldMu100_v",
      "HLT_TkMu100_v",
      "HLT_Ele35_WPTight_Gsf_v",
      "HLT_Photon200_v",
      "HLT_Ele115_CaloIdVT_GsfTrkIdT_v",
      "HLT_DoublePhoton70_v",
    };
  }
  else if(DataYear==2018){
    triggers = {
      "HLT_Mu50_v",
      "HLT_OldMu100_v",
      "HLT_TkMu100_v",
      "HLT_Ele32_WPTight_Gsf_v",
      "HLT_Photon200_v",
      "HLT_Ele115_CaloIdVT_GsfTrkIdT_v",
      "HLT_DoublePhoton70_v",
    };
  }
  else{
    cout << "[SkimTree_LRSMHighPt::initializeAnalyzer] DataYear is wrong : " << DataYear << endl;
  }

  cout << "[SkimTree_LRSMHighPt::initializeAnalyzer] triggers to skim = " << endl;
  for(unsigned int i=0; i<triggers.size(); i++){
    cout << "[SkimTree_LRSMHighPt::initializeAnalyzer]   " << triggers.at(i) << endl;
  }


}

void SkimTree_LRSMHighPt::executeEvent(){

  Event ev;
  ev.SetTrigger(*HLT_TriggerName);

  if( ev.PassTrigger(triggers) ){
    newtree->Fill();
  }

}

void SkimTree_LRSMHighPt::executeEventFromParameter(AnalyzerParameter param){

}

SkimTree_LRSMHighPt::SkimTree_LRSMHighPt(){

  newtree = NULL;

}

SkimTree_LRSMHighPt::~SkimTree_LRSMHighPt(){

}

void SkimTree_LRSMHighPt::WriteHist(){

  outfile->mkdir("recoTree");
  outfile->cd("recoTree");
  newtree->Write();
  outfile->cd();

}


