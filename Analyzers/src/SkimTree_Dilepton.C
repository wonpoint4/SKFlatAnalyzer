#include "SkimTree_Dilepton.h"

void SkimTree_Dilepton::initializeAnalyzer(){

  outfile->cd();
  cout << "[SkimTree_Dilepton::initializeAnalyzer()] gDirectory = " << gDirectory->GetName() << endl;
  newtree = fChain->CloneTree(0);

  double_triggers.clear();
  if(DataYear==2016){
    double_triggers = {
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v",
      "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
      "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
      "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",
      "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",     // B-G
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",    // B-G
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",  // H
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"  // H
    };
    single_triggers = {
      "HLT_IsoMu24_v",
      "HLT_IsoTkMu24_v",
      "HLT_Ele27_WPTight_Gsf_v",
    };
  }else if(DataYear==2017){
    double_triggers = {
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v",
      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v",
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"
    };
    single_triggers = {
      "HLT_IsoMu24_v",
      "HLT_IsoMu27_v",
      "HLT_Ele27_WPTight_Gsf_v",
      "HLT_Ele32_WPTight_Gsf_v",
    };
  }else if(DataYear==2018){
    double_triggers = {
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v",
      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v",
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"
    };
    single_triggers = {
      "HLT_IsoMu24_v",
      "HLT_Ele27_WPTight_Gsf_v",
      "HLT_Ele28_WPTight_Gsf_v",
      "HLT_Ele32_WPTight_Gsf_v",
    };
  }else{
    cout<<"[SkimTree_Dilepton::initializeAnalyzer] DataYear is wrong : " << DataYear << endl;
  }

  cout << "[SkimTree_Dilepton::initializeAnalyzer] triggers to skim = " << endl;
  for(unsigned int i=0; i<double_triggers.size(); i++){
    cout << "[SkimTree_Dilepton::initializeAnalyzer]   " << double_triggers.at(i) << endl;
  }
  for(unsigned int i=0; i<single_triggers.size(); i++){
    cout << "[SkimTree_Dilepton::initializeAnalyzer]   " << single_triggers.at(i) << endl;
  }

}

void SkimTree_Dilepton::executeEvent(){

  Event ev;
  ev.SetTrigger(*HLT_TriggerName);

  if( ev.PassTrigger(double_triggers) ){
    newtree->Fill();
  }else if(ev.PassTrigger(single_triggers)){
    bool save=false;
    vector<Muon> muons=GetAllMuons();
    std::sort(muons.begin(),muons.end(),PtComparing);
    if(muons.size()>=2)
      if(muons.at(0).Pt()>22 && muons.at(1).Pt()>7)
	save=true;

    vector<Electron> electrons=GetAllElectrons();
    std::sort(electrons.begin(),electrons.end(),PtComparing);
    if(electrons.size()>=2)
      if(electrons.at(0).Pt()>27 && electrons.at(1).Pt()>7)
	save=true;

    if(save) newtree->Fill();
  }

}

void SkimTree_Dilepton::executeEventFromParameter(AnalyzerParameter param){

}

SkimTree_Dilepton::SkimTree_Dilepton(){
  newtree=NULL;
}

SkimTree_Dilepton::~SkimTree_Dilepton(){

}

void SkimTree_Dilepton::WriteHist(){

  outfile->mkdir("recoTree");
  outfile->cd("recoTree");
  newtree->Write();
  outfile->cd();

}
