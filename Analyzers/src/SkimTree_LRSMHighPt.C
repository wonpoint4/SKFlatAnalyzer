#include "SkimTree_LRSMHighPt.h"

void SkimTree_LRSMHighPt::initializeAnalyzer(){

  outfile->cd();
  cout << "[SkimTree_LRSMHighPt::initializeAnalyzer()] gDirectory = " << gDirectory->GetName() << endl;
  newtree = fChain->CloneTree(0);

  if( MCSample.Contains("DY") ){
    cout << "[SkimTree_LRSMHighPt::initializeAnalyzer()] This is DY sample, so saving gen_* and LHE_* for Z-pt reweighting" << endl;
  }
  else{
    cout << "[SkimTree_LRSMHighPt::initializeAnalyzer()] Throwing away gen_* and LHE_*" << endl;
    if(!IsDATA){
      newtree->SetBranchStatus("gen_*",0);
      newtree->SetBranchStatus("LHE_*",0);
      newtree->SetBranchStatus("gen_weight",1); // for MCweight()
    }
  }

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

  //==== WR (Z') : lepton pt > 53 (65) GeV
  LeptonPtCut = 40.;
  //==== WR (Z') : jet pt > 40 GeV
  AK4JetPtCut = 30.;
  //==== WR (Z') : Jet pt > 200 (300) GeV
  AK8JetPtCut = 170.;

  cout << "[SkimTree_LRSMHighPt::initializeAnalyzer] LeptonPtCut = " << LeptonPtCut << endl;
  cout << "[SkimTree_LRSMHighPt::initializeAnalyzer] AK4JetPtCut = " << AK4JetPtCut << endl;
  cout << "[SkimTree_LRSMHighPt::initializeAnalyzer] AK8JetPtCut = " << AK8JetPtCut << endl;

}

void SkimTree_LRSMHighPt::executeEvent(){

  Event ev;
  ev.SetTrigger(*HLT_TriggerName);

  //==== Skim 1 ) trigger
  if(! (ev.PassTrigger(triggers)) ) return;

  //==== Skim 2) at least two leptons (e or mu) with pt > "LeptonPtCut"

  vector<Muon> allmuons = UseTunePMuon( GetAllMuons() );
  int Nmuon = 0;
  for(unsigned int i=0; i<allmuons.size(); i++){
    if( allmuons.at(i).Pt() > LeptonPtCut ) Nmuon++;
  }
  vector<Electron> electrons = GetElectrons("NOCUT", LeptonPtCut, 2.4);
  int Nelectron = electrons.size();
  if( Nmuon+Nelectron < 2 ) return;

  //==== Skim 3) Jets

  vector<Jet> jets = GetJets("tightLepVeto", AK4JetPtCut, 2.4);
  int Njet = jets.size();

  vector<FatJet> allfatjets = puppiCorr->Correct( GetFatJets("tight", AK8JetPtCut, 2.4) ); //==== corret SDMass
  int Nfatjet_SDMassCut = 0;
  for(unsigned int i=0; i<allfatjets.size(); i++){
    if(allfatjets.at(i).SDMass()>40.) Nfatjet_SDMassCut++;
  }

  //==== WR
  bool WRJetReq = (Njet>=2) || (Nfatjet_SDMassCut>=1);

  //==== Z'
  bool ZpJetReq = (Njet+Nfatjet_SDMassCut)>=2;

  if( !(WRJetReq||ZpJetReq) ) return;

  //=============================
  //==== If survived, fill tree
  //=============================

  newtree->Fill();

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


