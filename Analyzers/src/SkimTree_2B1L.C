#include "SkimTree_2B1L.h"

void SkimTree_2B1L::initializeAnalyzer(){

  outfile->cd();
  cout << "[SkimTree_2B1L::initializeAnalyzer()] gDirectory = " << gDirectory->GetName() << endl;
  newtree = fChain->CloneTree(0);

  single_electron_triggers.clear();
  single_muon_triggers.clear();
  if(DataYear==2016){
    single_muon_triggers = {
      "HLT_IsoMu24_v",
      "HLT_IsoTkMu24_v",
    };
    single_electron_triggers = {
      "HLT_Ele27_WPTight_Gsf_v",
    };
  }else if(DataYear==2017){
    single_muon_triggers = {
      "HLT_IsoMu24_v",
      "HLT_IsoMu27_v",
    };
    single_electron_triggers = {
      "HLT_Ele27_WPTight_Gsf_v",
      "HLT_Ele32_WPTight_Gsf_v",
      "HLT_Ele32_WPTight_Gsf_L1DoubleEG_v",
    };
  }else if(DataYear==2018){
    single_muon_triggers = {
      "HLT_IsoMu24_v",
    };
    single_electron_triggers = {
      "HLT_Ele27_WPTight_Gsf_v",
      "HLT_Ele28_WPTight_Gsf_v",
      "HLT_Ele32_WPTight_Gsf_v",
    };
  }else{
    cout<<"[SkimTree_2B1L::initializeAnalyzer] DataYear is wrong : " << DataYear << endl;
  }

  mcCorr->SetJetTaggingParameters({
    JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Tight, JetTagging::incl, JetTagging::comb),
  });

  cout << "[SkimTree_2B1L::initializeAnalyzer] triggers to skim = " << endl;
  for(unsigned int i=0; i<single_muon_triggers.size(); i++){
    cout << "[SkimTree_2B1L::initializeAnalyzer]   " << single_muon_triggers.at(i) << endl;
  }
  for(unsigned int i=0; i<single_electron_triggers.size(); i++){
    cout << "[SkimTree_2B1L::initializeAnalyzer]   " << single_electron_triggers.at(i) << endl;
  }

}

void SkimTree_2B1L::executeEvent(){

  Event ev;
  ev.SetTrigger(*HLT_TriggerName);

  if(!(ev.PassTrigger(single_muon_triggers) || ev.PassTrigger(single_electron_triggers))) return;

  vector<Jet> jets = SelectJets(GetAllJets(), "tightLepVeto", 20, (DataYear == 2016? 2.4: 2.5));
  std::sort(jets.begin(), jets.end(), PtComparing);
  JetTagging::Parameters DeepJet_Tight = JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Tight, JetTagging::incl, JetTagging::comb);
  std::vector<Jet> bjets = {};
  for(const auto& jet:jets){
    if(jet.GetTaggerResult(DeepJet_Tight.j_Tagger) > mcCorr->GetJetTaggingCutValue(DeepJet_Tight.j_Tagger, DeepJet_Tight.j_WP)) bjets.push_back(jet);
  }

  if(bjets.size() < 2) return;

  newtree->Fill();
}

void SkimTree_2B1L::executeEventFromParameter(AnalyzerParameter param){

}

SkimTree_2B1L::SkimTree_2B1L(){
  newtree=NULL;
}

SkimTree_2B1L::~SkimTree_2B1L(){

}

void SkimTree_2B1L::WriteHist(){

  outfile->mkdir("recoTree");
  outfile->cd("recoTree");
  newtree->Write();
  outfile->cd();

}
