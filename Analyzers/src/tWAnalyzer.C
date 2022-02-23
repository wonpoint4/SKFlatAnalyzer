#include "tWAnalyzer.h"

void tWAnalyzer::initializeAnalyzer(){

  MuonIDs = { "DeepJetM_MuM", "DeepJetM_MuT", "DeepJetT_MuM", "DeepJetT_MuT" };

  if(DataEra=="2018"){
    IsoMuTriggerName = "HLT_IsoMu24_v";
    TriggerSafePtCut = 26.;
  }
  else if(DataEra=="2017"){
    IsoMuTriggerName = "HLT_IsoMu27_v";
    TriggerSafePtCut = 29.;
  }

  cout << "[tWAnalyzer::initializeAnalyzer] IsoMuTriggerName = " << IsoMuTriggerName << endl;
  cout << "[tWAnalyzer::initializeAnalyzer TriggerSafePtCut = " << TriggerSafePtCut << endl;

  //==== B-Tagging
  //==== add taggers and WP that you want to use in analysis
  std::vector<JetTagging::Parameters> jtps;
  //==== If you want to use 1a or 2a method,
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Loose, JetTagging::incl, JetTagging::comb) );
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::comb) );
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Tight, JetTagging::incl, JetTagging::comb) );
  //==== set
  mcCorr->SetJetTaggingParameters(jtps);

}

void tWAnalyzer::executeEvent(){

  weight_Prefire = GetPrefireWeight(0);
  AnalyzerParameter param;

  for(unsigned int it_MuonID=0; it_MuonID<MuonIDs.size(); it_MuonID++){

    TString MuonID = MuonIDs.at(it_MuonID).Contains("MuM")? "POGMediumWithLooseTrkIso": "POGTightWithTightIso";
    TString MuonIDSFKey = "ID_SF_MediumID_trkIsoLoose_Q";
    TString MuonTrigSFKey = DataEra=="2017"? "IsoMu27_MediumID_trkIsoLoose_Q": "IsoMu24_MediumID_trkIsoLoose_Q";

    //==== clear parameter set
    param.Clear();

    param.syst_ = AnalyzerParameter::Central;
    param.Name = MuonIDs.at(it_MuonID)+GetEra();
    //==== You can define lepton ID string here
    param.Muon_Tight_ID = MuonID;
    param.Muon_ID_SF_Key = MuonIDSFKey;
    param.Muon_Trigger_SF_Key = MuonTrigSFKey;
    param.Jet_ID = "tightLepVeto";
    param.Muon_Veto_ID = "POGLoose";
    param.Electron_Veto_ID = "passLooseID";

    //==== Now, all parameters are set. Run executeEventFromParameter() with this parameter set
    executeEventFromParameter(param);
  }
}

void tWAnalyzer::executeEventFromParameter(AnalyzerParameter param){

  //=============
  //==== No Cut
  //=============
  TString prefix = param.Name+"/";
  TString hprefix = "Nocut";
  FillHist(prefix+hprefix+"", 0., 1., 1, 0., 1.);

  //========================
  //==== MET Filter
  //========================

  if(!PassMETFilter()) return;

  Event ev = GetEvent();
  Particle METv = ev.GetMETVector();
  TLorentzVector met = GetEvent().GetMETVector();
  double weight = 1.;
  if(!IsDATA){
    weight *= weight_norm_1invpb*ev.GetTriggerLumi("Full");
    weight *= ev.MCweight();
  }
  if(param.syst_ == AnalyzerParameter::Central){
  }
  else{
    cout << "[tWAnalyzer::executeEventFromParameter] Wrong syst" << endl;
    exit(EXIT_FAILURE);
  }

  //==================================================
  //==== Then, apply ID selections using this_AllXXX
  //==================================================

  JetTagging::Tagger btagger = JetTagging::DeepJet;
  JetTagging::WP btagLoose = JetTagging::Loose;
  JetTagging::WP btagWP = JetTagging::Medium;
  if(param.Name.Contains("DeepJetT")) btagWP = JetTagging::Tight;

  vector<Muon> muons = SMPGetMuons(param.Muon_Tight_ID, 20., 2.4);
  vector<Muon> vetomus = SMPGetMuons(param.Muon_Veto_ID, 20., 2.4);
  vector<Electron> vetoels = SMPGetElectrons(param.Electron_Veto_ID, 20., 2.4);
  vector<Jet> basicjets = GetJets(param.Jet_ID, 30., 2.4);

  std::sort(muons.begin(), muons.end(), PtComparing);
  std::sort(basicjets.begin(), basicjets.end(), PtComparing);

  vector<Jet> jets, bjets, loosebjets;
  jets.clear();
  bjets.clear();
  loosebjets.clear();
  for(unsigned int i=0; i<basicjets.size(); i++){
    double this_discr = basicjets.at(i).GetTaggerResult(btagger);
    if(this_discr > mcCorr->GetJetTaggingCutValue(btagger, btagWP)) bjets.push_back(basicjets.at(i));
    else{
      if(this_discr > mcCorr->GetJetTaggingCutValue(btagger, btagLoose)) loosebjets.push_back(basicjets.at(i));
      jets.push_back(basicjets.at(i));
    }
  }

  if(muons.size() != 1) return;
  FillCutflow(prefix+"cutflow","oneMu",weight);
  hprefix = "oneMu_";

  if(! (ev.PassTrigger(IsoMuTriggerName) )) return;
  if(muons.at(0).Pt() <= TriggerSafePtCut ) return;
  FillCutflow(prefix+"cutflow","trigger",weight);
  hprefix = "trigger_";

  FillHist(prefix+hprefix+"Nvetomu", vetomus.size(), weight, 15, 0, 15);
  FillHist(prefix+hprefix+"Nvetoel", vetoels.size(), weight, 15, 0, 15);
  FillHist(prefix+hprefix+"Nmuon", muons.size(), weight, 15, 0, 15);
  FillHist(prefix+hprefix+"MET", met.Pt(), weight, 60, 0., 300.);
  FillHist(prefix+hprefix+"mupT", muons.at(0).Pt(), weight, 60, 0., 300.);
  //FillHist(prefix+hprefix+"bpT", bjets.at(0).Pt(), weight, 60, 0., 300.);
  //FillHist(prefix+hprefix+"j0pT", jets.at(0).Pt(), weight, 60, 0., 300.);
  //FillHist(prefix+hprefix+"j1pT", jets.at(1).Pt(), weight, 60, 0., 300.);
  //FillHist(prefix+hprefix+"mWj", (jets.at(0)+jets.at(1)).M(), weight, 60, 0., 300.);
  //FillHist(prefix+hprefix+"mTop", (bjets.at(0)+jets.at(0)+jets.at(1)).M(), weight, 80, 0., 400.);
  //FillHist(prefix+hprefix+"mTopl", (bjets.at(0)+muons.at(0)+met).M(), weight, 80, 0., 400.);
  FillHist(prefix+hprefix+"mW", (muons.at(0)+met).M(), weight, 60, 0., 300.);
  if(muons.at(0).Charge() > 0) FillHist(prefix+hprefix+"mWp", (muons.at(0)+met).M(), weight, 60, 0., 300.);
  else FillHist(prefix+hprefix+"mWm", (muons.at(0)+met).M(), weight, 60, 0., 300.);

  if(vetomus.size()+vetoels.size() > 1) return;
  FillCutflow(prefix+"cutflow","lepVeto",weight);
  hprefix = "lepVeto_";

  FillHist(prefix+hprefix+"Nvetomu", vetomus.size(), weight, 15, 0, 15);
  FillHist(prefix+hprefix+"Nvetoel", vetoels.size(), weight, 15, 0, 15);
  FillHist(prefix+hprefix+"Nmuon", muons.size(), weight, 15, 0, 15);
  FillHist(prefix+hprefix+"Njet", jets.size(), weight, 15, 0, 15);
  FillHist(prefix+hprefix+"Nbjet", bjets.size(), weight, 10, 0, 10);
  FillHist(prefix+hprefix+"Nloosebjet", loosebjets.size(), weight, 10, 0, 10);
  FillHist(prefix+hprefix+"MET", met.Pt(), weight, 60, 0., 300.);
  FillHist(prefix+hprefix+"mupT", muons.at(0).Pt(), weight, 60, 0., 300.);
  //FillHist(prefix+hprefix+"bpT", bjets.at(0).Pt(), weight, 60, 0., 300.);
  //FillHist(prefix+hprefix+"j0pT", jets.at(0).Pt(), weight, 60, 0., 300.);
  //FillHist(prefix+hprefix+"j1pT", jets.at(1).Pt(), weight, 60, 0., 300.);
  //FillHist(prefix+hprefix+"mWj", (jets.at(0)+jets.at(1)).M(), weight, 60, 0., 300.);
  //FillHist(prefix+hprefix+"mTop", (bjets.at(0)+jets.at(0)+jets.at(1)).M(), weight, 80, 0., 400.);
  //FillHist(prefix+hprefix+"mTopl", (bjets.at(0)+muons.at(0)+met).M(), weight, 80, 0., 400.);
  FillHist(prefix+hprefix+"mW", (muons.at(0)+met).M(), weight, 60, 0., 300.);
  if(muons.at(0).Charge() > 0) FillHist(prefix+hprefix+"mWp", (muons.at(0)+met).M(), weight, 60, 0., 300.);
  else FillHist(prefix+hprefix+"mWm", (muons.at(0)+met).M(), weight, 60, 0., 300.);

  if(bjets.size() != 1) return;
  FillCutflow(prefix+"cutflow","1bjet",weight);
  hprefix = "1bjet_";

  FillHist(prefix+hprefix+"Njet", jets.size(), weight, 15, 0, 15);
  FillHist(prefix+hprefix+"Nbjet", bjets.size(), weight, 10, 0, 10);
  FillHist(prefix+hprefix+"Nloosebjet", loosebjets.size(), weight, 10, 0, 10);
  FillHist(prefix+hprefix+"MET", met.Pt(), weight, 60, 0., 300.);
  FillHist(prefix+hprefix+"mupT", muons.at(0).Pt(), weight, 60, 0., 300.);
  FillHist(prefix+hprefix+"bpT", bjets.at(0).Pt(), weight, 60, 0., 300.);
  //FillHist(prefix+hprefix+"j0pT", jets.at(0).Pt(), weight, 60, 0., 300.);
  //FillHist(prefix+hprefix+"j1pT", jets.at(1).Pt(), weight, 60, 0., 300.);
  //FillHist(prefix+hprefix+"mWj", (jets.at(0)+jets.at(1)).M(), weight, 60, 0., 300.);
  //FillHist(prefix+hprefix+"mTop", (bjets.at(0)+jets.at(0)+jets.at(1)).M(), weight, 80, 0., 400.);
  FillHist(prefix+hprefix+"mTopl", (bjets.at(0)+muons.at(0)+met).M(), weight, 80, 0., 400.);
  FillHist(prefix+hprefix+"mW", (muons.at(0)+met).M(), weight, 60, 0., 300.);
  if(muons.at(0).Charge() > 0) FillHist(prefix+hprefix+"mWp", (muons.at(0)+met).M(), weight, 60, 0., 300.);
  else FillHist(prefix+hprefix+"mWm", (muons.at(0)+met).M(), weight, 60, 0., 300.);

  if(loosebjets.size() != 0) return;
  FillCutflow(prefix+"cutflow","loosebVeto",weight);
  hprefix = "loosebVeto_";

  FillHist(prefix+hprefix+"Njet", jets.size(), weight, 15, 0, 15);
  FillHist(prefix+hprefix+"Nbjet", bjets.size(), weight, 10, 0, 10);
  FillHist(prefix+hprefix+"Nloosebjet", loosebjets.size(), weight, 10, 0, 10);
  FillHist(prefix+hprefix+"MET", met.Pt(), weight, 60, 0., 300.);
  FillHist(prefix+hprefix+"mupT", muons.at(0).Pt(), weight, 60, 0., 300.);
  FillHist(prefix+hprefix+"bpT", bjets.at(0).Pt(), weight, 60, 0., 300.);
  //FillHist(prefix+hprefix+"j0pT", jets.at(0).Pt(), weight, 60, 0., 300.);
  //FillHist(prefix+hprefix+"j1pT", jets.at(1).Pt(), weight, 60, 0., 300.);
  //FillHist(prefix+hprefix+"mWj", (jets.at(0)+jets.at(1)).M(), weight, 60, 0., 300.);
  //FillHist(prefix+hprefix+"mTop", (bjets.at(0)+jets.at(0)+jets.at(1)).M(), weight, 80, 0., 400.);
  FillHist(prefix+hprefix+"mTopl", (bjets.at(0)+muons.at(0)+met).M(), weight, 80, 0., 400.);
  FillHist(prefix+hprefix+"mW", (muons.at(0)+met).M(), weight, 60, 0., 300.);
  if(muons.at(0).Charge() > 0) FillHist(prefix+hprefix+"mWp", (muons.at(0)+met).M(), weight, 60, 0., 300.);
  else FillHist(prefix+hprefix+"mWm", (muons.at(0)+met).M(), weight, 60, 0., 300.);

  if(jets.size() != 2) return;
  FillCutflow(prefix+"cutflow","2jet",weight);
  hprefix = "2jet_";

  FillHist(prefix+hprefix+"Njet", jets.size(), weight, 15, 0, 15);
  FillHist(prefix+hprefix+"Nbjet", bjets.size(), weight, 10, 0, 10);
  FillHist(prefix+hprefix+"Nloosebjet", loosebjets.size(), weight, 10, 0, 10);
  FillHist(prefix+hprefix+"MET", met.Pt(), weight, 60, 0., 300.);
  FillHist(prefix+hprefix+"mupT", muons.at(0).Pt(), weight, 60, 0., 300.);
  FillHist(prefix+hprefix+"bpT", bjets.at(0).Pt(), weight, 60, 0., 300.);
  FillHist(prefix+hprefix+"j0pT", jets.at(0).Pt(), weight, 60, 0., 300.);
  FillHist(prefix+hprefix+"j1pT", jets.at(1).Pt(), weight, 60, 0., 300.);
  FillHist(prefix+hprefix+"mWj", (jets.at(0)+jets.at(1)).M(), weight, 60, 0., 300.);
  FillHist(prefix+hprefix+"mTop", (bjets.at(0)+jets.at(0)+jets.at(1)).M(), weight, 80, 0., 400.);
  FillHist(prefix+hprefix+"mTopl", (bjets.at(0)+muons.at(0)+met).M(), weight, 80, 0., 400.);
  FillHist(prefix+hprefix+"mW", (muons.at(0)+met).M(), weight, 60, 0., 300.);
  if(muons.at(0).Charge() > 0) FillHist(prefix+hprefix+"mWp", (muons.at(0)+met).M(), weight, 60, 0., 300.);
  else FillHist(prefix+hprefix+"mWm", (muons.at(0)+met).M(), weight, 60, 0., 300.);

  //===================
  //==== Event weight
  //===================

  //==== If MC
  if(!IsDATA){
    weight *= weight_Prefire;
    FillCutflow(prefix+"cutflow","prefire",weight);

    if(MCSample.Contains("SingleTop_tW")){
      vector<LHE> lhes = GetLHEs();
      //PrintLHEs(lhes);
      vector<Gen> gens = GetGens();
      for(unsigned int i=0;i<gens.size();i++){
        //gens.at(i).Print();
        if(gens.at(i).isHardProcess() && abs(gens.at(i).PID())==13){
          Gen MuMother = gens.at(gens.at(i).MotherIndex());
          if(abs(MuMother.PID()) != 24) cout<<"MuMother isn't W boson, Really? "<<endl;
          FillHist(prefix+"muMotherPID", MuMother.PID(), 1, 60, -30, 30);
          Gen MuGrandMother = gens.at(MuMother.MotherIndex());
          while(abs(MuGrandMother.PID()) == 24){
            MuGrandMother = gens.at(MuGrandMother.MotherIndex());
          }
          FillHist(prefix+"muGrandMotherPID", MuGrandMother.PID(), 1, 60, -30, 30);
        }
      }
    }

    //==== Example of applying Muon scale factors
    for(unsigned int i=0; i<muons.size(); i++){
      Lepton *l = (Lepton *)(&muons.at(i));
      double this_idsf  = Lepton_SF(param.Muon_ID_SF_Key, l, 0);
      double this_isosf = 1.;
      weight *= this_idsf*this_isosf;
    }
    FillCutflow(prefix+"cutflow","IDSF",weight);
    double this_trigsf = LeptonTrigger_SF(param.Muon_Trigger_SF_Key, MakeLeptonPointerVector(muons), 0);
    weight *= this_trigsf;
    FillCutflow(prefix+"cutflow","trigSF",weight);

    JetTagging::Parameters jtp = JetTagging::Parameters(btagger, btagWP, JetTagging::incl, JetTagging::comb);
    double btagweight = mcCorr->GetBTaggingReweight_1a(jets, jtp);
    weight *= btagweight;
    FillCutflow(prefix+"cutflow","btagSF",weight);
  }else{
    FillCutflow(prefix+"cutflow","prefire",weight);
    FillCutflow(prefix+"cutflow","IDSF",weight);
    FillCutflow(prefix+"cutflow","trigSF",weight);
    FillCutflow(prefix+"cutflow","btagSF",weight);
  }

  hprefix = "";
  FillHist(prefix+hprefix+"MET", met.Pt(), weight, 60, 0., 300.);
  FillHist(prefix+hprefix+"mupT", muons.at(0).Pt(), weight, 60, 0., 300.);
  FillHist(prefix+hprefix+"bpT", bjets.at(0).Pt(), weight, 60, 0., 300.);
  FillHist(prefix+hprefix+"j0pT", jets.at(0).Pt(), weight, 60, 0., 300.);
  FillHist(prefix+hprefix+"j1pT", jets.at(1).Pt(), weight, 60, 0., 300.);
  FillHist(prefix+hprefix+"mWj", (jets.at(0)+jets.at(1)).M(), weight, 60, 0., 300.);
  FillHist(prefix+hprefix+"mTop", (bjets.at(0)+jets.at(0)+jets.at(1)).M(), weight, 80, 0., 400.);
  FillHist(prefix+hprefix+"mTopl", (bjets.at(0)+muons.at(0)+met).M(), weight, 80, 0., 400.);
  FillHist(prefix+hprefix+"mW", (muons.at(0)+met).M(), weight, 60, 0., 300.);
  if(muons.at(0).Charge() > 0) FillHist(prefix+hprefix+"mWp", (muons.at(0)+met).M(), weight, 60, 0., 300.);
  else FillHist(prefix+hprefix+"mWm", (muons.at(0)+met).M(), weight, 60, 0., 300.);

  if(fabs(80.4 - (jets.at(0)+jets.at(1)).M()) > 15) return;
  FillCutflow(prefix+"cutflow","Whad15",weight);
  hprefix = "Whad15_";

  FillHist(prefix+hprefix+"MET", met.Pt(), weight, 60, 0., 300.);
  FillHist(prefix+hprefix+"mupT", muons.at(0).Pt(), weight, 60, 0., 300.);
  FillHist(prefix+hprefix+"bpT", bjets.at(0).Pt(), weight, 60, 0., 300.);
  FillHist(prefix+hprefix+"j0pT", jets.at(0).Pt(), weight, 60, 0., 300.);
  FillHist(prefix+hprefix+"j1pT", jets.at(1).Pt(), weight, 60, 0., 300.);
  FillHist(prefix+hprefix+"mWj", (jets.at(0)+jets.at(1)).M(), weight, 60, 0., 300.);
  FillHist(prefix+hprefix+"mTop", (bjets.at(0)+jets.at(0)+jets.at(1)).M(), weight, 80, 0., 400.);
  FillHist(prefix+hprefix+"mTopl", (bjets.at(0)+muons.at(0)+met).M(), weight, 80, 0., 400.);
  FillHist(prefix+hprefix+"mW", (muons.at(0)+met).M(), weight, 60, 0., 300.);
  if(muons.at(0).Charge() > 0) FillHist(prefix+hprefix+"mWp", (muons.at(0)+met).M(), weight, 60, 0., 300.);
  else FillHist(prefix+hprefix+"mWm", (muons.at(0)+met).M(), weight, 60, 0., 300.);

  if(fabs(172.5 - (bjets.at(0)+jets.at(0)+jets.at(1)).M()) > 30) return;
  FillCutflow(prefix+"cutflow","Top30",weight);
  hprefix = "Top30_";

  FillHist(prefix+hprefix+"MET", met.Pt(), weight, 60, 0., 300.);
  FillHist(prefix+hprefix+"mupT", muons.at(0).Pt(), weight, 60, 0., 300.);
  FillHist(prefix+hprefix+"bpT", bjets.at(0).Pt(), weight, 60, 0., 300.);
  FillHist(prefix+hprefix+"j0pT", jets.at(0).Pt(), weight, 60, 0., 300.);
  FillHist(prefix+hprefix+"j1pT", jets.at(1).Pt(), weight, 60, 0., 300.);
  FillHist(prefix+hprefix+"mWj", (jets.at(0)+jets.at(1)).M(), weight, 60, 0., 300.);
  FillHist(prefix+hprefix+"mTop", (bjets.at(0)+jets.at(0)+jets.at(1)).M(), weight, 80, 0., 400.);
  FillHist(prefix+hprefix+"mTopl", (bjets.at(0)+muons.at(0)+met).M(), weight, 80, 0., 400.);
  FillHist(prefix+hprefix+"mW", (muons.at(0)+met).M(), weight, 60, 0., 300.);
  if(muons.at(0).Charge() > 0) FillHist(prefix+hprefix+"mWp", (muons.at(0)+met).M(), weight, 60, 0., 300.);
  else FillHist(prefix+hprefix+"mWm", (muons.at(0)+met).M(), weight, 60, 0., 300.);

  if(fabs(172.5 - (bjets.at(0)+muons.at(0)+met).M()) < 30) return;
  FillCutflow(prefix+"cutflow","Topl30Veto",weight);
  hprefix = "Topl30Veto_";

  FillHist(prefix+hprefix+"MET", met.Pt(), weight, 60, 0., 300.);
  FillHist(prefix+hprefix+"mupT", muons.at(0).Pt(), weight, 60, 0., 300.);
  FillHist(prefix+hprefix+"bpT", bjets.at(0).Pt(), weight, 60, 0., 300.);
  FillHist(prefix+hprefix+"j0pT", jets.at(0).Pt(), weight, 60, 0., 300.);
  FillHist(prefix+hprefix+"j1pT", jets.at(1).Pt(), weight, 60, 0., 300.);
  FillHist(prefix+hprefix+"mWj", (jets.at(0)+jets.at(1)).M(), weight, 60, 0., 300.);
  FillHist(prefix+hprefix+"mTop", (bjets.at(0)+jets.at(0)+jets.at(1)).M(), weight, 80, 0., 400.);
  FillHist(prefix+hprefix+"mTopl", (bjets.at(0)+muons.at(0)+met).M(), weight, 80, 0., 400.);
  FillHist(prefix+hprefix+"mW", (muons.at(0)+met).M(), weight, 60, 0., 300.);
  if(muons.at(0).Charge() > 0) FillHist(prefix+hprefix+"mWp", (muons.at(0)+met).M(), weight, 60, 0., 300.);
  else FillHist(prefix+hprefix+"mWm", (muons.at(0)+met).M(), weight, 60, 0., 300.);

  if(fabs(80.4 - (muons.at(0)+met).M()) < 15) return;
  FillCutflow(prefix+"cutflow","Wlep15",weight);
  hprefix = "Wlep15_";

  FillHist(prefix+hprefix+"MET", met.Pt(), weight, 60, 0., 300.);
  FillHist(prefix+hprefix+"mupT", muons.at(0).Pt(), weight, 60, 0., 300.);
  FillHist(prefix+hprefix+"bpT", bjets.at(0).Pt(), weight, 60, 0., 300.);
  FillHist(prefix+hprefix+"j0pT", jets.at(0).Pt(), weight, 60, 0., 300.);
  FillHist(prefix+hprefix+"j1pT", jets.at(1).Pt(), weight, 60, 0., 300.);
  FillHist(prefix+hprefix+"mWj", (jets.at(0)+jets.at(1)).M(), weight, 60, 0., 300.);
  FillHist(prefix+hprefix+"mTop", (bjets.at(0)+jets.at(0)+jets.at(1)).M(), weight, 80, 0., 400.);
  FillHist(prefix+hprefix+"mTopl", (bjets.at(0)+muons.at(0)+met).M(), weight, 80, 0., 400.);
  FillHist(prefix+hprefix+"mW", (muons.at(0)+met).M(), weight, 60, 0., 300.);
  if(muons.at(0).Charge() > 0) FillHist(prefix+hprefix+"mWp", (muons.at(0)+met).M(), weight, 60, 0., 300.);
  else FillHist(prefix+hprefix+"mWm", (muons.at(0)+met).M(), weight, 60, 0., 300.);
}

tWAnalyzer::tWAnalyzer(){

}

tWAnalyzer::~tWAnalyzer(){

}

double tWAnalyzer::Lepton_SF(TString histkey,const Lepton* lep,int sys){
  if(IsDATA) return 1.;
  if(histkey=="") return 1.;
  if(histkey=="Default") return 1.;
  double this_pt,this_eta;
  TH2* this_hist=NULL;
  if(histkey.Contains(TRegexp("_Q$"))){
    if(lep->Charge()>0) histkey+="Plus";
    else histkey+="Minus";
  }else if(histkey.Contains("_Q_")){
    if(lep->Charge()>0) histkey.ReplaceAll("_Q_","_QPlus_");
    else histkey.ReplaceAll("_Q_","_QMinus_");;
  }
  if(lep->LeptonFlavour()==Lepton::MUON){
    this_pt=((Muon*)lep)->MiniAODPt();
    this_eta=lep->Eta();
    this_hist=mcCorr->map_hist_Muon[histkey];
    if(!this_hist && DataYear==2016 && !histkey.Contains("_BCDEF$") && !histkey.Contains("_GH$")){
      double lumi_periodB = 5750.490644035;
      double lumi_periodC = 2572.903488748;
      double lumi_periodD = 4242.291556970;
      double lumi_periodE = 4025.228136967;
      double lumi_periodF = 3104.509131800;
      double lumi_periodG = 7575.824256098;
      double lumi_periodH = 8650.628380028;
      double total_lumi = (lumi_periodB+lumi_periodC+lumi_periodD+lumi_periodE+lumi_periodF+lumi_periodG+lumi_periodH);

      double WeightBtoF = (lumi_periodB+lumi_periodC+lumi_periodD+lumi_periodE+lumi_periodF)/total_lumi;
      double WeightGtoH = (lumi_periodG+lumi_periodH)/total_lumi;

      if(histkey.Contains("_SF_")){
        TString histkey_data=histkey;
        histkey_data.ReplaceAll("_SF_","_Eff_DATA_");
        TString histkey_mc=histkey;
        histkey_mc.ReplaceAll("_SF_","_Eff_MC_");
        double data_eff=WeightBtoF*Lepton_SF(histkey_data+"_BCDEF",lep,sys)+WeightGtoH*Lepton_SF(histkey_data+"_GH",lep,sys);
        double mc_eff=WeightBtoF*Lepton_SF(histkey_mc+"_BCDEF",lep,-sys)+WeightGtoH*Lepton_SF(histkey_mc+"_GH",lep,-sys);
        if(mc_eff==0) return 1;
        else return data_eff/mc_eff;
      }else if(histkey.Contains("_Eff_")){
        return WeightBtoF*Lepton_SF(histkey+"_BCDEF",lep,sys)+WeightGtoH*Lepton_SF(histkey+"_GH",lep,sys);
      }
    }
  }else if(lep->LeptonFlavour()==Lepton::ELECTRON){
    this_pt=((Electron*)lep)->UncorrPt();
    this_eta=((Electron*)lep)->scEta();
    this_hist=mcCorr->map_hist_Electron[histkey];
  }else{
    cout <<"[tWAnalyzer::Lepton_SF] It is not lepton"<<endl;
    exit(EXIT_FAILURE);
  }
  if(!this_hist){
    cout <<"[tWAnalyzer::Lepton_SF] no hist "<<histkey<<endl;
    exit(EXIT_FAILURE);
  }
  double this_x,this_y;
  if(this_hist->GetXaxis()->GetXmax()>this_hist->GetYaxis()->GetXmax()){
    if(histkey.Contains("_Eff_") && this_pt<this_hist->GetXaxis()->GetXmin()) return 0;
    this_x=this_pt;
    this_y=this_eta;
  }else{
    if(histkey.Contains("_Eff_") && this_pt<this_hist->GetYaxis()->GetXmin()) return 0;
    this_x=this_eta;
    this_y=this_pt;
  }
  return GetBinContentUser(this_hist,this_x,this_y,sys);
}

double tWAnalyzer::LeptonTrigger_SF(TString triggerSF_key,const vector<Lepton*>& leps,int sys){
  if(IsDATA) return 1;
  if(triggerSF_key=="") return 1;
  if(triggerSF_key=="Default") return 1;

  double data_eff=1.,mc_eff=1.;
  for(const auto& lep:leps){
    data_eff*=1-Lepton_SF("Trigger_Eff_DATA_"+triggerSF_key,lep,sys);
    mc_eff*=1-Lepton_SF("Trigger_Eff_MC_"+triggerSF_key,lep,-sys);
  }
  data_eff=1-data_eff;
  mc_eff=1-mc_eff;
  if(mc_eff==0) return 1.;
  else return data_eff/mc_eff;
}
