#include "ExampleRun_kinFitter.h"

ExampleRun_kinFitter::ExampleRun_kinFitter(){

}

void ExampleRun_kinFitter::initializeAnalyzer(){

  MuonIDs = { "DeepCSV", "DeepJet", "DeepCSV_Tight", "DeepJet_Tight" };
  MuonIDSFKeys = { "IDISO_SF_MediumID_trkIsoLoose_Q", "IDISO_SF_MediumID_trkIsoLoose_Q", "IDISO_SF_MediumID_trkIsoLoose_Q", "IDISO_SF_MediumID_trkIsoLoose_Q" };
  MuonTrigSFKeys = { "IsoMu27_MediumID_trkIsoLoose_Q", "IsoMu27_MediumID_trkIsoLoose_Q", "IsoMu27_MediumID_trkIsoLoose_Q", "IsoMu27_MediumID_trkIsoLoose_Q" };

  if(DataYear==2016){
    IsoMuTriggerName = "HLT_IsoMu24_v";
    TriggerSafePtCut = 26.;
  }
  else if(DataYear==2017){
    IsoMuTriggerName = "HLT_IsoMu27_v";
    TriggerSafePtCut = 29.;
  }

  cout << "[ExampleRun_kinFitter::initializeAnalyzer] IsoMuTriggerName = " << IsoMuTriggerName << endl;
  cout << "[ExampleRun_kinFitter::initializeAnalyzer TriggerSafePtCut = " << TriggerSafePtCut << endl;

  //==== B-Tagging
  SetupMCJetTagEff();
  //==== add taggers and WP that you want to use in analysis
  std::vector<JetTagging::Parameters> jtps;
  //==== If you want to use 1a or 2a method,
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb) );
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::comb) );
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Tight, JetTagging::incl, JetTagging::comb) );
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Tight, JetTagging::incl, JetTagging::comb) );
  //==== set
  mcCorr->SetJetTaggingParameters(jtps);

  fitter = new TKinFitterDriver(DataYear);
}

ExampleRun_kinFitter::~ExampleRun_kinFitter(){

  //==== Destructor of this Analyzer
  delete fitter;
}

void ExampleRun_kinFitter::executeEvent(){

  AllMuons = GetAllMuons();
  weight_Prefire = GetPrefireWeight(0);
  AnalyzerParameter param;

  for(unsigned int it_MuonID=0; it_MuonID<MuonIDs.size(); it_MuonID++){

    TString MuonID = "POGMediumWithLooseTrkIso";
    TString MuonIDSFKey = MuonIDSFKeys.at(it_MuonID);
    TString MuonTrigSFKey = MuonTrigSFKeys.at(it_MuonID);

    //==== clear parameter set
    param.Clear();

    param.syst_ = AnalyzerParameter::Central;
    param.Name = MuonIDs.at(it_MuonID);
    //==== You can define lepton ID string here
    param.Muon_Tight_ID = MuonID;
    param.Muon_ID_SF_Key = MuonIDSFKey;
    param.Muon_Trigger_SF_Key = MuonTrigSFKey;
    param.Jet_ID = "tightLepVeto";

    //==== Now, all parameters are set. Run executeEventFromParameter() with this parameter set
    executeEventFromParameter(param);
  }
}

void ExampleRun_kinFitter::executeEventFromParameter(AnalyzerParameter param){

  //=============
  //==== No Cut
  //=============
  TString prefix = param.Name+"/";
  FillHist(prefix+"NoCut_", 0., 1., 1, 0., 1.);

  //========================
  //==== MET Filter
  //========================

  if(!PassMETFilter()) return;

  Event ev = GetEvent();
  Particle METv = ev.GetMETVector();
  TLorentzVector met = GetEvent().GetMETVector();
  if(met.Pt() < 20.) return;

  if(! (ev.PassTrigger(IsoMuTriggerName) )) return;

  if(param.syst_ == AnalyzerParameter::Central){
  }
  else{
    cout << "[ExampleRun_kinFitter::executeEventFromParameter] Wrong syst" << endl;
    exit(EXIT_FAILURE);
  }

  //==================================================
  //==== Then, apply ID selections using this_AllXXX
  //==================================================

  JetTagging::Tagger btagger = JetTagging::DeepCSV;
  JetTagging::WP btagWP = JetTagging::Medium;
  if(param.Name.Contains("DeepJet")) btagger = JetTagging::DeepJet;
  if(param.Name.Contains("Tight")) btagWP = JetTagging::Tight;

  vector<Muon> this_AllMuons = AllMuons;
  vector<Muon> muons = SMPGetMuons(param.Muon_Tight_ID, 20., 2.4);
  vector<Jet> basic_jets = GetJets(param.Jet_ID, 30., 2.4);
  //vector<Jet> jets = GetJets(param.Jet_ID, 30., 2.4);
  std::sort(muons.begin(), muons.end(), PtComparing);
  std::sort(basic_jets.begin(), basic_jets.end(), PtComparing);
  //std::sort(jets.begin(), jets.end(), PtComparing);

  vector<Jet> jets;
  jets.clear();
  unsigned int nmax_bjets = 15, nmax_jets = 15, nbjets=0, njets=0;
  for(unsigned int i=0; i<basic_jets.size(); i++){
    double this_discr = basic_jets.at(i).GetTaggerResult(btagger);
    if( this_discr > mcCorr->GetJetTaggingCutValue(btagger, btagWP)) nbjets++;
    else njets++;
    if( this_discr > mcCorr->GetJetTaggingCutValue(btagger, btagWP) && nbjets < nmax_bjets+1) jets.push_back(basic_jets.at(i));
    else if(this_discr < mcCorr->GetJetTaggingCutValue(btagger, btagWP) && njets < nmax_jets+1) jets.push_back(basic_jets.at(i));
  }
  // this jets set is the input of Kinematic fitter
  std::sort(jets.begin(), jets.end(), PtComparing);

  if(muons.size() != 1) return;
  if(muons.at(0).Pt() <= TriggerSafePtCut ) return;

  //===================
  //==== Event weight
  //===================

  double weight = 1.;
  //==== If MC
  if(!IsDATA){

    weight *= weight_norm_1invpb*ev.GetTriggerLumi("Full");
    weight *= ev.MCweight();
    weight *= weight_Prefire;

    //==== Example of applying Muon scale factors
    for(unsigned int i=0; i<muons.size(); i++){
      Lepton *l = (Lepton *)(&muons.at(i));
      double this_idsf  = Lepton_SF(param.Muon_ID_SF_Key, l, 0);
      double this_isosf = 1.;
      weight *= this_idsf*this_isosf;
    }
    double this_trigsf = LeptonTrigger_SF(param.Muon_Trigger_SF_Key, MakeLeptonPointerVector(muons), 0);
    weight *= this_trigsf;

    JetTagging::Parameters jtp = JetTagging::Parameters(btagger, btagWP, JetTagging::incl, JetTagging::comb);
    btagweight = SMPAnalyzerCore::GetBTaggingReweight_1a(basic_jets, jtp);
    weight *= btagweight;
    FillHist(prefix+"btagweight", btagweight, 1, 100, 0, 10);
  }

  FillHist(prefix+"basicjet_Num", basic_jets.size(), weight, 15, 0, 15);
  FillHist(prefix+"basicbjet_Num", nbjets, weight, 10, 0, 10);

  if(nbjets!=2) return;
  if(jets.size()<4) return;

  std::vector<bool> btag_vector{};
  for(unsigned int ij = 0 ; ij < jets.size(); ij++){
    double this_discr = jets.at(ij).GetTaggerResult(btagger);
    if( this_discr > mcCorr->GetJetTaggingCutValue(btagger, btagWP)) btag_vector.push_back(true);
    else btag_vector.push_back(false);
  }

  FillHist(prefix+"jet_Num", jets.size(), weight, 10, 0, 10);
  FillHist(prefix+"bjet_Num", std::accumulate(btag_vector.begin(), btag_vector.end(),0), weight, 10, 0, 10);

  //=======================
  //==== Kinematic Fitter
  //=======================
  std::vector<TLorentzVector> jet_vector{};
  TLorentzVector lepton{};
  for(auto& jet : jets){
    jet_vector.emplace_back(jet.Px(),jet.Py(),jet.Pz(),jet.E());
  }
  if(jet_vector.size() != btag_vector.size()){
    cout << " ExampleRun_kinFitter, jet_vector.size() != btag_vector.size()" << endl;
    exit(1);
  }
  lepton = (TLorentzVector)(muons.at(0));
  fitter->SetAllObjects(jet_vector,
        	        btag_vector,
        		lepton,
        		met
        	       );
  fitter->FindBestChi2Fit();
  auto fitter_results = fitter->GetResults();

  FillHist(prefix+"Fit_results_Num", fitter_results->size(), weight, 20, 0, 20);
  FillHist(prefix+"Fit_fail_rate", fitter_results->size()<1?0:1, weight, 4, 0, 4);
  if(fitter_results->size() < 1) return;

  // Before Fit variables
  int lept_b_idx = fitter_results->at(0).leptonic_top_b_jet_idx;
  int hadt_b_idx = fitter_results->at(0).hadronic_top_b_jet_idx;
  int W_up_jet_idx = fitter_results->at(0).w_ch_up_type_jet_idx;
  int W_down_jet_idx = fitter_results->at(0).w_ch_down_type_jet_idx;
  Jet W_up_jet = jets.at(W_up_jet_idx);
  Jet W_down_jet = jets.at(W_down_jet_idx);
  Jet lept_b = jets.at(lept_b_idx);
  Jet hadt_b = jets.at(hadt_b_idx);
  double lept_b_charge = bjetCharge(lept_b,1);
  double hadt_b_charge = bjetCharge(hadt_b,1);

  double hadronic_W_M = (W_up_jet+W_down_jet).M();
  double leptonic_W_M = fitter_results->at(0).leptonic_W_M;
  double hadronic_top_M = fitter_results->at(0).hadronic_top_M;
  double leptonic_top_M = fitter_results->at(0).leptonic_top_M;

  // After Fit variables
  TLorentzVector fitted_lept_b = fitter_results->at(0).fitted_lept_bjet;
  TLorentzVector fitted_hadt_b = fitter_results->at(0).fitted_hadt_bjet;
  TLorentzVector fitted_W_j1 = fitter_results->at(0).fitted_jet1;
  TLorentzVector fitted_W_j2 = fitter_results->at(0).fitted_jet2;
  TLorentzVector fitted_lep = fitter_results->at(0).fitted_lep;
  TLorentzVector fitted_neu = fitter_results->at(0).fitted_neu;

  double fitted_hadronic_W_M = (fitted_W_j1+fitted_W_j2).M();
  double fitted_leptonic_W_M = (fitted_lep+fitted_neu).M();
  double fitted_hadronic_top_M = (fitted_hadt_b+fitted_W_j1+fitted_W_j2).M();
  double fitted_leptonic_top_M = (fitted_lept_b+fitted_lep+fitted_neu).M();

  //============================================================
  //==== Gen - Jet Matching (To check the performance of Fitter)
  //============================================================
  if(!IsDATA && MCSample.Contains("TTLJ")){
    vector<Gen> gens=GetGens();
    Gen gen_parton0,gen_parton1, gen_b0,gen_b1, gen_l0,gen_l1, gen_j0,gen_j1;
    GetTTLJGenParticles(gens, gen_parton0,gen_parton1, gen_b0,gen_b1, gen_l0,gen_l1, gen_j0,gen_j1,3);

    //PrintGens(gens);
    if((gen_j0+gen_j1).M() > 200 || (gen_j0+gen_j1).M() ==0){
      cout<<"Whad wrong"<<endl;
      cout<<"(j0 index, pid, pt, E, px, py, pz) = "<<gen_j0.Index()<<", "<<gen_j0.PID()<<", "<<gen_j0.Pt()<<", "<<gen_j0.E()<<", "<<gen_j0.Px()<<", "<<gen_j0.Py()<<", "<<gen_j0.Pz()<<endl;
      cout<<"(j1 index, pid, pt, E, px, py, pz) = "<<gen_j1.Index()<<", "<<gen_j1.PID()<<", "<<gen_j1.Pt()<<", "<<gen_j1.E()<<", "<<gen_j1.Px()<<", "<<gen_j1.Py()<<", "<<gen_j1.Pz()<<endl;
      PrintGens(gens);
    }
    if((gen_l0+gen_l1).M() > 200 || (gen_l0+gen_l1).M() ==0){
      cout<<"Wlep wrong"<<endl;
      cout<<"(l0 index, pid, pt, E, px, py, pz) = "<<gen_l0.Index()<<", "<<gen_l0.PID()<<", "<<gen_l0.Pt()<<", "<<gen_l0.E()<<", "<<gen_l0.Px()<<", "<<gen_l0.Py()<<", "<<gen_l0.Pz()<<endl;
      cout<<"(l1 index, pid, pt, E, px, py, pz) = "<<gen_l1.Index()<<", "<<gen_l1.PID()<<", "<<gen_l1.Pt()<<", "<<gen_l1.E()<<", "<<gen_l1.Px()<<", "<<gen_l1.Py()<<", "<<gen_l1.Pz()<<endl;
      PrintGens(gens);}

    FillHist(prefix+"gen_W_had_Mass", (gen_j0+gen_j1).M(), weight, 50, 0., 200.);
    FillHist(prefix+"gen_W_lep_Mass", (gen_l0+gen_l1).M(), weight, 50, 0., 200.);
    FillHist(prefix+"gen_top_had1_Mass", (gen_b0+gen_j0+gen_j1).M(), weight, 75, 0., 300.);
    FillHist(prefix+"gen_top_had2_Mass", (gen_b1+gen_j0+gen_j1).M(), weight, 75, 0., 300.);
    FillHist(prefix+"gen_top_lep1_Mass", (gen_b0+gen_l0+gen_l1).M(), weight, 75, 0., 300.);
    FillHist(prefix+"gen_top_lep2_Mass", (gen_b1+gen_l0+gen_l1).M(), weight, 75, 0., 300.);
    if(gen_l0.PID()<0){ //gen_l0 = mu+
      FillHist(prefix+"gen_top_lep_Mass", (gen_b0+gen_l0+gen_l1).M(), weight, 75, 0., 300.);
      FillHist(prefix+"gen_top_had_Mass", (gen_b1+gen_j0+gen_j1).M(), weight, 75, 0., 300.);
    }else{
      FillHist(prefix+"gen_top_lep_Mass", (gen_b1+gen_l0+gen_l1).M(), weight, 75, 0., 300.);
      FillHist(prefix+"gen_top_had_Mass", (gen_b0+gen_j0+gen_j1).M(), weight, 75, 0., 300.);
    }

    double match_dR = 0.4;
    bool Gen_Fit_match = false;
    bool Gen_Fit_match_onlyb = false;
    //Easy charge means 0:negative, 1:positive

    //When lept_b is b, hadt_b is bbar
    if(gen_b0.DeltaR(lept_b) < match_dR && gen_b1.DeltaR(hadt_b) < match_dR){
      Gen_Fit_match_onlyb = true;
      //lept_b_charge = bjetCharge(lept_b,1,prefix+"Matched_b_",weight);
      //hadt_b_charge = bjetCharge(hadt_b,1,prefix+"Matched_bbar_",weight);
      if(muons.at(0).Charge() > 0) prefix += "Correct_";
      else prefix +="Wrong_";

      lept_b_charge = bjetCharge(lept_b,1,prefix+"Matched_b_",weight);
      hadt_b_charge = bjetCharge(hadt_b,1,prefix+"Matched_bbar_",weight);
      if(abs(lept_b_charge) > 0.2) FillHist(prefix+"lept_b_Charge_Matched_b", lept_b_charge, weight, 600, -6, 6);
      if(abs(lept_b_charge) > 0.2) FillHist(prefix+"Easy_lept_b_Charge_Matched_b", lept_b_charge<0?-1:1, weight, 8, -4., 4.);
      if(abs(hadt_b_charge) > 0.2) FillHist(prefix+"hadt_b_Charge_Matched_bbar", hadt_b_charge, weight, 600, -6, 6);
      if(abs(hadt_b_charge) > 0.2) FillHist(prefix+"Easy_hadt_b_Charge_Matched_bbar", hadt_b_charge<0?-1:1, weight, 8, -4., 4.);
      FillHist(prefix+"muon_Charge_lept_b_Matched_b", muons.at(0).Charge(), weight, 8, -4, 4);
      FillHist(prefix+"Easy_muon_Charge_lept_b_Matched_b", muons.at(0).Charge()<0?0:1, weight, 8, -4, 4);
      if(lept_b_charge < -0.2) FillHist(prefix+"muon_Charge_lept_b_negaive_Matched_b", muons.at(0).Charge(), weight, 8, -4, 4);
      if(lept_b_charge < -0.2) FillHist(prefix+"Easy_muon_Charge_lept_b_negaive_Matched_b", muons.at(0).Charge()<0?0:1, weight, 8, -4, 4);
      if(hadt_b_charge > 0.2) FillHist(prefix+"muon_Charge_hadt_b_positive_Matched_bbar", muons.at(0).Charge(), weight, 8, -4, 4);
      if(hadt_b_charge > 0.2) FillHist(prefix+"Easy_muon_Charge_hadt_b_positive_Matched_bbar", muons.at(0).Charge()<0?0:1, weight, 8, -4, 4);
      if(muons.at(0).Charge() > 0) FillHist(prefix+"lept_b_Charge_Muplus_Matched_b", lept_b_charge, weight, 600, -6, 6);
      if(abs(lept_b_charge) > 0.2 && muons.at(0).Charge() > 0) FillHist(prefix+"Easy_lept_b_Charge_Muplus_Matched_b", lept_b_charge<0?-1:1, weight, 8, -4, 4);
      if((gen_j0.DeltaR(W_up_jet) < match_dR && gen_j1.DeltaR(W_down_jet) < match_dR) || (gen_j0.DeltaR(W_down_jet) < match_dR && gen_j1.DeltaR(W_up_jet) < match_dR)) Gen_Fit_match = true;
      //When KinFit != gen
      if(muons.at(0).Charge() < 0) FillHist(prefix+"lept_b_Charge_Muminus_Matched_b", lept_b_charge, weight, 600, -6, 6);
      if(abs(lept_b_charge) > 0.2 && muons.at(0).Charge() < 0) FillHist(prefix+"Easy_lept_b_Charge_Muminus_Matched_b", lept_b_charge<0?-1:1, weight, 8, -4, 4);
    }
    //When lept_b is bbar, hadt_b is b
    else if(gen_b0.DeltaR(hadt_b) < match_dR && gen_b1.DeltaR(lept_b) < match_dR){
      Gen_Fit_match_onlyb = true;
      //lept_b_charge = bjetCharge(lept_b,1,prefix+"Matched_bbar_",weight);
      //hadt_b_charge = bjetCharge(hadt_b,1,prefix+"Matched_b_",weight);
      if(muons.at(0).Charge() < 0) prefix += "Correct_";
      else prefix +="Wrong_";

      lept_b_charge = bjetCharge(lept_b,1,prefix+"Matched_bbar_",weight);
      hadt_b_charge = bjetCharge(hadt_b,1,prefix+"Matched_b_",weight);
      if(abs(hadt_b_charge) > 0.2) FillHist(prefix+"hadt_b_Charge_Matched_b", hadt_b_charge, weight, 600, -6, 6);
      if(abs(hadt_b_charge) > 0.2) FillHist(prefix+"Easy_hadt_b_Charge_Matched_b", hadt_b_charge<0?-1:1, weight, 8, -4., 4.);
      if(abs(lept_b_charge) > 0.2) FillHist(prefix+"lept_b_Charge_Matched_bbar", lept_b_charge, weight, 600, -6, 6);
      if(abs(lept_b_charge) > 0.2) FillHist(prefix+"Easy_lept_b_Charge_Matched_bbar", lept_b_charge<0?-1:1, weight, 8, -4., 4.);
      FillHist(prefix+"muon_Charge_lept_b_Matched_bbar", muons.at(0).Charge(), weight, 8, -4, 4);
      FillHist(prefix+"Easy_muon_Charge_lept_b_Matched_bbar", muons.at(0).Charge()<0?0:1, weight, 8, -4, 4);
      if(lept_b_charge > 0.2) FillHist(prefix+"muon_Charge_lept_b_positive_Matched_bbar", muons.at(0).Charge(), weight, 8, -4, 4);
      if(lept_b_charge > 0.2) FillHist(prefix+"Easy_muon_Charge_lept_b_positive_Matched_bbar", muons.at(0).Charge()<0?0:1, weight, 8, -4, 4);
      if(hadt_b_charge < -0.2) FillHist(prefix+"muon_Charge_hadt_b_negative_Matched_b", muons.at(0).Charge(), weight, 8, -4, 4);
      if(hadt_b_charge < -0.2) FillHist(prefix+"Easy_muon_Charge_hadt_b_negative_Matched_b", muons.at(0).Charge()<0?0:1, weight, 8, -4, 4);
      if(muons.at(0).Charge() < 0) FillHist(prefix+"lept_b_Charge_Muminus_Matched_bbar", lept_b_charge, weight, 600, -6, 6);
      if(abs(lept_b_charge) > 0.2 && muons.at(0).Charge() < 0) FillHist(prefix+"Easy_lept_b_Charge_Muminus_Matched_bbar", lept_b_charge<0?-1:1, weight, 8, -4, 4);
      if((gen_j0.DeltaR(W_up_jet) < match_dR && gen_j1.DeltaR(W_down_jet) < match_dR) || (gen_j0.DeltaR(W_down_jet) < match_dR && gen_j1.DeltaR(W_up_jet) < match_dR)) Gen_Fit_match = true;
      //When KinFit != gen
      if(muons.at(0).Charge() > 0) FillHist(prefix+"lept_b_Charge_Muplus_Matched_bbar", lept_b_charge, weight, 600, -6, 6);
      if(abs(lept_b_charge) > 0.2 && muons.at(0).Charge() > 0) FillHist(prefix+"Easy_lept_b_Charge_Muplus_Matched_bbar", lept_b_charge<0?-1:1, weight, 8, -4, 4);
    }
    else prefix +="UnMatched_";

    FillHist(prefix+"Gen_b0_mindR", min(gen_b0.DeltaR(lept_b),gen_b0.DeltaR(hadt_b)), weight, 50, 0., 5.);
    FillHist(prefix+"Gen_b1_mindR", min(gen_b1.DeltaR(lept_b),gen_b1.DeltaR(hadt_b)), weight, 50, 0., 5.);
    FillHist(prefix+"Gen_j0_mindR", min(gen_j0.DeltaR(W_up_jet),gen_j0.DeltaR(W_down_jet)), weight, 50, 0., 5.);
    FillHist(prefix+"Gen_j1_mindR", min(gen_j1.DeltaR(W_up_jet),gen_j1.DeltaR(W_down_jet)), weight, 50, 0., 5.);
    FillHist(prefix+"Gen_l0_mindR", min(gen_l0.DeltaR(fitted_lep),gen_l0.DeltaR(fitted_neu)), weight, 50, 0., 5.);
    FillHist(prefix+"Gen_l1_mindR", min(gen_l1.DeltaR(fitted_neu),gen_l1.DeltaR(fitted_lep)), weight, 50, 0., 5.);
    FillHist(prefix+"Gen_b0b1_dR", gen_b0.DeltaR(gen_b1), weight, 50, 0., 5.);
    FillHist(prefix+"Gen_b0j0_dR", gen_b0.DeltaR(gen_j0), weight, 50, 0., 5.);
    FillHist(prefix+"Gen_b0j1_dR", gen_b0.DeltaR(gen_j1), weight, 50, 0., 5.);
    FillHist(prefix+"Gen_b1j0_dR", gen_b1.DeltaR(gen_j0), weight, 50, 0., 5.);
    FillHist(prefix+"Gen_b1j1_dR", gen_b1.DeltaR(gen_j1), weight, 50, 0., 5.);
    FillHist(prefix+"Gen_j0j1_dR", gen_j0.DeltaR(gen_j1), weight, 50, 0., 5.);
    FillHist(prefix+"Gen_b0_dR", gen_b0.DeltaR(lept_b), weight, 50, 0., 5.);
    FillHist(prefix+"Gen_b1_dR", gen_b1.DeltaR(hadt_b), weight, 50, 0., 5.);
    FillHist(prefix+"Gen_j0_dR", gen_j0.DeltaR(W_up_jet), weight, 50, 0., 5.);
    FillHist(prefix+"Gen_j1_dR", gen_j1.DeltaR(W_down_jet), weight, 50, 0., 5.);
    FillHist(prefix+"Gen_l0_dR", gen_l0.DeltaR(fitted_lep), weight, 50, 0., 5.);
    FillHist(prefix+"Gen_l1_dR", gen_l1.DeltaR(fitted_neu), weight, 50, 0., 5.);
    FillHist(prefix+"Gen_b0_pT", gen_b0.Pt(), weight, 100, 0, 300);
    FillHist(prefix+"Gen_b1_pT", gen_b1.Pt(), weight, 100, 0, 300);
    FillHist(prefix+"Gen_j0_pT", gen_j0.Pt(), weight, 100, 0, 300);
    FillHist(prefix+"Gen_j1_pT", gen_j1.Pt(), weight, 100, 0, 300);
    FillHist(prefix+"Gen_b0_pTFitdiff", gen_b0.Pt()-fitted_lept_b.Pt(), weight, 100, -50, 50);
    FillHist(prefix+"Gen_b1_pTFitdiff", gen_b1.Pt()-fitted_hadt_b.Pt(), weight, 100, -50, 50);
    FillHist(prefix+"Gen_j0_pTFitdiff", gen_j0.Pt()-fitted_W_j1.Pt(), weight, 100, -50, 50);
    FillHist(prefix+"Gen_j1_pTFitdiff", gen_j1.Pt()-fitted_W_j2.Pt(), weight, 100, -50, 50);
    FillHist(prefix+"Gen_l0_pTFitdiff", gen_l0.Pt()-fitted_lep.Pt(), weight, 100, -50, 50);
    FillHist(prefix+"Gen_l1_pTFitdiff", gen_l1.Pt()-fitted_neu.Pt(), weight, 100, -50, 50);
    FillHist(prefix+"Gen_b0_pTdiff", gen_b0.Pt()-lept_b.Pt(), weight, 100, -50, 50);
    FillHist(prefix+"Gen_b1_pTdiff", gen_b1.Pt()-hadt_b.Pt(), weight, 100, -50, 50);
    FillHist(prefix+"Gen_j0_pTdiff", gen_j0.Pt()-W_up_jet.Pt(), weight, 100, -50, 50);
    FillHist(prefix+"Gen_j1_pTdiff", gen_j1.Pt()-W_down_jet.Pt(), weight, 100, -50, 50);
    FillHist(prefix+"Gen_l0_pTdiff", gen_l0.Pt()-muons.at(0).Pt(), weight, 100, -50, 50);
    FillHist(prefix+"Gen_l1_pTdiff", gen_l1.Pt()-METv.Pt(), weight, 100, -50, 50);
    FillHist(prefix+"Gen_Fit_Match_onlyb", Gen_Fit_match_onlyb, weight, 2, 0, 2);
    FillHist(prefix+"Gen_Fit_Match", Gen_Fit_match, weight, 2, 0, 2);
  }

  //==========================
  //==== Now reco fill histograms
  //==========================

  FillHist(prefix+"W_had_Mass", hadronic_W_M, weight, 40, 0., 200.);
  FillHist(prefix+"W_had_FitMass", fitted_hadronic_W_M, weight, 80, 60., 100.);
  FillHist(prefix+"W_had_FitMassdiff", fitted_hadronic_W_M-hadronic_W_M, weight, 100, -50., 50.);
  FillHist(prefix+"W_lep_Mass", leptonic_W_M, weight, 40, 0., 200.);
  FillHist(prefix+"W_lep_FitMass", fitted_leptonic_W_M, weight, 150, 75., 150.);
  FillHist(prefix+"W_lep_FitMassdiff", fitted_leptonic_W_M-leptonic_W_M, weight, 100, -50., 50.);
  FillHist(prefix+"Top_had_Mass", hadronic_top_M, weight, 40, 100., 300.);
  FillHist(prefix+"Top_had_FitMass", fitted_hadronic_top_M, weight, 160, 160., 240.);
  FillHist(prefix+"Top_had_FitMassdiff", fitted_hadronic_top_M-hadronic_top_M, weight, 100, -50., 50.);
  FillHist(prefix+"Top_lep_Mass", leptonic_top_M, weight, 40, 100., 300.);
  FillHist(prefix+"Top_lep_FitMass", fitted_leptonic_top_M, weight, 80, 160., 200.);
  FillHist(prefix+"Top_lep_FitMassdiff", fitted_leptonic_top_M-leptonic_top_M, weight, 100, -50., 50.);

  FillHist(prefix+"lept_b_FitpT", fitted_lept_b.Pt(), weight, 100, 0., 300.);
  FillHist(prefix+"hadt_b_FitpT", fitted_hadt_b.Pt(), weight, 100, 0., 300.);
  FillHist(prefix+"W_had1_FitpT", fitted_W_j1.Pt(), weight, 100, 0., 300.);
  FillHist(prefix+"W_had2_FitpT", fitted_W_j2.Pt(), weight, 100, 0., 300.);
  FillHist(prefix+"lept_b_FitpTdiff", fitted_lept_b.Pt()-lept_b.Pt(), weight, 100, -50., 50.);
  FillHist(prefix+"hadt_b_FitpTdiff", fitted_hadt_b.Pt()-hadt_b.Pt(), weight, 100, -50., 50.);
  FillHist(prefix+"W_had1_FitpTdiff", fitted_W_j1.Pt()-W_up_jet.Pt(), weight, 100, -50., 50.);
  FillHist(prefix+"W_had2_FitpTdiff", fitted_W_j2.Pt()-W_down_jet.Pt(), weight, 100, -50., 50.);
  FillHist(prefix+"lept_b_hadt_b_FitdR", fitted_lept_b.DeltaR(fitted_hadt_b), weight, 100, 0., 5.);
  FillHist(prefix+"lept_b_hadt_b_FitdPhi", fitted_lept_b.DeltaPhi(fitted_hadt_b), weight, 200, -5., 5.);
  FillHist(prefix+"W_had1_had2_FitdR", fitted_W_j1.DeltaR(fitted_W_j2), weight, 100, 0., 5.);
  FillHist(prefix+"W_had1_had2_FitdPhi", fitted_W_j1.DeltaPhi(fitted_W_j2), weight, 200, -5., 5.);
  FillHist(prefix+"lept_b_FitP", fitted_lept_b.P(), weight, 100, 0., 300.);
  FillHist(prefix+"hadt_b_FitP", fitted_hadt_b.P(), weight, 100, 0., 300.);
  FillHist(prefix+"W_had1_FitP", fitted_W_j1.P(), weight, 100, 0., 300.);
  FillHist(prefix+"W_had2_FitP", fitted_W_j2.P(), weight, 100, 0., 300.);

  FillHist(prefix+"W_up_jet_idx",  W_up_jet_idx, weight, 10, 0., 10.);
  FillHist(prefix+"W_down_jet_idx",  W_down_jet_idx, weight, 10, 0., 10.);
  FillHist(prefix+"lept_b_idx", lept_b_idx, weight, 10, 0., 10.);
  FillHist(prefix+"hadt_b_idx", hadt_b_idx, weight, 10, 0., 10.);
  if(muons.at(0).Charge() < 0){
    if(abs(lept_b_charge) > 0.2) FillHist(prefix+"lept_b_Charge_Muminus", lept_b_charge, weight, 600, -6, 6);
    if(abs(hadt_b_charge) > 0.2) FillHist(prefix+"hadt_b_Charge_Muminus", hadt_b_charge, weight, 600, -6, 6);
    if(abs(lept_b_charge) > 0.2) FillHist(prefix+"Easy_lept_b_Charge_Muminus", lept_b_charge<0?-1:1, weight, 8, -4., 4.);
    if(abs(hadt_b_charge) > 0.2) FillHist(prefix+"Easy_hadt_b_Charge_Muminus", hadt_b_charge<0?-1:1, weight, 8, -4., 4.);
  }else{
    if(abs(lept_b_charge) > 0.2) FillHist(prefix+"lept_b_Charge_Muplus", lept_b_charge, weight, 600, -6, 6);
    if(abs(hadt_b_charge) > 0.2) FillHist(prefix+"hadt_b_Charge_Muplus", hadt_b_charge, weight, 600, -6, 6);
    if(abs(lept_b_charge) > 0.2) FillHist(prefix+"Easy_lept_b_Charge_Muplus", lept_b_charge<0?-1:1, weight, 8, -4., 4.);
    if(abs(hadt_b_charge) > 0.2) FillHist(prefix+"Easy_hadt_b_Charge_Muplus", hadt_b_charge<0?-1:1, weight, 8, -4., 4.);
  }
}

void ExampleRun_kinFitter::GetTTLJGenParticles(const vector<Gen>& gens,Gen& parton0,Gen& parton1,Gen& b0,Gen& b1,Gen& l0,Gen& l1,Gen& j0,Gen& j1, int mode){
  //mode 0:bare 1:dressed01 2:dressed04 3:beforeFSR
  vector<const Gen*> leptons;
  vector<const Gen*> photons;
  vector<const Gen*> jets;

  int ngen=gens.size();
  for(int i=0;i<ngen;i++){
    if(!gens.at(i).isPrompt()) continue;
    int genpid=gens.at(i).PID();
    if(gens.at(i).isHardProcess()){
      if(abs(genpid)<7||genpid==21){
        if(parton0.IsEmpty()) parton0=gens[i];
        else if(parton1.IsEmpty()) parton1=gens[i];
      }
    }
    if(gens.at(i).Status()==1){
      if(gens.at(i).PID()==22) photons.push_back(&gens[i]); //photon
    }
    if(!gens.at(i).isHardProcess()) continue;
    if((abs(genpid)>=11 && abs(genpid)<=18) && abs(gens.at(gens.at(i).MotherIndex()).PID()) == 24) leptons.push_back(&gens[i]); //leptons from W
    // b0 : b from t, b1 : bbar from tbar
    if(genpid==5 && gens.at(gens.at(i).MotherIndex()).PID() == 6) b0=gens[i];
    else if(genpid==-5 && gens.at(gens.at(i).MotherIndex()).PID() == -6) b1=gens[i];
    else if((abs(genpid)<=5 || genpid==21) && abs(gens.at(gens.at(i).MotherIndex()).PID()) == 24) jets.push_back(&gens[i]); //jets from W
  }

  // l0 : charged lepton, l1 : neutrino
  int nlepton=leptons.size();
  for(int i=0;i<nlepton;i++){
    for(int j=i+1;j<nlepton;j++){
      if(abs(leptons[i]->PID()+leptons[j]->PID()) != 1) continue;
      if((*leptons[i]+*leptons[j]).M()>(l0+l1).M()){ //leptonic W
        if(abs(leptons[i]->PID()) < abs(leptons[j]->PID())){
          l0=*leptons[i]; //charged lepton
          l1=*leptons[j]; //neutrino
        }else{
          l0=*leptons[j];
          l1=*leptons[i];
        }
      }
    }
  }

  // j0 : leading jet, j1 : subleading jet
  int njet=jets.size();
  for(int i=0;i<njet;i++){
    for(int j=i+1;j<njet;j++){
      if(!(abs(jets[i]->PID()+jets[j]->PID()) == 1 || abs(jets[i]->PID()+jets[j]->PID()) == 3)) continue;
      if((*jets[i]+*jets[j]).M()>(j0+j1).M()){ //hadronic W
	if(jets[i]->Pt() > jets[j]->Pt()){
	  j0=*jets[i];
	  j1=*jets[j];
	}else{
	  j0=*jets[j];
	  j1=*jets[i];
	}
      }
    }
  }

  //cout<<"(b0 index, pid, pt, p) = "<<b0.Index()<<", "<<b0.PID()<<", "<<b0.Pt()<<", "<<b0.P()<<endl;
  //cout<<"(b1 index, pid, pt, p) = "<<b1.Index()<<", "<<b1.PID()<<", "<<b1.Pt()<<", "<<b1.P()<<endl;
  //cout<<"(l0 index, pid, pt, p) = "<<l0.Index()<<", "<<l0.PID()<<", "<<l0.Pt()<<", "<<l0.P()<<endl;
  //cout<<"(l1 index, pid, pt, p) = "<<l1.Index()<<", "<<l1.PID()<<", "<<l1.Pt()<<", "<<l1.P()<<endl;
  //cout<<"(j0 index, pid, pt, p) = "<<j0.Index()<<", "<<j0.PID()<<", "<<j0.Pt()<<", "<<j0.P()<<endl;
  //cout<<"(j1 index, pid, pt, p) = "<<j1.Index()<<", "<<j1.PID()<<", "<<j1.Pt()<<", "<<j1.P()<<endl;

  if(mode>=3){
    if(nlepton>=4){
      for(int i=0;i<nlepton;i++){
        if(leptons[i]->Index()==l0.Index()||leptons[i]->Index()==l1.Index()) continue;
        for(int j=i+1;j<nlepton;j++){
          if(leptons[j]->Index()==l0.Index()||leptons[j]->Index()==l1.Index()) continue;
          if(!(leptons[i]->PID()+leptons[j]->PID()==0)) continue;
          vector<int> history_i=TrackGenSelfHistory(*leptons[i],gens);
          vector<int> history_j=TrackGenSelfHistory(*leptons[j],gens);
          if(history_i.at(1)==history_j.at(1)) photons.push_back(&gens[history_i.at(1)]);
        }
      }
    }
    for(const auto& photon:photons){
      vector<int> history=TrackGenSelfHistory(*photon,gens);
      if(gens[history.at(1)].PID()==l0.PID()) l0+=*photon;
      else if(gens[history.at(1)].PID()==l1.PID()) l1+=*photon;
    }
  }else if(mode>=1){
    double delr=mode==1?0.1:0.4;
    for(const auto& photon:photons){
      if(l0.DeltaR(*photon)>delr&&l1.DeltaR(*photon)>delr) continue;
      if(l0.DeltaR(*photon)<l1.DeltaR(*photon)) l0+=*photon;
      else l1+=*photon;
    }
  }
}
