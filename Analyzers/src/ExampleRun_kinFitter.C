#include "ExampleRun_kinFitter.h"

void ExampleRun_kinFitter::initializeAnalyzer(){
  SMPAnalyzerCore::initializeAnalyzer(); //setup eff zpt roc z0
  std::vector<JetTagging::Parameters> jtps = {
    JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::comb),
    JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Tight, JetTagging::incl, JetTagging::comb)
  };
  mcCorr->SetJetTaggingParameters(jtps);
  SetupPUJetWeight();

  fitter = new TKinFitterDriver(DataYear);
}

void ExampleRun_kinFitter::executeEvent(){
  ///////////////// GEN level /////////////////////
  //executeEventGen();

  ///////////////// RECO level /////////////////////
  if(!IsDATA || DataStream.Contains("SingleMuon")){
    executeEventWithParameter("m"+GetEraShort());
  }
  if(!IsDATA || DataStream.Contains("SingleElectron") || DataStream.Contains("EGamma")){
    executeEventWithParameter("e"+GetEraShort());
    executeEventWithParameter("E"+GetEraShort());
  }
}

void ExampleRun_kinFitter::executeEventWithParameter(TString channel){

  lepton0 = NULL;
  prefix = channel+"/"+gprefix, hprefix = "", suffix = "";

  // Weights Setup
  if(!IsDATA){
    lumiweight = reductionweight * MCweight()*_event.GetTriggerLumi("Full");
    PUweight = mcCorr->GetPileUpWeight(nPileUp,0);
    prefireweight = L1PrefireReweight_Central;
  }

  map_weight[""] = lumiweight;
  FillHist(prefix+hprefix+"weight_Lumi", lumiweight, map_weight[""], 200,-5,5);
  FillCutflow(prefix+hprefix+"cutflow"+suffix, "Lumi", map_weight[""]);

  // Trigger
  if(!IsFiredTriggers(channel)) return;
  FillCutflow(prefix+hprefix+"cutflow"+suffix, "PassTrig", map_weight[""]);

  // MET Filter
  if(!PassMETFilter()) return;
  FillCutflow(prefix+hprefix+"cutflow"+suffix, "METfilter", map_weight[""]);

  // Single Lepton + pT + MET
  if(!Hasleptons(channel)) return;

  // Jets
  vector<Jet> alljets = SelectJets(GetAllJets(), "tightLepVeto", 30, 2.4);
  vector<Jet> lepvetojets = {}, realjets = {}, bjets = {}, ajets = {};
  for(const auto jet:alljets){
    if(lepton0 && jet.DeltaR(*lepton0) < 0.4) continue;
    lepvetojets.push_back(jet);
  }
  for(const auto jet:lepvetojets){
    if(!PUJetIDPass(jet, "Loose")) continue;
    realjets.push_back(jet);
  }

  // B-tagging
  JetTagging::Parameters DeepJet_Tight = JetTagging::Parameters(JetTagging::DeepJet,JetTagging::Tight,JetTagging::incl,JetTagging::comb);
  JetTagging::Parameters DeepJet_Medium = JetTagging::Parameters(JetTagging::DeepJet,JetTagging::Medium,JetTagging::incl,JetTagging::comb);

  std::vector<bool> btag_vector = {};
  for(const auto& jet:realjets){
    if(jet.GetTaggerResult(DeepJet_Tight.j_Tagger) > mcCorr->GetJetTaggingCutValue(DeepJet_Tight.j_Tagger, DeepJet_Tight.j_WP)){
      btag_vector.push_back(true);
      bjets.push_back(jet);
    }else{
      btag_vector.push_back(false);
      ajets.push_back(jet);
    }
  }

  // Jet related weights
  if(!IsDATA){
    pujetSF = GetPUJetWeight(lepvetojets, "Loose", 0);
    btagSF = mcCorr->GetBTaggingReweight_1a(realjets, DeepJet_Tight, "central");
  }

  FillHist(prefix+"njets", realjets.size(), map_weight[""], 15,0,15);
  FillHist(prefix+"nbjets", bjets.size(), map_weight[""], 10,0,10);
  FillHist(prefix+"najets", ajets.size(), map_weight[""], 10,0,10);

  if(bjets.size() != 2) return;
  FillCutflow(prefix+hprefix+"cutflow"+suffix, "Tight2b", map_weight[""]);
  if(realjets.size() < 4) return;
  FillCutflow(prefix+hprefix+"cutflow"+suffix, "4Jets", map_weight[""]);
  bjet0 = &bjets.at(0); bcharge0 = jetCharge(*bjet0);
  bjet1 = &bjets.at(1); bcharge1 = jetCharge(*bjet1);
  ajet0 = &ajets.at(0); acharge0 = jetCharge(*ajet0);
  ajet1 = &ajets.at(1); acharge1 = jetCharge(*ajet1);

  map_weight["_noWts"] = map_weight[""];
  // Weights
  map_weight[""] *= PUweight;
  FillHist(prefix+hprefix+"weight_PU", PUweight, map_weight[""], 200,-5,5);
  FillCutflow(prefix+hprefix+"cutflow"+suffix, "PU", map_weight[""]);

  map_weight[""] *= prefireweight;
  FillHist(prefix+hprefix+"weight_Prefire", prefireweight, map_weight[""], 200,-5,5);
  FillCutflow(prefix+hprefix+"cutflow"+suffix, "Prefire", map_weight[""]);

  map_weight[""] *= topptweight;
  FillHist(prefix+hprefix+"weight_Toppt", topptweight, map_weight[""], 200,-5,5);
  FillCutflow(prefix+hprefix+"cutflow"+suffix, "Toppt", map_weight[""]);

  // Lepton Efficiency Correction
  map_weight["_noEffSF"] = map_weight[""];
  leptonTrackingSF = 1.;
  leptonRECOSF = 1.;
  leptonIDSF = 1.;
  leptonTriggerSF = 1.;

  if(!IsDATA){
    TString trigSFkey = "IsoMu24_MediumID_trkIsoLoose";
    if(channel.Contains("m")){
      leptonTrackingSF *= fEff->GetEfficiencySF("Muon_Tracking", lepton0, 0,0);
      leptonRECOSF *= fEff->GetEfficiencySF("Muon_RECO", lepton0, 0,0);
      leptonIDSF *= fEff->GetEfficiencySF("Muon_MediumID_trkIsoLoose", lepton0, 0,0);
      if(DataYear == 2017) trigSFkey = "IsoMu27_MediumID_trkIsoLoose";
      leptonTriggerSF *= GetLeptonTriggerSF(trigSFkey, leptons, 0,0);
    }else if(channel.Contains("e")){
      leptonRECOSF *= fEff->GetEfficiencySF("Electron_RECO", lepton0, 0,0);
      leptonIDSF *= fEff->GetEfficiencySF("Electron_MediumID", lepton0, 0,0);
      trigSFkey = "Ele27_MediumID";
      if(DataYear == 2017) trigSFkey = "Ele32_MediumID";
      leptonTriggerSF *= GetLeptonTriggerSF(trigSFkey, leptons, 0,0);
    }else if(channel.Contains("E")){
      leptonRECOSF *= fEff->GetEfficiencySF("Electron_RECO", lepton0, 0,0);
      leptonIDSF *= fEff->GetEfficiencySF("Electron_SelQ_MediumID", lepton0, 0,0);
      trigSFkey = "Ele27_SelQ_MediumID";
      if(DataYear == 2017) trigSFkey = "Ele32_SelQ_MediumID";
      leptonTriggerSF *= GetLeptonTriggerSF(trigSFkey, leptons, 0,0);
    }
  }

  map_weight[""] *= leptonTrackingSF;
  FillHist(prefix+hprefix+"weight_TrackingSF", leptonTrackingSF, map_weight[""], 200,-5,5);
  FillCutflow(prefix+hprefix+"cutflow"+suffix, "TrackingSF", map_weight[""]);

  map_weight[""] *= leptonRECOSF;
  FillHist(prefix+hprefix+"weight_RECOSF", leptonRECOSF, map_weight[""], 200,-5,5);
  FillCutflow(prefix+hprefix+"cutflow"+suffix, "RECOSF", map_weight[""]);

  map_weight[""] *= leptonIDSF;
  FillHist(prefix+hprefix+"weight_IDSF", leptonIDSF, map_weight[""], 200,-5,5);
  FillCutflow(prefix+hprefix+"cutflow"+suffix, "IDSF", map_weight[""]);

  map_weight[""] *= leptonTriggerSF;
  FillHist(prefix+hprefix+"weight_TriggerSF", leptonTriggerSF, map_weight[""], 200,-5,5);
  FillCutflow(prefix+hprefix+"cutflow"+suffix, "TriggerSF", map_weight[""]);

  map_weight["_noPUjetSF"] = map_weight[""];
  map_weight[""] *= pujetSF;
  FillHist(prefix+hprefix+"weight_PUjetSF", pujetSF, map_weight[""], 200,-5,5);
  FillCutflow(prefix+hprefix+"cutflow"+suffix, "PUjetSF", map_weight[""]);

  map_weight["_nobtagSF"] = map_weight[""];
  map_weight[""] *= btagSF;
  FillHist(prefix+hprefix+"weight_btagSF", btagSF, map_weight[""], 200,-5,5);
  FillCutflow(prefix+hprefix+"cutflow"+suffix, "btagSF", map_weight[""]);

  map_weight["_bChargeSF0"] = map_weight[""] * GetbChargeSFWeight(bjets, 0, 0);
  FillHist(prefix+hprefix+"weight_bChargeSF0", GetbChargeSFWeight(bjets, 0, 0), map_weight[""], 200,-5,5);
  map_weight["_bChargeSF1"] = map_weight[""] * GetbChargeSFWeight(bjets, 1, 0);
  FillHist(prefix+hprefix+"weight_bChargeSF1", GetbChargeSFWeight(bjets, 1, 0), map_weight[""], 200,-5,5);

  //=======================
  //==== Kinematic Fitter
  //=======================
  std::vector<TLorentzVector> jet_vector{};
  std::vector<double> jet_pt_resolution_vector{};
  TLorentzVector lepton{};
  for(auto& jet : realjets){
    jet_vector.emplace_back(jet.Px(),jet.Py(),jet.Pz(),jet.E());
    jet_pt_resolution_vector.push_back((jet.EnShift(1) - jet.EnShift(-1)) * 2);
  }
  if(jet_vector.size() != btag_vector.size()){
    cout << " ExampleRun_kinFitter, jet_vector.size() != btag_vector.size()" << endl;
    exit(1);
  }
  lepton = (TLorentzVector)(*lepton0);
  fitter->SetAllObjects(jet_vector, btag_vector, jet_pt_resolution_vector, lepton, met);
  fitter->FindBestChi2Fit();
  auto fitter_results = fitter->GetResults();

  FillHist(prefix+"kinFit_Nresults", fitter_results->size(), map_weight, 20,0,20);
  FillHist(prefix+"kinFit_efficiency", fitter_results->size()<1?0:1, map_weight, 2,0,2);
  if(fitter_results->size() < 1) return;
  FillCutflow(prefix+hprefix+"cutflow"+suffix, "kinFit", map_weight[""]);

  FillHist(prefix+"kinFit_chi2", fitter_results->at(0).chi2, map_weight, 500,0,250);
  FillHist(prefix+"kinFit_chi2_lep", fitter_results->at(0).chi2_lep, map_weight, 500,0,250);
  FillHist(prefix+"kinFit_chi2_had", fitter_results->at(0).chi2_had, map_weight, 500,0,250);
  FillHist(prefix+"kinFit_chi2_mass", fitter_results->at(0).chi2_mass, map_weight, 500,0,250);

  // Before Fit variables
  int lepb_idx = fitter_results->at(0).leptonic_top_b_jet_idx;
  int hadb_idx = fitter_results->at(0).hadronic_top_b_jet_idx;
  int W_up_jet_idx = fitter_results->at(0).w_ch_up_type_jet_idx;
  int W_down_jet_idx = fitter_results->at(0).w_ch_down_type_jet_idx;
  Jet W_up_jet = realjets.at(W_up_jet_idx);
  Jet W_down_jet = realjets.at(W_down_jet_idx);
  Jet lepb = realjets.at(lepb_idx);
  Jet hadb = realjets.at(hadb_idx);
  double lepb_charge = jetCharge(lepb);
  double hadb_charge = jetCharge(hadb);

  double hadronic_W_M = (W_up_jet+W_down_jet).M();
  double leptonic_W_M = fitter_results->at(0).leptonic_W_M;
  double hadronic_top_M = fitter_results->at(0).hadronic_top_M;
  double leptonic_top_M = fitter_results->at(0).leptonic_top_M;

  // After Fit variables
  TLorentzVector fitted_lepb = fitter_results->at(0).fitted_lept_bjet;
  TLorentzVector fitted_hadb = fitter_results->at(0).fitted_hadt_bjet;
  TLorentzVector fitted_W_j1 = fitter_results->at(0).fitted_jet1;
  TLorentzVector fitted_W_j2 = fitter_results->at(0).fitted_jet2;
  TLorentzVector fitted_lep = fitter_results->at(0).fitted_lep;
  TLorentzVector fitted_neu = fitter_results->at(0).fitted_neu;

  double fitted_hadronic_W_M = (fitted_W_j1+fitted_W_j2).M();
  double fitted_leptonic_W_M = (fitted_lep+fitted_neu).M();
  double fitted_hadronic_top_M = (fitted_hadb+fitted_W_j1+fitted_W_j2).M();
  double fitted_leptonic_top_M = (fitted_lepb+fitted_lep+fitted_neu).M();

  //============================================================
  //==== Gen - Jet Matching (To check the performance of Fitter)
  //============================================================
  if(!IsDATA && MCSample.Contains("TTLJ")){
    vector<Gen> gens=GetGens();
    Gen gen_parton0,gen_parton1, gen_b0,gen_b1, gen_l0,gen_l1, gen_j0,gen_j1;
    GetTTLJGenParticles(gens, gen_parton0,gen_parton1, gen_b0,gen_b1, gen_l0,gen_l1, gen_j0,gen_j1,3);

    FillHist(prefix+"gen_W_had_Mass", (gen_j0+gen_j1).M(), map_weight[""], 50,0,200);
    FillHist(prefix+"gen_W_lep_Mass", (gen_l0+gen_l1).M(), map_weight[""], 50,0,200);
    FillHist(prefix+"gen_top_had1_Mass", (gen_b0+gen_j0+gen_j1).M(), map_weight[""], 75,0,300);
    FillHist(prefix+"gen_top_had2_Mass", (gen_b1+gen_j0+gen_j1).M(), map_weight[""], 75,0,300);
    FillHist(prefix+"gen_top_lep1_Mass", (gen_b0+gen_l0+gen_l1).M(), map_weight[""], 75,0,300);
    FillHist(prefix+"gen_top_lep2_Mass", (gen_b1+gen_l0+gen_l1).M(), map_weight[""], 75,0,300);
    if(gen_l0.PID()<0){ //gen_l0 = mu+ el+
      FillHist(prefix+"gen_top_lep_Mass", (gen_b0+gen_l0+gen_l1).M(), map_weight[""], 75,0,300);
      FillHist(prefix+"gen_top_had_Mass", (gen_b1+gen_j0+gen_j1).M(), map_weight[""], 75,0,300);
    }else{
      FillHist(prefix+"gen_top_lep_Mass", (gen_b1+gen_l0+gen_l1).M(), map_weight[""], 75,0,300);
      FillHist(prefix+"gen_top_had_Mass", (gen_b0+gen_j0+gen_j1).M(), map_weight[""], 75,0,300);
    }

    double match_dR = 0.4;
    bool Gen_kinFit_match_onlyb = false;
    bool Gen_kinFit_match_full = false;
    // chargeEasy => 0:negative, 1:positive
    // purity => -1:UnMatched_, 0:Wrong_, 1:Correct_ in TTLJ
    // correct => 0:Wrong_, 1:Correct_ in TTLJ

    //When lepb = b, hadb = bbar -> lep+
    if(gen_b0.DeltaR(lepb) < match_dR && gen_b1.DeltaR(hadb) < match_dR){
      Gen_kinFit_match_onlyb = true;
      if((gen_j0.DeltaR(W_up_jet) < match_dR && gen_j1.DeltaR(W_down_jet) < match_dR) || (gen_j0.DeltaR(W_down_jet) < match_dR && gen_j1.DeltaR(W_up_jet) < match_dR)) Gen_kinFit_match_full = true;

      FillHist(prefix+"Matched_bjetCharge", lepb_charge, map_weight, 200,-5,5);
      FillHist(prefix+"Matched_bbarjetCharge", hadb_charge, map_weight, 200,-5,5);
      if(lepton0->Charge() > 0){
        FillHist(prefix+"kinFit_purity", 1, map_weight, 4,-2,2);
        FillHist(prefix+"kinFit_correct", 1, map_weight, 2,0,2);
        FillHist(prefix+"kinFit_purity_Lp", 1, map_weight, 4,-2,2);
        FillHist(prefix+"kinFit_correct_Lp", 1, map_weight, 2,0,2);
        prefix += "Correct_"; // lep+
      }else{
        FillHist(prefix+"kinFit_purity", 0, map_weight, 4,-2,2);
        FillHist(prefix+"kinFit_correct", 0, map_weight, 2,0,2);
        FillHist(prefix+"kinFit_purity_Lm", 0, map_weight, 4,-2,2);
        FillHist(prefix+"kinFit_correct_Lm", 0, map_weight, 2,0,2);
        prefix +="Wrong_"; // lep-
      }

      for(unsigned int i=1; i<afb_chbinnum+1; i++){
        if(afb_chbin[i-1] < abs(lepb_charge) && abs(lepb_charge) < afb_chbin[i]){
          FillHist(Form(prefix+"lepbjetCharge%d_Matched_b",i-1), lepb_charge, map_weight, 200,-5,5);
          FillHist(Form(prefix+"lepbjetCharge%dEasy_Matched_b",i-1), lepb_charge<0?0:1, map_weight, 2,0,2);
        }
        if(afb_chbin[i-1] < abs(hadb_charge) && abs(hadb_charge) < afb_chbin[i]){
          FillHist(Form(prefix+"hadbjetCharge%d_Matched_bbar",i-1), hadb_charge, map_weight, 200,-5,5);
          FillHist(Form(prefix+"hadbjetCharge%dEasy_Matched_bbar",i-1), hadb_charge<0?0:1, map_weight, 2,0,2);
        }
        if(-1 * afb_chbin[i] < lepb_charge && lepb_charge < -1 * afb_chbin[i-1]){
          FillHist(Form(prefix+"lcharge_lepb_negaive%d_Matched_b",i-1), lepton0->Charge(), map_weight, 4,-2,2);
          FillHist(Form(prefix+"lchargeEasy_lepb_negaive%d_Matched_b",i-1), lepton0->Charge()<0?0:1, map_weight, 2,0,2);
        }
        if(afb_chbin[i-1] < hadb_charge && hadb_charge < afb_chbin[i]){
          FillHist(Form(prefix+"lcharge_hadb_positive%d_Matched_bbar",i-1), lepton0->Charge(), map_weight, 4,-2,2);
          FillHist(Form(prefix+"lchargeEasy_hadb_positive%d_Matched_bbar",i-1), lepton0->Charge()<0?0:1, map_weight, 2,0,2);
        }
      }
      FillHist(prefix+"lepbjetCharge_Matched_b", lepb_charge, map_weight, 200,-5,5);
      FillHist(prefix+"lepbjetChargeEasy_Matched_b", lepb_charge<0?0:1, map_weight, 2,0,2);
      FillHist(prefix+"hadbjetCharge_Matched_bbar", hadb_charge, map_weight, 200,-5,5);
      FillHist(prefix+"hadbjetChargeEasy_Matched_bbar", hadb_charge<0?0:1, map_weight, 2,0,2);
      FillHist(prefix+"lcharge_lepb_Matched_b", lepton0->Charge(), map_weight, 4,-2,2);
      FillHist(prefix+"lchargeEasy_lepb_Matched_b", lepton0->Charge()<0?0:1, map_weight, 2,0,2);
      if(lepb_charge < 0){
	FillHist(prefix+"lcharge_lepb_negaive_Matched_b", lepton0->Charge(), map_weight, 4,-2,2);
	FillHist(prefix+"lchargeEasy_lepb_negaive_Matched_b", lepton0->Charge()<0?0:1, map_weight, 2,0,2);
      }
      if(0 < hadb_charge){
	FillHist(prefix+"lcharge_hadb_positive_Matched_bbar", lepton0->Charge(), map_weight, 4,-2,2);
	FillHist(prefix+"lchargeEasy_hadb_positive_Matched_bbar", lepton0->Charge()<0?0:1, map_weight, 2,0,2);
      }
    }
    //When lepb = bbar, hadb = b => lep-
    else if(gen_b0.DeltaR(hadb) < match_dR && gen_b1.DeltaR(lepb) < match_dR){
      Gen_kinFit_match_onlyb = true;
      if((gen_j0.DeltaR(W_up_jet) < match_dR && gen_j1.DeltaR(W_down_jet) < match_dR) || (gen_j0.DeltaR(W_down_jet) < match_dR && gen_j1.DeltaR(W_up_jet) < match_dR)) Gen_kinFit_match_full = true;

      FillHist(prefix+"Matched_bbarjetCharge", lepb_charge, map_weight, 200,-5,5);
      FillHist(prefix+"Matched_bjetCharge", hadb_charge, map_weight, 200,-5,5);
      if(lepton0->Charge() < 0){
        FillHist(prefix+"kinFit_purity", 1, map_weight, 4,-2,2);
        FillHist(prefix+"kinFit_correct", 1, map_weight, 2,0,2);
        FillHist(prefix+"kinFit_purity_Lm", 1, map_weight, 4,-2,2);
        FillHist(prefix+"kinFit_correct_Lm", 1, map_weight, 2,0,2);
        prefix += "Correct_"; // lep-
      }else{
        FillHist(prefix+"kinFit_purity", 0, map_weight, 4,-2,2);
        FillHist(prefix+"kinFit_correct", 0, map_weight, 2,0,2);
        FillHist(prefix+"kinFit_purity_Lp", 0, map_weight, 4,-2,2);
        FillHist(prefix+"kinFit_correct_Lp", 0, map_weight, 2,0,2);
        prefix +="Wrong_"; // lep+
      }

      for(unsigned int i=1; i<afb_chbinnum+1; i++){
        if(afb_chbin[i-1] < abs(hadb_charge) && abs(hadb_charge) < afb_chbin[i]){
          FillHist(Form(prefix+"hadbjetCharge%d_Matched_b",i-1), hadb_charge, map_weight, 200,-5,5);
          FillHist(Form(prefix+"hadbjetCharge%dEasy_Matched_b",i-1), hadb_charge<0?0:1, map_weight, 2,0,2);
        }
        if(afb_chbin[i-1] < abs(lepb_charge) && abs(lepb_charge) < afb_chbin[i]){
          FillHist(Form(prefix+"lepbjetCharge%d_Matched_bbar",i-1), lepb_charge, map_weight, 200,-5,5);
          FillHist(Form(prefix+"lepbjetCharge%dEasy_Matched_bbar",i-1), lepb_charge<0?0:1, map_weight, 2,0,2);
        }
        if(afb_chbin[i-1] < lepb_charge && lepb_charge < afb_chbin[i]){
          FillHist(Form(prefix+"lcharge_lepb_positive%d_Matched_bbar",i-1), lepton0->Charge(), map_weight, 4,-2,2);
          FillHist(Form(prefix+"lchargeEasy_lepb_positive%d_Matched_bbar",i-1), lepton0->Charge()<0?0:1, map_weight, 2,0,2);
        }
        if(-1 * afb_chbin[i] < hadb_charge && hadb_charge < -1 * afb_chbin[i-1]){
          FillHist(Form(prefix+"lcharge_hadb_negative%d_Matched_b",i-1), lepton0->Charge(), map_weight, 4,-2,2);
          FillHist(Form(prefix+"lchargeEasy_hadb_negative%d_Matched_b",i-1), lepton0->Charge()<0?0:1, map_weight, 2,0,2);
        }
      }
      FillHist(prefix+"hadbjetCharge_Matched_b", hadb_charge, map_weight, 200,-5,5);
      FillHist(prefix+"hadbjetChargeEasy_Matched_b", hadb_charge<0?0:1, map_weight, 2,0,2);
      FillHist(prefix+"lepbjetCharge_Matched_bbar", lepb_charge, map_weight, 200,-5,5);
      FillHist(prefix+"lepbjetChargeEasy_Matched_bbar", lepb_charge<0?0:1, map_weight, 2,0,2);
      FillHist(prefix+"lcharge_lepb_Matched_bbar", lepton0->Charge(), map_weight, 4,-2,2);
      FillHist(prefix+"lchargeEasy_lepb_Matched_bbar", lepton0->Charge()<0?0:1, map_weight, 2,0,2);
      if(0 < lepb_charge){
        FillHist(prefix+"lcharge_lepb_positive_Matched_bbar", lepton0->Charge(), map_weight, 4,-2,2);
        FillHist(prefix+"lchargeEasy_lepb_positive_Matched_bbar", lepton0->Charge()<0?0:1, map_weight, 2,0,2);
      }
      if(hadb_charge < 0){
        FillHist(prefix+"lcharge_hadb_negative_Matched_b", lepton0->Charge(), map_weight, 4,-2,2);
        FillHist(prefix+"lchargeEasy_hadb_negative_Matched_b", lepton0->Charge()<0?0:1, map_weight, 2,0,2);
      }
    }
    else{
      FillHist(prefix+"kinFit_purity", -1, map_weight, 4,-2,2);
      if(lepton0->Charge() > 0) FillHist(prefix+"kinFit_purity_Lp", -1, map_weight, 4,-2,2);
      if(lepton0->Charge() < 0) FillHist(prefix+"kinFit_purity_Lm", -1, map_weight, 4,-2,2);
      prefix +="UnMatched_";
    }

    FillHist(prefix+"Gen_b0_mindR", min(gen_b0.DeltaR(lepb),gen_b0.DeltaR(hadb)), map_weight[""], 50,0,5);
    FillHist(prefix+"Gen_b1_mindR", min(gen_b1.DeltaR(lepb),gen_b1.DeltaR(hadb)), map_weight[""], 50,0,5);
    FillHist(prefix+"Gen_j0_mindR", min(gen_j0.DeltaR(W_up_jet),gen_j0.DeltaR(W_down_jet)), map_weight[""], 50,0,5);
    FillHist(prefix+"Gen_j1_mindR", min(gen_j1.DeltaR(W_up_jet),gen_j1.DeltaR(W_down_jet)), map_weight[""], 50,0,5);
    FillHist(prefix+"Gen_l0_mindR", min(gen_l0.DeltaR(fitted_lep),gen_l0.DeltaR(fitted_neu)), map_weight[""], 50,0,5);
    FillHist(prefix+"Gen_l1_mindR", min(gen_l1.DeltaR(fitted_neu),gen_l1.DeltaR(fitted_lep)), map_weight[""], 50,0,5);
    FillHist(prefix+"Gen_b0b1_dR", gen_b0.DeltaR(gen_b1), map_weight[""], 50,0,5);
    FillHist(prefix+"Gen_b0j0_dR", gen_b0.DeltaR(gen_j0), map_weight[""], 50,0,5);
    FillHist(prefix+"Gen_b0j1_dR", gen_b0.DeltaR(gen_j1), map_weight[""], 50,0,5);
    FillHist(prefix+"Gen_b1j0_dR", gen_b1.DeltaR(gen_j0), map_weight[""], 50,0,5);
    FillHist(prefix+"Gen_b1j1_dR", gen_b1.DeltaR(gen_j1), map_weight[""], 50,0,5);
    FillHist(prefix+"Gen_j0j1_dR", gen_j0.DeltaR(gen_j1), map_weight[""], 50,0,5);
    FillHist(prefix+"Gen_b0_dR", gen_b0.DeltaR(lepb), map_weight[""], 50,0,5);
    FillHist(prefix+"Gen_b1_dR", gen_b1.DeltaR(hadb), map_weight[""], 50,0,5);
    FillHist(prefix+"Gen_j0_dR", gen_j0.DeltaR(W_up_jet), map_weight[""], 50,0,5);
    FillHist(prefix+"Gen_j1_dR", gen_j1.DeltaR(W_down_jet), map_weight[""], 50,0,5);
    FillHist(prefix+"Gen_l0_dR", gen_l0.DeltaR(fitted_lep), map_weight[""], 50,0,5);
    FillHist(prefix+"Gen_l1_dR", gen_l1.DeltaR(fitted_neu), map_weight[""], 50,0,5);
    FillHist(prefix+"Gen_b0_pT", gen_b0.Pt(), map_weight[""], 100,0,300);
    FillHist(prefix+"Gen_b1_pT", gen_b1.Pt(), map_weight[""], 100,0,300);
    FillHist(prefix+"Gen_j0_pT", gen_j0.Pt(), map_weight[""], 100,0,300);
    FillHist(prefix+"Gen_j1_pT", gen_j1.Pt(), map_weight[""], 100,0,300);
    FillHist(prefix+"Gen_b0_pTFitdiff", gen_b0.Pt()-fitted_lepb.Pt(), map_weight[""], 100,-50,50);
    FillHist(prefix+"Gen_b1_pTFitdiff", gen_b1.Pt()-fitted_hadb.Pt(), map_weight[""], 100,-50,50);
    FillHist(prefix+"Gen_j0_pTFitdiff", gen_j0.Pt()-fitted_W_j1.Pt(), map_weight[""], 100,-50,50);
    FillHist(prefix+"Gen_j1_pTFitdiff", gen_j1.Pt()-fitted_W_j2.Pt(), map_weight[""], 100,-50,50);
    FillHist(prefix+"Gen_l0_pTFitdiff", gen_l0.Pt()-fitted_lep.Pt(), map_weight[""], 100,-50,50);
    FillHist(prefix+"Gen_l1_pTFitdiff", gen_l1.Pt()-fitted_neu.Pt(), map_weight[""], 100,-50,50);
    FillHist(prefix+"Gen_b0_pTdiff", gen_b0.Pt()-lepb.Pt(), map_weight[""], 100,-50,50);
    FillHist(prefix+"Gen_b1_pTdiff", gen_b1.Pt()-hadb.Pt(), map_weight[""], 100,-50,50);
    FillHist(prefix+"Gen_j0_pTdiff", gen_j0.Pt()-W_up_jet.Pt(), map_weight[""], 100,-50,50);
    FillHist(prefix+"Gen_j1_pTdiff", gen_j1.Pt()-W_down_jet.Pt(), map_weight[""], 100,-50,50);
    FillHist(prefix+"Gen_l0_pTdiff", gen_l0.Pt()-lepton0->Pt(), map_weight[""], 100,-50,50);
    FillHist(prefix+"Gen_l1_pTdiff", gen_l1.Pt()-met.Pt(), map_weight[""], 100,-50,50);

    FillHist("Gen_kinFit_match_onlyb", Gen_kinFit_match_onlyb, map_weight[""], 2,0,2);
    FillHist("Gen_kinFit_match_full", Gen_kinFit_match_full, map_weight[""], 2,0,2);
    FillHist(prefix+"Gen_kinFit_match_onlyb", Gen_kinFit_match_onlyb, map_weight[""], 2,0,2);
    FillHist(prefix+"Gen_kinFit_match_full", Gen_kinFit_match_full, map_weight[""], 2,0,2);
    FillHist(channel+"/Gen_kinFit_match_onlyb", Gen_kinFit_match_onlyb, map_weight[""], 2,0,2);
    FillHist(channel+"/Gen_kinFit_match_full", Gen_kinFit_match_full, map_weight[""], 2,0,2);

    FillHist(prefix+"gen_W_had_Mass", (gen_j0+gen_j1).M(), map_weight[""], 50,0,200);
    FillHist(prefix+"gen_W_lep_Mass", (gen_l0+gen_l1).M(), map_weight[""], 50,0,200);
    FillHist(prefix+"gen_top_had1_Mass", (gen_b0+gen_j0+gen_j1).M(), map_weight[""], 75,0,300);
    FillHist(prefix+"gen_top_had2_Mass", (gen_b1+gen_j0+gen_j1).M(), map_weight[""], 75,0,300);
    FillHist(prefix+"gen_top_lep1_Mass", (gen_b0+gen_l0+gen_l1).M(), map_weight[""], 75,0,300);
    FillHist(prefix+"gen_top_lep2_Mass", (gen_b1+gen_l0+gen_l1).M(), map_weight[""], 75,0,300);
    if(gen_l0.PID()<0){ //gen_l0 = mu+ el+
      FillHist(prefix+"gen_top_lep_Mass", (gen_b0+gen_l0+gen_l1).M(), map_weight[""], 75,0,300);
      FillHist(prefix+"gen_top_had_Mass", (gen_b1+gen_j0+gen_j1).M(), map_weight[""], 75,0,300);
    }else{
      FillHist(prefix+"gen_top_lep_Mass", (gen_b1+gen_l0+gen_l1).M(), map_weight[""], 75,0,300);
      FillHist(prefix+"gen_top_had_Mass", (gen_b0+gen_j0+gen_j1).M(), map_weight[""], 75,0,300);
    }
  }

  //==========================
  //==== Now reco fill histograms
  //==========================

  FillHist(prefix+"W_had_Mass", hadronic_W_M, map_weight, 40,0,200);
  FillHist(prefix+"W_had_FitMass", fitted_hadronic_W_M, map_weight, 80,60,100);
  FillHist(prefix+"W_had_FitMassdiff", fitted_hadronic_W_M-hadronic_W_M, map_weight, 100,-50,50);
  FillHist(prefix+"W_lep_Mass", leptonic_W_M, map_weight, 40,0,200);
  FillHist(prefix+"W_lep_FitMass", fitted_leptonic_W_M, map_weight, 150,75,150);
  FillHist(prefix+"W_lep_FitMassdiff", fitted_leptonic_W_M-leptonic_W_M, map_weight, 100,-50,50);
  FillHist(prefix+"Top_had_Mass", hadronic_top_M, map_weight, 40,100,300);
  FillHist(prefix+"Top_had_FitMass", fitted_hadronic_top_M, map_weight, 160,160,240);
  FillHist(prefix+"Top_had_FitMassdiff", fitted_hadronic_top_M-hadronic_top_M, map_weight, 100,-50,50);
  FillHist(prefix+"Top_lep_Mass", leptonic_top_M, map_weight, 40,100,300);
  FillHist(prefix+"Top_lep_FitMass", fitted_leptonic_top_M, map_weight, 80,160,200);
  FillHist(prefix+"Top_lep_FitMassdiff", fitted_leptonic_top_M-leptonic_top_M, map_weight, 100,-50,50);

  FillHist(prefix+"lepb_FitpT", fitted_lepb.Pt(), map_weight, 100,0,300);
  FillHist(prefix+"hadb_FitpT", fitted_hadb.Pt(), map_weight, 100,0,300);
  FillHist(prefix+"W_had1_FitpT", fitted_W_j1.Pt(), map_weight, 100,0,300);
  FillHist(prefix+"W_had2_FitpT", fitted_W_j2.Pt(), map_weight, 100,0,300);
  FillHist(prefix+"lepb_FitpTdiff", fitted_lepb.Pt()-lepb.Pt(), map_weight, 100,-50,50);
  FillHist(prefix+"hadb_FitpTdiff", fitted_hadb.Pt()-hadb.Pt(), map_weight, 100,-50,50);
  FillHist(prefix+"W_had1_FitpTdiff", fitted_W_j1.Pt()-W_up_jet.Pt(), map_weight, 100,-50,50);
  FillHist(prefix+"W_had2_FitpTdiff", fitted_W_j2.Pt()-W_down_jet.Pt(), map_weight, 100,-50,50);
  FillHist(prefix+"lepb_hadb_FitdR", fitted_lepb.DeltaR(fitted_hadb), map_weight, 100,0,5);
  FillHist(prefix+"lepb_hadb_FitdPhi", fitted_lepb.DeltaPhi(fitted_hadb), map_weight, 200,-5,5);
  FillHist(prefix+"W_had1_had2_FitdR", fitted_W_j1.DeltaR(fitted_W_j2), map_weight, 100,0,5);
  FillHist(prefix+"W_had1_had2_FitdPhi", fitted_W_j1.DeltaPhi(fitted_W_j2), map_weight, 200,-5,5);
  FillHist(prefix+"lepb_FitP", fitted_lepb.P(), map_weight, 100,0,300);
  FillHist(prefix+"hadb_FitP", fitted_hadb.P(), map_weight, 100,0,300);
  FillHist(prefix+"W_had1_FitP", fitted_W_j1.P(), map_weight, 100,0,300);
  FillHist(prefix+"W_had2_FitP", fitted_W_j2.P(), map_weight, 100,0,300);

  FillHist(prefix+"W_up_jet_idx",  W_up_jet_idx, map_weight, 10,0,10);
  FillHist(prefix+"W_down_jet_idx",  W_down_jet_idx, map_weight, 10,0,10);
  FillHist(prefix+"lepb_idx", lepb_idx, map_weight, 10,0,10);
  FillHist(prefix+"hadb_idx", hadb_idx, map_weight, 10,0,10);
  for(unsigned int i=1; i<afb_chbinnum+1; i++){
    if(lepton0->Charge() < 0){
      if(afb_chbin[i-1] < abs(lepb_charge) && abs(lepb_charge) < afb_chbin[i]){
        FillHist(Form(prefix+"lepbjetcharge%d_Lm",i-1), lepb.Charge(), map_weight, 200,-2,2);
        FillHist(Form(prefix+"lepbjetcharge%dEasy_Lm",i-1), lepb.Charge()<0?0:1, map_weight, 2,0,2);
        FillHist(Form(prefix+"lepbjetCharge%d_Lm",i-1), lepb_charge, map_weight, 200,-5,5);
        FillHist(Form(prefix+"lepbjetCharge%dEasy_Lm",i-1), lepb_charge<0?0:1, map_weight, 2,0,2);
      }
      if(afb_chbin[i-1] < abs(hadb_charge) && abs(hadb_charge) < afb_chbin[i]){
        FillHist(Form(prefix+"hadbjetcharge%d_Lm",i-1), hadb.Charge(), map_weight, 200,-2,2);
        FillHist(Form(prefix+"hadbjetcharge%dEasy_Lm",i-1), hadb.Charge()<0?0:1, map_weight, 2,0,2);
        FillHist(Form(prefix+"hadbjetCharge%d_Lm",i-1), hadb_charge, map_weight, 200,-5,5);
        FillHist(Form(prefix+"hadbjetCharge%dEasy_Lm",i-1), hadb_charge<0?0:1, map_weight, 2,0,2);
      }
    }else{
      if(afb_chbin[i-1] < abs(lepb_charge) && abs(lepb_charge) < afb_chbin[i]){
        FillHist(Form(prefix+"lepbjetcharge%d_Lp",i-1), lepb.Charge(), map_weight, 200,-2,2);
        FillHist(Form(prefix+"lepbjetcharge%dEasy_Lp",i-1), lepb.Charge()<0?0:1, map_weight, 2,0,2);
        FillHist(Form(prefix+"lepbjetCharge%d_Lp",i-1), lepb_charge, map_weight, 200,-5,5);
        FillHist(Form(prefix+"lepbjetCharge%dEasy_Lp",i-1), lepb_charge<0?0:1, map_weight, 2,0,2);
      }
      if(afb_chbin[i-1] < abs(hadb_charge) && abs(hadb_charge) < afb_chbin[i]){
        FillHist(Form(prefix+"hadbjetcharge%d_Lp",i-1), hadb.Charge(), map_weight, 200,-2,2);
        FillHist(Form(prefix+"hadbjetcharge%dEasy_Lp",i-1), hadb.Charge()<0?0:1, map_weight, 2,0,2);
        FillHist(Form(prefix+"hadbjetCharge%d_Lp",i-1), hadb_charge, map_weight, 200,-5,5);
        FillHist(Form(prefix+"hadbjetCharge%dEasy_Lp",i-1), hadb_charge<0?0:1, map_weight, 2,0,2);
      }
    }
  }
  if(lepton0->Charge() < 0){
    FillHist(prefix+"lepbjetcharge_Lm", lepb.Charge(), map_weight, 200,-2,2);
    FillHist(prefix+"lepbjetchargeEasy_Lm", lepb.Charge()<0?0:1, map_weight, 2,0,2);
    FillHist(prefix+"hadbjetcharge_Lm", hadb.Charge(), map_weight, 200,-2,2);
    FillHist(prefix+"hadbjetchargeEasy_Lm", hadb.Charge()<0?0:1, map_weight, 2,0,2);
    FillHist(prefix+"lepbjetCharge_Lm", lepb_charge, map_weight, 200,-5,5);
    FillHist(prefix+"lepbjetChargeEasy_Lm", lepb_charge<0?0:1, map_weight, 2,0,2);
    FillHist(prefix+"hadbjetCharge_Lm", hadb_charge, map_weight, 200,-5,5);
    FillHist(prefix+"hadbjetChargeEasy_Lm", hadb_charge<0?0:1, map_weight, 2,0,2);
  }else{
    FillHist(prefix+"lepbjetcharge_Lp", lepb.Charge(), map_weight, 200,-2,2);
    FillHist(prefix+"lepbjetchargeEasy_Lp", lepb.Charge()<0?0:1, map_weight, 2,0,2);
    FillHist(prefix+"hadbjetcharge_Lp", hadb.Charge(), map_weight, 200,-2,2);
    FillHist(prefix+"hadbjetchargeEasy_Lp", hadb.Charge()<0?0:1, map_weight, 2,0,2);
    FillHist(prefix+"lepbjetCharge_Lp", lepb_charge, map_weight, 200,-5,5);
    FillHist(prefix+"lepbjetChargeEasy_Lp", lepb_charge<0?0:1, map_weight, 2,0,2);
    FillHist(prefix+"hadbjetCharge_Lp", hadb_charge, map_weight, 200,-5,5);
    FillHist(prefix+"hadbjetChargeEasy_Lp", hadb_charge<0?0:1, map_weight, 2,0,2);
  }
  FillHist(prefix+"lepbjetcharge", lepb.Charge(), map_weight, 200,-2,2);
  FillHist(prefix+"lepbjetchargeEasy", lepb.Charge()<0?0:1, map_weight, 2,0,2);
  FillHist(prefix+"hadbjetcharge", hadb.Charge(), map_weight, 200,-2,2);
  FillHist(prefix+"hadbjetchargeEasy", hadb.Charge()<0?0:1, map_weight, 2,0,2);
  FillHist(prefix+"lepbjetCharge", lepb_charge, map_weight, 200,-5,5);
  FillHist(prefix+"lepbjetChargeEasy", lepb_charge<0?0:1, map_weight, 2,0,2);
  FillHist(prefix+"hadbjetCharge", hadb_charge, map_weight, 200,-5,5);
  FillHist(prefix+"hadbjetChargeEasy", hadb_charge<0?0:1, map_weight, 2,0,2);

  FillHist(prefix+hprefix+"ptl", lepton0->Pt(), map_weight, AFBAnalyzer::lptbinnum,AFBAnalyzer::lptbin);
  FillHist(prefix+hprefix+"etal", lepton0->Eta(), map_weight, 50,-2.5,2.5);
  FillHist(prefix+hprefix+"bjetcharge0"+suffix, bjets.at(0).Charge(), map_weight, 200,-2,2);
  FillHist(prefix+hprefix+"bjetcharge1"+suffix, bjets.at(1).Charge(), map_weight, 200,-2,2);
  FillHist(prefix+hprefix+"bjetschargeSum"+suffix, bjets.at(0).Charge() + bjets.at(1).Charge(), map_weight, 400,-4,4);
  FillHist(prefix+hprefix+"bjetschargeAbsSum"+suffix, (bjets.at(0).Charge()<0?-1:1) + (bjets.at(1).Charge()<0?-1:1), map_weight, 8,-4,4);
  FillHist(prefix+hprefix+"bjetCharge0"+suffix, bcharge0, map_weight, 200,-5,5);
  FillHist(prefix+hprefix+"bjetCharge1"+suffix, bcharge1, map_weight, 200,-5,5);
  FillHist(prefix+hprefix+"bjetsChargeSum"+suffix, bcharge0 + bcharge1, map_weight, 400,-10,10);
  FillHist(prefix+hprefix+"bjetsChargeAbsSum"+suffix, (bcharge0<0?-1:1) + (bcharge1<0?-1:1), map_weight, 8,-4,4);

  FillHist(prefix+hprefix+"ajetcharge0"+suffix, ajets.at(0).Charge(), map_weight, 200,-2,2);
  FillHist(prefix+hprefix+"ajetcharge1"+suffix, ajets.at(1).Charge(), map_weight, 200,-2,2);
  FillHist(prefix+hprefix+"ajetschargeSum"+suffix, ajets.at(0).Charge() + ajets.at(1).Charge(), map_weight, 400,-4,4);
  FillHist(prefix+hprefix+"ajetschargeAbsSum"+suffix, (ajets.at(0).Charge()<0?-1:1) + (ajets.at(1).Charge()<0?-1:1), map_weight, 8,-4,4);
  FillHist(prefix+hprefix+"ajetCharge0"+suffix, acharge0, map_weight, 200,-5,5);
  FillHist(prefix+hprefix+"ajetCharge1"+suffix, acharge1, map_weight, 200,-5,5);
  FillHist(prefix+hprefix+"ajetsChargeSum"+suffix, acharge0 + acharge1, map_weight, 400,-10,10);
  FillHist(prefix+hprefix+"ajetsChargeAbsSum"+suffix, (acharge0<0?-1:1) + (acharge1<0?-1:1), map_weight, 8,-4,4);

  FillHist(prefix+hprefix+"jetschargeSum"+suffix, bjets.at(0).Charge() + bjets.at(1).Charge() + ajets.at(0).Charge() + ajets.at(1).Charge(), map_weight, 800,-8,8);
  FillHist(prefix+hprefix+"jetschargeAbsSum"+suffix, (bjets.at(0).Charge()<0?-1:1) + (bjets.at(1).Charge()<0?-1:1) + (ajets.at(0).Charge()<0?-1:1) + (ajets.at(1).Charge()<0?-1:1), map_weight, 16,-8,8);
  FillHist(prefix+hprefix+"jetsChargeSum"+suffix, bcharge0 + bcharge1 + acharge0 + acharge1, map_weight, 800,-20,20);
  FillHist(prefix+hprefix+"jetsChargeAbsSum"+suffix, (bcharge0<0?-1:1) + (bcharge1<0?-1:1) + (acharge0<0?-1:1) + (acharge1<0?-1:1), map_weight, 16,-8,8);
}

bool ExampleRun_kinFitter::Hasleptons(TString channel){
  bool moreleptons = false;
  double l0pt = 26.;
  if(channel.Contains("m")){
    if(DataYear == 2017) l0pt = 29.;
    muons = MuonMomentumCorrection(SMPGetMuons("POGMediumWithLooseTrkIso", 8.0,2.4), 0,0,0);
    if(muons.size() > 0) lepton0 = &muons.at(0);
    if(muons.size() > 1) moreleptons = true;
  }else if(channel.Contains("e")){
    l0pt = 30.;
    if(DataYear > 2016) l0pt = 35.;
    electrons = ElectronEnergyCorrection(SMPGetElectrons("passMediumID", 8.0,2.5), 0,0);
    if(electrons.size() > 0) lepton0 = &electrons.at(0);
    if(electrons.size() > 1) moreleptons = true;
  }else if(channel.Contains("E")){
    l0pt = 30.;
    if(DataYear > 2016) l0pt = 35.;
    electrons = ElectronEnergyCorrection(SMPGetElectrons("passMediumID_SelQ", 8.0,2.5), 0,0);
    if(electrons.size() > 0) lepton0 = &electrons.at(0);
    if(electrons.size() > 1) moreleptons = true;
  }

  if(!lepton0) return false;
  FillCutflow(prefix+hprefix+"cutflow"+suffix, "1Leptons", map_weight[""]);
  if(lepton0->Pt() < l0pt) return false;
  FillCutflow(prefix+hprefix+"cutflow"+suffix, "LepPt", map_weight[""]);
  if(moreleptons) return false;
  FillCutflow(prefix+hprefix+"cutflow"+suffix, "1Lepton", map_weight[""]);
  met = GetEvent().GetMETVector();
  if(met.Pt() < 20.) return false;
  FillCutflow(prefix+hprefix+"cutflow"+suffix, "MET20", map_weight[""]);

  leptons ={};
  leptons.push_back(lepton0);

  return true;
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

ExampleRun_kinFitter::ExampleRun_kinFitter(){}
ExampleRun_kinFitter::~ExampleRun_kinFitter(){
  //==== Destructor of this Analyzer
  delete fitter;
}
