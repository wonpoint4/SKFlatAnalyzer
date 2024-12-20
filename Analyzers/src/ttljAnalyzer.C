#include "ttljAnalyzer.h"

void ttljAnalyzer::initializeAnalyzer(){
  SMPAnalyzerCore::initializeAnalyzer(); //setup eff zpt roc z0
  std::vector<JetTagging::Parameters> jtps = {
    JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::comb),
    JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Tight, JetTagging::incl, JetTagging::comb)
  };
  mcCorr->SetJetTaggingParameters(jtps);
  SetupPUJetWeight();
  SetupLikelihoods(0, 0);
}

void ttljAnalyzer::executeEvent(){
  ///////////////// GEN level /////////////////////
  if(!IsDATA && MCSample.Contains("TTLJ")) executeEventGen();

  ///////////////// RECO level /////////////////////
  if(!IsDATA || DataStream.Contains("SingleMuon")){
    executeEventWithParameter("m"+GetEraShort());
    executeEventWithParameter("m"+GetEraShort()+"_2mb");
    executeEventWithParameter("m"+GetEraShort()+"_3b");
    executeEventWithParameter("m"+GetEraShort()+"_PUjet");
    executeEventWithParameter("m"+GetEraShort()+"_LR_mode1"); // likelihood_mass_bl
    executeEventWithParameter("m"+GetEraShort()+"_LR_mode2"); // likelihood_mass_bl, mass_bjj
    executeEventWithParameter("m"+GetEraShort()+"_LR_mode3"); // likelihood_mass_bl, mass_bjj, mass_jj
    executeEventWithParameter("m"+GetEraShort()+"_LR_mode4"); // likelihood_mass_bl, mass_bjj, mass_blMET
    executeEventWithParameter("m"+GetEraShort()+"_LR_mode5"); // likelihood_mass_bl, mass_bjj, mass_jj, mass_blMET, dPhi_tt
    executeEventWithParameter("m"+GetEraShort()+"_LR_mode6"); // likelihood_mass_bl, mass_bjj, mass_jj, mass_blMET, dPhi_tt, dR_tt
    //executeEventWithParameter("m"+GetEraShort()+"_LR_mode7");
    //executeEventWithParameter("m"+GetEraShort()+"_LR_mode8");
    //executeEventWithParameter("m"+GetEraShort()+"_LR_mode9");
  }
  if(!IsDATA || DataStream.Contains("SingleElectron") || DataStream.Contains("EGamma")){
    executeEventWithParameter("e"+GetEraShort());
    executeEventWithParameter("E"+GetEraShort());
    executeEventWithParameter("E"+GetEraShort()+"_2mb");
    executeEventWithParameter("E"+GetEraShort()+"_3b");
    executeEventWithParameter("E"+GetEraShort()+"_PUjet");
    executeEventWithParameter("E"+GetEraShort()+"_LR_mode1"); // likelihood_mass_bl
    executeEventWithParameter("E"+GetEraShort()+"_LR_mode2"); // likelihood_mass_bl, mass_bjj
    executeEventWithParameter("E"+GetEraShort()+"_LR_mode3"); // likelihood_mass_bl, mass_bjj, mass_jj
    executeEventWithParameter("E"+GetEraShort()+"_LR_mode4"); // likelihood_mass_bl, mass_bjj, mass_blMET
    executeEventWithParameter("E"+GetEraShort()+"_LR_mode5"); // likelihood_mass_bl, mass_bjj, mass_jj, mass_blMET, dPhi_tt
    executeEventWithParameter("E"+GetEraShort()+"_LR_mode6"); // likelihood_mass_bl, mass_bjj, mass_jj, mass_blMET, dPhi_tt, dR_tt
    //executeEventWithParameter("E"+GetEraShort()+"_LR_mode7");
    //executeEventWithParameter("E"+GetEraShort()+"_LR_mode8");
    //executeEventWithParameter("E"+GetEraShort()+"_LR_mode9");
  }
}

void ttljAnalyzer::executeEventWithParameter(TString channel){

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
    if(channel.Contains("PUjet") && !PUJetIDPass(jet, "Loose")) continue;
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

  if(!channel.Contains("b")){
    if(bjets.size() != 2) return;
    FillCutflow(prefix+hprefix+"cutflow"+suffix, "Tight2b", map_weight[""]);
  }else if(channel.Contains("2mb")){
    if(bjets.size() < 2) return;
    FillCutflow(prefix+hprefix+"cutflow"+suffix, "Tight2mb", map_weight[""]);
  }else if(channel.Contains("3b")){
    if(bjets.size() < 3) return;
    FillCutflow(prefix+hprefix+"cutflow"+suffix, "Tight3b", map_weight[""]);
  }
  if(ajets.size() < 2) return;
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
    if(channel.Contains("m"+GetEraShort())){
      leptonTrackingSF *= fEff->GetEfficiencySF("Muon_Tracking", lepton0, 0,0);
      leptonRECOSF *= fEff->GetEfficiencySF("Muon_RECO", lepton0, 0,0);
      leptonIDSF *= fEff->GetEfficiencySF("Muon_MediumID_trkIsoLoose", lepton0, 0,0);
      if(DataYear == 2017) trigSFkey = "IsoMu27_MediumID_trkIsoLoose";
      leptonTriggerSF *= GetLeptonTriggerSF(trigSFkey, leptons, 0,0);
    }else if(channel.Contains("e"+GetEraShort())){
      leptonRECOSF *= fEff->GetEfficiencySF("Electron_RECO", lepton0, 0,0);
      leptonIDSF *= fEff->GetEfficiencySF("Electron_MediumID", lepton0, 0,0);
      trigSFkey = "Ele27_MediumID";
      if(DataYear == 2017) trigSFkey = "Ele32_MediumID";
      leptonTriggerSF *= GetLeptonTriggerSF(trigSFkey, leptons, 0,0);
    }else if(channel.Contains("E"+GetEraShort())){
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

  map_weight["_nobtagSF"] = map_weight[""];
  map_weight[""] *= btagSF;
  FillHist(prefix+hprefix+"weight_btagSF", btagSF, map_weight[""], 200,-5,5);
  FillCutflow(prefix+hprefix+"cutflow"+suffix, "btagSF", map_weight[""]);

  map_weight["_PUjetSF"] = map_weight[""] * pujetSF;
  FillHist(prefix+hprefix+"weight_PUjetSF", pujetSF, map_weight[""], 200,-5,5);

  map_weight["_bChargeSF0"] = map_weight[""] * GetbChargeSFWeight(bjets, 0, 0);
  FillHist(prefix+hprefix+"weight_bChargeSF0", GetbChargeSFWeight(bjets, 0, 0), map_weight[""], 200,-5,5);
  map_weight["_bChargeSF1"] = map_weight[""] * GetbChargeSFWeight(bjets, 1, 0);
  FillHist(prefix+hprefix+"weight_bChargeSF1", GetbChargeSFWeight(bjets, 1, 0), map_weight[""], 200,-5,5);

  //==== Making Likelihood
  if(!IsDATA && MCSample.Contains("TTLJ")) FillingLikelihood(channel, bjets, ajets, map_weight[""], 0, 0); // ByungHun Oh's method - drop events with ambiguity
  if(!IsDATA && MCSample.Contains("TTLJ")) FillingLikelihood(channel, bjets, ajets, map_weight[""], 1, 0); // Charmonium guy's method - match smaller dR < 0.3

  bool goodKinematic = false;
  goodKinematic = Kinematic_Cut(bjets, ajets);
  FillHist(prefix+"LR_efficiency", goodKinematic?1:0, map_weight[""], 2,0,2);
  if(!goodKinematic) return;
  FillCutflow(prefix+hprefix+"cutflow"+suffix, "LR", map_weight[""]);

  //==== Finding the correct bbjj combination
  vector<unsigned int> idx_bbjj = {0, 0, 0, 0};
  vector<double> LRs = {0, 0, 0, 0, 0, 0};
  if(channel.Contains("LR_mode1")) idx_bbjj = Finding_bbjj_byLikelihood(channel, bjets, ajets, LRs, 1);
  else if(channel.Contains("LR_mode2")) idx_bbjj = Finding_bbjj_byLikelihood(channel, bjets, ajets, LRs, 2);
  else if(channel.Contains("LR_mode3")) idx_bbjj = Finding_bbjj_byLikelihood(channel, bjets, ajets, LRs, 3);
  else if(channel.Contains("LR_mode4")) idx_bbjj = Finding_bbjj_byLikelihood(channel, bjets, ajets, LRs, 4);
  else if(channel.Contains("LR_mode5")) idx_bbjj = Finding_bbjj_byLikelihood(channel, bjets, ajets, LRs, 5);
  else if(channel.Contains("LR_mode6")) idx_bbjj = Finding_bbjj_byLikelihood(channel, bjets, ajets, LRs, 6);
  else if(channel.Contains("LR_mode7")) idx_bbjj = Finding_bbjj_byLikelihood(channel, bjets, ajets, LRs, 7);
  else if(channel.Contains("LR_mode8")) idx_bbjj = Finding_bbjj_byLikelihood(channel, bjets, ajets, LRs, 8);
  else if(channel.Contains("LR_mode9")) idx_bbjj = Finding_bbjj_byLikelihood(channel, bjets, ajets, LRs, 9);
  else idx_bbjj = Finding_bbjj_byLikelihood(channel, bjets, ajets, LRs);

  Jet lepb = bjets.at(idx_bbjj.at(0));
  Jet hadb = bjets.at(idx_bbjj.at(1));
  Jet Wj0 = ajets.at(idx_bbjj.at(2));
  Jet Wj1 = ajets.at(idx_bbjj.at(3));
  double lepb_charge = jetCharge(lepb);
  double hadb_charge = jetCharge(hadb);

  if(!IsDATA && MCSample.Contains("TTLJ")){
    double match_dR = 0.4;
    bool Gen_LR_match_onlyb = false;
    bool Gen_LR_match_Whad = false;
    bool Gen_LR_match_full = false;
    // chargeEasy => 0:negative, 1:positive
    // purity => -1:UnMatched, 0:Wrong, 1:Correct in TTLJ
    // correct => 0:Wrong, 1:Correct in TTLJ

    if((gen_j0.DeltaR(Wj0) < match_dR && gen_j1.DeltaR(Wj1) < match_dR) || (gen_j0.DeltaR(Wj1) < match_dR && gen_j1.DeltaR(Wj0) < match_dR)) Gen_LR_match_Whad = true;
    //When lepb = b, hadb = bbar -> lep+
    if(gen_b0.DeltaR(lepb) < match_dR && gen_b1.DeltaR(hadb) < match_dR){
      Gen_LR_match_onlyb = true;
      if(Gen_LR_match_Whad) Gen_LR_match_full = true;
      if(lepton0->Charge() > 0){
        FillHist(prefix+"LR_purity", 1, map_weight[""], 4,-2,2);
        FillHist(prefix+"LR_correct", 1, map_weight[""], 2,0,2);
        FillHist(prefix+"LR_purity_Lp", 1, map_weight[""], 4,-2,2);
        FillHist(prefix+"LR_correct_Lp", 1, map_weight[""], 2,0,2);
        prefix += "Correct_"; // lep+
      }else{
        FillHist(prefix+"LR_purity", 0, map_weight[""], 4,-2,2);
        FillHist(prefix+"LR_correct", 0, map_weight[""], 2,0,2);
        FillHist(prefix+"LR_purity_Lm", 0, map_weight[""], 4,-2,2);
        FillHist(prefix+"LR_correct_Lm", 0, map_weight[""], 2,0,2);
        prefix +="Wrong_"; // lep-
      }
      FillHist(prefix+"lepbjetCharge_Matched_b", lepb_charge, map_weight[""], 200,-5,5);
      FillHist(prefix+"lepbjetChargeEasy_Matched_b", lepb_charge<0?0:1, map_weight[""], 2,0,2);
      FillHist(prefix+"hadbjetCharge_Matched_bbar", hadb_charge, map_weight[""], 200,-5,5);
      FillHist(prefix+"hadbjetChargeEasy_Matched_bbar", hadb_charge<0?0:1, map_weight[""], 2,0,2);
      FillHist(prefix+"lcharge_lepb_Matched_b", lepton0->Charge(), map_weight[""], 4,-2,2);
      FillHist(prefix+"lchargeEasy_lepb_Matched_b", lepton0->Charge()<0?0:1, map_weight[""], 2,0,2);
      if(lepb_charge < 0){
        FillHist(prefix+"lcharge_lepb_negaive_Matched_b", lepton0->Charge(), map_weight[""], 4,-2,2);
        FillHist(prefix+"lchargeEasy_lepb_negaive_Matched_b", lepton0->Charge()<0?0:1, map_weight[""], 2,0,2);
      }
      if(0 < hadb_charge){
        FillHist(prefix+"lcharge_hadb_positive_Matched_bbar", lepton0->Charge(), map_weight[""], 4,-2,2);
        FillHist(prefix+"lchargeEasy_hadb_positive_Matched_bbar", lepton0->Charge()<0?0:1, map_weight[""], 2,0,2);
      }
    }
    //When lepb = bbar, hadb = b => lep-
    else if(gen_b0.DeltaR(hadb) < match_dR && gen_b1.DeltaR(lepb) < match_dR){
      Gen_LR_match_onlyb = true;
      if(Gen_LR_match_Whad) Gen_LR_match_full = true;
      if(lepton0->Charge() < 0){
        FillHist(prefix+"LR_purity", 1, map_weight[""], 4,-2,2);
        FillHist(prefix+"LR_correct", 1, map_weight[""], 2,0,2);
        FillHist(prefix+"LR_purity_Lm", 1, map_weight[""], 4,-2,2);
        FillHist(prefix+"LR_correct_Lm", 1, map_weight[""], 2,0,2);
        prefix += "Correct_"; // lep-
      }else{
        FillHist(prefix+"LR_purity", 0, map_weight[""], 4,-2,2);
        FillHist(prefix+"LR_correct", 0, map_weight[""], 2,0,2);
        FillHist(prefix+"LR_purity_Lp", 0, map_weight[""], 4,-2,2);
        FillHist(prefix+"LR_correct_Lp", 0, map_weight[""], 2,0,2);
        prefix +="Wrong_"; // lep+
      }
      FillHist(prefix+"hadbjetCharge_Matched_b", hadb_charge, map_weight[""], 200,-5,5);
      FillHist(prefix+"hadbjetChargeEasy_Matched_b", hadb_charge<0?0:1, map_weight[""], 2,0,2);
      FillHist(prefix+"lepbjetCharge_Matched_bbar", lepb_charge, map_weight[""], 200,-5,5);
      FillHist(prefix+"lepbjetChargeEasy_Matched_bbar", lepb_charge<0?0:1, map_weight[""], 2,0,2);
      FillHist(prefix+"lcharge_lepb_Matched_bbar", lepton0->Charge(), map_weight[""], 4,-2,2);
      FillHist(prefix+"lchargeEasy_lepb_Matched_bbar", lepton0->Charge()<0?0:1, map_weight[""], 2,0,2);
      if(0 < lepb_charge){
        FillHist(prefix+"lcharge_lepb_positive_Matched_bbar", lepton0->Charge(), map_weight[""], 4,-2,2);
        FillHist(prefix+"lchargeEasy_lepb_positive_Matched_bbar", lepton0->Charge()<0?0:1, map_weight[""], 2,0,2);
      }
      if(hadb_charge < 0){
        FillHist(prefix+"lcharge_hadb_negative_Matched_b", lepton0->Charge(), map_weight[""], 4,-2,2);
        FillHist(prefix+"lchargeEasy_hadb_negative_Matched_b", lepton0->Charge()<0?0:1, map_weight[""], 2,0,2);
      }
    }
    else{
      FillHist(prefix+"LR_purity", -1, map_weight[""], 4,-2,2);
      if(lepton0->Charge() > 0) FillHist(prefix+"LR_purity_Lp", -1, map_weight[""], 4,-2,2);
      if(lepton0->Charge() < 0) FillHist(prefix+"LR_purity_Lm", -1, map_weight[""], 4,-2,2);
      prefix +="UnMatched_";
    }

    FillHist("Gen_LR_match_onlyb", Gen_LR_match_onlyb, map_weight[""], 2,0,2);
    FillHist("Gen_LR_match_Whad", Gen_LR_match_Whad, map_weight[""], 2,0,2);
    FillHist("Gen_LR_match_full", Gen_LR_match_full, map_weight[""], 2,0,2);
    FillHist(prefix+"Gen_LR_match_onlyb", Gen_LR_match_onlyb, map_weight[""], 2,0,2);
    FillHist(prefix+"Gen_LR_match_Whad", Gen_LR_match_Whad, map_weight[""], 2,0,2);
    FillHist(prefix+"Gen_LR_match_full", Gen_LR_match_full, map_weight[""], 2,0,2);
    FillHist(channel+"/Gen_LR_match_onlyb", Gen_LR_match_onlyb, map_weight[""], 2,0,2);
    FillHist(channel+"/Gen_LR_match_Whad", Gen_LR_match_Whad, map_weight[""], 2,0,2);
    FillHist(channel+"/Gen_LR_match_full", Gen_LR_match_full, map_weight[""], 2,0,2);
  }

  //==========================
  //==== Now reco fill histograms
  //==========================

  FillHist(prefix+"mass_jj", (Wj0 + Wj1).M(), map_weight, 40,0,200);
  FillHist(prefix+"mass_bjj", (hadb + Wj0 + Wj1).M(), map_weight, 40,100,300);
  //FillHist(prefix+"mass_Toplep", (lepb + *lepton0 + neutrino).M(), map_weight, 40,100,300);
  FillHist(prefix+"mass_blMET", (lepb + *lepton0 + met).M(), map_weight, 40,100,300);
  FillHist(prefix+"mass_bl", (lepb + *lepton0).M(), map_weight, 40,0,200);

  FillHist(prefix+"dR_bjj_blMET", (hadb + Wj0 + Wj1).DeltaR(lepb + *lepton0 + met), map_weight, 200,0,10);
  FillHist(prefix+"dR_bjj_lMET", (hadb + Wj0 + Wj1).DeltaR(*lepton0 + met), map_weight, 200,0,10);
  FillHist(prefix+"dR_bjj_jj", (hadb + Wj0 + Wj1).DeltaR(Wj0 + Wj1), map_weight, 200,0,10);
  FillHist(prefix+"dR_blMET_lMET", (lepb + *lepton0 + met).DeltaR(*lepton0 + met), map_weight, 200,0,10);
  FillHist(prefix+"dR_blMET_jj", (lepb + *lepton0 + met).DeltaR(Wj0 + Wj1), map_weight, 200,0,10);
  FillHist(prefix+"dR_lMET_jj", (*lepton0 + met).DeltaR(Wj0 + Wj1), map_weight, 200,0,10);

  FillHist(prefix+"dPhi_bjj_blMET", fabs((hadb + Wj0 + Wj1).DeltaPhi(lepb + *lepton0 + met)), map_weight, 200,0,10);
  FillHist(prefix+"dPhi_bjj_lMET", fabs((hadb + Wj0 + Wj1).DeltaPhi(*lepton0 + met)), map_weight, 200,0,10);
  FillHist(prefix+"dPhi_bjj_jj", fabs((hadb + Wj0 + Wj1).DeltaPhi(Wj0 + Wj1)), map_weight, 200,0,10);
  FillHist(prefix+"dPhi_blMET_lMET", fabs((lepb + *lepton0 + met).DeltaPhi(*lepton0 + met)), map_weight, 200,0,10);
  FillHist(prefix+"dPhi_blMET_jj", fabs((lepb + *lepton0 + met).DeltaPhi(Wj0 + Wj1)), map_weight, 200,0,10);
  FillHist(prefix+"dPhi_lMET_jj", fabs((*lepton0 + met).DeltaPhi(Wj0 + Wj1)), map_weight, 200,0,10);

  //FillHist(prefix+"likelihood_ratio", LRs.at(0) * LRs.at(1) * LRs.at(2) * LRs.at(3), map_weight, 100,0,1);
  FillHist(prefix+"likelihood_ratio", LRs.at(0) * LRs.at(1) * LRs.at(2) * LRs.at(3) * LRs.at(4) * LRs.at(5), map_weight, 100,0,1);
  FillHist(prefix+"likelihood_ratio_Mbl", LRs.at(0), map_weight, 100,0,1);
  FillHist(prefix+"likelihood_ratio_MblMET", LRs.at(1), map_weight, 100,0,1);
  FillHist(prefix+"likelihood_ratio_Mbjj", LRs.at(2), map_weight, 100,0,1);
  FillHist(prefix+"likelihood_ratio_Mjj", LRs.at(3), map_weight, 100,0,1);
  FillHist(prefix+"likelihood_ratio_dRtt", LRs.at(4), map_weight, 100,0,1);
  FillHist(prefix+"likelihood_ratio_dPhitt", LRs.at(5), map_weight, 100,0,1);

  FillHist(prefix+"lepbjetcharge", lepb.Charge(), map_weight, 200,-2,2);
  FillHist(prefix+"lepbjetchargeEasy", lepb.Charge()<0?0:1, map_weight, 2,0,2);
  FillHist(prefix+"hadbjetcharge", hadb.Charge(), map_weight, 200,-2,2);
  FillHist(prefix+"hadbjetchargeEasy", hadb.Charge()<0?0:1, map_weight, 2,0,2);
  FillHist(prefix+"lepbjetCharge", lepb_charge, map_weight, 200,-5,5);
  FillHist(prefix+"lepbjetChargeEasy", lepb_charge<0?0:1, map_weight, 2,0,2);
  FillHist(prefix+"hadbjetCharge", hadb_charge, map_weight, 200,-5,5);
  FillHist(prefix+"hadbjetChargeEasy", hadb_charge<0?0:1, map_weight, 2,0,2);
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
}

bool ttljAnalyzer::Hasleptons(TString channel){
  bool moreleptons = false;
  double l0pt = 26.;
  if(channel.Contains("m"+GetEraShort())){
    if(DataYear == 2017) l0pt = 29.;
    muons = MuonMomentumCorrection(SMPGetMuons("POGMediumWithLooseTrkIso", 8.0,2.4), 0,0,0);
    if(muons.size() > 0) lepton0 = &muons.at(0);
    if(muons.size() > 1) moreleptons = true;
  }else if(channel.Contains("e"+GetEraShort())){
    l0pt = 30.;
    if(DataYear > 2016) l0pt = 35.;
    electrons = ElectronEnergyCorrection(SMPGetElectrons("passMediumID", 8.0,2.5), 0,0);
    if(electrons.size() > 0) lepton0 = &electrons.at(0);
    if(electrons.size() > 1) moreleptons = true;
  }else if(channel.Contains("E"+GetEraShort())){
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

void ttljAnalyzer::executeEventGen(){
  FillHist("gen/executeEventGen", 1, 1, 2,0,2);
  vector<Gen> gens=GetGens();
  GetTTLJGenParticles(gens, gen_parton0,gen_parton1, gen_b0,gen_b1, gen_l0,gen_l1, gen_j0,gen_j1, 3);

  FillHist("gen/Idx_parton0", gen_parton0.Index(), 1, 100,0,100);
  FillHist("gen/Idx_parton1", gen_parton1.Index(), 1, 100,0,100);
  FillHist("gen/Idx_b0", gen_b0.Index(), 1, 100,0,100);
  FillHist("gen/Idx_b1", gen_b1.Index(), 1, 100,0,100);
  FillHist("gen/Idx_l0", gen_l0.Index(), 1, 100,0,100);
  FillHist("gen/Idx_l1", gen_l1.Index(), 1, 100,0,100);
  FillHist("gen/Idx_j0", gen_j0.Index(), 1, 100,0,100);
  FillHist("gen/Idx_j1", gen_j1.Index(), 1, 100,0,100);

  FillHist("gen/PID_parton0", gen_parton0.PID(), 1, 60,-30,30);
  FillHist("gen/PID_parton1", gen_parton1.PID(), 1, 60,-30,30);
  FillHist("gen/PID_b0", gen_b0.PID(), 1, 60,-30,30);
  FillHist("gen/PID_b1", gen_b1.PID(), 1, 60,-30,30);
  FillHist("gen/PID_l0", gen_l0.PID(), 1, 60,-30,30);
  FillHist("gen/PID_l1", gen_l1.PID(), 1, 60,-30,30);
  FillHist("gen/PID_j0", gen_j0.PID(), 1, 60,-30,30);
  FillHist("gen/PID_j1", gen_j1.PID(), 1, 60,-30,30);

  FillHist("gen/Pt_parton0", gen_parton0.Pt(), 1, 100,0,400);
  FillHist("gen/Pt_parton1", gen_parton1.Pt(), 1, 100,0,400);
  FillHist("gen/Pt_b0", gen_b0.Pt(), 1, 100,0,400);
  FillHist("gen/Pt_b1", gen_b1.Pt(), 1, 100,0,400);
  FillHist("gen/Pt_l0", gen_l0.Pt(), 1, 100,0,400);
  FillHist("gen/Pt_l1", gen_l1.Pt(), 1, 100,0,400);
  FillHist("gen/Pt_j0", gen_j0.Pt(), 1, 100,0,400);
  FillHist("gen/Pt_j1", gen_j1.Pt(), 1, 100,0,400);

  FillHist("gen/Eta_parton0", gen_parton0.Eta(), 1, 200,-5,5);
  FillHist("gen/Eta_parton1", gen_parton1.Eta(), 1, 200,-5,5);
  FillHist("gen/Eta_b0", gen_b0.Eta(), 1, 200,-5,5);
  FillHist("gen/Eta_b1", gen_b1.Eta(), 1, 200,-5,5);
  FillHist("gen/Eta_l0", gen_l0.Eta(), 1, 200,-5,5);
  FillHist("gen/Eta_l1", gen_l1.Eta(), 1, 200,-5,5);
  FillHist("gen/Eta_j0", gen_j0.Eta(), 1, 200,-5,5);
  FillHist("gen/Eta_j1", gen_j1.Eta(), 1, 200,-5,5);

  FillHist("gen/Phi_parton0", gen_parton0.Phi(), 1, 64,-3.2,3.2);
  FillHist("gen/Phi_parton1", gen_parton1.Phi(), 1, 64,-3.2,3.2);
  FillHist("gen/Phi_b0", gen_b0.Phi(), 1, 64,-3.2,3.2);
  FillHist("gen/Phi_b1", gen_b1.Phi(), 1, 64,-3.2,3.2);
  FillHist("gen/Phi_l0", gen_l0.Phi(), 1, 64,-3.2,3.2);
  FillHist("gen/Phi_l1", gen_l1.Phi(), 1, 64,-3.2,3.2);
  FillHist("gen/Phi_j0", gen_j0.Phi(), 1, 64,-3.2,3.2);
  FillHist("gen/Phi_j1", gen_j1.Phi(), 1, 64,-3.2,3.2);

  FillHist("gen/Mass_Whad", (gen_j0 + gen_j1).M(), 1, 50,0,200);
  FillHist("gen/Mass_Wlep", (gen_l0 + gen_l1).M(), 1, 50,0,200);
  FillHist("gen/Mass_Tophad1", (gen_b0 + gen_j0 + gen_j1).M(), 1, 75,0,300);
  FillHist("gen/Mass_Tophad2", (gen_b1 + gen_j0 + gen_j1).M(), 1, 75,0,300);
  FillHist("gen/Mass_Toplep1", (gen_b0 + gen_l0 + gen_l1).M(), 1, 75,0,300);
  FillHist("gen/Mass_Toplep2", (gen_b1 + gen_l0 + gen_l1).M(), 1, 75,0,300);
  if(gen_l0.PID() < 0){ // lepton is l+, thus lepb, l+ from t
    gen_lepb = &gen_b0;
    gen_hadb = &gen_b1;
  }else{
    gen_lepb = &gen_b1;
    gen_hadb = &gen_b0;
  }
  FillHist("gen/Mass_Toplep", (*gen_lepb + gen_l0 + gen_l1).M(), 1, 75,0,300);
  FillHist("gen/Mass_Tophad", (*gen_hadb + gen_j0 + gen_j1).M(), 1, 75,0,300);

  FillHist("gen/dR_TophadToplep", (*gen_hadb + gen_j0 + gen_j1).DeltaR(*gen_lepb + gen_l0 + gen_l1), 1, 100,0,10);
  FillHist("gen/dR_TophadWlep", (*gen_hadb + gen_j0 + gen_j1).DeltaR(gen_l0 + gen_l1), 1, 100,0,10);
  FillHist("gen/dR_TophadWhad", (*gen_hadb + gen_j0 + gen_j1).DeltaR(gen_j0 + gen_j1), 1, 100,0,10);
  FillHist("gen/dR_ToplepWlep", (*gen_lepb + gen_l0 + gen_l1).DeltaR(gen_l0 + gen_l1), 1, 100,0,10);
  FillHist("gen/dR_ToplepWhad", (*gen_lepb + gen_l0 + gen_l1).DeltaR(gen_j0 + gen_j1), 1, 100,0,10);
  FillHist("gen/dR_WlepWhad", (gen_l0 + gen_l1).DeltaR(gen_j0 + gen_j1), 1, 100,0,10);

  FillHist("gen/dPhi_TophadToplep", (*gen_hadb + gen_j0 + gen_j1).DeltaPhi(*gen_lepb + gen_l0 + gen_l1), 1, 64,-3.2,3.2);
  FillHist("gen/dPhi_TophadWlep", (*gen_hadb + gen_j0 + gen_j1).DeltaPhi(gen_l0 + gen_l1), 1, 64,-3.2,3.2);
  FillHist("gen/dPhi_TophadWhad", (*gen_hadb + gen_j0 + gen_j1).DeltaPhi(gen_j0 + gen_j1), 1, 64,-3.2,3.2);
  FillHist("gen/dPhi_ToplepWlep", (*gen_lepb + gen_l0 + gen_l1).DeltaPhi(gen_l0 + gen_l1), 1, 64,-3.2,3.2);
  FillHist("gen/dPhi_ToplepWhad", (*gen_lepb + gen_l0 + gen_l1).DeltaPhi(gen_j0 + gen_j1), 1, 64,-3.2,3.2);
  FillHist("gen/dPhi_WlepWhad", (gen_l0 + gen_l1).DeltaPhi(gen_j0 + gen_j1), 1, 64,-3.2,3.2);
}

void ttljAnalyzer::GetTTLJGenParticles(const vector<Gen>& gens, Gen& parton0, Gen& parton1, Gen& b0, Gen& b1, Gen& l0, Gen& l1, Gen& j0, Gen& j1, int mode){
  //mode 0:bare 1:dressed01 2:dressed04 3:beforeFSR
  vector<const Gen*> leptons;
  vector<const Gen*> photons;
  vector<const Gen*> jets;

  parton0=Gen();
  parton1=Gen();
  b0=Gen();
  b1=Gen();
  l0=Gen();
  l1=Gen();
  j0=Gen();
  j1=Gen();
  int ngen=gens.size();
  for(int i=0;i<ngen;i++){
    if(!gens.at(i).isPrompt()) continue;
    int genpid=gens.at(i).PID();
    if(gens.at(i).isHardProcess()){
      if(abs(genpid)<7||genpid==21||genpid==22){
        if(parton0.IsEmpty()) parton0=gens[i];
        else if(parton1.IsEmpty()) parton1=gens[i];
      }
    }
    if(gens.at(i).Status()==1){
      if(gens.at(i).PID()==22) photons.push_back(&gens[i]); //photon
    }
    if(!gens.at(i).isHardProcess()) continue;
    if((abs(genpid)>=11 && abs(genpid)<=18) && abs(gens.at(gens.at(i).MotherIndex()).PID()) == 24) leptons.push_back(&gens[i]); //leptons from W
    // b0 : b from t, b1 : b~ from t~
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

  // j0 : leading jet, j1 : subleading jet from hadronic W
  int njet=jets.size();
  for(int i=0;i<njet;i++){
    for(int j=i+1;j<njet;j++){
      if(!(abs(jets[i]->PID()+jets[j]->PID()) == 1 || abs(jets[i]->PID()+jets[j]->PID()) == 3)) continue;
      if((*jets[i]+*jets[j]).M()>(j0+j1).M()){
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

  if(mode>=3){
    if(nlepton>=3){
      for(int i=0;i<nlepton;i++){
        if(leptons[i]->Index()==l0.Index()||leptons[i]->Index()==l1.Index()) continue;
        for(int j=i+1;j<nlepton;j++){
          if(leptons[j]->Index()==l0.Index()||leptons[j]->Index()==l1.Index()) continue;
          if(!(leptons[i]->PID()+leptons[j]->PID()==0)) continue;
          vector<int> history_i=TrackGenSelfHistory(*leptons[i],gens);
          vector<int> history_j=TrackGenSelfHistory(*leptons[j],gens);
          if(history_i.at(1)==history_j.at(1)){
            photons.push_back(leptons[i]);
            photons.push_back(leptons[j]);
          }
        }
      }
    }
    for(const auto& photon:photons){
      vector<int> history=TrackGenSelfHistory(*photon,gens);
      if(gens[history.at(1)].PID()==l0.PID()) l0+=*photon;
      else if(gens[history.at(1)].PID()==23) l0+=*photon; // for minnlo+photos
    }
  }else if(mode>=1){
    double delr=mode==1?0.1:0.4;
    for(const auto& photon:photons){
      if(l0.DeltaR(*photon)>delr&&l1.DeltaR(*photon)>delr) continue;
      if(l0.DeltaR(*photon)<l1.DeltaR(*photon)) l0+=*photon;
    }
  }
}

void ttljAnalyzer::FillingLikelihood(TString channel, vector<Jet> bjets, vector<Jet> ajets, double weight, unsigned int mode1, unsigned int mode2){
  FillCutflow(Form(channel+"/cutflow_goodMatching_mode1_%d_ForLikelihood", mode1), "Preselection", weight);
  double match_dR = 0.4;
  if(gen_l0.DeltaR(*lepton0) > match_dR) return;
  if(mode1 == 1) match_dR = 0.3;
  FillCutflow(Form(channel+"/cutflow_goodMatching_mode1_%d_ForLikelihood", mode1), "1lep", weight);

  Jet* lepb = NULL;
  Jet* hadb = NULL;
  Jet* Wj0 = NULL;
  Jet* Wj1 = NULL;

  for(auto& bjet : bjets){
    if(gen_lepb->DeltaR(bjet) < match_dR){
      if(!lepb) lepb = &bjet;
      else if(mode1 == 1 && (gen_lepb->DeltaR(bjet) < gen_lepb->DeltaR(*lepb))) lepb = &bjet; // Do matching with mindR (when mode1 == 1) like Charmonium guys
      else if(mode1 == 0) return; // Not use events with ambiguity like ByungHun Oh
    }
    if(gen_hadb->DeltaR(bjet) < match_dR){
      if(!hadb) hadb = &bjet;
      else if(mode1 == 1 && (gen_hadb->DeltaR(bjet) < gen_hadb->DeltaR(*hadb))) hadb = &bjet; // Do matching with mindR (when mode1 == 1) like Charmonium guys
      else if(mode1 == 0) return; // Not use events with ambiguity like ByungHun Oh
    }
  }
  FillCutflow(Form(channel+"/cutflow_goodMatching_mode1_%d_ForLikelihood", mode1), "2b-1gen", weight);

  for(auto& jet : ajets){
    if(gen_j0.DeltaR(jet) < match_dR){
      if(!Wj0) Wj0 = &jet;
      else if(mode1 == 1 && (gen_j0.DeltaR(jet) < gen_j0.DeltaR(*Wj0))) Wj0 = &jet; // Do matching with mindR (when mode1 == 1) like Charmonium guys
      else if(mode1 == 0) return; // Not use events with ambiguity like ByungHun Oh
    }
    if(gen_j1.DeltaR(jet) < match_dR){
      if(!Wj1) Wj1 = &jet;
      else if(mode1 == 1 && (gen_j1.DeltaR(jet) < gen_j1.DeltaR(*Wj1))) Wj1 = &jet; // Do matching with mindR (when mode1 == 1) like Charmonium guys
      else if(mode1 == 0) return; // Not use events with ambiguity like ByungHun Oh
    }
  }
  FillCutflow(Form(channel+"/cutflow_goodMatching_mode1_%d_ForLikelihood", mode1), "2j-1gen", weight);

  if(!lepb) return;
  FillCutflow(Form(channel+"/cutflow_goodMatching_mode1_%d_ForLikelihood", mode1), "nolepb", weight);
  if(!hadb) return;
  FillCutflow(Form(channel+"/cutflow_goodMatching_mode1_%d_ForLikelihood", mode1), "nohadb", weight);
  if(!Wj0) return;
  FillCutflow(Form(channel+"/cutflow_goodMatching_mode1_%d_ForLikelihood", mode1), "noWj0", weight);
  if(!Wj1) return;
  FillCutflow(Form(channel+"/cutflow_goodMatching_mode1_%d_ForLikelihood", mode1), "noWj1", weight);

  if(!lepb || !hadb || !Wj0 || !Wj1) return;
  FillCutflow(Form(channel+"/cutflow_goodMatching_mode1_%d_ForLikelihood", mode1), "4jets", weight);

  if(lepb == hadb) return;
  FillCutflow(Form(channel+"/cutflow_goodMatching_mode1_%d_ForLikelihood", mode1), "1b-2gens", weight);
  if(Wj0 == Wj1) return;
  FillCutflow(Form(channel+"/cutflow_goodMatching_mode1_%d_ForLikelihood", mode1), "1j-2gens", weight);

  // Filling Likelihood Histograms
  // mode2 = 0 (Wrong : wrong && lepb,hadb switched)
  for(unsigned int a=0; a<bjets.size(); a++){
    for(unsigned int b=a+1; b<bjets.size(); b++){
      if(&bjets.at(a) == lepb && &bjets.at(b) == hadb){
        FillHist(Form(channel+"/likelihood_mode1_%d_Mbl_Correct_mode2_0", mode1), (bjets.at(a) + *lepton0).M(), weight, 50,0,200);
        FillHist(Form(channel+"/likelihood_mode1_%d_MblMET_Correct_mode2_0", mode1), (bjets.at(a) + *lepton0 + met).M(), weight, 50,0,300);
        FillHist(Form(channel+"/likelihood_mode1_%d_Mbl_Wrong_mode2_0", mode1), (bjets.at(b) + *lepton0).M(), weight, 50,0,200);
        FillHist(Form(channel+"/likelihood_mode1_%d_MblMET_Wrong_mode2_0", mode1), (bjets.at(b) + *lepton0 + met).M(), weight, 50,0,300);
        for(unsigned int c=0; c<ajets.size(); c++){
          for(unsigned int d=c+1; d<ajets.size(); d++){
            if((&ajets.at(c) == Wj0 && &ajets.at(d) == Wj1) || (&ajets.at(c) == Wj1 && &ajets.at(d) == Wj0)){
              FillHist(Form(channel+"/likelihood_mode1_%d_Mbjj_Correct_mode2_0", mode1), (bjets.at(b) + ajets.at(c) + ajets.at(d)).M(), weight, 50,0,300);
              FillHist(Form(channel+"/likelihood_mode1_%d_Mjj_Correct_mode2_0", mode1), (ajets.at(c) + ajets.at(d)).M(), weight, 50,0,200);
              FillHist(Form(channel+"/likelihood_mode1_%d_dRtt_Correct_mode2_0", mode1), (bjets.at(b) + ajets.at(c) + ajets.at(d)).DeltaR(bjets.at(a) + *lepton0 + met), weight, 60,0,6);
              FillHist(Form(channel+"/likelihood_mode1_%d_dPhitt_Correct_mode2_0", mode1), fabs((bjets.at(b) + ajets.at(c) + ajets.at(d)).DeltaPhi(bjets.at(a) + *lepton0 + met)), weight, 30,0,3);
              FillHist(Form(channel+"/likelihood_mode1_%d_Mbjj_Wrong_mode2_0", mode1), (bjets.at(a) + ajets.at(c) + ajets.at(d)).M(), weight, 50,0,300);
              FillHist(Form(channel+"/likelihood_mode1_%d_dRtt_Wrong_mode2_0", mode1), (bjets.at(a) + ajets.at(c) + ajets.at(d)).DeltaR(bjets.at(b) + *lepton0 + met), weight, 60,0,6);
              FillHist(Form(channel+"/likelihood_mode1_%d_dPhitt_Wrong_mode2_0", mode1), fabs((bjets.at(a) + ajets.at(c) + ajets.at(d)).DeltaPhi(bjets.at(b) + *lepton0 + met)), weight, 30,0,3);
            }
            else{
              FillHist(Form(channel+"/likelihood_mode1_%d_Mbjj_Wrong_mode2_0", mode1), (bjets.at(a) + ajets.at(c) + ajets.at(d)).M(), weight, 50,0,300);
              FillHist(Form(channel+"/likelihood_mode1_%d_Mjj_Wrong_mode2_0", mode1), (ajets.at(c) + ajets.at(d)).M(), weight, 50,0,200);
              FillHist(Form(channel+"/likelihood_mode1_%d_dRtt_Wrong_mode2_0", mode1), (bjets.at(a) + ajets.at(c) + ajets.at(d)).DeltaR(bjets.at(b) + *lepton0 + met), weight, 60,0,6);
              FillHist(Form(channel+"/likelihood_mode1_%d_dPhitt_Wrong_mode2_0", mode1), fabs((bjets.at(a) + ajets.at(c) + ajets.at(d)).DeltaPhi(bjets.at(b) + *lepton0 + met)), weight, 30,0,3);
            }
          }
        }
      }else if(&bjets.at(a) == hadb && &bjets.at(b) == lepb){
        FillHist(Form(channel+"/likelihood_mode1_%d_Mbl_Correct_mode2_0", mode1), (bjets.at(b) + *lepton0).M(), weight, 50,0,200);
        FillHist(Form(channel+"/likelihood_mode1_%d_MblMET_Correct_mode2_0", mode1), (bjets.at(b) + *lepton0 + met).M(), weight, 50,0,300);
        FillHist(Form(channel+"/likelihood_mode1_%d_Mbl_Wrong_mode2_0", mode1), (bjets.at(a) + *lepton0).M(), weight, 50,0,200);
        FillHist(Form(channel+"/likelihood_mode1_%d_MblMET_Wrong_mode2_0", mode1), (bjets.at(a) + *lepton0 + met).M(), weight, 50,0,300);
        for(unsigned int c=0; c<ajets.size(); c++){
          for(unsigned int d=c+1; d<ajets.size(); d++){
            if((&ajets.at(c) == Wj0 && &ajets.at(d) == Wj1) || (&ajets.at(c) == Wj1 && &ajets.at(d) == Wj0)){
              FillHist(Form(channel+"/likelihood_mode1_%d_Mbjj_Correct_mode2_0", mode1), (bjets.at(a) + ajets.at(c) + ajets.at(d)).M(), weight, 50,0,300);
              FillHist(Form(channel+"/likelihood_mode1_%d_Mjj_Correct_mode2_0", mode1), (ajets.at(c) + ajets.at(d)).M(), weight, 50,0,200);
              FillHist(Form(channel+"/likelihood_mode1_%d_dRtt_Correct_mode2_0", mode1), (bjets.at(a) + ajets.at(c) + ajets.at(d)).DeltaR(bjets.at(b) + *lepton0 + met), weight, 60,0,6);
              FillHist(Form(channel+"/likelihood_mode1_%d_dPhitt_Correct_mode2_0", mode1), fabs((bjets.at(a) + ajets.at(c) + ajets.at(d)).DeltaPhi(bjets.at(b) + *lepton0 + met)), weight, 30,0,3);
              FillHist(Form(channel+"/likelihood_mode1_%d_Mbjj_Wrong_mode2_0", mode1), (bjets.at(b) + ajets.at(c) + ajets.at(d)).M(), weight, 50,0,300);
              FillHist(Form(channel+"/likelihood_mode1_%d_dRtt_Wrong_mode2_0", mode1), (bjets.at(b) + ajets.at(c) + ajets.at(d)).DeltaR(bjets.at(a) + *lepton0 + met), weight, 60,0,6);
              FillHist(Form(channel+"/likelihood_mode1_%d_dPhitt_Wrong_mode2_0", mode1), fabs((bjets.at(b) + ajets.at(c) + ajets.at(d)).DeltaPhi(bjets.at(a) + *lepton0 + met)), weight, 30,0,3);
            }
            else{
              FillHist(Form(channel+"/likelihood_mode1_%d_Mbjj_Wrong_mode2_0", mode1), (bjets.at(b) + ajets.at(c) + ajets.at(d)).M(), weight, 50,0,300);
              FillHist(Form(channel+"/likelihood_mode1_%d_Mjj_Wrong_mode2_0", mode1), (ajets.at(c) + ajets.at(d)).M(), weight, 50,0,200);
              FillHist(Form(channel+"/likelihood_mode1_%d_dRtt_Wrong_mode2_0", mode1), (bjets.at(b) + ajets.at(c) + ajets.at(d)).DeltaR(bjets.at(a) + *lepton0 + met), weight, 60,0,6);
              FillHist(Form(channel+"/likelihood_mode1_%d_dPhitt_Wrong_mode2_0", mode1), fabs((bjets.at(b) + ajets.at(c) + ajets.at(d)).DeltaPhi(bjets.at(a) + *lepton0 + met)), weight, 30,0,3);
            }
          }
        }
      }
    }

    // mode2 = 1 (Wrong : simply wrong)
    if(&bjets.at(a) == lepb){
      FillHist(Form(channel+"/likelihood_mode1_%d_Mbl_Correct_mode2_1", mode1), (bjets.at(a) + *lepton0).M(), weight, 50,0,200);
      FillHist(Form(channel+"/likelihood_mode1_%d_MblMET_Correct_mode2_1", mode1), (bjets.at(a) + *lepton0 + met).M(), weight, 50,0,300);
    }else{
      FillHist(Form(channel+"/likelihood_mode1_%d_Mbl_Wrong_mode2_1", mode1), (bjets.at(a) + *lepton0).M(), weight, 50,0,200);
      FillHist(Form(channel+"/likelihood_mode1_%d_MblMET_Wrong_mode2_1", mode1), (bjets.at(a) + *lepton0 + met).M(), weight, 50,0,300);
    }

    if(&bjets.at(a) == hadb){
      for(unsigned int c=0; c<ajets.size(); c++){
        for(unsigned int d=c+1; d<ajets.size(); d++){
          if((&ajets.at(c) == Wj0 && &ajets.at(d) == Wj1) || (&ajets.at(c) == Wj1 && &ajets.at(d) == Wj0)){
            FillHist(Form(channel+"/likelihood_mode1_%d_Mbjj_Correct_mode2_1", mode1), (bjets.at(a) + ajets.at(c) + ajets.at(d)).M(), weight, 50,0,300);
            FillHist(Form(channel+"/likelihood_mode1_%d_Mjj_Correct_mode2_1", mode1), (ajets.at(c) + ajets.at(d)).M(), weight, 50,0,200);
          }
          else{
            FillHist(Form(channel+"/likelihood_mode1_%d_Mbjj_Wrong_mode2_1", mode1), (bjets.at(a) + ajets.at(c) + ajets.at(d)).M(), weight, 50,0,300);
            FillHist(Form(channel+"/likelihood_mode1_%d_Mjj_Wrong_mode2_1", mode1), (ajets.at(c) + ajets.at(d)).M(), weight, 50,0,200);
          }
        }
      }
      for(unsigned int b=a+1; b<bjets.size(); b++){
        for(unsigned int c=0; c<ajets.size(); c++){
          for(unsigned int d=c+1; d<ajets.size(); d++){
            if(&bjets.at(b) == lepb){
              if((&ajets.at(c) == Wj0 && &ajets.at(d) == Wj1) || (&ajets.at(c) == Wj1 && &ajets.at(d) == Wj0)){
                FillHist(Form(channel+"/likelihood_mode1_%d_dRtt_Correct_mode2_1", mode1), (bjets.at(a) + ajets.at(c) + ajets.at(d)).DeltaR(bjets.at(b) + *lepton0 + met), weight, 60,0,6);
                FillHist(Form(channel+"/likelihood_mode1_%d_dPhitt_Correct_mode2_1", mode1), fabs((bjets.at(a) + ajets.at(c) + ajets.at(d)).DeltaPhi(bjets.at(b) + *lepton0 + met)), weight, 30,0,3);
              }
            }else{
              FillHist(Form(channel+"/likelihood_mode1_%d_dRtt_Wrong_mode2_1", mode1), (bjets.at(a) + ajets.at(c) + ajets.at(d)).DeltaR(bjets.at(b) + *lepton0 + met), weight, 60,0,6);
              FillHist(Form(channel+"/likelihood_mode1_%d_dPhitt_Wrong_mode2_1", mode1), fabs((bjets.at(a) + ajets.at(c) + ajets.at(d)).DeltaPhi(bjets.at(b) + *lepton0 + met)), weight, 30,0,3);
            }
          }
        }
      }
    }else{
      for(unsigned int c=0; c<ajets.size(); c++){
        for(unsigned int d=c+1; d<ajets.size(); d++){
          if((&ajets.at(c) == Wj0 && &ajets.at(d) == Wj1) || (&ajets.at(c) == Wj1 && &ajets.at(d) == Wj0)){
            FillHist(Form(channel+"/likelihood_mode1_%d_Mbjj_Wrong_mode2_1", mode1), (bjets.at(a) + ajets.at(c) + ajets.at(d)).M(), weight, 50,0,300);
            FillHist(Form(channel+"/likelihood_mode1_%d_Mjj_Correct_mode2_1", mode1), (ajets.at(c) + ajets.at(d)).M(), weight, 50,0,200);
          }
          else{
            FillHist(Form(channel+"/likelihood_mode1_%d_Mbjj_Wrong_mode2_1", mode1), (bjets.at(a) + ajets.at(c) + ajets.at(d)).M(), weight, 50,0,300);
            FillHist(Form(channel+"/likelihood_mode1_%d_Mjj_Wrong_mode2_1", mode1), (ajets.at(c) + ajets.at(d)).M(), weight, 50,0,200);
          }
        }
      }
      for(unsigned int b=a+1; b<bjets.size(); b++){
        for(unsigned int c=0; c<ajets.size(); c++){
          for(unsigned int d=c+1; d<ajets.size(); d++){
            FillHist(Form(channel+"/likelihood_mode1_%d_dRtt_Wrong_mode2_1", mode1), (bjets.at(a) + ajets.at(c) + ajets.at(d)).DeltaR(bjets.at(b) + *lepton0 + met), weight, 60,0,6);
            FillHist(Form(channel+"/likelihood_mode1_%d_dPhitt_Wrong_mode2_1", mode1), fabs((bjets.at(a) + ajets.at(c) + ajets.at(d)).DeltaPhi(bjets.at(b) + *lepton0 + met)), weight, 30,0,3);
          }
        }
      }
    }
  }
}

void ttljAnalyzer::SetupLikelihoods(unsigned int mode1, unsigned int mode2){
  TString path = getenv("DATA_DIR")+TString("/")+GetEra()+"/JME/Likelihoods.root";
  if(IsExists(path)){
    cout<<"[ttljAnalyzer::SetupLikelihoods] using file "+path<<" for mode1 = "<<mode1<<", mode2 = "<<mode2<<endl;
  }else{
    cout<<"[ttljAnalyzer::SetupLikelihoods] no "+path<<endl;
    return;
  }

  TFile f(path);
  hMbl_correct_E = (TH1*)f.Get(Form("E"+GetEra()+"/likelihood_mode1_%d_Mbl_Correct_mode2_%d", mode1, mode2));
  hMbl_wrong_E = (TH1*)f.Get(Form("E"+GetEra()+"/likelihood_mode1_%d_Mbl_Wrong_mode2_%d", mode1, mode2));
  hMblMET_correct_E = (TH1*)f.Get(Form("E"+GetEra()+"/likelihood_mode1_%d_MblMET_Correct_mode2_%d", mode1, mode2));
  hMblMET_wrong_E = (TH1*)f.Get(Form("E"+GetEra()+"/likelihood_mode1_%d_MblMET_Wrong_mode2_%d", mode1, mode2));
  hMbjj_correct_E = (TH1*)f.Get(Form("E"+GetEra()+"/likelihood_mode1_%d_Mbjj_Correct_mode2_%d", mode1, mode2));
  hMbjj_wrong_E = (TH1*)f.Get(Form("E"+GetEra()+"/likelihood_mode1_%d_Mbjj_Wrong_mode2_%d", mode1, mode2));
  hMjj_correct_E = (TH1*)f.Get(Form("E"+GetEra()+"/likelihood_mode1_%d_Mjj_Correct_mode2_%d", mode1, mode2));
  hMjj_wrong_E = (TH1*)f.Get(Form("E"+GetEra()+"/likelihood_mode1_%d_Mjj_Wrong_mode2_%d", mode1, mode2));
  hdRtt_correct_E = (TH1*)f.Get(Form("E"+GetEra()+"/likelihood_mode1_%d_dRtt_Correct_mode2_%d", mode1, mode2));
  hdRtt_wrong_E = (TH1*)f.Get(Form("E"+GetEra()+"/likelihood_mode1_%d_dRtt_Wrong_mode2_%d", mode1, mode2));
  hdPhitt_correct_E = (TH1*)f.Get(Form("E"+GetEra()+"/likelihood_mode1_%d_dPhitt_Correct_mode2_%d", mode1, mode2));
  hdPhitt_wrong_E = (TH1*)f.Get(Form("E"+GetEra()+"/likelihood_mode1_%d_dPhitt_Wrong_mode2_%d", mode1, mode2));

  hMbl_correct_m = (TH1*)f.Get(Form("m"+GetEra()+"/likelihood_mode1_%d_Mbl_Correct_mode2_%d", mode1, mode2));
  hMbl_wrong_m = (TH1*)f.Get(Form("m"+GetEra()+"/likelihood_mode1_%d_Mbl_Wrong_mode2_%d", mode1, mode2));
  hMblMET_correct_m = (TH1*)f.Get(Form("m"+GetEra()+"/likelihood_mode1_%d_MblMET_Correct_mode2_%d", mode1, mode2));
  hMblMET_wrong_m = (TH1*)f.Get(Form("m"+GetEra()+"/likelihood_mode1_%d_MblMET_Wrong_mode2_%d", mode1, mode2));
  hMbjj_correct_m = (TH1*)f.Get(Form("m"+GetEra()+"/likelihood_mode1_%d_Mbjj_Correct_mode2_%d", mode1, mode2));
  hMbjj_wrong_m = (TH1*)f.Get(Form("m"+GetEra()+"/likelihood_mode1_%d_Mbjj_Wrong_mode2_%d", mode1, mode2));
  hMjj_correct_m = (TH1*)f.Get(Form("m"+GetEra()+"/likelihood_mode1_%d_Mjj_Correct_mode2_%d", mode1, mode2));
  hMjj_wrong_m = (TH1*)f.Get(Form("m"+GetEra()+"/likelihood_mode1_%d_Mjj_Wrong_mode2_%d", mode1, mode2));
  hdRtt_correct_m = (TH1*)f.Get(Form("m"+GetEra()+"/likelihood_mode1_%d_dRtt_Correct_mode2_%d", mode1, mode2));
  hdRtt_wrong_m = (TH1*)f.Get(Form("m"+GetEra()+"/likelihood_mode1_%d_dRtt_Wrong_mode2_%d", mode1, mode2));
  hdPhitt_correct_m = (TH1*)f.Get(Form("m"+GetEra()+"/likelihood_mode1_%d_dPhitt_Correct_mode2_%d", mode1, mode2));
  hdPhitt_wrong_m = (TH1*)f.Get(Form("m"+GetEra()+"/likelihood_mode1_%d_dPhitt_Wrong_mode2_%d", mode1, mode2));

  if(hMbl_correct_E) hMbl_correct_E->SetDirectory(0);
  if(hMbl_wrong_E) hMbl_wrong_E->SetDirectory(0);
  if(hMblMET_correct_E) hMblMET_correct_E->SetDirectory(0);
  if(hMblMET_wrong_E) hMblMET_wrong_E->SetDirectory(0);
  if(hMbjj_correct_E) hMbjj_correct_E->SetDirectory(0);
  if(hMbjj_wrong_E) hMbjj_wrong_E->SetDirectory(0);
  if(hMjj_correct_E) hMjj_correct_E->SetDirectory(0);
  if(hMjj_wrong_E) hMjj_wrong_E->SetDirectory(0);
  if(hdRtt_correct_E) hdRtt_correct_E->SetDirectory(0);
  if(hdRtt_wrong_E) hdRtt_wrong_E->SetDirectory(0);
  if(hdPhitt_correct_E) hdPhitt_correct_E->SetDirectory(0);
  if(hdPhitt_wrong_E) hdPhitt_wrong_E->SetDirectory(0);

  if(hMbl_correct_m) hMbl_correct_m->SetDirectory(0);
  if(hMbl_wrong_m) hMbl_wrong_m->SetDirectory(0);
  if(hMblMET_correct_m) hMblMET_correct_m->SetDirectory(0);
  if(hMblMET_wrong_m) hMblMET_wrong_m->SetDirectory(0);
  if(hMbjj_correct_m) hMbjj_correct_m->SetDirectory(0);
  if(hMbjj_wrong_m) hMbjj_wrong_m->SetDirectory(0);
  if(hMjj_correct_m) hMjj_correct_m->SetDirectory(0);
  if(hMjj_wrong_m) hMjj_wrong_m->SetDirectory(0);
  if(hdRtt_correct_m) hdRtt_correct_m->SetDirectory(0);
  if(hdRtt_wrong_m) hdRtt_wrong_m->SetDirectory(0);
  if(hdPhitt_correct_m) hdPhitt_correct_m->SetDirectory(0);
  if(hdPhitt_wrong_m) hdPhitt_wrong_m->SetDirectory(0);

  cout<<"[ttljAnalyzer::SetupLikelihoods] All Likelihoods are set "<<endl;
  f.Close();
}

vector<unsigned int> ttljAnalyzer::Finding_bbjj_byLikelihood(TString channel, vector<Jet> bjets, vector<Jet> ajets, vector<double>& Likelihood_ratios, unsigned int mode3){
  unsigned int lb = 0, hb = 0, j0 = 0, j1 = 0;
  TH1* hMbl_correct = NULL;
  TH1* hMbl_wrong = NULL;
  TH1* hMblMET_correct = NULL;
  TH1* hMblMET_wrong = NULL;
  TH1* hMbjj_correct = NULL;
  TH1* hMbjj_wrong = NULL;
  TH1* hMjj_correct = NULL;
  TH1* hMjj_wrong = NULL;
  TH1* hdRtt_correct = NULL;
  TH1* hdRtt_wrong = NULL;
  TH1* hdPhitt_correct = NULL;
  TH1* hdPhitt_wrong = NULL;

  if(channel.Contains("E"+GetEraShort()) || channel.Contains("e"+GetEraShort())){
    hMbl_correct = hMbl_correct_E;
    hMbl_wrong = hMbl_wrong_E;
    hMblMET_correct = hMblMET_correct_E;
    hMblMET_wrong = hMblMET_wrong_E;
    hMbjj_correct = hMbjj_correct_E;
    hMbjj_wrong = hMbjj_wrong_E;
    hMjj_correct = hMjj_correct_E;
    hMjj_wrong = hMjj_wrong_E;
    hdRtt_correct = hdRtt_correct_E;
    hdRtt_wrong = hdRtt_wrong_E;
    hdPhitt_correct = hdPhitt_correct_E;
    hdPhitt_wrong = hdPhitt_wrong_E;
  }else if(channel.Contains("m"+GetEraShort())){
    hMbl_correct = hMbl_correct_m;
    hMbl_wrong = hMbl_wrong_m;
    hMblMET_correct = hMblMET_correct_m;
    hMblMET_wrong = hMblMET_wrong_m;
    hMbjj_correct = hMbjj_correct_m;
    hMbjj_wrong = hMbjj_wrong_m;
    hMjj_correct = hMjj_correct_m;
    hMjj_wrong = hMjj_wrong_m;
    hdRtt_correct = hdRtt_correct_m;
    hdRtt_wrong = hdRtt_wrong_m;
    hdPhitt_correct = hdPhitt_correct_m;
    hdPhitt_wrong = hdPhitt_wrong_m;
  }else{
  cout<<"[ttljAnalyzer::Finding_bbjj_byLikelihood] channel seems weird, channel = "+channel<<endl;
  return {lb, hb, j0, j1};
  }

  double max_likelihood_ratio = -1;
  for(unsigned int a=0; a<bjets.size(); a++){
    for(unsigned int b=0; b<bjets.size(); b++){
      if(a == b) continue;
      for(unsigned int c=0; c<ajets.size(); c++){
        for(unsigned int d=c+1; d<ajets.size(); d++){

          // Let bjets.at(a) = lepb, bjets.at(b) = hadb
          double Mbl = (bjets.at(a) + *lepton0).M();
          double MblMET = (bjets.at(a) + *lepton0 + met).M();
          double Mbjj = (bjets.at(b) + ajets.at(c) + ajets.at(d)).M();
          double Mjj = (ajets.at(c) + ajets.at(d)).M();
          double dRtt = (bjets.at(b) + ajets.at(c) + ajets.at(d)).DeltaR(bjets.at(a) + *lepton0 + met);
          double dPhitt = fabs((bjets.at(b) + ajets.at(c) + ajets.at(d)).DeltaPhi(bjets.at(a) + *lepton0 + met));

          double Likelihood_Mbl_correct = hMbl_correct->GetBinContent(hMbl_correct->FindBin(Mbl)) / hMbl_correct->Integral();
          double Likelihood_Mbl_wrong = hMbl_wrong->GetBinContent(hMbl_wrong->FindBin(Mbl)) / hMbl_wrong->Integral();
          double Likelihood_MblMET_correct = hMblMET_correct->GetBinContent(hMblMET_correct->FindBin(MblMET)) / hMblMET_correct->Integral();
          double Likelihood_MblMET_wrong = hMblMET_wrong->GetBinContent(hMblMET_wrong->FindBin(MblMET)) / hMblMET_wrong->Integral();
          double Likelihood_Mbjj_correct = hMbjj_correct->GetBinContent(hMbjj_correct->FindBin(Mbjj)) / hMbjj_correct->Integral();
          double Likelihood_Mbjj_wrong = hMbjj_wrong->GetBinContent(hMbjj_wrong->FindBin(Mbjj)) / hMbjj_wrong->Integral();
          double Likelihood_Mjj_correct = hMjj_correct->GetBinContent(hMjj_correct->FindBin(Mjj)) / hMjj_correct->Integral();
          double Likelihood_Mjj_wrong = hMjj_wrong->GetBinContent(hMjj_wrong->FindBin(Mjj)) / hMjj_wrong->Integral();
          double Likelihood_dRtt_correct = hdRtt_correct->GetBinContent(hdRtt_correct->FindBin(dRtt)) / hdRtt_correct->Integral();
          double Likelihood_dRtt_wrong = hdRtt_wrong->GetBinContent(hdRtt_wrong->FindBin(dRtt)) / hdRtt_wrong->Integral();
          double Likelihood_dPhitt_correct = hdPhitt_correct->GetBinContent(hdPhitt_correct->FindBin(dPhitt)) / hdPhitt_correct->Integral();
          double Likelihood_dPhitt_wrong = hdPhitt_wrong->GetBinContent(hdPhitt_wrong->FindBin(dPhitt)) / hdPhitt_wrong->Integral();

          double Likelihood_ratio_Mbl = Likelihood_Mbl_correct / (Likelihood_Mbl_correct + Likelihood_Mbl_wrong);
          double Likelihood_ratio_MblMET = Likelihood_MblMET_correct / (Likelihood_MblMET_correct + Likelihood_MblMET_wrong);
          double Likelihood_ratio_Mbjj = Likelihood_Mbjj_correct / (Likelihood_Mbjj_correct + Likelihood_Mbjj_wrong);
          double Likelihood_ratio_Mjj = Likelihood_Mjj_correct / (Likelihood_Mjj_correct + Likelihood_Mjj_wrong);
          double Likelihood_ratio_dRtt = Likelihood_dRtt_correct / (Likelihood_dRtt_correct + Likelihood_dRtt_wrong);
          double Likelihood_ratio_dPhitt = Likelihood_dPhitt_correct / (Likelihood_dPhitt_correct + Likelihood_dPhitt_wrong);
          double Likelihood_ratio = Likelihood_ratio_Mbl * Likelihood_ratio_MblMET * Likelihood_ratio_Mbjj * Likelihood_ratio_Mjj;
          if(mode3 == 1) Likelihood_ratio = Likelihood_ratio_Mbl;
          else if(mode3 == 2) Likelihood_ratio = Likelihood_ratio_Mbl * Likelihood_ratio_Mbjj;
          else if(mode3 == 3) Likelihood_ratio = Likelihood_ratio_Mbl * Likelihood_ratio_Mbjj * Likelihood_ratio_Mjj;
          else if(mode3 == 4) Likelihood_ratio = Likelihood_ratio_Mbl * Likelihood_ratio_Mbjj * Likelihood_ratio_MblMET;
          else if(mode3 == 5) Likelihood_ratio = Likelihood_ratio_Mbl * Likelihood_ratio_Mbjj * Likelihood_ratio_Mjj * Likelihood_ratio_MblMET * Likelihood_ratio_dPhitt;
          else if(mode3 == 6) Likelihood_ratio = Likelihood_ratio_Mbl * Likelihood_ratio_Mbjj * Likelihood_ratio_Mjj * Likelihood_ratio_MblMET * Likelihood_ratio_dPhitt * Likelihood_ratio_dRtt;
          //else if(mode3 == 2) Likelihood_ratio = Likelihood_ratio_Mbjj;
          //else if(mode3 == 3) Likelihood_ratio = Likelihood_ratio_Mjj;
          //else if(mode3 == 4) Likelihood_ratio = Likelihood_ratio_MblMET;
          //else if(mode3 == 5) Likelihood_ratio = Likelihood_ratio_dPhitt;
          //else if(mode3 == 6) Likelihood_ratio = Likelihood_ratio_dRtt;
          //else if(mode3 == 7) Likelihood_ratio = Likelihood_ratio_Mbjj * Likelihood_ratio_Mjj;
          //else if(mode3 == 8) Likelihood_ratio = Likelihood_ratio_Mbl * Likelihood_ratio_Mbjj * Likelihood_ratio_Mjj;
          //else if(mode3 == 9) Likelihood_ratio *= Likelihood_ratio_Mbjj * Likelihood_ratio_Mjj;

          if(Likelihood_ratio > max_likelihood_ratio){
            max_likelihood_ratio = Likelihood_ratio;
            Likelihood_ratios = {Likelihood_ratio_Mbl, Likelihood_ratio_MblMET, Likelihood_ratio_Mbjj, Likelihood_ratio_Mjj, Likelihood_ratio_dRtt, Likelihood_ratio_dPhitt};
            lb = a;
            hb = b;
            j0 = c;
            j1 = d;
          }
        }
      }
    }
  }

  return {lb, hb, j0, j1};
}

bool ttljAnalyzer::Kinematic_Cut(vector<Jet> bjets, vector<Jet> ajets){

  bool kinematic_cut = false;
  for(unsigned int a=0; a<bjets.size(); a++){
    for(unsigned int b=0; b<bjets.size(); b++){
      if(a == b) continue;
      for(unsigned int c=0; c<ajets.size(); c++){
        for(unsigned int d=c+1; d<ajets.size(); d++){

          TLorentzVector leptonic_top = bjets.at(a) + *lepton0 + met;
          TLorentzVector hadronic_top = bjets.at(b) + ajets.at(c) + ajets.at(d);

          //if((100 < hadronic_top.M() && hadronic_top.M() < 240) && ((bjets.at(a) + *lepton0).M() < 170) && (fabs(hadronic_top.DeltaPhi(leptonic_top)) > 1.5)){
          if(((bjets.at(a) + *lepton0).M() < 170) && (fabs(hadronic_top.DeltaPhi(leptonic_top)) > 2.0)){
            kinematic_cut = true;
            break;
          }
        }
      }
    }
  }
  return kinematic_cut;
}

ttljAnalyzer::ttljAnalyzer(){}
ttljAnalyzer::~ttljAnalyzer(){
  //==== Destructor of this Analyzer
}
