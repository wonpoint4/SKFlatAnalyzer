#include "ttljAnalyzer.h"

void ttljAnalyzer::initializeAnalyzer(){
  SMPAnalyzerCore::initializeAnalyzer(); //setup eff zpt roc z0
  IsNominalRun = !HasFlag("SYS") && !HasFlag("PDFSYS");

  PDFbase = LHAPDF::mkPDF(306000);
  PDFnf4 = LHAPDF::mkPDF(325500);
  mcCorr->SetJetTaggingParameters({
    JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Tight, JetTagging::incl, JetTagging::comb),
  });
  SetupPUJetWeight();
  SetupLikelihoods(0, 1);
}

void ttljAnalyzer::executeEvent(){
  ///////////////// GEN level /////////////////////
  if(!IsDATA) executeEventGen();

  ///////////////// RECO level /////////////////////
  if(!IsDATA || DataStream.Contains("SingleMuon")){
    executeEventWithParameter("m"+GetEraShort());
    if(nPV <= 10) executeEventWithParameter("m"+GetEraShort()+"F");
    else if(nPV <= 20) executeEventWithParameter("m"+GetEraShort()+"S");
    else if(nPV <= 30) executeEventWithParameter("m"+GetEraShort()+"L");
    else if(nPV <= 40) executeEventWithParameter("m"+GetEraShort()+"M");
    else if(nPV <= 50) executeEventWithParameter("m"+GetEraShort()+"H");
    else executeEventWithParameter("m"+GetEraShort()+"V");
    if(HasFlag("SYS")){
      for(TString syst:{"jet_scale_up", "jet_scale_down", "jet_smear_up", "jet_smear_down"}){
        if(syst.Contains("scale") || !IsDATA){
          executeEventWithParameter("m"+GetEraShort(), syst);
          if(nPV <= 10) executeEventWithParameter("m"+GetEraShort()+"F", syst);
          else if(nPV <= 20) executeEventWithParameter("m"+GetEraShort()+"S", syst);
          else if(nPV <= 30) executeEventWithParameter("m"+GetEraShort()+"L", syst);
          else if(nPV <= 40) executeEventWithParameter("m"+GetEraShort()+"M", syst);
          else if(nPV <= 50) executeEventWithParameter("m"+GetEraShort()+"H", syst);
          else executeEventWithParameter("m"+GetEraShort()+"V", syst);
        }
      }
    }
  }
  if(!IsDATA || DataStream.Contains("SingleElectron") || DataStream.Contains("EGamma")){
    executeEventWithParameter("e"+GetEraShort());
    executeEventWithParameter("E"+GetEraShort());
    if(nPV <= 10) executeEventWithParameter("E"+GetEraShort()+"F");
    else if(nPV <= 20) executeEventWithParameter("E"+GetEraShort()+"S");
    else if(nPV <= 30) executeEventWithParameter("E"+GetEraShort()+"L");
    else if(nPV <= 40) executeEventWithParameter("E"+GetEraShort()+"M");
    else if(nPV <= 50) executeEventWithParameter("E"+GetEraShort()+"H");
    else executeEventWithParameter("E"+GetEraShort()+"V");
    if(HasFlag("SYS")){
      for(TString syst:{"jet_scale_up", "jet_scale_down", "jet_smear_up", "jet_smear_down"}){
        if(syst.Contains("scale") || !IsDATA) executeEventWithParameter("e"+GetEraShort(), syst);
        if(syst.Contains("scale") || !IsDATA){
          executeEventWithParameter("E"+GetEraShort(), syst);
          if(nPV <= 10) executeEventWithParameter("E"+GetEraShort()+"F", syst);
          else if(nPV <= 20) executeEventWithParameter("E"+GetEraShort()+"S", syst);
          else if(nPV <= 30) executeEventWithParameter("E"+GetEraShort()+"L", syst);
          else if(nPV <= 40) executeEventWithParameter("E"+GetEraShort()+"M", syst);
          else if(nPV <= 50) executeEventWithParameter("E"+GetEraShort()+"H", syst);
          else executeEventWithParameter("E"+GetEraShort()+"V", syst);
        }
      }
    }
  }
}

void ttljAnalyzer::executeEventWithParameter(TString channel, TString option){

  lepton0 = NULL;
  prefix = channel+"/"+gprefix, hprefix = "", suffix = "";
  if(option.Contains("jet_scale_up")) suffix += "_jet_scale_up";
  else if(option.Contains("jet_scale_down")) suffix += "_jet_scale_down";
  else if(option.Contains("jet_smear_up")) suffix += "_jet_smear_up";
  else if(option.Contains("jet_smear_down")) suffix += "_jet_smear_down";
  if(IsNominalRun || option != "") IsNominalLike = true;
  else IsNominalLike = false;

  // Weights Setup
  map_weight.clear();
  lumiweight = 1., PUweight = 1., prefireweight = 1.;
  if(!IsDATA){
    lumiweight = reductionweight * MCweight() * _event.GetTriggerLumi("Full");
    PUweight = mcCorr->GetPileUpWeight(nPileUp, 0);
    prefireweight = L1PrefireReweight_Central;
  }

  map_weight[""] = lumiweight;
  if(IsNominalLike){
    FillHist(prefix+hprefix+"weight_Lumi"+suffix, lumiweight, map_weight[""], 200,-5,5);
    FillCutflow(prefix+hprefix+"cutflow"+suffix, "Lumi", map_weight[""]);
  }
  // Trigger
  if(!IsFiredTriggers(channel)) return;
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "PassTrig", map_weight[""]);

  // MET Filter
  if(!PassMETFilter()) return;
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "METfilter", map_weight[""]);

  // Single Lepton + pT + MET
  if(!Hasleptons(channel)) return;

  // Jets
  vector<Jet> alljets = {};
  if(option.Contains("jet_scale_up")) alljets = SelectJets(ScaleJets(GetAllJets(), 1), "tightLepVeto", 30, (DataYear == 2016? 2.4: 2.5));
  else if(option.Contains("jet_scale_down")) alljets = SelectJets(ScaleJets(GetAllJets(), -1), "tightLepVeto", 30, (DataYear == 2016? 2.4: 2.5));
  else if(option.Contains("jet_smear_up")) alljets = SelectJets(SmearJets(GetAllJets(), 1), "tightLepVeto", 30, (DataYear == 2016? 2.4: 2.5));
  else if(option.Contains("jet_smear_down")) alljets = SelectJets(SmearJets(GetAllJets(), -1), "tightLepVeto", 30, (DataYear == 2016? 2.4: 2.5));
  else alljets = SelectJets(GetAllJets(), "tightLepVeto", 30, (DataYear == 2016? 2.4: 2.5));
  std::sort(alljets.begin(), alljets.end(), PtComparing);

  vector<Jet> lepvetojets = {}, realjets = {}, bjets = {}, ajets = {};
  for(const auto& jet:alljets){
    if(lepton0 && jet.DeltaR(*lepton0) < 0.4) continue;
    lepvetojets.push_back(jet);
  }
  for(const auto& jet:lepvetojets){
    if(!PUJetIDPass(jet, "Loose")) continue;
    realjets.push_back(jet);
  }

  // b-tagging
  JetTagging::Parameters DeepJet_Tight = JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Tight, JetTagging::incl, JetTagging::comb);
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
  pujetSF = 1., btagSF = 1.;
  if(!IsDATA){
    pujetSF = GetPUJetWeight(lepvetojets, "Loose", 0);
    btagSF = mcCorr->GetBTaggingReweight_1a(realjets, DeepJet_Tight);
  }

  FillHist(prefix+"njets"+suffix, realjets.size(), map_weight[""], 15,0,15);
  FillHist(prefix+"nbjets"+suffix, bjets.size(), map_weight[""], 10,0,10);
  FillHist(prefix+"najets"+suffix, ajets.size(), map_weight[""], 10,0,10);

  if(bjets.size() != 2) return;
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "Tight2b", map_weight[""]);
  if(ajets.size() < 2) return;
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, ">=2j", map_weight[""]);

  // Weights
  if(!IsDATA && IsNominalLike) map_weight["_noWts"] = map_weight[""];
  map_weight[""] *= PUweight;
  if(IsNominalLike){
    FillHist(prefix+hprefix+"weight_PU"+suffix, PUweight, map_weight[""], 200,-5,5);
    FillCutflow(prefix+hprefix+"cutflow"+suffix, "PU", map_weight[""]);
  }
  map_weight[""] *= prefireweight;
  if(IsNominalLike){
    FillHist(prefix+hprefix+"weight_Prefire"+suffix, prefireweight, map_weight[""], 200,-5,5);
    FillCutflow(prefix+hprefix+"cutflow"+suffix, "Prefire", map_weight[""]);
  }
  map_weight[""] *= topptweight;
  if(IsNominalLike){
    FillHist(prefix+hprefix+"weight_Toppt"+suffix, topptweight, map_weight[""], 200,-5,5);
    FillCutflow(prefix+hprefix+"cutflow"+suffix, "Toppt", map_weight[""]);
  }

  // Lepton Efficiency Correction
  if(!IsDATA && IsNominalLike) map_weight["_noEffSF"] = map_weight[""];
  leptonTrackingSF = 1.;
  leptonRECOSF = 1.;
  leptonIDSF = 1.;
  leptonTriggerSF = 1.;

  if(!IsDATA){
    if(channel.Contains("m"+GetEraShort())){
      leptonTrackingSF *= fEff->GetEfficiencySF("Muon_Tracking", lepton0, 0,0);
      leptonRECOSF *= fEff->GetEfficiencySF("Muon_RECO", lepton0, 0,0);
      leptonIDSF *= fEff->GetEfficiencySF("Muon_MediumID_trkIsoLoose", lepton0, 0,0);
      TString trigSFkey = "IsoMu24_MediumID_trkIsoLoose";
      if(DataYear != 2017) leptonTriggerSF *= GetLeptonTriggerSF(trigSFkey, leptons, 0,0);
      else{
        vector<TString> trigSFkeys = {"IsoMu24_MediumID_trkIsoLoose", "IsoMu27_MediumID_trkIsoLoose"};
        leptonTriggerSF *= GetLeptonTriggerORSF(trigSFkeys, leptons, 0,0);
      }
    }else if(channel.Contains("e"+GetEraShort())){
      leptonRECOSF *= fEff->GetEfficiencySF("Electron_RECO", lepton0, 0,0);
      leptonIDSF *= fEff->GetEfficiencySF("Electron_MediumID", lepton0, 0,0);
      vector<TString> trigSFkeys = {"Ele27_MediumID"};
      if(DataYear == 2017) trigSFkeys = {"Ele27_MediumID", "Ele32_MediumID"};
      if(DataYear == 2018) trigSFkeys = {"Ele28_MediumID", "Ele32_MediumID"};
      if(DataYear == 2016) leptonTriggerSF *= GetLeptonTriggerSF(trigSFkeys[0], leptons, 0,0);
      else leptonTriggerSF *= GetLeptonTriggerORSF(trigSFkeys, leptons, 0,0);
    }else if(channel.Contains("E"+GetEraShort())){
      leptonRECOSF *= fEff->GetEfficiencySF("Electron_RECO", lepton0, 0,0);
      leptonIDSF *= fEff->GetEfficiencySF("Electron_SelQ_MediumID", lepton0, 0,0);
      vector<TString> trigSFkeys = {"Ele27_SelQ_MediumID"};
      if(DataYear == 2017) trigSFkeys = {"Ele27_SelQ_MediumID", "Ele32_SelQ_MediumID"};
      if(DataYear == 2018) trigSFkeys = {"Ele28_SelQ_MediumID", "Ele32_SelQ_MediumID"};
      if(DataYear == 2016) leptonTriggerSF *= GetLeptonTriggerSF(trigSFkeys[0], leptons, 0,0);
      else leptonTriggerSF *= GetLeptonTriggerORSF(trigSFkeys, leptons, 0,0);
    }
  }

  map_weight[""] *= leptonTrackingSF;
  if(IsNominalLike){
    FillHist(prefix+hprefix+"weight_TrackingSF"+suffix, leptonTrackingSF, map_weight[""], 200,-5,5);
    FillCutflow(prefix+hprefix+"cutflow"+suffix, "TrackingSF", map_weight[""]);
  }
  map_weight[""] *= leptonRECOSF;
  if(IsNominalLike){
    FillHist(prefix+hprefix+"weight_RECOSF"+suffix, leptonRECOSF, map_weight[""], 200,-5,5);
    FillCutflow(prefix+hprefix+"cutflow"+suffix, "RECOSF", map_weight[""]);
  }
  map_weight[""] *= leptonIDSF;
  if(IsNominalLike){
    FillHist(prefix+hprefix+"weight_IDSF"+suffix, leptonIDSF, map_weight[""], 200,-5,5);
    FillCutflow(prefix+hprefix+"cutflow"+suffix, "IDSF", map_weight[""]);
  }
  map_weight[""] *= leptonTriggerSF;
  if(IsNominalLike){
    FillHist(prefix+hprefix+"weight_TriggerSF"+suffix, leptonTriggerSF, map_weight[""], 200,-5,5);
    FillCutflow(prefix+hprefix+"cutflow"+suffix, "TriggerSF", map_weight[""]);
  }
  map_weight[""] *= pujetSF;
  if(IsNominalLike){
    FillHist(prefix+hprefix+"weight_PUjetSF"+suffix, pujetSF, map_weight[""], 200,-5,5);
    FillCutflow(prefix+hprefix+"cutflow"+suffix, "PUjetSF", map_weight[""]);
  }
  map_weight[""] *= btagSF;
  if(IsNominalLike){
    FillHist(prefix+hprefix+"weight_btagSF"+suffix, btagSF, map_weight[""], 200,-5,5);
    FillCutflow(prefix+hprefix+"cutflow"+suffix, "btagSF", map_weight[""]);
  }

  if(IsNominalLike){
    if(!IsDATA) map_weight["_bChargeSF0"] = map_weight[""] * GetbChargeSFWeight(bjets, 0, 0);
    FillHist(prefix+hprefix+"weight_bChargeSF0"+suffix, GetbChargeSFWeight(bjets, 0, 0), map_weight[""], 200,-5,5);
    if(!IsDATA) map_weight["_bChargeSF1"] = map_weight[""] * GetbChargeSFWeight(bjets, 1, 0);
    FillHist(prefix+hprefix+"weight_bChargeSF1"+suffix, GetbChargeSFWeight(bjets, 1, 0), map_weight[""], 200,-5,5);
  }

  //==== Making Likelihood
  if(IsNominalLike){
    if(!IsDATA && MCSample.Contains("TTLJ") && !channel.Contains("e"+GetEraShort())) FillingLikelihood(bjets, ajets, map_weight[""], suffix, 0, 0); // ByungHun Oh's method - drop events with ambiguity
    if(!IsDATA && MCSample.Contains("TTLJ") && !channel.Contains("e"+GetEraShort())) FillingLikelihood(bjets, ajets, map_weight[""], suffix, 1, 0); // Charmonium guy's method - match smaller dR < 0.3
  }

  //==== Finding the correct bbjj combination
  vector<unsigned int> idx_bbjj = {0, 0, 0, 0};
  vector<double> LRs = {0, 0, 0, 0, 0, 0};
  idx_bbjj = Finding_bbjj_byLikelihood(channel, bjets, ajets, LRs);

  bool goodKinematic = true;
  if((idx_bbjj.at(0) == idx_bbjj.at(1)) || (idx_bbjj.at(2) == idx_bbjj.at(3))) goodKinematic = false;
  if(IsNominalLike) FillHist(prefix+"LR_efficiency"+suffix, (goodKinematic? 1: 0), map_weight[""], 2,0,2);
  if(!goodKinematic) return;
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "Kin. cuts", map_weight[""]);

  Jet lepb = bjets.at(idx_bbjj.at(0));
  Jet hadb = bjets.at(idx_bbjj.at(1));
  Jet Wj0 = ajets.at(idx_bbjj.at(2));
  Jet Wj1 = ajets.at(idx_bbjj.at(3));
  double lepb_charge = jetCharge(lepb);
  double hadb_charge = jetCharge(hadb);
  double a0charge = jetCharge(Wj0);
  double a1charge = jetCharge(Wj1);

  //==== Weights of Systematics
  if(!IsDATA && HasFlag("SYS") && option == ""){
    // Prefiring weight
    map_weight["_noprefireweight"] =  map_weight[""] / prefireweight;
    map_weight["_prefireweight_up"] =  map_weight[""] / prefireweight * L1PrefireReweight_Up;
    map_weight["_prefireweight_down"] = map_weight[""] / prefireweight * L1PrefireReweight_Down;

    // PU reweight
    map_weight["_noPUweight"] =  map_weight[""] / PUweight;
    map_weight["_PUweight_up"] =  map_weight[""] / PUweight * mcCorr->GetPileUpWeight(nPileUp, 1);
    map_weight["_PUweight_down"] = map_weight[""] / PUweight * mcCorr->GetPileUpWeight(nPileUp, -1);

    // b-tagging SF
    map_weight["_nobtagSF"] =  map_weight[""] / btagSF;
    map_weight["_btagSF_hup"] = map_weight[""] / btagSF * mcCorr->GetBTaggingReweight_1a(realjets, DeepJet_Tight, "SystUpHTag");
    map_weight["_btagSF_hdown"] = map_weight[""] / btagSF * mcCorr->GetBTaggingReweight_1a(realjets, DeepJet_Tight, "SystDownHTag");
    map_weight["_btagSF_hcorr"] = map_weight[""] / btagSF * mcCorr->GetBTaggingReweight_1a(realjets, DeepJet_Tight, "SystUpHTagCorr");
    map_weight["_btagSF_huncorr2016a"] = map_weight[""];
    map_weight["_btagSF_huncorr2016b"] = map_weight[""];
    map_weight["_btagSF_huncorr2017"] = map_weight[""];
    map_weight["_btagSF_huncorr2018"] = map_weight[""];
    map_weight["_btagSF_huncorr"+GetEraShort()] = map_weight[""] / btagSF * mcCorr->GetBTaggingReweight_1a(realjets, DeepJet_Tight, "SystUpHTagUnCorr");
    map_weight["_btagSF_lup"] = map_weight[""] / btagSF * mcCorr->GetBTaggingReweight_1a(realjets, DeepJet_Tight, "SystUpLTag");
    map_weight["_btagSF_ldown"] = map_weight[""] / btagSF * mcCorr->GetBTaggingReweight_1a(realjets, DeepJet_Tight, "SystDownLTag");
    map_weight["_btagSF_lcorr"] = map_weight[""] / btagSF * mcCorr->GetBTaggingReweight_1a(realjets, DeepJet_Tight, "SystUpLTagCorr");;
    map_weight["_btagSF_luncorr2016a"] = map_weight[""];
    map_weight["_btagSF_luncorr2016b"] = map_weight[""];
    map_weight["_btagSF_luncorr2017"] = map_weight[""];
    map_weight["_btagSF_luncorr2018"] = map_weight[""];
    map_weight["_btagSF_luncorr"+GetEraShort()] = map_weight[""] / btagSF * mcCorr->GetBTaggingReweight_1a(realjets, DeepJet_Tight, "SystUpLTagUnCorr");

    // PUjetID SF
    map_weight["_noPUjetSF"] =  map_weight[""] / pujetSF;
    map_weight["_PUjetSF_up"] =  map_weight[""] / pujetSF * GetPUJetWeight(lepvetojets, "Loose", 1);
    map_weight["_PUjetSF_down"] = map_weight[""] / pujetSF * GetPUJetWeight(lepvetojets, "Loose", -1);

    // bChargeID SF
    map_weight["_bChargeSF0_up"] = map_weight[""] * GetbChargeSFWeight(bjets, 0, 1);
    map_weight["_bChargeSF0_down"] = map_weight[""] * GetbChargeSFWeight(bjets, 0, -1);
    map_weight["_bChargeSF1_up"] = map_weight[""] * GetbChargeSFWeight(bjets, 1, 1);
    map_weight["_bChargeSF1_down"] = map_weight[""] * GetbChargeSFWeight(bjets, 1, -1);
  }
  if(HasFlag("SYS") && option == "") map_weight.erase("");

  if(!IsDATA && MCSample.Contains("TTLJ")){
    bool isTTLJinAcceptance = true;
    for(Gen gen:{gen_b0, gen_b1, gen_j0, gen_j1}){
      if(gen.Pt() < 30) isTTLJinAcceptance = false;
      if(fabs(gen.Eta()) > (DataYear == 2016? 2.4: 2.5)) isTTLJinAcceptance = false;
    }
    FillHist(prefix+"LR_acceptance"+suffix, (isTTLJinAcceptance? 1: 0), map_weight, 2,0,2);

    double match_dR = 0.4;
    bool Gen_LR_match_onlyb = false;
    bool Gen_LR_match_Whad = false;
    bool Gen_LR_match_full = false;
    if((gen_j0.DeltaR(Wj0) < match_dR && gen_j1.DeltaR(Wj1) < match_dR) || (gen_j0.DeltaR(Wj1) < match_dR && gen_j1.DeltaR(Wj0) < match_dR)) Gen_LR_match_Whad = true;

    // ChargeEasy => 0:negative, 1:positive
    // purity => -1:UnMatched, 0:Wrong, 1:Correct in TTLJ
    // correct => 0:Wrong, 1:Correct in TTLJ

    //When lepb = b, hadb = bbar -> lep+
    if(gen_b0.DeltaR(lepb) < match_dR && gen_b1.DeltaR(hadb) < match_dR){
      Gen_LR_match_onlyb = true;
      if(Gen_LR_match_Whad) Gen_LR_match_full = true;
      if(lepton0->Charge() > 0){
        FillHist(prefix+"LR_purity"+suffix, 1, map_weight, 4,-2,2);
        FillHist(prefix+"LR_correct"+suffix, 1, map_weight, 2,0,2);
        if(isTTLJinAcceptance) FillHist(prefix+"LR_purity_inAcceptance"+suffix, 1, map_weight, 2,0,2);
        if(isTTLJinAcceptance) FillHist(prefix+"LR_correct_inAcceptance"+suffix, 1, map_weight, 2,0,2);
        FillHist(prefix+"LR_purity_Lp"+suffix, 1, map_weight, 4,-2,2);
        FillHist(prefix+"LR_correct_Lp"+suffix, 1, map_weight, 2,0,2);

        // For bCharge Accuracy (gen-info)
        FillHist(prefix+"genbCharge_Lp"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        FillHist(prefix+"genbbarCharge_Lp"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        for(unsigned int i=1; i<afb_chbinnum+1; i++){
          if(afb_chbin[i-1] < abs(lepb_charge) && abs(lepb_charge) < afb_chbin[i]) FillHist(Form(prefix+"genbCharge%d_Lp"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
          if(afb_chbin[i-1] < abs(hadb_charge) && abs(hadb_charge) < afb_chbin[i]) FillHist(Form(prefix+"genbbarCharge%d_Lp"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        }

        prefix += "Correct"; // lep+
      }else{
        FillHist(prefix+"LR_purity"+suffix, 0, map_weight, 4,-2,2);
        FillHist(prefix+"LR_correct"+suffix, 0, map_weight, 2,0,2);
        if(isTTLJinAcceptance) FillHist(prefix+"LR_purity_inAcceptance"+suffix, 0, map_weight, 2,0,2);
        if(isTTLJinAcceptance) FillHist(prefix+"LR_correct_inAcceptance"+suffix, 0, map_weight, 2,0,2);
        FillHist(prefix+"LR_purity_Lm"+suffix, 0, map_weight, 4,-2,2);
        FillHist(prefix+"LR_correct_Lm"+suffix, 0, map_weight, 2,0,2);

        // For bCharge Accuracy (gen-info)
        FillHist(prefix+"genbCharge_Lm"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        FillHist(prefix+"genbbarCharge_Lm"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        for(unsigned int i=1; i<afb_chbinnum+1; i++){
          if(afb_chbin[i-1] < abs(lepb_charge) && abs(lepb_charge) < afb_chbin[i]) FillHist(Form(prefix+"genbCharge%d_Lm"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
          if(afb_chbin[i-1] < abs(hadb_charge) && abs(hadb_charge) < afb_chbin[i]) FillHist(Form(prefix+"genbbarCharge%d_Lm"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        }

        prefix += "Wrong"; // lep-
      }
    }
    //When lepb = bbar, hadb = b => lep-
    else if(gen_b0.DeltaR(hadb) < match_dR && gen_b1.DeltaR(lepb) < match_dR){
      Gen_LR_match_onlyb = true;
      if(Gen_LR_match_Whad) Gen_LR_match_full = true;
      if(lepton0->Charge() < 0){
        FillHist(prefix+"LR_purity"+suffix, 1, map_weight, 4,-2,2);
        FillHist(prefix+"LR_correct"+suffix, 1, map_weight, 2,0,2);
        if(isTTLJinAcceptance) FillHist(prefix+"LR_purity_inAcceptance"+suffix, 1, map_weight, 2,0,2);
        if(isTTLJinAcceptance) FillHist(prefix+"LR_correct_inAcceptance"+suffix, 1, map_weight, 2,0,2);
        FillHist(prefix+"LR_purity_Lm"+suffix, 1, map_weight, 4,-2,2);
        FillHist(prefix+"LR_correct_Lm"+suffix, 1, map_weight, 2,0,2);

        // For bCharge Accuracy (gen-info)
        FillHist(prefix+"genbCharge_Lm"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        FillHist(prefix+"genbbarCharge_Lm"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        for(unsigned int i=1; i<afb_chbinnum+1; i++){
          if(afb_chbin[i-1] < abs(hadb_charge) && abs(hadb_charge) < afb_chbin[i]) FillHist(Form(prefix+"genbCharge%d_Lm"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
          if(afb_chbin[i-1] < abs(lepb_charge) && abs(lepb_charge) < afb_chbin[i]) FillHist(Form(prefix+"genbbarCharge%d_Lm"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        }

        prefix += "Correct"; // lep-
      }else{
        FillHist(prefix+"LR_purity"+suffix, 0, map_weight, 4,-2,2);
        FillHist(prefix+"LR_correct"+suffix, 0, map_weight, 2,0,2);
        if(isTTLJinAcceptance) FillHist(prefix+"LR_purity_inAcceptance"+suffix, 0, map_weight, 2,0,2);
        if(isTTLJinAcceptance) FillHist(prefix+"LR_correct_inAcceptance"+suffix, 0, map_weight, 2,0,2);
        FillHist(prefix+"LR_purity_Lp"+suffix, 0, map_weight, 4,-2,2);
        FillHist(prefix+"LR_correct_Lp"+suffix, 0, map_weight, 2,0,2);

        // For bCharge Accuracy (gen-info)
        FillHist(prefix+"genbCharge_Lp"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        FillHist(prefix+"genbbarCharge_Lp"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        for(unsigned int i=1; i<afb_chbinnum+1; i++){
          if(afb_chbin[i-1] < abs(hadb_charge) && abs(hadb_charge) < afb_chbin[i]) FillHist(Form(prefix+"genbCharge%d_Lp"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
          if(afb_chbin[i-1] < abs(lepb_charge) && abs(lepb_charge) < afb_chbin[i]) FillHist(Form(prefix+"genbbarCharge%d_Lp"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        }

        prefix += "Wrong"; // lep+
      }
    }
    else{
      FillHist(prefix+"LR_purity"+suffix, -1, map_weight, 4,-2,2);
      if(isTTLJinAcceptance) FillHist(prefix+"LR_purity_inAcceptance"+suffix, -1, map_weight, 2,0,2);
      if(lepton0->Charge() > 0) FillHist(prefix+"LR_purity_Lp"+suffix, -1, map_weight, 4,-2,2);
      if(lepton0->Charge() < 0) FillHist(prefix+"LR_purity_Lm"+suffix, -1, map_weight, 4,-2,2);
      prefix += "UnMatched";
    }
    prefix += (isTTLJinAcceptance? "0_": "1_");

    FillHist("Gen_LR_match_onlyb"+suffix, Gen_LR_match_onlyb, map_weight, 2,0,2);
    FillHist("Gen_LR_match_Whad"+suffix, Gen_LR_match_Whad, map_weight, 2,0,2);
    FillHist("Gen_LR_match_full"+suffix, Gen_LR_match_full, map_weight, 2,0,2);
    FillHist(prefix+"Gen_LR_match_onlyb"+suffix, Gen_LR_match_onlyb, map_weight, 2,0,2);
    FillHist(prefix+"Gen_LR_match_Whad"+suffix, Gen_LR_match_Whad, map_weight, 2,0,2);
    FillHist(prefix+"Gen_LR_match_full"+suffix, Gen_LR_match_full, map_weight, 2,0,2);
    FillHist(channel+"/Gen_LR_match_onlyb"+suffix, Gen_LR_match_onlyb, map_weight, 2,0,2);
    FillHist(channel+"/Gen_LR_match_Whad"+suffix, Gen_LR_match_Whad, map_weight, 2,0,2);
    FillHist(channel+"/Gen_LR_match_full"+suffix, Gen_LR_match_full, map_weight, 2,0,2);
    if(isTTLJinAcceptance){
      FillHist("Gen_LR_match_onlyb_inAcceptance"+suffix, Gen_LR_match_onlyb, map_weight, 2,0,2);
      FillHist("Gen_LR_match_Whad_inAcceptance"+suffix, Gen_LR_match_Whad, map_weight, 2,0,2);
      FillHist("Gen_LR_match_full_inAcceptance"+suffix, Gen_LR_match_full, map_weight, 2,0,2);
      FillHist(prefix+"Gen_LR_match_onlyb_inAcceptance"+suffix, Gen_LR_match_onlyb, map_weight, 2,0,2);
      FillHist(prefix+"Gen_LR_match_Whad_inAcceptance"+suffix, Gen_LR_match_Whad, map_weight, 2,0,2);
      FillHist(prefix+"Gen_LR_match_full_inAcceptance"+suffix, Gen_LR_match_full, map_weight, 2,0,2);
      FillHist(channel+"/Gen_LR_match_onlyb_inAcceptance"+suffix, Gen_LR_match_onlyb, map_weight, 2,0,2);
      FillHist(channel+"/Gen_LR_match_Whad_inAcceptance"+suffix, Gen_LR_match_Whad, map_weight, 2,0,2);
      FillHist(channel+"/Gen_LR_match_full_inAcceptance"+suffix, Gen_LR_match_full, map_weight, 2,0,2);
    }
  }

  //==========================
  //==== Now reco fill histograms
  //==========================
  FillHist(prefix+"mass_jj"+suffix, (Wj0 + Wj1).M(), map_weight, 40,0,200);
  FillHist(prefix+"mass_bjj"+suffix, (hadb + Wj0 + Wj1).M(), map_weight, 40,100,300);
  //FillHist(prefix+"mass_Toplep"+suffix, (lepb + *lepton0 + neutrino).M(), map_weight, 40,100,300);
  FillHist(prefix+"mass_blMET"+suffix, (lepb + *lepton0 + met).M(), map_weight, 40,100,300);
  FillHist(prefix+"mass_bl"+suffix, (lepb + *lepton0).M(), map_weight, 40,0,200);

  FillHist(prefix+"dR_bjj_blMET"+suffix, (hadb + Wj0 + Wj1).DeltaR(lepb + *lepton0 + met), map_weight, 200,0,10);
  FillHist(prefix+"dR_bjj_lMET"+suffix, (hadb + Wj0 + Wj1).DeltaR(*lepton0 + met), map_weight, 200,0,10);
  FillHist(prefix+"dR_bjj_jj"+suffix, (hadb + Wj0 + Wj1).DeltaR(Wj0 + Wj1), map_weight, 200,0,10);
  FillHist(prefix+"dR_blMET_lMET"+suffix, (lepb + *lepton0 + met).DeltaR(*lepton0 + met), map_weight, 200,0,10);
  FillHist(prefix+"dR_blMET_jj"+suffix, (lepb + *lepton0 + met).DeltaR(Wj0 + Wj1), map_weight, 200,0,10);
  FillHist(prefix+"dR_lMET_jj"+suffix, (*lepton0 + met).DeltaR(Wj0 + Wj1), map_weight, 200,0,10);

  FillHist(prefix+"dPhi_bjj_blMET"+suffix, fabs((hadb + Wj0 + Wj1).DeltaPhi(lepb + *lepton0 + met)), map_weight, 200,0,10);
  FillHist(prefix+"dPhi_bjj_lMET"+suffix, fabs((hadb + Wj0 + Wj1).DeltaPhi(*lepton0 + met)), map_weight, 200,0,10);
  FillHist(prefix+"dPhi_bjj_jj"+suffix, fabs((hadb + Wj0 + Wj1).DeltaPhi(Wj0 + Wj1)), map_weight, 200,0,10);
  FillHist(prefix+"dPhi_blMET_lMET"+suffix, fabs((lepb + *lepton0 + met).DeltaPhi(*lepton0 + met)), map_weight, 200,0,10);
  FillHist(prefix+"dPhi_blMET_jj"+suffix, fabs((lepb + *lepton0 + met).DeltaPhi(Wj0 + Wj1)), map_weight, 200,0,10);
  FillHist(prefix+"dPhi_lMET_jj"+suffix, fabs((*lepton0 + met).DeltaPhi(Wj0 + Wj1)), map_weight, 200,0,10);

  FillHist(prefix+"likelihood_ratio"+suffix, LRs.at(0) * LRs.at(1) * LRs.at(2) * LRs.at(3), map_weight, 1000,0,50);
  FillHist(prefix+"likelihood_ratio_Mbl"+suffix, LRs.at(0), map_weight, 100,0,5);
  FillHist(prefix+"likelihood_ratio_MblMET"+suffix, LRs.at(1), map_weight, 100,0,5);
  FillHist(prefix+"likelihood_ratio_Mbjj"+suffix, LRs.at(2), map_weight, 100,0,5);
  FillHist(prefix+"likelihood_ratio_Mjj"+suffix, LRs.at(3), map_weight, 100,0,5);
  FillHist(prefix+"likelihood_ratio_dRtt"+suffix, LRs.at(4), map_weight, 100,0,5);
  FillHist(prefix+"likelihood_ratio_dPhitt"+suffix, LRs.at(5), map_weight, 100,0,5);

  FillHist(prefix+"lepbChargeRaw"+suffix, lepb_charge, map_weight, 200,-5,5);
  FillHist(prefix+"lepbCharge"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
  FillHist(prefix+"hadbChargeRaw"+suffix, hadb_charge, map_weight, 200,-5,5);
  FillHist(prefix+"hadbCharge"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
  FillHist(prefix+"a0ChargeRaw"+suffix, a0charge, map_weight, 200,-5,5);
  FillHist(prefix+"a0Charge"+suffix, (a0charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
  FillHist(prefix+"a1ChargeRaw"+suffix, a1charge, map_weight, 200,-5,5);
  FillHist(prefix+"a1Charge"+suffix, (a1charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
  if(lepton0->Charge() < 0){
    FillHist(prefix+"lepbChargeRaw_Lm"+suffix, lepb_charge, map_weight, 200,-5,5);
    FillHist(prefix+"lepbCharge_Lm"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
    FillHist(prefix+"hadbChargeRaw_Lm"+suffix, hadb_charge, map_weight, 200,-5,5);
    FillHist(prefix+"hadbCharge_Lm"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
    FillHist(prefix+"recobbarChargeRaw"+suffix, lepb_charge, map_weight, 200,-5,5);
    FillHist(prefix+"recobbarCharge"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
    FillHist(prefix+"recobChargeRaw"+suffix, hadb_charge, map_weight, 200,-5,5);
    FillHist(prefix+"recobCharge"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
  }else{
    FillHist(prefix+"lepbChargeRaw_Lp"+suffix, lepb_charge, map_weight, 200,-5,5);
    FillHist(prefix+"lepbCharge_Lp"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
    FillHist(prefix+"hadbChargeRaw_Lp"+suffix, hadb_charge, map_weight, 200,-5,5);
    FillHist(prefix+"hadbCharge_Lp"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
    FillHist(prefix+"recobChargeRaw"+suffix, lepb_charge, map_weight, 200,-5,5);
    FillHist(prefix+"recobCharge"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
    FillHist(prefix+"recobbarChargeRaw"+suffix, hadb_charge, map_weight, 200,-5,5);
    FillHist(prefix+"recobbarCharge"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
  }
  for(unsigned int i=1; i<afb_chbinnum+1; i++){
    if(lepton0->Charge() < 0){
      if(afb_chbin[i-1] < abs(lepb_charge) && abs(lepb_charge) < afb_chbin[i]){
        FillHist(Form(prefix+"lepbCharge%dRaw_Lm"+suffix, i-1), lepb_charge, map_weight, 200,-5,5);
        FillHist(Form(prefix+"lepbCharge%d_Lm"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        FillHist(Form(prefix+"recobbarCharge%dRaw"+suffix, i-1), lepb_charge, map_weight, 200,-5,5);
        FillHist(Form(prefix+"recobbarCharge%d"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
      }
      if(afb_chbin[i-1] < abs(hadb_charge) && abs(hadb_charge) < afb_chbin[i]){
        FillHist(Form(prefix+"hadbCharge%dRaw_Lm"+suffix, i-1), hadb_charge, map_weight, 200,-5,5);
        FillHist(Form(prefix+"hadbCharge%d_Lm"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        FillHist(Form(prefix+"recobCharge%dRaw"+suffix, i-1), hadb_charge, map_weight, 200,-5,5);
        FillHist(Form(prefix+"recobCharge%d"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
      }
    }else{
      if(afb_chbin[i-1] < abs(lepb_charge) && abs(lepb_charge) < afb_chbin[i]){
        FillHist(Form(prefix+"lepbCharge%dRaw_Lp"+suffix, i-1), lepb_charge, map_weight, 200,-5,5);
        FillHist(Form(prefix+"lepbCharge%d_Lp"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        FillHist(Form(prefix+"recobCharge%dRaw"+suffix, i-1), lepb_charge, map_weight, 200,-5,5);
        FillHist(Form(prefix+"recobCharge%d"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
      }
      if(afb_chbin[i-1] < abs(hadb_charge) && abs(hadb_charge) < afb_chbin[i]){
        FillHist(Form(prefix+"hadbCharge%dRaw_Lp"+suffix, i-1), hadb_charge, map_weight, 200,-5,5);
        FillHist(Form(prefix+"hadbCharge%d_Lp"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        FillHist(Form(prefix+"recobbarCharge%dRaw"+suffix, i-1), hadb_charge, map_weight, 200,-5,5);
        FillHist(Form(prefix+"recobbarCharge%d"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
      }
    }
  }

  FillHist(prefix+hprefix+"lpt"+suffix, lepton0->Pt(), map_weight, AFBAnalyzer::lptbinnum,AFBAnalyzer::lptbin);
  FillHist(prefix+hprefix+"bpt"+suffix, lepb.Pt(), map_weight, AFBAnalyzer::lptbinnum,AFBAnalyzer::lptbin);
  FillHist(prefix+hprefix+"bpt"+suffix, hadb.Pt(), map_weight, AFBAnalyzer::lptbinnum,AFBAnalyzer::lptbin);
  FillHist(prefix+hprefix+"lepbpt"+suffix, lepb.Pt(), map_weight, AFBAnalyzer::lptbinnum,AFBAnalyzer::lptbin);
  FillHist(prefix+hprefix+"hadbpt"+suffix, hadb.Pt(), map_weight, AFBAnalyzer::lptbinnum,AFBAnalyzer::lptbin);
  FillHist(prefix+hprefix+"jpt"+suffix, Wj0.Pt(), map_weight, AFBAnalyzer::lptbinnum,AFBAnalyzer::lptbin);
  FillHist(prefix+hprefix+"jpt"+suffix, Wj1.Pt(), map_weight, AFBAnalyzer::lptbinnum,AFBAnalyzer::lptbin);
  FillHist(prefix+hprefix+"j0pt"+suffix, Wj0.Pt(), map_weight, AFBAnalyzer::lptbinnum,AFBAnalyzer::lptbin);
  FillHist(prefix+hprefix+"j1pt"+suffix, Wj1.Pt(), map_weight, AFBAnalyzer::lptbinnum,AFBAnalyzer::lptbin);
  FillHist(prefix+hprefix+"j0idx"+suffix, idx_bbjj.at(2), map_weight, 15,0,5);
  FillHist(prefix+hprefix+"j1idx"+suffix, idx_bbjj.at(3), map_weight, 15,0,5);

  FillHist(prefix+hprefix+"leta"+suffix, lepton0->Eta(), map_weight, 50,-2.5,2.5);
  FillHist(prefix+hprefix+"beta"+suffix, lepb.Eta(), map_weight, 100,-5,5);
  FillHist(prefix+hprefix+"beta"+suffix, hadb.Eta(), map_weight, 100,-5,5);
  FillHist(prefix+hprefix+"lepbeta"+suffix, lepb.Eta(), map_weight, 100,-5,5);
  FillHist(prefix+hprefix+"hadbeta"+suffix, hadb.Eta(), map_weight, 100,-5,5);
  FillHist(prefix+hprefix+"jeta"+suffix, Wj0.Eta(), map_weight, 100,-5,5);
  FillHist(prefix+hprefix+"jeta"+suffix, Wj1.Eta(), map_weight, 100,-5,5);
  FillHist(prefix+hprefix+"j0eta"+suffix, Wj0.Eta(), map_weight, 100,-5,5);
  FillHist(prefix+hprefix+"j1eta"+suffix, Wj1.Eta(), map_weight, 100,-5,5);

  FillHist(prefix+hprefix+"asChargeSum"+suffix, a0charge + a1charge, map_weight, 400,-10,10);
  FillHist(prefix+hprefix+"asChargeAbsSum"+suffix, (a0charge < 0? -1: 1) + (a1charge < 0? -1: 1), map_weight, 8,-4,4);
  FillHist(prefix+hprefix+"bsChargeSum"+suffix, lepb_charge + hadb_charge, map_weight, 400,-10,10);
  FillHist(prefix+hprefix+"bsChargeAbsSum"+suffix, (lepb_charge < 0? -1: 1) + (hadb_charge < 0? -1: 1), map_weight, 8,-4,4);

  FillHist(prefix+hprefix+"nPV"+suffix, nPV, map_weight, 100,0,100);
}

bool ttljAnalyzer::Hasleptons(TString channel){
  bool moreleptons = false;
  double l0pt = 26.;
  if(channel.Contains("m"+GetEraShort())){
    muons = MuonMomentumCorrection(SMPGetMuons("POGMediumWithLooseTrkIso", 8.0, 2.4), 0,0,0);
    if(muons.size() > 0) lepton0 = &muons.at(0);
    if(muons.size() > 1) moreleptons = true;
  }else if(channel.Contains("e"+GetEraShort())){
    l0pt = 30.;
    electrons = ElectronEnergyCorrection(SMPGetElectrons("passMediumID", 8.0, 2.5), 0,0);
    if(electrons.size() > 0) lepton0 = &electrons.at(0);
    if(electrons.size() > 1) moreleptons = true;
  }else if(channel.Contains("E"+GetEraShort())){
    l0pt = 30.;
    electrons = ElectronEnergyCorrection(SMPGetElectrons("passMediumID_SelQ", 8.0, 2.5), 0,0);
    if(electrons.size() > 0) lepton0 = &electrons.at(0);
    if(electrons.size() > 1) moreleptons = true;
  }else{
    cout<<"[ttljAnalyzer::Hasleptons] channel="<<channel<<" is weird"<<endl;
  }

  if(!lepton0) return false;
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "1Leptons", map_weight[""]);
  if(lepton0->Pt() < l0pt) return false;
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "LepPt", map_weight[""]);
  if(moreleptons) return false;
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "1Lepton", map_weight[""]);
  met = GetEvent().GetMETVector();
  if(met.Pt() < 20.) return false;
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "MET20", map_weight[""]);

  leptons ={};
  leptons.push_back(lepton0);

  return true;
}

void ttljAnalyzer::executeEventGen(){
  FillHist("gen/executeEventGen", 1, 1, 2,0,2);
  vector<Gen> gens=GetGens();
  if(IsTTSample) topptweight=mcCorr->GetTopPtReweight(gens);
  if(!MCSample.Contains("TTLJ")) return;

  GetTTLJGenParticles(gens, gen_parton0,gen_parton1, gen_b0,gen_b1, gen_l0,gen_l1, gen_j0,gen_j1, 3);
  if(gen_l0.PID() < 0){ // lepton is l+, thus lepb, l+ from t
    gen_lepb = &gen_b0;
    gen_hadb = &gen_b1;
  }else{
    gen_lepb = &gen_b1;
    gen_hadb = &gen_b0;
  }

  if(!IsNominalRun) return;

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
  FillHist("gen/Pt_j", gen_j0.Pt(), 1, 100,0,400);
  FillHist("gen/Pt_j", gen_j1.Pt(), 1, 100,0,400);
  FillHist("gen/Pt_j0", gen_j0.Pt(), 1, 100,0,400);
  FillHist("gen/Pt_j1", gen_j1.Pt(), 1, 100,0,400);

  FillHist("gen/Eta_parton0", gen_parton0.Eta(), 1, 200,-5,5);
  FillHist("gen/Eta_parton1", gen_parton1.Eta(), 1, 200,-5,5);
  FillHist("gen/Eta_b0", gen_b0.Eta(), 1, 200,-5,5);
  FillHist("gen/Eta_b1", gen_b1.Eta(), 1, 200,-5,5);
  FillHist("gen/Eta_l0", gen_l0.Eta(), 1, 200,-5,5);
  FillHist("gen/Eta_l1", gen_l1.Eta(), 1, 200,-5,5);
  FillHist("gen/Eta_j", gen_j0.Eta(), 1, 200,-5,5);
  FillHist("gen/Eta_j", gen_j1.Eta(), 1, 200,-5,5);
  FillHist("gen/Eta_j0", gen_j0.Eta(), 1, 200,-5,5);
  FillHist("gen/Eta_j1", gen_j1.Eta(), 1, 200,-5,5);

  FillHist("gen/2D_pteta_b", gen_b0.Pt(), gen_b0.Eta(), 1, 100,0,400, 100,-5,5);
  FillHist("gen/2D_pteta_b", gen_b1.Pt(), gen_b1.Eta(), 1, 100,0,400, 100,-5,5);
  FillHist("gen/2D_pteta_j", gen_j0.Pt(), gen_j0.Eta(), 1, 100,0,400, 100,-5,5);
  FillHist("gen/2D_pteta_j", gen_j1.Pt(), gen_j1.Eta(), 1, 100,0,400, 100,-5,5);
  FillHist("gen/2D_pteta_j0", gen_j0.Pt(), gen_j0.Eta(), 1, 100,0,400, 100,-5,5);
  FillHist("gen/2D_pteta_j1", gen_j1.Pt(), gen_j1.Eta(), 1, 100,0,400, 100,-5,5);

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

  // Checks for various Acceptances
  bool acceptance_2b = false;
  bool acceptance_2j = false;
  bool acceptance_2l = false;
  if(gen_b0.Pt() > 30 && fabs(gen_b0.Eta()) < 2.4 && gen_b1.Pt() > 30 && fabs(gen_b1.Eta()) < 2.4) acceptance_2b = true;
  if(gen_j0.Pt() > 30 && fabs(gen_j0.Eta()) < 2.4 && gen_j1.Pt() > 30 && fabs(gen_j1.Eta()) < 2.4) acceptance_2j = true;
  if(gen_l0.Pt() > 30 && fabs(gen_l0.Eta()) < 2.4 && gen_l1.Pt() > 20) acceptance_2l = true;

  FillHist("gen/acceptance_pteta30_2b", (acceptance_2b? 1: 0), 1, 2,0,2);
  FillHist("gen/acceptance_pteta30_2j", (acceptance_2j? 1: 0), 1, 2,0,2);
  FillHist("gen/acceptance_pteta30_2l", (acceptance_2l? 1: 0), 1, 2,0,2);
  FillHist("gen/acceptance_pteta30_2b2j", ((acceptance_2b && acceptance_2j)? 1: 0), 1, 2,0,2);
  FillHist("gen/acceptance_pteta30_2b2j2l", ((acceptance_2b && acceptance_2j && acceptance_2l)? 1: 0), 1, 2,0,2);

  acceptance_2b = false;
  acceptance_2j = false;
  acceptance_2l = false;
  if(gen_b0.Pt() > 25 && fabs(gen_b0.Eta()) < 2.4 && gen_b1.Pt() > 25 && fabs(gen_b1.Eta()) < 2.4) acceptance_2b = true;
  if(gen_j0.Pt() > 25 && fabs(gen_j0.Eta()) < 2.4 && gen_j1.Pt() > 25 && fabs(gen_j1.Eta()) < 2.4) acceptance_2j = true;
  if(gen_l0.Pt() > 25 && fabs(gen_l0.Eta()) < 2.4 && gen_l1.Pt() > 20) acceptance_2l = true;

  FillHist("gen/acceptance_pteta_2b", (acceptance_2b? 1: 0), 1, 2,0,2);
  FillHist("gen/acceptance_pteta_2j", (acceptance_2j? 1: 0), 1, 2,0,2);
  FillHist("gen/acceptance_pteta_2l", (acceptance_2l? 1: 0), 1, 2,0,2);
  FillHist("gen/acceptance_pteta_2b2j", ((acceptance_2b && acceptance_2j)? 1: 0), 1, 2,0,2);
  FillHist("gen/acceptance_pteta_2b2j2l", ((acceptance_2b && acceptance_2j && acceptance_2l)? 1: 0), 1, 2,0,2);

  acceptance_2b = false;
  acceptance_2j = false;
  acceptance_2l = false;
  if(gen_b0.Pt() > 20 && fabs(gen_b0.Eta()) < 2.4 && gen_b1.Pt() > 20 && fabs(gen_b1.Eta()) < 2.4) acceptance_2b = true;
  if(gen_j0.Pt() > 20 && fabs(gen_j0.Eta()) < 2.4 && gen_j1.Pt() > 20 && fabs(gen_j1.Eta()) < 2.4) acceptance_2j = true;
  if(gen_l0.Pt() > 20 && fabs(gen_l0.Eta()) < 2.4 && gen_l1.Pt() > 20) acceptance_2l = true;

  FillHist("gen/acceptance_pteta20_2b", (acceptance_2b? 1: 0), 1, 2,0,2);
  FillHist("gen/acceptance_pteta20_2j", (acceptance_2j? 1: 0), 1, 2,0,2);
  FillHist("gen/acceptance_pteta20_2l", (acceptance_2l? 1: 0), 1, 2,0,2);
  FillHist("gen/acceptance_pteta20_2b2j", ((acceptance_2b && acceptance_2j)? 1: 0), 1, 2,0,2);
  FillHist("gen/acceptance_pteta20_2b2j2l", ((acceptance_2b && acceptance_2j && acceptance_2l)? 1: 0), 1, 2,0,2);

  acceptance_2b = false;
  acceptance_2j = false;
  acceptance_2l = false;
  if(gen_b0.Pt() > 25 && gen_b1.Pt() > 25) acceptance_2b = true;
  if(gen_j0.Pt() > 25 && gen_j1.Pt() > 25) acceptance_2j = true;
  if(gen_l0.Pt() > 25 && gen_l1.Pt() > 20) acceptance_2l = true;

  FillHist("gen/acceptance_pt_2b", (acceptance_2b? 1: 0), 1, 2,0,2);
  FillHist("gen/acceptance_pt_2j", (acceptance_2j? 1: 0), 1, 2,0,2);
  FillHist("gen/acceptance_pt_2l", (acceptance_2l? 1: 0), 1, 2,0,2);
  FillHist("gen/acceptance_pt_2b2j", ((acceptance_2b && acceptance_2j)? 1: 0), 1, 2,0,2);
  FillHist("gen/acceptance_pt_2b2j2l", ((acceptance_2b && acceptance_2j && acceptance_2l)? 1: 0), 1, 2,0,2);

  acceptance_2b = false;
  acceptance_2j = false;
  acceptance_2l = false;
  if(fabs(gen_b0.Eta()) < 2.4 && fabs(gen_b1.Eta()) < 2.4) acceptance_2b = true;
  if(fabs(gen_j0.Eta()) < 2.4 && fabs(gen_j1.Eta()) < 2.4) acceptance_2j = true;
  if(fabs(gen_l0.Eta()) < 2.4) acceptance_2l = true;

  FillHist("gen/acceptance_eta_2b", (acceptance_2b? 1: 0), 1, 2,0,2);
  FillHist("gen/acceptance_eta_2j", (acceptance_2j? 1: 0), 1, 2,0,2);
  FillHist("gen/acceptance_eta_2l", (acceptance_2l? 1: 0), 1, 2,0,2);
  FillHist("gen/acceptance_eta_2b2j", ((acceptance_2b && acceptance_2j)? 1: 0), 1, 2,0,2);
  FillHist("gen/acceptance_eta_2b2j2l", ((acceptance_2b && acceptance_2j && acceptance_2l)? 1: 0), 1, 2,0,2);
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

void ttljAnalyzer::FillingLikelihood(vector<Jet> bjets, vector<Jet> ajets, double weight, TString suffix, unsigned int mode1, unsigned int mode2){
  FillCutflow(Form("likelihood"+GetEraShort()+"/cutflow_goodMatching_mode1_%d_ForLikelihood"+suffix, mode1), "Preselection", weight);
  double match_dR = 0.4;
  if(gen_l0.DeltaR(*lepton0) > match_dR) return;
  if(mode1 == 1) match_dR = 0.3;
  FillCutflow(Form("likelihood"+GetEraShort()+"/cutflow_goodMatching_mode1_%d_ForLikelihood"+suffix, mode1), "1lep", weight);

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
  FillCutflow(Form("likelihood"+GetEraShort()+"/cutflow_goodMatching_mode1_%d_ForLikelihood"+suffix, mode1), "2b-1gen", weight);

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
  FillCutflow(Form("likelihood"+GetEraShort()+"/cutflow_goodMatching_mode1_%d_ForLikelihood"+suffix, mode1), "2j-1gen", weight);

  if(!lepb) return;
  FillCutflow(Form("likelihood"+GetEraShort()+"/cutflow_goodMatching_mode1_%d_ForLikelihood"+suffix, mode1), "nolepb", weight);
  if(!hadb) return;
  FillCutflow(Form("likelihood"+GetEraShort()+"/cutflow_goodMatching_mode1_%d_ForLikelihood"+suffix, mode1), "nohadb", weight);
  if(!Wj0) return;
  FillCutflow(Form("likelihood"+GetEraShort()+"/cutflow_goodMatching_mode1_%d_ForLikelihood"+suffix, mode1), "noWj0", weight);
  if(!Wj1) return;
  FillCutflow(Form("likelihood"+GetEraShort()+"/cutflow_goodMatching_mode1_%d_ForLikelihood"+suffix, mode1), "noWj1", weight);

  if(!lepb || !hadb || !Wj0 || !Wj1) return;
  FillCutflow(Form("likelihood"+GetEraShort()+"/cutflow_goodMatching_mode1_%d_ForLikelihood"+suffix, mode1), "4jets", weight);

  if(lepb == hadb) return;
  FillCutflow(Form("likelihood"+GetEraShort()+"/cutflow_goodMatching_mode1_%d_ForLikelihood"+suffix, mode1), "1b-2gens", weight);
  if(Wj0 == Wj1) return;
  FillCutflow(Form("likelihood"+GetEraShort()+"/cutflow_goodMatching_mode1_%d_ForLikelihood"+suffix, mode1), "1j-2gens", weight);

  // Filling Likelihood Histograms
  // mode2 = 0 (Wrong : wrong && lepb,hadb switched)
  for(unsigned int a=0; a<bjets.size(); a++){
    for(unsigned int b=a+1; b<bjets.size(); b++){
      if(&bjets.at(a) == lepb && &bjets.at(b) == hadb){
        FillHist(Form("likelihood"+GetEraShort()+"/Correct_Mbl_mode1_%d_mode2_0"+suffix, mode1), (bjets.at(a) + *lepton0).M(), weight, 200,0,200);
        FillHist(Form("likelihood"+GetEraShort()+"/Correct_MblMET_mode1_%d_mode2_0"+suffix, mode1), (bjets.at(a) + *lepton0 + met).M(), weight, 200,0,300);
        FillHist(Form("likelihood"+GetEraShort()+"/Wrong_Mbl_mode1_%d_mode2_0"+suffix, mode1), (bjets.at(b) + *lepton0).M(), weight, 200,0,200);
        FillHist(Form("likelihood"+GetEraShort()+"/Wrong_MblMET_mode1_%d_mode2_0"+suffix, mode1), (bjets.at(b) + *lepton0 + met).M(), weight, 200,0,300);
        for(unsigned int c=0; c<ajets.size(); c++){
          for(unsigned int d=c+1; d<ajets.size(); d++){
            if((&ajets.at(c) == Wj0 && &ajets.at(d) == Wj1) || (&ajets.at(c) == Wj1 && &ajets.at(d) == Wj0)){
              FillHist(Form("likelihood"+GetEraShort()+"/Correct_Mbjj_mode1_%d_mode2_0"+suffix, mode1), (bjets.at(b) + ajets.at(c) + ajets.at(d)).M(), weight, 200,0,300);
              FillHist(Form("likelihood"+GetEraShort()+"/Correct_Mjj_mode1_%d_mode2_0"+suffix, mode1), (ajets.at(c) + ajets.at(d)).M(), weight, 200,0,200);
              FillHist(Form("likelihood"+GetEraShort()+"/Correct_dRtt_mode1_%d_mode2_0"+suffix, mode1), (bjets.at(b) + ajets.at(c) + ajets.at(d)).DeltaR(bjets.at(a) + *lepton0 + met), weight, 120,0,6);
              FillHist(Form("likelihood"+GetEraShort()+"/Correct_mode1_%d_dPhitt_mode2_0"+suffix, mode1), fabs((bjets.at(b) + ajets.at(c) + ajets.at(d)).DeltaPhi(bjets.at(a) + *lepton0 + met)), weight, 60,0,3);
              FillHist(Form("likelihood"+GetEraShort()+"/Wrong_Mbjj_mode1_%d_mode2_0"+suffix, mode1), (bjets.at(a) + ajets.at(c) + ajets.at(d)).M(), weight, 200,0,300);
              FillHist(Form("likelihood"+GetEraShort()+"/Wrong_dRtt_mode1_%d_mode2_0"+suffix, mode1), (bjets.at(a) + ajets.at(c) + ajets.at(d)).DeltaR(bjets.at(b) + *lepton0 + met), weight, 120,0,6);
              FillHist(Form("likelihood"+GetEraShort()+"/Wrong_mode1_%d_dPhitt_mode2_0"+suffix, mode1), fabs((bjets.at(a) + ajets.at(c) + ajets.at(d)).DeltaPhi(bjets.at(b) + *lepton0 + met)), weight, 60,0,3);
            }
            else{
              FillHist(Form("likelihood"+GetEraShort()+"/Wrong_Mbjj_mode1_%d_mode2_0"+suffix, mode1), (bjets.at(a) + ajets.at(c) + ajets.at(d)).M(), weight, 200,0,300);
              FillHist(Form("likelihood"+GetEraShort()+"/Wrong_Mjj_mode1_%d_mode2_0"+suffix, mode1), (ajets.at(c) + ajets.at(d)).M(), weight, 200,0,200);
              FillHist(Form("likelihood"+GetEraShort()+"/Wrong_dRtt_mode1_%d_mode2_0"+suffix, mode1), (bjets.at(a) + ajets.at(c) + ajets.at(d)).DeltaR(bjets.at(b) + *lepton0 + met), weight, 120,0,6);
              FillHist(Form("likelihood"+GetEraShort()+"/Wrong_mode1_%d_dPhitt_mode2_0"+suffix, mode1), fabs((bjets.at(a) + ajets.at(c) + ajets.at(d)).DeltaPhi(bjets.at(b) + *lepton0 + met)), weight, 60,0,3);
            }
          }
        }
      }else if(&bjets.at(a) == hadb && &bjets.at(b) == lepb){
        FillHist(Form("likelihood"+GetEraShort()+"/Correct_Mbl_mode1_%d_mode2_0"+suffix, mode1), (bjets.at(b) + *lepton0).M(), weight, 200,0,200);
        FillHist(Form("likelihood"+GetEraShort()+"/Correct_MblMET_mode1_%d_mode2_0"+suffix, mode1), (bjets.at(b) + *lepton0 + met).M(), weight, 200,0,300);
        FillHist(Form("likelihood"+GetEraShort()+"/Wrong_Mbl_mode1_%d_mode2_0"+suffix, mode1), (bjets.at(a) + *lepton0).M(), weight, 200,0,200);
        FillHist(Form("likelihood"+GetEraShort()+"/Wrong_MblMET_mode1_%d_mode2_0"+suffix, mode1), (bjets.at(a) + *lepton0 + met).M(), weight, 200,0,300);
        for(unsigned int c=0; c<ajets.size(); c++){
          for(unsigned int d=c+1; d<ajets.size(); d++){
            if((&ajets.at(c) == Wj0 && &ajets.at(d) == Wj1) || (&ajets.at(c) == Wj1 && &ajets.at(d) == Wj0)){
              FillHist(Form("likelihood"+GetEraShort()+"/Correct_Mbjj_mode1_%d_mode2_0"+suffix, mode1), (bjets.at(a) + ajets.at(c) + ajets.at(d)).M(), weight, 200,0,300);
              FillHist(Form("likelihood"+GetEraShort()+"/Correct_Mjj_mode1_%d_mode2_0"+suffix, mode1), (ajets.at(c) + ajets.at(d)).M(), weight, 200,0,200);
              FillHist(Form("likelihood"+GetEraShort()+"/Correct_dRtt_mode1_%d_mode2_0"+suffix, mode1), (bjets.at(a) + ajets.at(c) + ajets.at(d)).DeltaR(bjets.at(b) + *lepton0 + met), weight, 120,0,6);
              FillHist(Form("likelihood"+GetEraShort()+"/Correct_dPhitt_mode1_%d_mode2_0"+suffix, mode1), fabs((bjets.at(a) + ajets.at(c) + ajets.at(d)).DeltaPhi(bjets.at(b) + *lepton0 + met)), weight, 60,0,3);
              FillHist(Form("likelihood"+GetEraShort()+"/Wrong_Mbjj_mode1_%d_mode2_0"+suffix, mode1), (bjets.at(b) + ajets.at(c) + ajets.at(d)).M(), weight, 200,0,300);
              FillHist(Form("likelihood"+GetEraShort()+"/Wrong_dRtt_mode1_%d_mode2_0"+suffix, mode1), (bjets.at(b) + ajets.at(c) + ajets.at(d)).DeltaR(bjets.at(a) + *lepton0 + met), weight, 120,0,6);
              FillHist(Form("likelihood"+GetEraShort()+"/Wrong_dPhitt_mode1_%d_mode2_0"+suffix, mode1), fabs((bjets.at(b) + ajets.at(c) + ajets.at(d)).DeltaPhi(bjets.at(a) + *lepton0 + met)), weight, 60,0,3);
            }
            else{
              FillHist(Form("likelihood"+GetEraShort()+"/Wrong_Mbjj_mode1_%d_mode2_0"+suffix, mode1), (bjets.at(b) + ajets.at(c) + ajets.at(d)).M(), weight, 200,0,300);
              FillHist(Form("likelihood"+GetEraShort()+"/Wrong_Mjj_mode1_%d_mode2_0"+suffix, mode1), (ajets.at(c) + ajets.at(d)).M(), weight, 200,0,200);
              FillHist(Form("likelihood"+GetEraShort()+"/Wrong_dRtt_mode1_%d_mode2_0"+suffix, mode1), (bjets.at(b) + ajets.at(c) + ajets.at(d)).DeltaR(bjets.at(a) + *lepton0 + met), weight, 120,0,6);
              FillHist(Form("likelihood"+GetEraShort()+"/Wrong_dPhitt_mode1_%d_mode2_0"+suffix, mode1), fabs((bjets.at(b) + ajets.at(c) + ajets.at(d)).DeltaPhi(bjets.at(a) + *lepton0 + met)), weight, 60,0,3);
            }
          }
        }
      }
    }

    // mode2 = 1 (Wrong : simply wrong)
    if(&bjets.at(a) == lepb){
      FillHist(Form("likelihood"+GetEraShort()+"/Correct_Mbl_mode1_%d_mode2_1"+suffix, mode1), (bjets.at(a) + *lepton0).M(), weight, 200,0,200);
      FillHist(Form("likelihood"+GetEraShort()+"/Correct_MblMET_mode1_%d_mode2_1"+suffix, mode1), (bjets.at(a) + *lepton0 + met).M(), weight, 200,0,300);
    }else{
      FillHist(Form("likelihood"+GetEraShort()+"/Wrong_Mbl_mode1_%d_mode2_1"+suffix, mode1), (bjets.at(a) + *lepton0).M(), weight, 200,0,200);
      FillHist(Form("likelihood"+GetEraShort()+"/Wrong_MblMET_mode1_%d_mode2_1"+suffix, mode1), (bjets.at(a) + *lepton0 + met).M(), weight, 200,0,300);
    }

    if(&bjets.at(a) == hadb){
      for(unsigned int c=0; c<ajets.size(); c++){
        for(unsigned int d=c+1; d<ajets.size(); d++){
          if((&ajets.at(c) == Wj0 && &ajets.at(d) == Wj1) || (&ajets.at(c) == Wj1 && &ajets.at(d) == Wj0)){
            FillHist(Form("likelihood"+GetEraShort()+"/Correct_Mbjj_mode1_%d_mode2_1_idxWj0"+suffix, mode1), (c < d? c: d), weight, 30,0,30);
            FillHist(Form("likelihood"+GetEraShort()+"/Correct_Mbjj_mode1_%d_mode2_1_idxWj1"+suffix, mode1), (c > d? c: d), weight, 30,0,30);
            double c_ratio = c, d_ratio = d;
            c_ratio /= ajets.size(); d_ratio /= ajets.size();
            FillHist(Form("likelihood"+GetEraShort()+"/Correct_Mbjj_mode1_%d_mode2_1_idxratioWj0"+suffix, mode1), (c < d? c_ratio: d_ratio), weight, 11,0,1.1);
            FillHist(Form("likelihood"+GetEraShort()+"/Correct_Mbjj_mode1_%d_mode2_1_idxratioWj1"+suffix, mode1), (c > d? c_ratio: d_ratio), weight, 11,0,1.1);
            FillHist(Form("likelihood"+GetEraShort()+"/Correct_Mbjj_mode1_%d_mode2_1_pWj0"+suffix, mode1), (c < d? ajets.at(c).P(): ajets.at(d).P()), weight, 50,0,500);
            FillHist(Form("likelihood"+GetEraShort()+"/Correct_Mbjj_mode1_%d_mode2_1_pWj1"+suffix, mode1), (c > d? ajets.at(c).P(): ajets.at(d).P()), weight, 50,0,500);
            FillHist(Form("likelihood"+GetEraShort()+"/Correct_Mbjj_mode1_%d_mode2_1_ptWj0"+suffix, mode1), (c < d? ajets.at(c).Pt(): ajets.at(d).Pt()), weight, 50,0,500);
            FillHist(Form("likelihood"+GetEraShort()+"/Correct_Mbjj_mode1_%d_mode2_1_ptWj1"+suffix, mode1), (c > d? ajets.at(c).Pt(): ajets.at(d).Pt()), weight, 50,0,500);
            FillHist(Form("likelihood"+GetEraShort()+"/Correct_Mbjj_mode1_%d_mode2_1_etaWj0"+suffix, mode1), (c < d? ajets.at(c).Eta(): ajets.at(d).Eta()), weight, 100,-5,5);
            FillHist(Form("likelihood"+GetEraShort()+"/Correct_Mbjj_mode1_%d_mode2_1_etaWj1"+suffix, mode1), (c > d? ajets.at(c).Eta(): ajets.at(d).Eta()), weight, 100,-5,5);

            FillHist(Form("likelihood"+GetEraShort()+"/Correct_Mbjj_mode1_%d_mode2_1"+suffix, mode1), (bjets.at(a) + ajets.at(c) + ajets.at(d)).M(), weight, 200,0,300);
            FillHist(Form("likelihood"+GetEraShort()+"/Correct_Mjj_mode1_%d_mode2_1"+suffix, mode1), (ajets.at(c) + ajets.at(d)).M(), weight, 200,0,200);
          }
          else{
            FillHist(Form("likelihood"+GetEraShort()+"/Wrong_Mbjj_mode1_%d_mode2_1"+suffix, mode1), (bjets.at(a) + ajets.at(c) + ajets.at(d)).M(), weight, 200,0,300);
            FillHist(Form("likelihood"+GetEraShort()+"/Wrong_Mjj_mode1_%d_mode2_1"+suffix, mode1), (ajets.at(c) + ajets.at(d)).M(), weight, 200,0,200);
          }
        }
      }
      for(unsigned int b=a+1; b<bjets.size(); b++){
        for(unsigned int c=0; c<ajets.size(); c++){
          for(unsigned int d=c+1; d<ajets.size(); d++){
            if(&bjets.at(b) == lepb){
              if((&ajets.at(c) == Wj0 && &ajets.at(d) == Wj1) || (&ajets.at(c) == Wj1 && &ajets.at(d) == Wj0)){
                FillHist(Form("likelihood"+GetEraShort()+"/Correct_dRtt_mode1_%d_mode2_1"+suffix, mode1), (bjets.at(a) + ajets.at(c) + ajets.at(d)).DeltaR(bjets.at(b) + *lepton0 + met), weight, 120,0,6);
                FillHist(Form("likelihood"+GetEraShort()+"/Correct_dPhitt_mode1_%d_mode2_1"+suffix, mode1), fabs((bjets.at(a) + ajets.at(c) + ajets.at(d)).DeltaPhi(bjets.at(b) + *lepton0 + met)), weight, 60,0,3);
              }
            }else{
              FillHist(Form("likelihood"+GetEraShort()+"/Wrong_dRtt_mode1_%d_mode2_1"+suffix, mode1), (bjets.at(a) + ajets.at(c) + ajets.at(d)).DeltaR(bjets.at(b) + *lepton0 + met), weight, 120,0,6);
              FillHist(Form("likelihood"+GetEraShort()+"/Wrong_dPhitt_mode1_%d_mode2_1"+suffix, mode1), fabs((bjets.at(a) + ajets.at(c) + ajets.at(d)).DeltaPhi(bjets.at(b) + *lepton0 + met)), weight, 60,0,3);
            }
          }
        }
      }
    }else{
      for(unsigned int c=0; c<ajets.size(); c++){
        for(unsigned int d=c+1; d<ajets.size(); d++){
          if((&ajets.at(c) == Wj0 && &ajets.at(d) == Wj1) || (&ajets.at(c) == Wj1 && &ajets.at(d) == Wj0)){
            FillHist(Form("likelihood"+GetEraShort()+"/Wrong_Mbjj_mode1_%d_mode2_1"+suffix, mode1), (bjets.at(a) + ajets.at(c) + ajets.at(d)).M(), weight, 200,0,300);
            FillHist(Form("likelihood"+GetEraShort()+"/Correct_Mjj_mode1_%d_mode2_1"+suffix, mode1), (ajets.at(c) + ajets.at(d)).M(), weight, 200,0,200);
          }
          else{
            FillHist(Form("likelihood"+GetEraShort()+"/Wrong_Mbjj_mode1_%d_mode2_1"+suffix, mode1), (bjets.at(a) + ajets.at(c) + ajets.at(d)).M(), weight, 200,0,300);
            FillHist(Form("likelihood"+GetEraShort()+"/Wrong_Mjj_mode1_%d_mode2_1"+suffix, mode1), (ajets.at(c) + ajets.at(d)).M(), weight, 200,0,200);
          }
        }
      }
      for(unsigned int b=a+1; b<bjets.size(); b++){
        for(unsigned int c=0; c<ajets.size(); c++){
          for(unsigned int d=c+1; d<ajets.size(); d++){
            FillHist(Form("likelihood"+GetEraShort()+"/Wrong_dRtt_mode1_%d_mode2_1"+suffix, mode1), (bjets.at(a) + ajets.at(c) + ajets.at(d)).DeltaR(bjets.at(b) + *lepton0 + met), weight, 120,0,6);
            FillHist(Form("likelihood"+GetEraShort()+"/Wrong_dPhitt_mode1_%d_mode2_1"+suffix, mode1), fabs((bjets.at(a) + ajets.at(c) + ajets.at(d)).DeltaPhi(bjets.at(b) + *lepton0 + met)), weight, 60,0,3);
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
  hMbl_correct = (TH1*)f.Get(Form("likelihood"+GetEraShort()+"/Correct_Mbl_mode1_%d_mode2_%d", mode1, mode2));
  hMbl_wrong = (TH1*)f.Get(Form("likelihood"+GetEraShort()+"/Wrong_Mbl_mode1_%d_mode2_%d", mode1, mode2));
  hMblMET_correct = (TH1*)f.Get(Form("likelihood"+GetEraShort()+"/Correct_MblMET_mode1_%d_mode2_%d", mode1, mode2));
  hMblMET_wrong = (TH1*)f.Get(Form("likelihood"+GetEraShort()+"/Wrong_MblMET_mode1_%d_mode2_%d", mode1, mode2));
  hMbjj_correct = (TH1*)f.Get(Form("likelihood"+GetEraShort()+"/Correct_Mbjj_mode1_%d_mode2_%d", mode1, mode2));
  hMbjj_wrong = (TH1*)f.Get(Form("likelihood"+GetEraShort()+"/Wrong_Mbjj_mode1_%d_mode2_%d", mode1, mode2));
  hMjj_correct = (TH1*)f.Get(Form("likelihood"+GetEraShort()+"/Correct_Mjj_mode1_%d_mode2_%d", mode1, mode2));
  hMjj_wrong = (TH1*)f.Get(Form("likelihood"+GetEraShort()+"/Wrong_Mjj_mode1_%d_mode2_%d", mode1, mode2));
  hdRtt_correct = (TH1*)f.Get(Form("likelihood"+GetEraShort()+"/Correct_dRtt_mode1_%d_mode2_%d", mode1, mode2));
  hdRtt_wrong = (TH1*)f.Get(Form("likelihood"+GetEraShort()+"/Wrong_dRtt_mode1_%d_mode2_%d", mode1, mode2));
  hdPhitt_correct = (TH1*)f.Get(Form("likelihood"+GetEraShort()+"/Correct_dPhitt_mode1_%d_mode2_%d", mode1, mode2));
  hdPhitt_wrong = (TH1*)f.Get(Form("likelihood"+GetEraShort()+"/Wrong_dPhitt_mode1_%d_mode2_%d", mode1, mode2));

  if(hMbl_correct) hMbl_correct->SetDirectory(0);
  if(hMbl_wrong) hMbl_wrong->SetDirectory(0);
  if(hMblMET_correct) hMblMET_correct->SetDirectory(0);
  if(hMblMET_wrong) hMblMET_wrong->SetDirectory(0);
  if(hMbjj_correct) hMbjj_correct->SetDirectory(0);
  if(hMbjj_wrong) hMbjj_wrong->SetDirectory(0);
  if(hMjj_correct) hMjj_correct->SetDirectory(0);
  if(hMjj_wrong) hMjj_wrong->SetDirectory(0);
  if(hdRtt_correct) hdRtt_correct->SetDirectory(0);
  if(hdRtt_wrong) hdRtt_wrong->SetDirectory(0);
  if(hdPhitt_correct) hdPhitt_correct->SetDirectory(0);
  if(hdPhitt_wrong) hdPhitt_wrong->SetDirectory(0);

  cout<<"[ttljAnalyzer::SetupLikelihoods] All Likelihoods are set "<<endl;
  f.Close();
}

vector<unsigned int> ttljAnalyzer::Finding_bbjj_byLikelihood(TString channel, vector<Jet> bjets, vector<Jet> ajets, vector<double>& Likelihood_ratios, unsigned int mode3){
  unsigned int lb = 0, hb = 0, j0 = 0, j1 = 0;
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

          // Kinematic Cuts
          //if(c > 4 || d > 4) break;
          if(Mbjj < 100 || 240 < Mbjj) continue;
          if(Mbl > 170) continue;
          if(fabs(Mjj - 80.4) > 30) continue;

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

          double Likelihood_ratio_Mbl = Likelihood_Mbl_correct / Likelihood_Mbl_wrong;
          double Likelihood_ratio_MblMET = Likelihood_MblMET_correct / Likelihood_MblMET_wrong;
          double Likelihood_ratio_Mbjj = Likelihood_Mbjj_correct / Likelihood_Mbjj_wrong;
          double Likelihood_ratio_Mjj = Likelihood_Mjj_correct / Likelihood_Mjj_wrong;
          double Likelihood_ratio_dRtt = Likelihood_dRtt_correct / Likelihood_dRtt_wrong;
          double Likelihood_ratio_dPhitt = Likelihood_dPhitt_correct / Likelihood_dPhitt_wrong;

          double Likelihood_ratio = Likelihood_ratio_Mbl * Likelihood_ratio_MblMET * Likelihood_ratio_Mbjj * Likelihood_ratio_Mjj;
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

// From Hyonsan's [SMPAnalyzerCore::GetLeptonTriggerORSF]
double ttljAnalyzer::GetLeptonTriggerORSF(const vector<TString> trigkeys, const vector<Lepton*>& leps, int set, int mem, TString option){
  if(IsDATA) return 1;
  if(trigkeys.size() != 2){
    cout<<"[ttljAnalyzer::LeptonTriggerOR_SF] trigkeys.size()= "<<trigkeys.size()<<endl;
    exit(EXIT_FAILURE);
  }

  double lumi = _event.GetTriggerLumi("Full");
  double lumi0, lumi1, lumi01;
  TString trig0, trig1;
  if(DataYear == 2017 && trigkeys[0].Contains("IsoMu24") && trigkeys[1].Contains("IsoMu27")){
    trig0 = "HLT_IsoMu24_v";
    trig1 = "HLT_IsoMu27_v";
    lumi0 = _event.GetTriggerLumi(trig0);
    lumi1 = _event.GetTriggerLumi(trig1);
    lumi01 = lumi0;
  }else if(DataYear == 2017 && trigkeys[0].Contains("Ele27") && trigkeys[1].Contains("Ele32")){
    trig0 = "HLT_Ele27_WPTight_Gsf_v";
    trig1 = "HLT_Ele32_WPTight_Gsf_v";
    lumi0 = _event.GetTriggerLumi(trig0);
    lumi1 = _event.GetTriggerLumi(trig1);
    lumi01 = 17599.732185;
  }else if(DataYear == 2018 && trigkeys[0].Contains("Ele28") && trigkeys[1].Contains("Ele32")){
    trig0 = "HLT_Ele28_WPTight_Gsf_v";
    trig1 = "HLT_Ele32_WPTight_Gsf_v";
    lumi0 = _event.GetTriggerLumi(trig0);
    lumi1 = _event.GetTriggerLumi(trig1);
    lumi01 = lumi0;
  }else{
    cout<<"[ttljAnalyzer::GetLeptonTriggerORSF] not available combination '"<<trigkeys[0]<<"'||'"<<trigkeys[1]<<"' for "<<DataEra<<endl;
    exit(EXIT_FAILURE);
  }

  double data_eff0 = 1., sim_eff0 = 1.;
  double data_eff1 = 1., sim_eff1 = 1.;
  for(const auto& lep:leps){
    if(!lep) continue;
    data_eff0 *= 1 - fEff->GetDataEfficiency(trigkeys[0], lep, set, mem, option);
    sim_eff0 *= 1 - fEff->GetSimEfficiency(trigkeys[0], lep, set, mem, option);
    data_eff1 *= 1 - fEff->GetDataEfficiency(trigkeys[1], lep, set, mem, option);
    sim_eff1 *= 1 - fEff->GetSimEfficiency(trigkeys[1], lep, set, mem, option);
  }
  data_eff0 = 1 - data_eff0;
  sim_eff0 = 1 - sim_eff0;
  data_eff1 = 1 - data_eff1;
  sim_eff1 = 1 - sim_eff1;

  double sf=0.;
  if(_event.PassTrigger(trig1)){
    double this_sf = (lumi1 - lumi01) / lumi;
    if(sim_eff1) this_sf *= data_eff1 / sim_eff1;
    sf += this_sf;
  }
  if(_event.PassTrigger(trig0)){
    double this_sf = (lumi0 - lumi01) / lumi;
    if(sim_eff0) this_sf *= data_eff0 / sim_eff0;
    sf += this_sf;
  }
  //overlap region
  if(_event.PassTrigger(trig0)){
    double this_sf = lumi01 / lumi / 2;
    if(sim_eff0) this_sf *= data_eff0 / sim_eff0;
    sf += this_sf;
  }
  if(_event.PassTrigger(trig1)){
    double this_sf = lumi01 / lumi / 2;
    if(sim_eff1) this_sf *= data_eff1 / sim_eff1;
    sf += this_sf;
  }else if(_event.PassTrigger(trig0)){
    double this_sf = lumi01 / lumi / 2;
    if(sim_eff0) this_sf *= data_eff0 / sim_eff0;
    sf += this_sf;
  }

  return sf;
}

ttljAnalyzer::ttljAnalyzer(){}
ttljAnalyzer::~ttljAnalyzer(){
  //==== Destructor of this Analyzer
}
