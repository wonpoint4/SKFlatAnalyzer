#include "dybAnalyzer.h"

void dybAnalyzer::initializeAnalyzer(){
  SMPAnalyzerCore::initializeAnalyzer(); //setup eff zpt roc z0
  IsSkimmed = GetSkimName() != ""? true: false;
  IsNominalRun = !HasFlag("SYS") && !HasFlag("PDFSYS") && IsSkimmed;

  PDFbase = LHAPDF::mkPDF(306000);
  PDFnf4 = LHAPDF::mkPDF(325500);
  mcCorr->SetJetTaggingParameters({
    JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Tight, JetTagging::incl, JetTagging::comb),
    JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Loose, JetTagging::incl, JetTagging::comb),
  });
  SetupPUJetWeight();
}

void dybAnalyzer::executeEvent(){
  //// FIXME some events of DYJets has nan PDF weights. I don't know why...
  if(MCSample == "DYJets" && !isnormal(weight_Scale->at(0))) return;

  ///////////////// GEN level /////////////////////
  executeEventGen();

  ///////////////// RECO level /////////////////////
  if(!IsDATA || DataStream.Contains("DoubleMuon")){
    executeEventWithParameter("mm"+GetEraShort());
    if(HasFlag("SYS")){
      for(TString syst:{"jet_scale_up", "jet_scale_down", "jet_smear_up", "jet_smear_down"}){
        if(syst.Contains("scale") || !IsDATA) executeEventWithParameter("mm"+GetEraShort(), syst);
      }
    }
  }
  if(!IsDATA || DataStream.Contains("DoubleEG") || DataStream.Contains("EGamma")){
    executeEventWithParameter("ee"+GetEraShort());
    if(HasFlag("SYS")){
      for(TString syst:{"jet_scale_up", "jet_scale_down", "jet_smear_up", "jet_smear_down"}){
        if(syst.Contains("scale") || !IsDATA) executeEventWithParameter("ee"+GetEraShort(), syst);
      }
    }
  }
}

void dybAnalyzer::executeEventWithParameter(TString channel, TString option){

  lepton0 = NULL;
  lepton1 = NULL;
  jet0 = NULL;
  bcharge = 0;
  prefix = channel+"/"+gprefix, hprefix = "", suffix = "";

  if(option.Contains("jet_scale_up")) suffix += "_jet_scale_up";
  else if(option.Contains("jet_scale_down")) suffix += "_jet_scale_down";
  else if(option.Contains("jet_smear_up")) suffix += "_jet_smear_up";
  else if(option.Contains("jet_smear_down")) suffix += "_jet_smear_down";
  if(IsNominalRun || option != "") IsNominalLike = true;
  else IsNominalLike = false;
  map_weight.clear();

  // Weights Setup
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

  // Dilepton + pT + OS + Mass
  if(!HasDileptons(channel)) return;

  // Jets
  vector<Jet> alljets = {};
  if(option.Contains("jet_scale_up")) alljets = SelectJets(ScaleJets(GetAllJets(), 1), "tightLepVeto", 20, 2.4);
  else if(option.Contains("jet_scale_down")) alljets = SelectJets(ScaleJets(GetAllJets(), -1), "tightLepVeto", 20, 2.4);
  else if(option.Contains("jet_smear_up")) alljets = SelectJets(SmearJets(GetAllJets(), 1), "tightLepVeto", 20, 2.4);
  else if(option.Contains("jet_smear_down")) alljets = SelectJets(SmearJets(GetAllJets(), -1), "tightLepVeto", 20, 2.4);
  else alljets = SelectJets(GetAllJets(), "tightLepVeto", 20, 2.4);
  std::sort(alljets.begin(), alljets.end(), PtComparing);

  vector<Jet> lepvetojets = {}, realjets = {}, bjets = {}, ajets = {};
  for(const auto jet:alljets){
    if(lepton0 && jet.DeltaR(*lepton0) < 0.4) continue;
    if(lepton1 && jet.DeltaR(*lepton1) < 0.4) continue;
    lepvetojets.push_back(jet);
  }
  for(const auto jet:lepvetojets){
    if(!PUJetIDPass(jet, "Loose")) continue;
    realjets.push_back(jet);
  }

  // b-tagging
  JetTagging::Parameters DeepJet_Tight = JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Tight, JetTagging::incl, JetTagging::comb);
  JetTagging::Parameters DeepJet_Loose = JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Loose, JetTagging::incl, JetTagging::comb);

  for(const auto& jet:realjets){
    if(jet.Pt() > 30 && jet.GetTaggerResult(DeepJet_Tight.j_Tagger) > mcCorr->GetJetTaggingCutValue(DeepJet_Tight.j_Tagger, DeepJet_Tight.j_WP)) bjets.push_back(jet);
    else if(jet.Pt() > 20 && jet.GetTaggerResult(DeepJet_Loose.j_Tagger) > mcCorr->GetJetTaggingCutValue(DeepJet_Loose.j_Tagger, DeepJet_Loose.j_WP)) ajets.push_back(jet);
  }

  // Jet related weights
  if(!IsDATA){
    pujetSF = GetPUJetWeight(lepvetojets, "Loose", 0);
    btagSF = GetBTaggingReweight_1a_2WP(realjets, DeepJet_Tight, DeepJet_Loose);
  }

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
  map_weight[""] *= zptweight;
  if(IsNominalLike){
    FillHist(prefix+hprefix+"weight_Zpt"+suffix, zptweight, map_weight[""], 200,-5,5);
    FillCutflow(prefix+hprefix+"cutflow"+suffix, "Zpt", map_weight[""]);
  }
  map_weight[""] *= weakweight;
  if(IsNominalLike){
    FillHist(prefix+hprefix+"weight_Weak"+suffix, weakweight, map_weight[""], 200,-5,5);
    FillCutflow(prefix+hprefix+"cutflow"+suffix, "Weak", map_weight[""]);
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
    TString DZSF = "";
    if(channel.Contains("mm"+GetEraShort())){
      for(const Lepton* lepton:leptons){
        leptonTrackingSF *= fEff->GetEfficiencySF("Muon_Tracking", lepton, 0,0);
        leptonRECOSF *= fEff->GetEfficiencySF("Muon_RECO", lepton, 0,0);
        leptonIDSF *= fEff->GetEfficiencySF("Muon_MediumID_trkIsoLoose", lepton, 0,0);
      }
      if(GetEraShort() !="2016a") DZSF = "DZ_MediumID_trkIsoLoose";
      leptonTriggerSF *= GetDileptonTriggerSF("Mu17Leg1_MediumID_trkIsoLoose", "Mu8Leg2_MediumID_trkIsoLoose", DZSF, leptons, 0,0);
    }else if(channel.Contains("ee"+GetEraShort())){
      for(const Lepton* lepton:leptons){
        leptonRECOSF *= fEff->GetEfficiencySF("Electron_RECO", lepton, 0,0);
        leptonIDSF *= fEff->GetEfficiencySF("Electron_MediumID", lepton, 0,0);
      }
      if(GetEraShort().Contains("2016")) DZSF = "DZ_MediumID";
      leptonTriggerSF *= GetDileptonTriggerSF("Ele23Leg1_MediumID", "Ele12Leg2_MediumID", DZSF, leptons, 0,0);
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

  //==== Inclusive DY
  double dimass = (*lepton0 + *lepton1).M();
  double dirap = (*lepton0 + *lepton1).Rapidity();
  double dipt = (*lepton0 + *lepton1).Pt();
  double costhetaCS = GetCosThetaCS(lepton0, lepton1);

  FillHist(prefix+hprefix+"mll_IncDY"+suffix, dimass, map_weight[""], 80,70,110);
  FillHist(prefix+hprefix+"yll_IncDY"+suffix, dirap, map_weight[""], 96,-2.4,2.4);
  FillHist(prefix+hprefix+"ptll_IncDY"+suffix, dipt, map_weight[""], AFBAnalyzer::unfold_0bjet_ptbinnum_reco,AFBAnalyzer::unfold_0bjet_ptbin_reco);
  FillHist(prefix+hprefix+"ptl_IncDY"+suffix, lepton0->Pt(), map_weight[""], AFBAnalyzer::lptbinnum,AFBAnalyzer::lptbin);
  FillHist(prefix+hprefix+"etal_IncDY"+suffix, lepton0->Eta(), map_weight[""], 50,-2.5,2.5);
  FillHist(prefix+hprefix+"ptl_IncDY"+suffix, lepton1->Pt(), map_weight[""], AFBAnalyzer::lptbinnum,AFBAnalyzer::lptbin);
  FillHist(prefix+hprefix+"etal_IncDY"+suffix, lepton1->Eta(), map_weight[""], 50,-2.5,2.5);
  FillHist(prefix+hprefix+"costhetaCS_IncDY"+suffix, dimass, dirap, dipt, costhetaCS, map_weight, afb_mbinnum,(double*)afb_mbin, afb_ybinnum,(double*)afb_ybin, afb_ptbinnum,(double*)afb_ptbin, 20,-1,1);

  FillHist(prefix+hprefix+"nalljets_IncDY"+suffix, alljets.size(), map_weight[""], 15,0,15);
  FillHist(prefix+hprefix+"nlepvetojets_IncDY"+suffix, lepvetojets.size(), map_weight[""], 15,0,15);
  FillHist(prefix+hprefix+"nrealjets_IncDY"+suffix, realjets.size(), map_weight[""], 15,0,15);
  FillHist(prefix+hprefix+"nbjets_IncDY"+suffix, bjets.size(), map_weight[""], 10,0,10);

  if(bjets.size() != 1) return;
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "Tight1b", map_weight[""]);
  jet0 = &bjets.at(0);
  bcharge = jetCharge(*jet0);
  double costhetaRecoil = GetCosThetaRecoil(lepton0, lepton1, jet0);

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

  FillHist(prefix+hprefix+"mll_Tight1b"+suffix, dimass, map_weight[""], 80,70,110);
  FillHist(prefix+hprefix+"yll_Tight1b"+suffix, dirap, map_weight[""], 96,-2.4,2.4);
  FillHist(prefix+hprefix+"ptll_Tight1b"+suffix, dipt, map_weight[""], AFBAnalyzer::unfold_0bjet_ptbinnum_reco,AFBAnalyzer::unfold_0bjet_ptbin_reco);
  FillHist(prefix+hprefix+"ptl_Tight1b"+suffix, lepton0->Pt(), map_weight[""], AFBAnalyzer::lptbinnum,AFBAnalyzer::lptbin);
  FillHist(prefix+hprefix+"etal_Tight1b"+suffix, lepton0->Eta(), map_weight[""], 50,-2.5,2.5);
  FillHist(prefix+hprefix+"ptl_Tight1b"+suffix, lepton1->Pt(), map_weight[""], AFBAnalyzer::lptbinnum,AFBAnalyzer::lptbin);
  FillHist(prefix+hprefix+"etal_Tight1b"+suffix, lepton1->Eta(), map_weight[""], 50,-2.5,2.5);
  FillHist(prefix+hprefix+"bjetCharge_Tight1b"+suffix, bcharge, map_weight[""], 200,-5,5);
  FillHist(prefix+hprefix+"costhetaRecoil_Tight1b"+suffix, dimass, fabs(bcharge), jet0->Pt(), costhetaRecoil, map_weight, afb_mbinnum,(double*)afb_mbin, afb_chbinnum,(double*)afb_chbin, afb_ptbinnum,(double*)afb_ptbin, 20,-1,1);

  if(ajets.size() > 0) return;
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "Veto2b", map_weight[""]);

  int n_30jet = 0;
  for(unsigned int k=0; k<realjets.size(); k++){
    if(realjets.at(k).Pt() > 30) n_30jet += 1;
  }
  if(n_30jet > 1) return;
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "Veto2j", map_weight[""]);

  if(PuppiMET_Type1_pt > 75) return;
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "MET75", map_weight[""]);

  if(abs((*lepton0 + *lepton1).DeltaPhi(*jet0)) < 1.6) return;
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "ZbdPhi1p6", map_weight[""]);

  if((*lepton0 + *lepton1 + *jet0).Pt() > 60) return;
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "ZbpT60", map_weight[""]);

  if((*lepton0 + *lepton1).Pt() < 15) return;
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "ZpT15", map_weight[""]);

  //==== Weights of Systematics
  if(!IsDATA && HasFlag("SYS") && option == ""){
    // PU reweight
    map_weight["_noPUweight"] =  map_weight[""] / PUweight;
    map_weight["_PUweight_up"] =  map_weight[""] / PUweight * GetPileUpWeight(nPileUp, 1);
    map_weight["_PUweight_down"] = map_weight[""] / PUweight * GetPileUpWeight(nPileUp, -1);

    // b-tagging SF
    map_weight["_nobtagSF"] =  map_weight[""] / btagSF;
    map_weight["_btagSF_hup"] = map_weight[""] / btagSF * GetBTaggingReweight_1a_2WP(realjets, DeepJet_Tight, DeepJet_Loose, "SystUpHTag");
    map_weight["_btagSF_hdown"] = map_weight[""] / btagSF * GetBTaggingReweight_1a_2WP(realjets, DeepJet_Tight, DeepJet_Loose, "SystDownHTag");
    map_weight["_btagSF_hcorr"] = map_weight[""] / btagSF * GetBTaggingReweight_1a_2WP(realjets, DeepJet_Tight, DeepJet_Loose, "SystUpHTagCorr");
    map_weight["_btagSF_huncorr"] = map_weight[""] / btagSF * GetBTaggingReweight_1a_2WP(realjets, DeepJet_Tight, DeepJet_Loose, "SystUpHTagUnCorr");
    map_weight["_btagSF_lup"] = map_weight[""] / btagSF * GetBTaggingReweight_1a_2WP(realjets, DeepJet_Tight, DeepJet_Loose, "SystUpLTag");
    map_weight["_btagSF_ldown"] = map_weight[""] / btagSF * GetBTaggingReweight_1a_2WP(realjets, DeepJet_Tight, DeepJet_Loose, "SystDownLTag");
    map_weight["_btagSF_lcorr"] = map_weight[""] / btagSF * GetBTaggingReweight_1a_2WP(realjets, DeepJet_Tight, DeepJet_Loose, "SystUpLTagCorr");;
    map_weight["_btagSF_luncorr"] = map_weight[""] / btagSF * GetBTaggingReweight_1a_2WP(realjets, DeepJet_Tight, DeepJet_Loose, "SystUpLTagUnCorr");
  }else if(!IsDATA && MCSample.Contains("MiNNLO") && IsNominalRun){
    for(unsigned int i=0;i<weight_sthw2->size();i++) map_weight[Form("_sthw2_%d",i)] = map_weight[""] * weight_sthw2->at(i);
  }
  if(HasFlag("SYS") && option == "") map_weight.erase("");

  FillHist(prefix+hprefix+"mll"+suffix, dimass, map_weight, 80,70,110);
  FillHist(prefix+hprefix+"yll"+suffix, dirap, map_weight, 96,-2.4,2.4);
  FillHist(prefix+hprefix+"ptll"+suffix, dipt, map_weight, AFBAnalyzer::unfold_nbjet_ptbinnum_reco,AFBAnalyzer::unfold_nbjet_ptbin_reco);
  FillHist(prefix+hprefix+"ptl"+suffix, lepton0->Pt(), map_weight, AFBAnalyzer::lptbinnum,AFBAnalyzer::lptbin);
  FillHist(prefix+hprefix+"etal"+suffix, lepton0->Eta(), map_weight, 50,-2.5,2.5);
  FillHist(prefix+hprefix+"ptl"+suffix, lepton1->Pt(), map_weight, AFBAnalyzer::lptbinnum,AFBAnalyzer::lptbin);
  FillHist(prefix+hprefix+"etal"+suffix, lepton1->Eta(), map_weight, 50,-2.5,2.5);
  FillHist(prefix+hprefix+"bjetCharge"+suffix, bcharge, map_weight, 200,-5,5);
  FillHist(prefix+hprefix+"bjetChargeEasy"+suffix, bcharge<0?0:1, map_weight, 2,0,2);
  for(unsigned int i=1; i<afb_chbinnum+1; i++){
    if(afb_chbin[i-1] < abs(bcharge) && abs(bcharge) < afb_chbin[i]){
      FillHist(Form(prefix+"bjetCharge%d"+suffix,i-1), bcharge, map_weight, 200,-5,5);
      FillHist(Form(prefix+"bjetCharge%dEasy"+suffix,i-1), bcharge<0?0:1, map_weight, 2,0,2);
    }
  }
  FillHist(prefix+hprefix+"costhetaRecoil"+suffix, dimass, fabs(bcharge), jet0->Pt(), costhetaRecoil, map_weight, afb_mbinnum,(double*)afb_mbin, afb_chbinnum,(double*)afb_chbin, afb_ptbinnum,(double*)afb_ptbin, 20,-1,1);
}

bool dybAnalyzer::HasDileptons(TString channel){
  double l0pt = 20., l1pt = 10.;
  if(channel.Contains("mm"+GetEraShort())){
    muons = MuonMomentumCorrection(SMPGetMuons("POGMediumWithLooseTrkIso", 8.0,2.4), 0,0,0);
    //muons = SMPGetMuons("POGMediumWithLooseTrkIso", 8.0,2.4);
    if(muons.size() > 0) lepton0 = &muons.at(0);
    if(muons.size() > 1) lepton1 = &muons.at(1);
  }else if(channel.Contains("ee"+GetEraShort())){
    l0pt = 25.;
    l1pt = 15.;
    electrons = ElectronEnergyCorrection(SMPGetElectrons("passMediumID", 8.0,2.5), 0,0);
    //electrons = SMPGetElectrons("passMediumID", 8.0,2.5);
    if(electrons.size() > 0) lepton0 = &electrons.at(0);
    if(electrons.size() > 1) lepton1 = &electrons.at(1);
  }else{
    cout<<"[dybAnalyzer::Hasleptons] channel="<<channel<<" is weird"<<endl;
  }

  if(!lepton0 || !lepton1) return false;                         // Dilepton
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "Dilepton", map_weight[""]);
  if(lepton0->Pt() < l0pt || lepton1->Pt() < l1pt) return false; // pT cut
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "LepPt", map_weight[""]);
  if(lepton0->Charge() * lepton1->Charge() > 0) prefix += "ss_"; // Opposite charge
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "Charge", map_weight[""]);
  if((*lepton0 + *lepton1).M() < 52) return false;               // Mass 52
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "Mass52", map_weight[""]);

  leptons ={};
  leptons.push_back(lepton0);
  leptons.push_back(lepton1);

  return true;
}
bool dybAnalyzer::IsFiredTriggers(TString channel){
  vector<TString> triggers = {};
  // Dilepton for DY
  if(DataYear == 2016 && channel.Contains("mm")){
    triggers = {
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v",
      "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",
      "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
      "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
      "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"
    };
  }else if(DataYear == 2017 && channel.Contains("mm")) triggers = {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v"};
  else if(DataYear == 2018 && channel.Contains("mm")) triggers = {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"};
  else if(DataYear == 2016 && channel.Contains("ee")) triggers = {"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"};
  else if(DataYear >= 2017 && channel.Contains("ee")) triggers = {"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"};

  // Single lepton for TTLJ
  else if(DataYear == 2016 && channel.Contains("m"+GetEraShort())){
    triggers = {
      "HLT_IsoMu24_v",
      "HLT_IsoTkMu24_v"
    };
  }else if(DataYear == 2017 && channel.Contains("m"+GetEraShort())) triggers = {"HLT_IsoMu27_v"};
  else if(DataYear == 2018 && channel.Contains("m"+GetEraShort())) triggers = {"HLT_IsoMu24_v"};
  else if(DataYear == 2016 && (channel.Contains("e"+GetEraShort()) || channel.Contains("E")+GetEraShort())) triggers = {"HLT_Ele27_WPTight_Gsf_v"};
  else if(DataYear >= 2017 && (channel.Contains("e"+GetEraShort()) || channel.Contains("E")+GetEraShort())){
    triggers = {
      "HLT_Ele27_WPTight_Gsf_v",
      "HLT_Ele32_WPTight_Gsf_v"
    };
  }else{
    cout<<"[dybAnalyzer::IsFiredTriggers] something channel is wrong"<<endl;
    exit(EXIT_FAILURE);
  }

  return _event.PassTrigger(triggers);
}
double dybAnalyzer::jetCharge(const Jet& jet){
  double jetCharge = jet.Charge();
  vector<Muon> allmus = GetAllMuons(); // Need RoccoR too?
  std::sort(allmus.begin(),allmus.end(),PtComparing);
  vector<Electron> allels = GetAllElectrons(); // Need Aepcor too?
  std::sort(allels.begin(),allels.end(),PtComparing);

  // Selections from Hyonsan's AN (AN-20-216 v3)
  vector<Muon> bmuon;
  for(const auto& mu:allmus){
    if(jet.DeltaR(mu) > 0.3) continue;
    if(mu.P() * sin(mu.Angle(jet.Vect())) < 1.0) continue;
    if(mu.IsType(5)) bmuon.push_back(mu); // PFMuon
  }

  vector<Electron> belectron;
  for(const auto& el:allels){
    if(jet.DeltaR(el) > 0.3) continue;
    if(el.P() * sin(el.Angle(jet.Vect())) < 1.0) continue;
    if(fabs(el.scEta()) <= 1.479){
      if(el.Full5x5_sigmaIetaIeta() > 0.0126) continue ;
      if(fabs(el.dEtaSeed()) > 0.00463) continue ;
      if(fabs(el.dPhiIn()) > 0.148) continue;
    }else{
      if(el.Full5x5_sigmaIetaIeta() > 0.0457) continue ;
      if(fabs(el.dEtaSeed()) > 0.00814) continue ;
      if(fabs(el.dPhiIn()) > 0.19) continue;
    }
    if(!el.PassConversionVeto()) continue;
    if(el.IsGsfCtfScPixChargeConsistent()) belectron.push_back(el);
  }

  //The jet has soft muon inside, and its charge will determine the jet charge
  if(bmuon.size() > 0) jetCharge += 2 * bmuon.at(0).Charge();
  else if(belectron.size() > 0) jetCharge += 4 * belectron.at(0).Charge();

  return jetCharge;
}

dybAnalyzer::dybAnalyzer(){}
dybAnalyzer::~dybAnalyzer(){}

// From Hyonsan's functions in SMPAnalyzerCore
void dybAnalyzer::executeEventGen(){
  gprefix = "";
  if(IsData) return;

  if(IsDYSample){
    // LHE Setting
    vector<LHE> lhes = GetLHEs();
    LHE lhe_l0, lhe_l1;
    for(int i=0; i<(int)lhes.size(); i++){
      if(lhe_l0.ID() == 0 && (abs(lhes[i].ID()) == 11 || abs(lhes[i].ID()) == 13 || abs(lhes[i].ID()) == 15)) lhe_l0 = lhes[i];
      if(lhe_l0.ID()      && (abs(lhes[i].ID()) == 11 || abs(lhes[i].ID()) == 13 || abs(lhes[i].ID()) == 15)) lhe_l1 = lhes[i];
    }
    if(lhe_l0.ID() == 0 || lhe_l1.ID() == 0){
      cout <<"[dybAnalyzer::executeEventGen] something lhe is wrong"<<endl;
      exit(EXIT_FAILURE);
    }
    if(lhe_l0.Pt() < lhe_l1.Pt()){
      LHE temp = lhe_l0;
      lhe_l0 = lhe_l1;
      lhe_l1 = temp;
    }

    // GEN Setting
    vector<Gen> gens = GetGens();
    Gen gen_l0, gen_l1;
    vector<const Gen*> leptons;
    vector<const Gen*> photons;
    for(int i=0; i<(int)gens.size(); i++){
      if(!gens.at(i).isPrompt()) continue;
      if(gens.at(i).Status() == 1){
        if(abs(gens.at(i).PID()) == 11 || abs(gens.at(i).PID()) == 13) leptons.push_back(&gens[i]);
        else if(gens.at(i).PID() == 22) photons.push_back(&gens[i]);
      }
    }
    const double maxdr=0.4;
    for(int i=0; i<(int)leptons.size(); i++){
      if(leptons[i]->PID() != lhe_l0.ID()) continue;
      if(leptons[i]->DeltaR(lhe_l0) > maxdr) continue;
      if(fabs(leptons[i]->E()-lhe_l0.E()) < fabs(gen_l0.E()-lhe_l0.E())) gen_l0 = *leptons[i];
    }
    if(gen_l0.PID() == 0){
      for(int i=0; i<(int)leptons.size(); i++){
        if(leptons[i]->PID() != lhe_l0.ID()) continue;
        if(leptons[i]->DeltaR(lhe_l0) < gen_l0.DeltaR(lhe_l0)) gen_l0 = *leptons[i];
      }
    }
    for(int i=0; i<(int)leptons.size(); i++){
      if(leptons[i]->PID() != lhe_l1.ID()) continue;
      if(leptons[i]->DeltaR(lhe_l1) > maxdr) continue;
      if(fabs(leptons[i]->E()-lhe_l1.E()) < fabs(gen_l1.E()-lhe_l1.E())) gen_l1 = *leptons[i];
    }
    if(gen_l1.PID() == 0){
      for(int i=0; i<(int)leptons.size(); i++){
        if(leptons[i]->PID() != lhe_l1.ID()) continue;
        if(leptons[i]->DeltaR(lhe_l1) < gen_l1.DeltaR(lhe_l1)) gen_l1 = *leptons[i];
      }
    }
    if(gen_l0.Pt() < gen_l1.Pt()){
      Gen tmp = gen_l0;
      gen_l0 = gen_l1;
      gen_l1 = tmp;
    }
    // dressing Gens
    if(leptons.size() >= 4){
      for(int i=0; i<(int)leptons.size(); i++){
        if(leptons[i]->Index() == gen_l0.Index() || leptons[i]->Index() == gen_l1.Index()) continue;
        for(int j=i+1; j<(int)leptons.size(); j++){
          if(leptons[j]->Index() == gen_l0.Index() || leptons[j]->Index() == gen_l1.Index()) continue;
          if(leptons[i]->PID() + leptons[j]->PID() != 0) continue;
          vector<int> history_i = TrackGenSelfHistory(*leptons[i],gens);
          vector<int> history_j = TrackGenSelfHistory(*leptons[j],gens);
          if(history_i.at(1) == history_j.at(1)){ // from the same mother(Z or lep)
            photons.push_back(leptons[i]);
            photons.push_back(leptons[j]);
          }
        }
      }
      for(const auto& photon:photons){
        vector<int> history = TrackGenSelfHistory(*photon,gens);
        if(gens[history.at(1)].PID() == gen_l0.PID()) gen_l0 += *photon;
        else if(gens[history.at(1)].PID() == gen_l1.PID()) gen_l1 += *photon;
        else if(gens[history.at(1)].PID() == 23){ // for minnlo+photos
          if(photon->DeltaR(gen_l0) < photon->DeltaR(gen_l1)) gen_l0 += *photon;
          else gen_l1 += *photon;
        }
      }    
    }

    if(abs(lhe_l0.ID()) == 11 || abs(lhe_l0.ID()) == 13){
      TLorentzVector genZ = (gen_l0 + gen_l1);
      zptweight = fZptCorrection->GetZptWeight(genZ.Pt(),genZ.Rapidity(),genZ.M());
      weakweight = GetDYWeakWeight(genZ.M());

      // Only qqbar collisions (LO DY)
      if(lhe_p0.ID() + lhe_p1.ID() == 0) gprefix += "";
      // Only qG collisions (NLO DY)
      else if((abs(lhe_p0.ID()) <= 5 && lhe_p1.ID() == 21) || (lhe_p0.ID() == 21 && abs(lhe_p1.ID()) <= 5)){
        if(lhe_p0.ID() == 5 || lhe_p1.ID() == 5) gprefix += "dyb_";
        else if(lhe_p0.ID() == -5 || lhe_p1.ID() == -5) gprefix += "dybbar_";
        else if(lhe_p0.ID() == 4 || lhe_p1.ID() == 4) gprefix += "dyc_";
        else if(lhe_p0.ID() == -4 || lhe_p1.ID() == -4) gprefix += "dycbar_";
        else gprefix += "dyudsg_";
      } // GG collisions or bq, cq collisions (NNLO DY) - find the heavy flavor parton with highest-pt within accptance
      else if((lhe_p0.ID()==21 && lhe_p1.ID()==21) || (abs(lhe_p0.ID())==4 || abs(lhe_p0.ID())==5 || abs(lhe_p1.ID())==4 || abs(lhe_p1.ID())==5)){
        Gen heavyparton = gens.at(0);
        int nheavyparton = 0;
        for(unsigned int i=0; i<gens.size(); i++){
          if(!gens.at(i).isHardProcess()) continue;
          if(abs(gens.at(i).PID()) >= 11 && abs(gens.at(i).PID()) <= 16) continue; // No Lepton
          if(gens.at(i).PID() == 22 || gens.at(i).PID() == 23) continue; // No Gamma, Z
          if(gens.at(i).Pt() < 30 || abs(gens.at(i).Eta()) > 2.4) continue; // In the acceptance

          if(nheavyparton == 0 && (abs(gens.at(i).PID()) == 4 || abs(gens.at(i).PID()) == 5)){
            heavyparton = gens.at(i);
            nheavyparton++;
            continue;
          }else if(nheavyparton > 0 && (abs(gens.at(i).PID()) == 4 || abs(gens.at(i).PID()) == 5)){
            heavyparton = (heavyparton.Pt() > gens.at(i).Pt()? heavyparton: gens.at(i));
            nheavyparton++;
            continue;//break;
          }
        }
        if(nheavyparton > 0 && heavyparton.PID() == 5) gprefix += "dyb_";
        else if(nheavyparton > 0 && heavyparton.PID() == -5) gprefix += "dybbar_";
        else if(nheavyparton > 0 && heavyparton.PID() == 4) gprefix += "dyc_";
        else if(nheavyparton > 0 && heavyparton.PID() == -4) gprefix += "dycbar_";
        else gprefix += "";
      }
    }else gprefix += "tau_";
  }
  if(IsTTSample) topptweight=mcCorr->GetTopPtReweight(gens);
}

double dybAnalyzer::GetCosThetaCS(const Particle *p0,const Particle *p1,int direction){
  if(!p0||!p1) return 0.;
  const TLorentzVector *l0,*l1;
  if(p0->Charge()<0&&p1->Charge()>0){
    l0=p0;
    l1=p1;
  }else if(p0->Charge()>0&&p1->Charge()<0){
    l0=p1;
    l1=p0;
  }else if(strcmp(p0->ClassName(),"LHE")==0){ 
    if(((LHE*)p0)->ID()>0&&((LHE*)p1)->ID()<0){
      l0=p0;
      l1=p1;
    }else if(((LHE*)p0)->ID()<0&&((LHE*)p1)->ID()>0){
      l0=p1;
      l1=p0;
    }else{
      gRandom->SetSeed((run<<15)+(lumi<<10)+(event<<5)+p0->Eta()*100);
      if(gRandom->Rndm()<0.5){
	l0=p0;
	l1=p1;
      }else{
	l0=p1;
	l1=p0;
      }      
    } 
  }else{
    gRandom->SetSeed((run<<15)+(lumi<<10)+(event<<5)+p0->Eta()*100);
    if(gRandom->Rndm()<0.5){
      l0=p0;
      l1=p1;
    }else{
      l0=p1;
      l1=p0;
    }      
  }

  TLorentzVector dilepton=*l0+*l1;
  double l0pp=(l0->E()+l0->Pz())/sqrt(2);
  double l0pm=(l0->E()-l0->Pz())/sqrt(2);
  double l1pp=(l1->E()+l1->Pz())/sqrt(2);
  double l1pm=(l1->E()-l1->Pz())/sqrt(2);
  double dimass=dilepton.M();
  double dipt=dilepton.Pt();
  if(direction==0) direction=dilepton.Pz()>0?1:-1;
  return direction*2*(l0pp*l1pm-l0pm*l1pp)/sqrt(dimass*dimass*(dimass*dimass+dipt*dipt));
}
double dybAnalyzer::GetCosThetaRecoil(const Particle *p0, const Particle *p1, Particle *b, int mode){
  if(!p0||!p1) return 0.;
  const TLorentzVector *lm,*lp;
  if(p0->Charge()<0&&p1->Charge()>0){
    lm=p0;
    lp=p1;
  }else if(p0->Charge()>0&&p1->Charge()<0){
    lm=p1;
    lp=p0;
  }else if(strcmp(p0->ClassName(),"LHE")==0){ 
    if(((LHE*)p0)->ID()>0&&((LHE*)p1)->ID()<0){
      lm=p0;
      lp=p1;
    }else if(((LHE*)p0)->ID()<0&&((LHE*)p1)->ID()>0){
      lm=p1;
      lp=p0;
    }else{
      if(gRandom->Rndm()<0.5){
        lm=p0;
        lp=p1;
      }else{
	lm=p1;
	lp=p0;
      }      
    } 
  }else{
    if(gRandom->Rndm()<0.5){
      lm=p0;
      lp=p1;
    }else{
      lm=p1;
      lp=p0;
    }      
  }
  int direction=0;
  if(bcharge>0) direction = -1;
  else direction=1;
  if(mode==0){
    return direction*((*lm-*lp)*(*b))/((*lm+*lp)*(*b));
  }else{
    TLorentzVector b_m0;
    b_m0.SetPtEtaPhiM(b->Pt(),b->Eta(),b->Phi(),0);
    b_m0*=b->E()/b_m0.E();
    return direction*((*lm-*lp)*(b_m0))/((*lm+*lp)*(b_m0));
  }
}

// Functions in MY SMPAnalyzerCore
double dybAnalyzer::GetBTaggingReweight_1a_2WP(const vector<Jet>& jets, JetTagging::Parameters jtpT, JetTagging::Parameters jtpL, string Syst){
  //Syst. usage ex.: "SystUpHTag"(all component variation for heavy flav(b,c).),
  //                 "SystUpHTagCorr"(variation of heavy flav(b,c) sf only for yearly correlated components)
  //change H->L for light flav., Up->Down for downward variation, Corr->UnCorr for yearly independent components

  if(IsDATA) return 1.;

  TString SystStr(Syst);
  double Prob_MC(1.), Prob_DATA(1.), SF(1.);
  bool Syst_HTag=false, Syst_LTag=false; int SystDir=0, CorrType=0;
  string SystKey;
  if(SystStr.Contains("Syst")){
    if     (SystStr.Contains("HTag")) Syst_HTag=true;
    else if(SystStr.Contains("LTag")) Syst_LTag=true;
    if     (SystStr.Contains("Up")  ) SystDir= 1;
    else if(SystStr.Contains("Down")) SystDir=-1;
    if     (SystStr.Contains("UnCorr")) CorrType=-1;
    else if(SystStr.Contains("Corr"))   CorrType= 1;
    if(SystDir==0){ cout<<"SystStr in not correct form"<<endl; exit(ENODATA); }
    if(!(Syst_HTag or Syst_LTag)){ cout<<"SystMode but no H/L mode assigned"<<endl; exit(ENODATA); }
  }

  for(unsigned int i=0; i<jets.size(); i++){
    int JetHadFlav = jets.at(i).hadronFlavour();
    bool ApplySyst=false;
    if     (Syst_HTag && (JetHadFlav==4 or JetHadFlav==5)){ ApplySyst=true; }
    else if(Syst_LTag && (JetHadFlav==0                 )){ ApplySyst=true; }

    if     (ApplySyst && CorrType==0) SystKey=SystDir>0? "up":"down";
    else if(ApplySyst && CorrType >0) SystKey=SystDir>0? "up_correlated":"down_correlated";
    else if(ApplySyst && CorrType <0) SystKey=SystDir>0? "up_uncorrelated":"down_uncorrelated";
    else                              SystKey="central";

    double this_MC_EffT = mcCorr->GetMCJetTagEff(jtpT.j_Tagger, jtpT.j_WP, jets.at(i).hadronFlavour(), jets.at(i).Pt(), jets.at(i).Eta());
    double this_MC_EffL = mcCorr->GetMCJetTagEff(jtpL.j_Tagger, jtpL.j_WP, jets.at(i).hadronFlavour(), jets.at(i).Pt(), jets.at(i).Eta());
    double this_SFT = mcCorr->GetJetTaggingSF(jtpT,
                                              jets.at(i).hadronFlavour(),
                                              jets.at(i).Pt(),
                                              jets.at(i).Eta(),
                                              jets.at(i).GetTaggerResult(jtpT.j_Tagger),
                                              SystKey );
    double this_SFL = mcCorr->GetJetTaggingSF(jtpL,
                                              jets.at(i).hadronFlavour(),
                                              jets.at(i).Pt(),
                                              jets.at(i).Eta(),
                                              jets.at(i).GetTaggerResult(jtpL.j_Tagger),
                                              SystKey );
    double this_DATA_EffT = this_MC_EffT*this_SFT;
    double this_DATA_EffL = this_MC_EffL*this_SFL;

    bool isTaggedT = jets.at(i).GetTaggerResult(jtpT.j_Tagger) > mcCorr->GetJetTaggingCutValue(jtpT.j_Tagger, jtpT.j_WP);
    bool isTaggedL = jets.at(i).GetTaggerResult(jtpL.j_Tagger) > mcCorr->GetJetTaggingCutValue(jtpL.j_Tagger, jtpL.j_WP);
    if(isTaggedT){
      Prob_MC *= this_MC_EffT;
      Prob_DATA *= this_DATA_EffT;
    }
    else if(isTaggedL){
      if(this_MC_EffL == this_MC_EffT) this_MC_EffL += 1E-10;
      Prob_MC *= this_MC_EffL - this_MC_EffT;
      Prob_DATA *= this_DATA_EffL - this_DATA_EffT;
    }
    else{
      Prob_MC *= 1.-this_MC_EffL;
      Prob_DATA *= 1.-this_DATA_EffL;
    }
  }

  if(Prob_MC>0. && Prob_DATA>0.) SF=Prob_DATA/Prob_MC;
  else SF=0.;

  return SF;
}

void dybAnalyzer::SetupPUJetWeight(){
  TString datapath = getenv("DATA_DIR");
  TFile fPUID(datapath+"/"+GetEra()+"/ID/PUJet/PUID.root");
  vector<TString> IDs = {"T", "M", "L"};
  for(unsigned int i=0; i<IDs.size(); i++){
    cout<<"[dybAnalyzer::SetupPUJetWeight] setting PUJetWeight with ID : "+IDs.at(i)<<endl;

    TString era = GetEra();
    if(era == "2016postVFP") era = "2016";
    else if(era == "2016preVFP") era = "2016APV";

    heff_data = (TH2F*)fPUID.Get("h2_eff_dataUL"+era+"_"+IDs.at(i));
    heff_mc   = (TH2F*)fPUID.Get("h2_eff_mcUL"+era+"_"+IDs.at(i));
    hmistag_data = (TH2F*)fPUID.Get("h2_mistag_dataUL"+era+"_"+IDs.at(i));
    hmistag_mc   = (TH2F*)fPUID.Get("h2_mistag_mcUL"+era+"_"+IDs.at(i));

    heff_data->SetDirectory(0);
    heff_mc->SetDirectory(0);
    hmistag_data->SetDirectory(0);
    hmistag_mc->SetDirectory(0);
  }

  fPUID.Close();
}

bool dybAnalyzer::PUJetIDPass(Jet jet, TString ID){
  if(jet.Pt() >= 50) return true;

  if(DataEra.Contains("2016")){
    if(ID == "Tight"){
      if(jet.Pt() >= 40){
        if(jet.PileupJetId() > 0.97) return true;
      }else if(jet.Pt() >= 30){
        if(jet.PileupJetId() > 0.94) return true;
      }else if(jet.Pt() >= 20){
        if(jet.PileupJetId() > 0.87) return true;
      }
    }
    else if(ID == "Medium"){
      if(jet.Pt() >= 40){
        if(jet.PileupJetId() > 0.93) return true;
      }else if(jet.Pt() >= 30){
        if(jet.PileupJetId() > 0.86) return true;
      }else if(jet.Pt() >= 20){
        if(jet.PileupJetId() > 0.62) return true;
      }
    }
    else{
      if(jet.Pt() >= 40){
        if(jet.PileupJetId() > -0.42) return true;
      }else if(jet.Pt() >= 30){
        if(jet.PileupJetId() > -0.71) return true;
      }else if(jet.Pt() >= 20){
        if(jet.PileupJetId() > -0.90) return true;
      }
    }
  }

  if(DataEra == "2017"){
    if(ID == "Tight"){
      if(jet.Pt() >= 40){
        if(jet.PileupJetId() > 0.98) return true;
      }else if(jet.Pt() >= 30){
        if(jet.PileupJetId() > 0.96) return true;
      }else if(jet.Pt() >= 20){
        if(jet.PileupJetId() > 0.90) return true;
      }
    }
    else if(ID == "Medium"){
      if(jet.Pt() >= 40){
        if(jet.PileupJetId() > 0.96) return true;
      }else if(jet.Pt() >= 30){
        if(jet.PileupJetId() > 0.90) return true;
      }else if(jet.Pt() >= 20){
        if(jet.PileupJetId() > 0.68) return true;
      }
    }
    else{
      if(fabs(jet.Eta()) < 2.5){
        if(jet.Pt() >= 40){
          if(jet.PileupJetId() > -0.19) return true;
        }else if(jet.Pt() >= 30){
          if(jet.PileupJetId() > -0.63) return true;
        }else if(jet.Pt() >= 20){
          if(jet.PileupJetId() > -0.88) return true;
        }
      }else if(fabs(jet.Eta()) < 2.75){
        if(jet.Pt() >= 40){
          if(jet.PileupJetId() > 0.22) return true;
        }else if(jet.Pt() >= 30){
          if(jet.PileupJetId() > -0.18) return true;
        }else if(jet.Pt() >= 20){
          if(jet.PileupJetId() > -0.55) return true;
        }
      }else if(fabs(jet.Eta()) < 3.0){
        if(jet.Pt() >= 40){
          if(jet.PileupJetId() > -0.13) return true;
        }else if(jet.Pt() >= 30){
          if(jet.PileupJetId() > -0.43) return true;
        }else if(jet.Pt() >= 20){
          if(jet.PileupJetId() > -0.60) return true;
        }
      }else{
        if(jet.Pt() >= 40){
          if(jet.PileupJetId() > -0.03) return true;
        }else if(jet.Pt() >= 30){
          if(jet.PileupJetId() > -0.24) return true;
        }else if(jet.Pt() >= 20){
          if(jet.PileupJetId() > -0.43) return true;
        }
      }
    }
  }

  if(DataEra == "2018"){
    if(ID == "Tight"){
      if(jet.Pt() >= 40){
        if(jet.PileupJetId() > 0.98) return true;
      }else if(jet.Pt() >= 30){
        if(jet.PileupJetId() > 0.96) return true;
      }else if(jet.Pt() >= 20){
        if(jet.PileupJetId() > 0.90) return true;
      }
    }
    else if(ID == "Medium"){
      if(jet.Pt() >= 40){
        if(jet.PileupJetId() > 0.96) return true;
      }else if(jet.Pt() >= 30){
        if(jet.PileupJetId() > 0.90) return true;
      }else if(jet.Pt() >= 20){
        if(jet.PileupJetId() > 0.68) return true;
      }
    }
    else{
      if(jet.Pt() >= 40){
        if(jet.PileupJetId() > -0.19) return true;
      }else if(jet.Pt() >= 30){
        if(jet.PileupJetId() > -0.63) return true;
      }else if(jet.Pt() >= 20){
        if(jet.PileupJetId() > -0.88) return true;
      }
    }
  }
  return false;
}

double dybAnalyzer::GetPUJetWeight(const vector<Jet>& jets, TString ID, int sys){
  sys = 0;
  if(IsDATA) return 1.;

  vector<Gen> gens=GetGens();

  double Prob_MC(1.), Prob_DATA(1.);
  for(unsigned int i=0; i<jets.size(); i++){
    double jetpt = jets.at(i).Pt();
    double jeteta = jets.at(i).Eta();
    if(jets.at(i).Pt() < 20) cout<<"jet pt < 20GeV, something wrong"<<endl;;
    if(jets.at(i).Pt() > 50) continue;
    if(abs(jets.at(i).Eta()) > 2.5) continue;

    double this_DATA_eff = heff_data->GetBinContent(heff_data->FindBin(jetpt, jeteta));
    double this_MC_eff = heff_mc->GetBinContent(heff_mc->FindBin(jetpt, jeteta));
    double this_DATA_mistag = hmistag_data->GetBinContent(hmistag_data->FindBin(jetpt, jeteta));
    double this_MC_mistag = hmistag_mc->GetBinContent(hmistag_mc->FindBin(jetpt, jeteta));
    if(this_DATA_eff * this_MC_eff * this_DATA_mistag * this_MC_mistag == 0.) continue;

    bool isRealJet = false;
    isRealJet = (jets.at(i).GenHFHadronMatcherFlavour() >= 0.);
    bool isPassID = PUJetIDPass(jets.at(i), ID);

    if(isRealJet){
      if(isPassID){
        if(this_MC_eff == 0) this_MC_eff += 1E-4;
        Prob_DATA *= this_DATA_eff;
        Prob_MC *= this_MC_eff;
      }else{
        if(this_MC_eff == 1) this_MC_eff -= 1E-4;
        Prob_DATA *= 1.-this_DATA_eff;
        Prob_MC *= 1.-this_MC_eff;
      }
    }else{
      if(isPassID){
        if(this_MC_mistag == 0) this_MC_mistag += 1E-4;
        Prob_DATA *= this_DATA_mistag;
        Prob_MC *= this_MC_mistag;
      }else{
        if(this_MC_mistag == 1) this_MC_mistag -= 1E-4;
        Prob_DATA *= 1.-this_DATA_mistag;
        Prob_MC *= 1.-this_MC_mistag;
      }
    }
  }

  return Prob_DATA/Prob_MC;
}

double dybAnalyzer::GetbChargeSFWeight(const vector<Jet>& jets, int mode, int sys){
  sys = 0;
  double weight = 1.;
  if(IsDATA) return weight;

  vector<Gen> gens=GetGens();
  for(const auto& jet:jets){
    int genpid = 0;
    double dR = 99.;
    for(int i=0; i<gens.size(); i++){
      if(!gens.at(i).isPrompt()) continue;
      if(!gens.at(i).isHardProcess()) continue;
      if(fabs(gens.at(i).PID()) != 5) continue;
      if(gens.at(i).DeltaR(jet) > dR) continue;
      dR = gens.at(i).DeltaR(jet);
      if(dR < 0.4) genpid = gens.at(i).PID();
    }
    if(genpid == 0) continue; // No matched b-partons

    double Charge = jetCharge(jet);
    double alpha_plus_DATA_eff = 0.631064;
    double alpha_minus_DATA_eff = 0.611111;
    double alpha_plus_MC_eff = 0.655499;
    double alpha_minus_MC_eff = 0.637514;
    if(mode == 1){ // 1D bChargeSF
      if(fabs(Charge) < afb_chbin[1]){// 0.1
        alpha_plus_DATA_eff = 0.526248;
        alpha_minus_DATA_eff = 0.519203;
        alpha_plus_MC_eff = 0.53133;
        alpha_minus_MC_eff = 0.525911;
      }else if(fabs(Charge) < afb_chbin[2]){// 0.3
        alpha_plus_DATA_eff = 0.577542;
        alpha_minus_DATA_eff = 0.560728;
        alpha_plus_MC_eff = 0.590962;
        alpha_minus_MC_eff = 0.57714;
      }else if(fabs(Charge) < afb_chbin[3]){// 0.6
        alpha_plus_DATA_eff = 0.67162;
        alpha_minus_DATA_eff = 0.645446;
        alpha_plus_MC_eff = 0.707146;
        alpha_minus_MC_eff = 0.679753;
      }else if(fabs(Charge) < afb_chbin[4]){// 1.0
        alpha_plus_DATA_eff = 0.784729;
        alpha_minus_DATA_eff = 0.73104;
        alpha_plus_MC_eff = 0.849374;
        alpha_minus_MC_eff = 0.812045;
      }else if(fabs(Charge) < afb_chbin[5]){// 3.0, soft muons
        alpha_plus_DATA_eff = 0.756391;
        alpha_minus_DATA_eff = 0.742828;
        alpha_plus_MC_eff = 0.768852;
        alpha_minus_MC_eff = 0.766965;
      }else{// 5.0, soft electrons
        alpha_plus_DATA_eff = 0.745411;
        alpha_minus_DATA_eff = 0.737282;
        alpha_plus_MC_eff = 0.750156;
        alpha_minus_MC_eff = 0.751799;
      }
    }

    bool isChargeCorrect = false;
    isChargeCorrect = ((Charge * genpid) <= 0? true: false);
    if(isChargeCorrect){
      if(genpid > 0 ) weight *= alpha_minus_DATA_eff / alpha_minus_MC_eff;
      else weight *= alpha_plus_DATA_eff / alpha_plus_MC_eff;
    }else{
      if(genpid > 0 ) weight *= (1. - alpha_minus_DATA_eff) / (1. - alpha_minus_MC_eff);
      else weight *= (1. - alpha_plus_DATA_eff) / (1. - alpha_plus_MC_eff);
    }
  }

  return weight;
}
