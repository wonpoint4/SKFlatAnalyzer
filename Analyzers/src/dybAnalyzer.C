#include "dybAnalyzer.h"

void dybAnalyzer::initializeAnalyzer(){
  SMPAnalyzerCore::initializeAnalyzer(); //setup eff zpt roc z0
  IsSkimmed = (GetSkimName() != ""? true: false);
  IsNominalRun = !HasFlag("SYS") && !HasFlag("PDFSYS") && !HasFlag("LEPSYS") && IsSkimmed;

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
  if(IsDYSample || IsTTSample) executeEventGen();

  ///////////////// RECO level /////////////////////
  jets_raw = GetAllJets(); // This can make event loops much slower in case of running over Unskimmed samples
  if(!IsDATA || DataStream.Contains("DoubleMuon")){
    muons = MuonMomentumCorrection(SMPGetMuons("POGMediumWithLooseTrkIso", 8.0, 2.4), 0,0, true);
    executeEventWithParameter("mm"+GetEraShort());
    if(HasFlag("SYS")){
      for(TString syst:{"_jet_scale_up", "_jet_scale_down", "_jet_smear_up", "_jet_smear_down"}){
        if(syst.Contains("scale") || !IsDATA) executeEventWithParameter("mm"+GetEraShort(), syst);
      }
    }else if(HasFlag("LEPSYS")){
      for(unsigned int s=0; s<nmem_muon.size(); s++){
        for(unsigned int m=0; m<nmem_muon.at(s); m++){
          executeEventWithParameter("mm"+GetEraShort(), Form("_MuonMomentum_s%dm%d", s, m), s,m);
        }
      }
    }
  }
  if(!IsDATA || DataStream.Contains("DoubleEG") || DataStream.Contains("EGamma")){
    electrons = ElectronEnergyCorrection(SMPGetElectrons("passMediumID", 8.0, 2.5), 0,0, true);
    executeEventWithParameter("ee"+GetEraShort());
    if(HasFlag("SYS")){
      for(TString syst:{"_jet_scale_up", "_jet_scale_down", "_jet_smear_up", "_jet_smear_down"}){
        if(syst.Contains("scale") || !IsDATA) executeEventWithParameter("ee"+GetEraShort(), syst);
      }
    }else if(HasFlag("LEPSYS")){
      for(unsigned int s=0; s<nmem_electron.size(); s++){
        for(unsigned int m=0; m<nmem_electron.at(s); m++){
          executeEventWithParameter("ee"+GetEraShort(), Form("_ElectronEnergy_s%dm%d", s, m), s,m);
        }
      }
    }
  }
}

void dybAnalyzer::executeEventWithParameter(TString channel, TString option, unsigned int set, unsigned int mem){

  lepton0 = NULL;
  lepton1 = NULL;
  jet0 = NULL;
  bcharge = 0;
  prefix = channel+"/"+gprefix, hprefix = "", suffix = "";
  suffix += option;
  if((IsNominalRun || option != "") && set != 1) IsNominalLike = true;
  else IsNominalLike = false;

  // Weights Setup
  map_weight.clear();
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
  if(!HasDileptons(channel, set, mem)) return;

  // Jets
  vector<Jet> alljets = {};
  if(option.Contains("jet_scale_up")) alljets = SelectJets(ScaleJets(jets_raw, 1), "tightLepVeto", 20, (DataYear == 2016? 2.4: 2.5));
  else if(option.Contains("jet_scale_down")) alljets = SelectJets(ScaleJets(jets_raw, -1), "tightLepVeto", 20, (DataYear == 2016? 2.4: 2.5));
  else if(option.Contains("jet_smear_up")) alljets = SelectJets(SmearJets(jets_raw, 1), "tightLepVeto", 20, (DataYear == 2016? 2.4: 2.5));
  else if(option.Contains("jet_smear_down")) alljets = SelectJets(SmearJets(jets_raw, -1), "tightLepVeto", 20, (DataYear == 2016? 2.4: 2.5));
  else alljets = SelectJets(jets_raw, "tightLepVeto", 20, (DataYear == 2016? 2.4: 2.5));
  std::sort(alljets.begin(), alljets.end(), PtComparing);

  vector<Jet> lepvetojets = {}, realjets = {}, bjets = {}, ajets = {};
  for(const auto& jet:alljets){
    if(lepton0 && jet.DeltaR(*lepton0) < 0.4) continue;
    if(lepton1 && jet.DeltaR(*lepton1) < 0.4) continue;
    lepvetojets.push_back(jet);
  }
  for(const auto& jet:lepvetojets){
    if(!PUJetIDPass(jet, "Loose")) continue;
    realjets.push_back(jet);
  }

  // b-tagging
  JetTagging::Parameters DeepJet_Tight = JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Tight, JetTagging::incl, JetTagging::comb);
  JetTagging::Parameters DeepJet_Loose = JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Loose, JetTagging::incl, JetTagging::comb);

  for(auto& jet:realjets){
    //jet *= jet.BJetNNCorrection(); // full bJetEnergyCorrection (BBjetRegression)?
    if(jet.Pt() > 25 && jet.GetTaggerResult(DeepJet_Tight.j_Tagger) > mcCorr->GetJetTaggingCutValue(DeepJet_Tight.j_Tagger, DeepJet_Tight.j_WP)) bjets.push_back(jet);
    else if(jet.GetTaggerResult(DeepJet_Loose.j_Tagger) > mcCorr->GetJetTaggingCutValue(DeepJet_Loose.j_Tagger, DeepJet_Loose.j_WP)) ajets.push_back(jet);
  }

  // Jet related weights
  pujetSF = 1., btagSF = 1., bchargeSF = 1.;
  if(!IsDATA){
    pujetSF = GetPUJetWeight(lepvetojets, "Loose", 0);
    btagSF = GetBTaggingReweight_1a_2WP(realjets, DeepJet_Tight, DeepJet_Loose);
  }

  // Weights
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
  muonTrackingSF = 1.;
  muonRECOSF = 1.;
  muonIDSF = 1.;
  muonTriggerSF = 1.;
  electronRECOSF = 1.;
  electronIDSF = 1.;
  electronTriggerSF = 1.;

  if(!IsDATA){
    // Efficiency SF keys
    TString muonTrackingSF_key = "Muon_Tracking";
    TString muonRECOSF_key = "Muon_RECO";
    TString muonIDSF_key = "Muon_MediumID_trkIsoLoose";
    TString muonTriggerLeg1SF_key = "Mu17Leg1_MediumID_trkIsoLoose";
    TString muonTriggerLeg2SF_key = "Mu8Leg2_MediumID_trkIsoLoose";
    TString muonTriggerDZSF_key = GetEraShort() !="2016a"? "DZ_MediumID_trkIsoLoose": "";
    TString electronRECOSF_key = "Electron_RECO";
    TString electronIDSF_key = "Electron_MediumID";
    TString electronTriggerLeg1SF_key = "Ele23Leg1_MediumID";
    TString electronTriggerLeg2SF_key = "Ele12Leg2_MediumID";
    TString electronTriggerDZSF_key = GetEraShort().Contains("2016")? "DZ_MediumID": "";

    if(HasFlag("LEPSYS") && option == "" && !hprefix.Contains("ss_")){
      muonTrackingSF_sys = Make2DWeights(fEff->GetStructure(muonTrackingSF_key));
      muonRECOSF_sys = Make2DWeights(fEff->GetStructure(muonRECOSF_key));
      muonIDSF_sys = Make2DWeights(fEff->GetStructure(muonIDSF_key));
      muonTriggerSF_sys = Make2DWeights(fEff->GetStructure(muonTriggerLeg1SF_key));
      electronRECOSF_sys = Make2DWeights(fEff->GetStructure(electronRECOSF_key));
      electronIDSF_sys = Make2DWeights(fEff->GetStructure(electronIDSF_key));
      electronTriggerSF_sys = Make2DWeights(fEff->GetStructure(electronTriggerLeg1SF_key));
    }

    // Tracking, RECO, ID SF per lepton
    for(const Lepton* lepton:leptons){
      if(lepton->LeptonFlavour() == Lepton::MUON){
        muonTrackingSF *= fEff->GetEfficiencySF(muonTrackingSF_key, lepton, 0,0);
        muonRECOSF *= fEff->GetEfficiencySF(muonRECOSF_key, lepton, 0,0);
        muonIDSF *= fEff->GetEfficiencySF(muonIDSF_key, lepton, 0,0);

        if(HasFlag("LEPSYS") && option == "" && !hprefix.Contains("ss_")){
          for(unsigned int s=0; s<muonTrackingSF_sys.size(); s++){
            for(unsigned int m=0; m<muonTrackingSF_sys[s].size(); m++){
              muonTrackingSF_sys[s][m] *= fEff->GetEfficiencySF(muonTrackingSF_key, lepton, s,m);
            }
          }
          for(unsigned int s=0; s<muonRECOSF_sys.size(); s++){
            for(unsigned int m=0; m<muonRECOSF_sys[s].size(); m++){
              muonRECOSF_sys[s][m] *= fEff->GetEfficiencySF(muonRECOSF_key, lepton, s,m);
            }
          }
          for(unsigned int s=0; s<muonIDSF_sys.size(); s++){
            for(unsigned int m=0; m<muonIDSF_sys[s].size(); m++){
              muonIDSF_sys[s][m] *= fEff->GetEfficiencySF(muonIDSF_key, lepton, s,m);
            }
          }
        }
      }else if(lepton->LeptonFlavour() == Lepton::ELECTRON){
        electronRECOSF *= fEff->GetEfficiencySF(electronRECOSF_key, lepton, 0,0);
        electronIDSF *= fEff->GetEfficiencySF(electronIDSF_key, lepton, 0,0);

        if(HasFlag("LEPSYS") && option == "" && !hprefix.Contains("ss_")){
          for(unsigned int s=0; s<electronRECOSF_sys.size(); s++){
            for(unsigned int m=0; m<electronRECOSF_sys[s].size(); m++){
              electronRECOSF_sys[s][m] *= fEff->GetEfficiencySF(electronRECOSF_key, lepton, s,m);
            }
          }
          for(unsigned int s=0; s<electronIDSF_sys.size(); s++){
            for(unsigned int m=0; m<electronIDSF_sys[s].size(); m++){
              electronIDSF_sys[s][m] *= fEff->GetEfficiencySF(electronIDSF_key, lepton, s,m);
            }
          }
        }
      }
    }

    // Trigger SF
    if(channel.Contains("mm"+GetEraShort())){
      muonTriggerSF *= GetDileptonTriggerSF(muonTriggerLeg1SF_key, muonTriggerLeg2SF_key, muonTriggerDZSF_key, leptons, 0,0);
      if(HasFlag("LEPSYS") && option == "" && !hprefix.Contains("ss_")){
        for(unsigned int s=0; s<muonTriggerSF_sys.size(); s++){
          for(unsigned int m=0; m<muonTriggerSF_sys[s].size(); m++){
            muonTriggerSF_sys[s][m] *= GetDileptonTriggerSF(muonTriggerLeg1SF_key, muonTriggerLeg2SF_key, muonTriggerDZSF_key, leptons, s,m);
          }
        }
      }
    }
    if(channel.Contains("ee"+GetEraShort())){
      electronTriggerSF *= GetDileptonTriggerSF(electronTriggerLeg1SF_key, electronTriggerLeg2SF_key, electronTriggerDZSF_key, leptons, 0,0);
      if(HasFlag("LEPSYS") && option == "" && !hprefix.Contains("ss_")){
        for(unsigned int s=0; s<electronTriggerSF_sys.size(); s++){
          for(unsigned int m=0; m<electronTriggerSF_sys[s].size(); m++){
            electronTriggerSF_sys[s][m] *= GetDileptonTriggerSF(electronTriggerLeg1SF_key, electronTriggerLeg2SF_key, electronTriggerDZSF_key, leptons, s,m);
          }
        }
      }
    }
  }

  map_weight[""] *= muonTrackingSF;
  if(IsNominalLike){
    FillHist(prefix+hprefix+"weight_muonTrackingSF"+suffix, muonTrackingSF, map_weight[""], 200,-5,5);
    FillCutflow(prefix+hprefix+"cutflow"+suffix, "muonTrackingSF", map_weight[""]);
  }
  map_weight[""] *= muonRECOSF;
  if(IsNominalLike){
    FillHist(prefix+hprefix+"weight_muonRECOSF"+suffix, muonRECOSF, map_weight[""], 200,-5,5);
    FillCutflow(prefix+hprefix+"cutflow"+suffix, "muonRECOSF", map_weight[""]);
  }
  map_weight[""] *= muonIDSF;
  if(IsNominalLike){
    FillHist(prefix+hprefix+"weight_muonIDSF"+suffix, muonIDSF, map_weight[""], 200,-5,5);
    FillCutflow(prefix+hprefix+"cutflow"+suffix, "muonIDSF", map_weight[""]);
  }
  map_weight[""] *= muonTriggerSF;
  if(IsNominalLike){
    FillHist(prefix+hprefix+"weight_muonTriggerSF"+suffix, muonTriggerSF, map_weight[""], 200,-5,5);
    FillCutflow(prefix+hprefix+"cutflow"+suffix, "muonTriggerSF", map_weight[""]);
  }
  map_weight[""] *= electronRECOSF;
  if(IsNominalLike){
    FillHist(prefix+hprefix+"weight_electronRECOSF"+suffix, electronRECOSF, map_weight[""], 200,-5,5);
    FillCutflow(prefix+hprefix+"cutflow"+suffix, "electronRECOSF", map_weight[""]);
  }
  map_weight[""] *= electronIDSF;
  if(IsNominalLike){
    FillHist(prefix+hprefix+"weight_electronIDSF"+suffix, electronIDSF, map_weight[""], 200,-5,5);
    FillCutflow(prefix+hprefix+"cutflow"+suffix, "electronIDSF", map_weight[""]);
  }
  map_weight[""] *= electronTriggerSF;
  if(IsNominalLike){
    FillHist(prefix+hprefix+"weight_electronTriggerSF"+suffix, electronTriggerSF, map_weight[""], 200,-5,5);
    FillCutflow(prefix+hprefix+"cutflow"+suffix, "electronTriggerSF", map_weight[""]);
  }

  //==== Inclusive DY
  double dimass = (*lepton0 + *lepton1).M();
  double dirap = (*lepton0 + *lepton1).Rapidity();
  double dipt = (*lepton0 + *lepton1).Pt();
  double costhetaCS = GetCosThetaCS(lepton0, lepton1);

  if(IsNominalLike){
    FillHist(prefix+hprefix+"mll_IncDY"+suffix, dimass, map_weight[""], 80,70,110);
    FillHist(prefix+hprefix+"yll_IncDY"+suffix, dirap, map_weight[""], 96,-2.4,2.4);
    FillHist(prefix+hprefix+"ptll_IncDY"+suffix, dipt, map_weight[""], AFBAnalyzer::unfold_0bjet_ptbinnum_reco,AFBAnalyzer::unfold_0bjet_ptbin_reco);
    FillHist(prefix+hprefix+"lpt_IncDY"+suffix, lepton0->Pt(), map_weight[""], AFBAnalyzer::lptbinnum,AFBAnalyzer::lptbin);
    FillHist(prefix+hprefix+"leta_IncDY"+suffix, lepton0->Eta(), map_weight[""], 50,-2.5,2.5);
    FillHist(prefix+hprefix+"lpt_IncDY"+suffix, lepton1->Pt(), map_weight[""], AFBAnalyzer::lptbinnum,AFBAnalyzer::lptbin);
    FillHist(prefix+hprefix+"leta_IncDY"+suffix, lepton1->Eta(), map_weight[""], 50,-2.5,2.5);
    FillHist(prefix+hprefix+"costhetaCS_IncDY"+suffix, dimass, dirap, dipt, costhetaCS, map_weight, afb_mbinnum,(double*)afb_mbin, afb_ybinnum,(double*)afb_ybin, afb_ptbinnum,(double*)afb_ptbin, 20,-1,1);

    FillHist(prefix+hprefix+"nalljets_IncDY"+suffix, alljets.size(), map_weight[""], 15,0,15);
    FillHist(prefix+hprefix+"nlepvetojets_IncDY"+suffix, lepvetojets.size(), map_weight[""], 15,0,15);
    FillHist(prefix+hprefix+"nrealjets_IncDY"+suffix, realjets.size(), map_weight[""], 15,0,15);
    FillHist(prefix+hprefix+"nbjets_IncDY"+suffix, bjets.size(), map_weight[""], 10,0,10);
  }

  if(bjets.size() == 0) return;
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "Tightnb", map_weight[""]);
  //bjets.at(0) *= bjets.at(0).BJetNNCorrection(); // bJetEnergyCorrection (BjetRegression)?
  jet0 = &bjets.at(0);
  bcharge = jetCharge(*jet0);
  bchargeSF = GetbChargeSFWeight(bjets, 1, 0);
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

  map_weight[""] *= bchargeSF;
  if(IsNominalLike){
    FillHist(prefix+hprefix+"weight_bChargeSF1"+suffix, bchargeSF, map_weight[""], 200,-5,5);
    FillCutflow(prefix+hprefix+"cutflow"+suffix, "bChargeSF1", map_weight[""]);
  }

  //==== Weights of Systematics
  if(!IsDATA && HasFlag("SYS") && option == "" && !hprefix.Contains("ss_")){
    // Prefiring weight
    map_weight["_noprefireweight"] =  map_weight[""] / prefireweight;
    map_weight["_prefireweight_up"] =  map_weight[""] / prefireweight * L1PrefireReweight_Up;
    map_weight["_prefireweight_down"] = map_weight[""] / prefireweight * L1PrefireReweight_Down;

    // PU reweight
    map_weight["_noPUweight"] = map_weight[""] / PUweight;
    map_weight["_PUweight_up"] = map_weight[""] / PUweight * GetPileUpWeight(nPileUp, 1);
    map_weight["_PUweight_down"] = map_weight[""] / PUweight * GetPileUpWeight(nPileUp, -1);

    // b-tagging SF
    map_weight["_nobtagSF"] =  map_weight[""] / btagSF;
    map_weight["_btagSF_hup"] = map_weight[""] / btagSF * GetBTaggingReweight_1a_2WP(realjets, DeepJet_Tight, DeepJet_Loose, "SystUpHTag");
    map_weight["_btagSF_hdown"] = map_weight[""] / btagSF * GetBTaggingReweight_1a_2WP(realjets, DeepJet_Tight, DeepJet_Loose, "SystDownHTag");
    map_weight["_btagSF_hcorr"] = map_weight[""] / btagSF * GetBTaggingReweight_1a_2WP(realjets, DeepJet_Tight, DeepJet_Loose, "SystUpHTagCorr");
    map_weight["_btagSF_huncorr"+GetEraShort()] = map_weight[""] / btagSF * GetBTaggingReweight_1a_2WP(realjets, DeepJet_Tight, DeepJet_Loose, "SystUpHTagUnCorr");
    map_weight["_btagSF_lup"] = map_weight[""] / btagSF * GetBTaggingReweight_1a_2WP(realjets, DeepJet_Tight, DeepJet_Loose, "SystUpLTag");
    map_weight["_btagSF_ldown"] = map_weight[""] / btagSF * GetBTaggingReweight_1a_2WP(realjets, DeepJet_Tight, DeepJet_Loose, "SystDownLTag");
    map_weight["_btagSF_lcorr"] = map_weight[""] / btagSF * GetBTaggingReweight_1a_2WP(realjets, DeepJet_Tight, DeepJet_Loose, "SystUpLTagCorr");;
    map_weight["_btagSF_luncorr"+GetEraShort()] = map_weight[""] / btagSF * GetBTaggingReweight_1a_2WP(realjets, DeepJet_Tight, DeepJet_Loose, "SystUpLTagUnCorr");

    // PUjetID SF
    map_weight["_noPUjetSF"] =  map_weight[""] / pujetSF;
    map_weight["_PUjetSF_up"] =  map_weight[""] / pujetSF * GetPUJetWeight(lepvetojets, "Loose", 1);
    map_weight["_PUjetSF_down"] = map_weight[""] / pujetSF * GetPUJetWeight(lepvetojets, "Loose", -1);

    // bChargeID SF
    map_weight["_bChargeSF0"] = map_weight[""] / bchargeSF * GetbChargeSFWeight(bjets, 0, 0);
    map_weight["_bChargeSF0_up"] = map_weight[""] / bchargeSF * GetbChargeSFWeight(bjets, 0, 1);
    map_weight["_bChargeSF0_down"] = map_weight[""] / bchargeSF * GetbChargeSFWeight(bjets, 0, -1);
    map_weight["_nobChargeSF1"] = map_weight[""] / bchargeSF;
    for(TString bCh:{"0", "1", "2", "3", "4", "5"}){
      map_weight["_bChargeSF1_up"+bCh] = map_weight[""] / bchargeSF * GetbChargeSFWeight(bjets, 1, 1, bCh);
      map_weight["_bChargeSF1_down"+bCh] = map_weight[""] / bchargeSF * GetbChargeSFWeight(bjets, 1, -1, bCh);
    }
    map_weight["_bChargeSFHS"] = map_weight[""] / bchargeSF * GetbChargeSFWeight(bjets, 2, 0);
    map_weight["_bChargeSFHS_up"] = map_weight[""] / bchargeSF * GetbChargeSFWeight(bjets, 2, 1);
    map_weight["_bChargeSFHS_down"] = map_weight[""] / bchargeSF * GetbChargeSFWeight(bjets, 2, -1);
  }else if(!IsDATA && HasFlag("PDFSYS") && option == "" && !hprefix.Contains("ss_")){
    // Zpt Reweight
    map_weight["_noZpt"] =  map_weight[""] / zptweight;
    map_weight["_Zpt_gym"] =  map_weight[""] / zptweight * zptweight_gym;

    // Weak Reweight
    map_weight["_noWeak"] =  map_weight[""] / weakweight;

    // Top pt Reweight
    map_weight["_noToppt"] =  map_weight[""] / topptweight;

    // AlphaS
    if(weight_AlphaS->size() == 2){
      map_weight["_alphaS_up"] = map_weight[""] * weight_AlphaS->at(1);
      map_weight["_alphaS_down"] = map_weight[""] * weight_AlphaS->at(0);
    }

    // FSR, ISR
    if(weight_PSSyst->size()){
      map_weight["_FSR_up"] = map_weight[""] * TMath::Range(-5., 5., weight_PSSyst->at(1));
      map_weight["_FSR_down"] = map_weight[""] * TMath::Range(-5., 5., weight_PSSyst->at(0));
      map_weight["_ISR_up"] = map_weight[""] * TMath::Range(-5., 5., weight_PSSyst->at(3));
      map_weight["_ISR_down"] = map_weight[""] * TMath::Range(-5., 5., weight_PSSyst->at(2));
    }

    // Scale Variations
    for(unsigned int i=0; i<weight_Scale->size(); i++) map_weight[Form("_scalevariation%d", i)] = map_weight[""] * TMath::Range(-10., 10., weight_Scale->at(i));

    // PDF
    for(unsigned int i=0; i<weight_PDF->size(); i++) map_weight[Form("_pdf%d", i)] = map_weight[""] * TMath::Range(-10., 10., weight_PDF->at(i));
  }else if(!IsDATA && HasFlag("LEPSYS") && option == "" && !hprefix.Contains("ss_")){
    // EfficiencySF - stat
    for(int j=0; j<fEff->nreplica; j++){
      double SF_stat = 1. / muonTrackingSF / muonRECOSF / muonIDSF / muonTriggerSF / electronRECOSF / electronIDSF / electronTriggerSF;
      SF_stat *= muonTrackingSF_sys.size()? (double)muonTrackingSF_sys[0][j]: 1.;
      SF_stat *= muonRECOSF_sys.size()? (double)muonRECOSF_sys[0][j]: 1.;
      SF_stat *= muonIDSF_sys.size()? (double)muonIDSF_sys[0][j]: 1.;
      SF_stat *= muonTriggerSF_sys.size()? (double)muonTriggerSF_sys[0][j]: 1.;
      SF_stat *= electronRECOSF_sys.size()? (double)electronRECOSF_sys[0][j]: 1.;
      SF_stat *= electronIDSF_sys.size()? (double)electronIDSF_sys[0][j]: 1.;
      SF_stat *= electronTriggerSF_sys.size()? (double)electronTriggerSF_sys[0][j]: 1.;
      map_weight[Form("_lepeffSF_stat%d", j)] = map_weight[""] * SF_stat;
    }

    // EfficiencySF - syst
    for(unsigned int i=0; i<muonTrackingSF_sys.size(); i++){
      for(unsigned int j=0; j<muonTrackingSF_sys[i].size(); j++){
        map_weight[Form("_muonTrackingeffSF_s%dm%d", i, j)] = map_weight[""] / muonTrackingSF * muonTrackingSF_sys[i][j];
      }
    }
    for(unsigned int i=0; i<muonRECOSF_sys.size(); i++){
      for(unsigned int j=0; j<muonRECOSF_sys[i].size(); j++){
        map_weight[Form("_muonRECOeffSF_s%dm%d", i, j)] = map_weight[""] / muonRECOSF * muonRECOSF_sys[i][j];
      }
    }
    for(unsigned int i=0; i<muonIDSF_sys.size(); i++){
      for(unsigned int j=0; j<muonIDSF_sys[i].size(); j++){
        map_weight[Form("_muonIDeffSF_s%dm%d", i, j)] = map_weight[""] / muonIDSF * muonIDSF_sys[i][j];
      }
    }
    for(unsigned int i=0; i<muonTriggerSF_sys.size(); i++){
      for(unsigned int j=0; j<muonTriggerSF_sys[i].size(); j++){
        map_weight[Form("_muonTriggereffSF_s%dm%d", i, j)] = map_weight[""] / muonTriggerSF * muonTriggerSF_sys[i][j];
      }
    }
    for(unsigned int i=0; i<electronRECOSF_sys.size(); i++){
      for(unsigned int j=0; j<electronRECOSF_sys[i].size(); j++){
        map_weight[Form("_electronRECOeffSF_s%dm%d", i, j)] = map_weight[""] / electronRECOSF * electronRECOSF_sys[i][j];
      }
    }
    for(unsigned int i=0; i<electronIDSF_sys.size(); i++){
      for(unsigned int j=0; j<electronIDSF_sys[i].size(); j++){
        map_weight[Form("_electronIDeffSF_s%dm%d", i, j)] = map_weight[""] / electronIDSF * electronIDSF_sys[i][j];
      }
    }
    for(unsigned int i=0; i<electronTriggerSF_sys.size(); i++){
      for(unsigned int j=0; j<electronTriggerSF_sys[i].size(); j++){
        map_weight[Form("_electronTriggereffSF_s%dm%d", i, j)] = map_weight[""] / electronTriggerSF * electronTriggerSF_sys[i][j];
      }
    }

    // Dummy for Aepcor, RoccoR
    if(channel.Contains("mm"+GetEraShort())){
      for(unsigned int i=0; i<nmem_electron.size(); i++){
        for(unsigned int j=0; j<nmem_electron.at(i); j++){
          map_weight[Form("_ElectronEnergy_s%dm%d", i, j)] = map_weight[""];
        }
      }
    }
    if(channel.Contains("ee"+GetEraShort())){
      for(unsigned int i=0; i<nmem_muon.size(); i++){
        for(unsigned int j=0; j<nmem_muon.at(i); j++){
          map_weight[Form("_MuonMomentum_s%dm%d", i, j)] = map_weight[""];
        }
      }
    }
  }else if(!IsDATA && MCSample.Contains("MiNNLO") && IsNominalRun && !hprefix.Contains("ss_")){
    for(unsigned int i=0; i<weight_sthw2->size(); i++) map_weight[Form("_sthw2_%d", i)] = map_weight[""] * weight_sthw2->at(i);
  }
  if((HasFlag("SYS") || HasFlag("PDFSYS") || HasFlag("LEPSYS")) && option == "") map_weight.erase("");

  // For comparison with Hyonsan's results (8 bins)
  if(jet0->Pt() > 40){
    FillHist(prefix+hprefix+"mll_nbjets"+suffix, dimass, map_weight, 80,70,110);
    FillHist(prefix+hprefix+"yll_nbjets"+suffix, dirap, map_weight, 96,-2.4,2.4);
    FillHist(prefix+hprefix+"lpt_nbjets"+suffix, lepton0->Pt(), map_weight, 200,0,200);
    FillHist(prefix+hprefix+"leta_nbjets"+suffix, lepton0->Eta(), map_weight, 50,-2.5,2.5);
    FillHist(prefix+hprefix+"lpt_nbjets"+suffix, lepton1->Pt(), map_weight, 200,0,200);
    FillHist(prefix+hprefix+"leta_nbjets"+suffix, lepton1->Eta(), map_weight, 50,-2.5,2.5);
    FillHist(prefix+hprefix+"bpt_nbjets"+suffix, jet0->Pt(), map_weight, 200,0,200);
    FillHist(prefix+hprefix+"beta_nbjets"+suffix, jet0->Eta(), map_weight, 50,-2.5,2.5);
    FillHist(prefix+hprefix+"bChargeRaw_nbjets"+suffix, bcharge, map_weight, 200,-5,5);
    FillHist(prefix+hprefix+"bCharge_nbjets"+suffix, (bcharge < 0? -0.5: 0.5), map_weight, 2,-1,1);
    FillHist(prefix+hprefix+"met_nbjets"+suffix, PuppiMET_Type1_pt, map_weight, 200,0,200);
    FillHist(prefix+hprefix+"ZbdPhi_nbjets"+suffix, abs((*lepton0 + *lepton1).DeltaPhi(*jet0)), map_weight, 200,0,10);
    FillHist(prefix+hprefix+"Zbpt_nbjets"+suffix, (*lepton0 + *lepton1 + *jet0).Pt(), map_weight, 200,0,200);
    FillHist(prefix+hprefix+"Zpt_nbjets"+suffix, (*lepton0 + *lepton1).Pt(), map_weight, 200,0,200);
    FillHist(prefix+hprefix+"costhetaRecoil_nbjets"+suffix, dimass, fabs(bcharge), jet0->Pt(), costhetaRecoil, map_weight, afb_mbinnum,(double*)afb_mbin, afb_chbinnum,(double*)afb_chbin, afb_ptbinnum,(double*)afb_ptbin, 20,-1,1);
    FillHist(prefix+hprefix+"costhetaRecoil2_nbjets"+suffix, dimass, fabs(bcharge), jet0->Pt(), costhetaRecoil, map_weight, afb_mbinnum_HS,(double*)afb_mbin_HS, afb_chbinnum,(double*)afb_chbin, afb_ptbinnum,(double*)afb_ptbin, 20,-1,1);
  }

  if(bjets.size() != 1) return;
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "Tight1b", map_weight[""]);

  FillHist(prefix+hprefix+"mll_Tight1b"+suffix, dimass, map_weight, 80,70,110);
  FillHist(prefix+hprefix+"yll_Tight1b"+suffix, dirap, map_weight, 96,-2.4,2.4);
  FillHist(prefix+hprefix+"lpt_Tight1b"+suffix, lepton0->Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"leta_Tight1b"+suffix, lepton0->Eta(), map_weight, 50,-2.5,2.5);
  FillHist(prefix+hprefix+"lpt_Tight1b"+suffix, lepton1->Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"leta_Tight1b"+suffix, lepton1->Eta(), map_weight, 50,-2.5,2.5);
  FillHist(prefix+hprefix+"bpt_Tight1b"+suffix, jet0->Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"beta_Tight1b"+suffix, jet0->Eta(), map_weight, 50,-2.5,2.5);
  FillHist(prefix+hprefix+"bChargeRaw_Tight1b"+suffix, bcharge, map_weight, 200,-5,5);
  FillHist(prefix+hprefix+"bCharge_Tight1b"+suffix, (bcharge < 0? -0.5: 0.5), map_weight, 2,-1,1);
  FillHist(prefix+hprefix+"met_Tight1b"+suffix, PuppiMET_Type1_pt, map_weight, 200,0,200);
  FillHist(prefix+hprefix+"ZbdPhi_Tight1b"+suffix, abs((*lepton0 + *lepton1).DeltaPhi(*jet0)), map_weight, 200,0,10);
  FillHist(prefix+hprefix+"Zbpt_Tight1b"+suffix, (*lepton0 + *lepton1 + *jet0).Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"Zpt_Tight1b"+suffix, (*lepton0 + *lepton1).Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"costhetaRecoil_Tight1b"+suffix, dimass, fabs(bcharge), jet0->Pt(), costhetaRecoil, map_weight, afb_mbinnum,(double*)afb_mbin, afb_chbinnum,(double*)afb_chbin, afb_ptbinnum,(double*)afb_ptbin, 20,-1,1);

  FillHist(prefix+hprefix+"na20_Tight1b"+suffix, ajets.size(), map_weight, 10,0,10);
  unsigned int najets_25 = 0, najets_30 = 0, najets_35 = 0, najets_40 = 0;
  for(const auto& jet:ajets){
    if(jet.Pt() > 25) najets_25++;
    if(jet.Pt() > 30) najets_30++;
    if(jet.Pt() > 35) najets_35++;
    if(jet.Pt() > 40) najets_40++;
  }
  FillHist(prefix+hprefix+"na25_Tight1b"+suffix, najets_25, map_weight, 10,0,10);
  FillHist(prefix+hprefix+"na30_Tight1b"+suffix, najets_30, map_weight, 10,0,10);
  FillHist(prefix+hprefix+"na35_Tight1b"+suffix, najets_35, map_weight, 10,0,10);
  FillHist(prefix+hprefix+"na40_Tight1b"+suffix, najets_40, map_weight, 10,0,10);

  if(ajets.size() > 0){
    FillHist(prefix+hprefix+"apt_Tight1b"+suffix, ajets.at(0).Pt(), map_weight, 200,0,200);
    if(ajets.at(0).Pt() > 25) return;
  }
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "Veto2b", map_weight[""]);

  FillHist(prefix+hprefix+"mll_Veto2b"+suffix, dimass, map_weight, 80,70,110);
  FillHist(prefix+hprefix+"yll_Veto2b"+suffix, dirap, map_weight, 96,-2.4,2.4);
  FillHist(prefix+hprefix+"lpt_Veto2b"+suffix, lepton0->Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"leta_Veto2b"+suffix, lepton0->Eta(), map_weight, 50,-2.5,2.5);
  FillHist(prefix+hprefix+"lpt_Veto2b"+suffix, lepton1->Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"leta_Veto2b"+suffix, lepton1->Eta(), map_weight, 50,-2.5,2.5);
  FillHist(prefix+hprefix+"bpt_Veto2b"+suffix, jet0->Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"beta_Veto2b"+suffix, jet0->Eta(), map_weight, 50,-2.5,2.5);
  FillHist(prefix+hprefix+"bChargeRaw_Veto2b"+suffix, bcharge, map_weight, 200,-5,5);
  FillHist(prefix+hprefix+"bCharge_Veto2b"+suffix, (bcharge < 0? -0.5: 0.5), map_weight, 2,-1,1);
  FillHist(prefix+hprefix+"met_Veto2b"+suffix, PuppiMET_Type1_pt, map_weight, 200,0,200);
  FillHist(prefix+hprefix+"ZbdPhi_Veto2b"+suffix, abs((*lepton0 + *lepton1).DeltaPhi(*jet0)), map_weight, 200,0,10);
  FillHist(prefix+hprefix+"Zbpt_Veto2b"+suffix, (*lepton0 + *lepton1 + *jet0).Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"Zpt_Veto2b"+suffix, (*lepton0 + *lepton1).Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"costhetaRecoil_Veto2b"+suffix, dimass, fabs(bcharge), jet0->Pt(), costhetaRecoil, map_weight, afb_mbinnum,(double*)afb_mbin, afb_chbinnum,(double*)afb_chbin, afb_ptbinnum,(double*)afb_ptbin, 20,-1,1);

  if(PuppiMET_Type1_pt > 60) return;
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "MET60", map_weight[""]);

  FillHist(prefix+hprefix+"mll_met60"+suffix, dimass, map_weight, 80,70,110);
  FillHist(prefix+hprefix+"yll_met60"+suffix, dirap, map_weight, 96,-2.4,2.4);
  FillHist(prefix+hprefix+"lpt_met60"+suffix, lepton0->Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"leta_met60"+suffix, lepton0->Eta(), map_weight, 50,-2.5,2.5);
  FillHist(prefix+hprefix+"lpt_met60"+suffix, lepton1->Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"leta_met60"+suffix, lepton1->Eta(), map_weight, 50,-2.5,2.5);
  FillHist(prefix+hprefix+"bpt_met60"+suffix, jet0->Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"beta_met60"+suffix, jet0->Eta(), map_weight, 50,-2.5,2.5);
  FillHist(prefix+hprefix+"bChargeRaw_met60"+suffix, bcharge, map_weight, 200,-5,5);
  FillHist(prefix+hprefix+"bCharge_met60"+suffix, (bcharge < 0? -0.5: 0.5), map_weight, 2,-1,1);
  FillHist(prefix+hprefix+"met_met60"+suffix, PuppiMET_Type1_pt, map_weight, 200,0,200);
  FillHist(prefix+hprefix+"ZbdPhi_met60"+suffix, abs((*lepton0 + *lepton1).DeltaPhi(*jet0)), map_weight, 200,0,10);
  FillHist(prefix+hprefix+"Zbpt_met60"+suffix, (*lepton0 + *lepton1 + *jet0).Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"Zpt_met60"+suffix, (*lepton0 + *lepton1).Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"costhetaRecoil_met60"+suffix, dimass, fabs(bcharge), jet0->Pt(), costhetaRecoil, map_weight, afb_mbinnum,(double*)afb_mbin, afb_chbinnum,(double*)afb_chbin, afb_ptbinnum,(double*)afb_ptbin, 20,-1,1);

  if((*lepton0 + *lepton1).Pt() < 15) return;
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "ZpT15", map_weight[""]);

  FillHist(prefix+hprefix+"mll"+suffix, dimass, map_weight, 80,70,110);
  FillHist(prefix+hprefix+"yll"+suffix, dirap, map_weight, 96,-2.4,2.4);
  FillHist(prefix+hprefix+"lpt"+suffix, lepton0->Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"leta"+suffix, lepton0->Eta(), map_weight, 50,-2.5,2.5);
  FillHist(prefix+hprefix+"lpt"+suffix, lepton1->Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"leta"+suffix, lepton1->Eta(), map_weight, 50,-2.5,2.5);
  FillHist(prefix+hprefix+"bpt"+suffix, jet0->Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"beta"+suffix, jet0->Eta(), map_weight, 50,-2.5,2.5);
  FillHist(prefix+hprefix+"bChargeRaw"+suffix, bcharge, map_weight, 200,-5,5);
  FillHist(prefix+hprefix+"bCharge"+suffix, (bcharge < 0? -0.5: 0.5), map_weight, 2,-1,1);
  FillHist(prefix+hprefix+"met"+suffix, PuppiMET_Type1_pt, map_weight, 200,0,200);
  FillHist(prefix+hprefix+"ZbdPhi"+suffix, abs((*lepton0 + *lepton1).DeltaPhi(*jet0)), map_weight, 200,0,10);
  FillHist(prefix+hprefix+"Zbpt"+suffix, (*lepton0 + *lepton1 + *jet0).Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"Zpt"+suffix, (*lepton0 + *lepton1).Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"costhetaRecoil"+suffix, dimass, fabs(bcharge), jet0->Pt(), costhetaRecoil, map_weight, afb_mbinnum,(double*)afb_mbin, afb_chbinnum,(double*)afb_chbin, afb_ptbinnum,(double*)afb_ptbin, 20,-1,1);
  FillHist(prefix+hprefix+"costhetaRecoil2"+suffix, dimass, fabs(bcharge), jet0->Pt(), costhetaRecoil, map_weight, afb_mbinnum_HS,(double*)afb_mbin_HS, afb_chbinnum,(double*)afb_chbin, afb_ptbinnum,(double*)afb_ptbin, 20,-1,1);
  FillHist(prefix+hprefix+"costhetaRecoil3"+suffix, dimass, fabs(bcharge), jet0->Pt(), costhetaRecoil, map_weight, afb_mbinnum_original,(double*)afb_mbin_original, afb_chbinnum,(double*)afb_chbin, afb_ptbinnum,(double*)afb_ptbin, 20,-1,1);
  FillHist(prefix+hprefix+"costhetaRecoil4"+suffix, dimass, fabs(bcharge), jet0->Pt(), costhetaRecoil, map_weight, afb_mbinnum_original2,(double*)afb_mbin_original2, afb_chbinnum,(double*)afb_chbin, afb_ptbinnum,(double*)afb_ptbin, 20,-1,1);
  FillHist(prefix+hprefix+"AbscosthetaRecoil"+suffix, dimass, fabs(bcharge), jet0->Pt(), fabs(costhetaRecoil), map_weight, afb_mbinnum,(double*)afb_mbin, afb_chbinnum,(double*)afb_chbin, afb_ptbinnum,(double*)afb_ptbin, 10,0,1);
  FillHist(prefix+hprefix+"costhetaCS"+suffix, dimass, dirap, dipt, costhetaCS, map_weight, afb_mbinnum,(double*)afb_mbin, afb_ybinnum,(double*)afb_ybin, afb_ptbinnum,(double*)afb_ptbin, 20,-1,1);

  // bCharges vs. nPV
  FillHist(prefix+hprefix+"nPV"+suffix, nPV, map_weight, 100,0,100);
  TString prefix_nPV = prefix+hprefix;
  if(nPV <= 10) prefix_nPV.ReplaceAll(GetEraShort(), GetEraShort()+"F");      // few
  else if(nPV <= 20) prefix_nPV.ReplaceAll(GetEraShort(), GetEraShort()+"S"); // some
  else if(nPV <= 30) prefix_nPV.ReplaceAll(GetEraShort(), GetEraShort()+"L"); // little
  else if(nPV <= 40) prefix_nPV.ReplaceAll(GetEraShort(), GetEraShort()+"M"); // middle
  else if(nPV <= 50) prefix_nPV.ReplaceAll(GetEraShort(), GetEraShort()+"H"); // high
  else prefix_nPV.ReplaceAll(GetEraShort(), GetEraShort()+"V");               // very high
  for(int i=1; i<afb_chbinnum+1; i++){
    if(afb_chbin[i-1] < abs(bcharge) && abs(bcharge) < afb_chbin[i]){
      FillHist(Form(prefix+hprefix+"bCharge%dRaw"+suffix, i-1), bcharge, map_weight, 200,-5,5);
      FillHist(Form(prefix+hprefix+"bCharge%d"+suffix, i-1), (bcharge < 0? -0.5: 0.5), map_weight, 2,-1,1);
      FillHist(prefix+hprefix+"bChargebin"+suffix, LHAPDF::sgn(bcharge) * (i - 0.5), map_weight, 12,-6,6);
      FillHist(prefix+hprefix+"bChargeAbsbin"+suffix, i - 0.5, map_weight, 6,0,6);
      FillHist(Form(prefix_nPV+"bCharge%dRaw"+suffix, i-1), bcharge, map_weight, 200,-5,5);
      FillHist(Form(prefix_nPV+"bCharge%d"+suffix, i-1), (bcharge < 0? -0.5: 0.5), map_weight, 2,-1,1);
      FillHist(prefix_nPV+"bChargebin"+suffix, LHAPDF::sgn(bcharge) * (i - 0.5), map_weight, 12,-6,6);
      FillHist(prefix_nPV+"bChargeAbsbin"+suffix, i - 0.5, map_weight, 6,0,6);
    }
  }

  FillHist(prefix+hprefix+"nj20"+suffix, realjets.size(), map_weight, 10,0,10);
  if(realjets.size() > 1) FillHist(prefix+hprefix+"j1pt"+suffix, realjets.at(1).Pt(), map_weight, 200,0,200);
  unsigned int njets_25 = 0, njets_30 = 0, njets_35 = 0, njets_40 = 0;
  for(const auto& jet:realjets){
    if(jet.Pt() > 25) njets_25++;
    if(jet.Pt() > 30) njets_30++;
    if(jet.Pt() > 35) njets_35++;
    if(jet.Pt() > 40) njets_40++;
  }
  FillHist(prefix+hprefix+"nj25"+suffix, njets_25, map_weight, 10,0,10);
  FillHist(prefix+hprefix+"nj30"+suffix, njets_30, map_weight, 10,0,10);
  FillHist(prefix+hprefix+"nj35"+suffix, njets_35, map_weight, 10,0,10);
  FillHist(prefix+hprefix+"nj40"+suffix, njets_40, map_weight, 10,0,10);
}

bool dybAnalyzer::HasDileptons(TString channel, unsigned int s, unsigned int m){
  double l0pt = 20., l1pt = 10.;
  if(channel.Contains("mm"+GetEraShort())){
    muons = MuonMomentumCorrection(muons, s,m, false);
    if(muons.size() > 0) lepton0 = &muons.at(0);
    if(muons.size() > 1) lepton1 = &muons.at(1);
  }else if(channel.Contains("ee"+GetEraShort())){
    l0pt = 25.;
    l1pt = 15.;
    electrons = ElectronEnergyCorrection(electrons, s,m, false);
    if(electrons.size() > 0) lepton0 = &electrons.at(0);
    if(electrons.size() > 1) lepton1 = &electrons.at(1);
  }else{
    cout<<"[dybAnalyzer::Hasleptons] channel="<<channel<<" is weird"<<endl;
  }

  if(!lepton0 || !lepton1) return false;                         // Dilepton
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "Dilepton", map_weight[""]);
  if(lepton0->Pt() < l0pt || lepton1->Pt() < l1pt) return false; // pT cut
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "LepPt", map_weight[""]);
  if(lepton0->Charge() * lepton1->Charge() > 0) hprefix += "ss_"; // Opposite charge
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "Charge", map_weight[""]);
  if((*lepton0 + *lepton1).M() < 52) return false;               // Mass 52
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "Mass52", map_weight[""]);

  leptons ={};
  leptons.push_back(lepton0);
  leptons.push_back(lepton1);

  if(hprefix.Contains("ss_") && suffix != "") return false; // To reduce SS histograms

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
  }else if(DataYear == 2017 && channel.Contains("m"+GetEraShort())){
    triggers = {
      "HLT_IsoMu24_v",
      "HLT_IsoMu27_v"
    };
  }
  else if(DataYear == 2018 && channel.Contains("m"+GetEraShort())) triggers = {"HLT_IsoMu24_v"};
  else if(DataYear == 2016 && (channel.Contains("e"+GetEraShort()) || channel.Contains("E")+GetEraShort())) triggers = {"HLT_Ele27_WPTight_Gsf_v"};
  else if(DataYear == 2017 && (channel.Contains("e"+GetEraShort()) || channel.Contains("E")+GetEraShort())){
    triggers = {
      "HLT_Ele27_WPTight_Gsf_v",
      "HLT_Ele32_WPTight_Gsf_v"
    };
  }else if(DataYear == 2018 && (channel.Contains("e"+GetEraShort()) || channel.Contains("E")+GetEraShort())){
    triggers = {
      "HLT_Ele28_WPTight_Gsf_v",
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
    if(mu.IsType(Muon::Type::PFMuon)) bmuon.push_back(mu);
  }

  vector<Electron> belectron;
  for(const auto& el:allels){
    if(jet.DeltaR(el) > 0.3) continue;
    if(el.P() * sin(el.Angle(jet.Vect())) < 1.0) continue;
    if(!PassID(&el, "SoftElectronID")) continue;
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

  vector<Gen> gens = GetGens();
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
      zptweight = fZptCorrection->GetZptWeight(genZ.Pt(), genZ.Rapidity());
      zptweight_gym = fZptCorrection->GetZptWeight(genZ.Pt(), genZ.Rapidity(), genZ.M());
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
  if(IsTTSample) topptweight = mcCorr->GetTopPtReweight(gens);
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
      gRandom->SetSeed((run<<15)+(lumi<<10)+(event<<5)+p0->Eta()*100);
      if(gRandom->Rndm()<0.5){
        lm=p0;
        lp=p1;
      }else{
	lm=p1;
	lp=p0;
      }      
    } 
  }else{
    gRandom->SetSeed((run<<15)+(lumi<<10)+(event<<5)+p0->Eta()*100);
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

vector<vector<double>> dybAnalyzer::Make2DWeights(const vector<int>& structure){
  vector<vector<double>> rt;
  for(const int& nmem:structure){
    rt.push_back(vector<double>(nmem,1.0));
  }
  return rt;
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
  TFile fPUID(datapath+"/"+GetEra()+"/ID/PUJet/PUID_106XTraining_ULRun2_EffSFandUncties_v1.root");

  TString era = GetEra();
  if(era == "2016postVFP") era = "2016";
  else if(era == "2016preVFP") era = "2016APV";

  // Currently, Medium and Tight are not used due to discontinuity at pt = 50 GeV
  heff_sf = (TH2F*)fPUID.Get("h2_eff_sfUL"+era+"_L");
  heff_sf_unc = (TH2F*)fPUID.Get("h2_eff_sfUL"+era+"_L_Systuncty");

  heff_sf->SetDirectory(0);
  heff_sf_unc->SetDirectory(0);

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
      }else{
        if(jet.PileupJetId() > 0.71) return true;
      }
    }
    else if(ID == "Medium"){
      if(jet.Pt() >= 40){
        if(jet.PileupJetId() > 0.93) return true;
      }else if(jet.Pt() >= 30){
        if(jet.PileupJetId() > 0.86) return true;
      }else if(jet.Pt() >= 20){
        if(jet.PileupJetId() > 0.62) return true;
      }else{
        if(jet.PileupJetId() > 0.20) return true;
      }
    }
    else{
      if(jet.Pt() >= 40){
        if(jet.PileupJetId() > -0.42) return true;
      }else if(jet.Pt() >= 30){
        if(jet.PileupJetId() > -0.71) return true;
      }else if(jet.Pt() >= 20){
        if(jet.PileupJetId() > -0.90) return true;
      }else{
        if(jet.PileupJetId() > -0.95) return true;
      }
    }
  }

  else if(DataEra == "2017" || DataEra == "2018"){
    if(ID == "Tight"){
      if(jet.Pt() >= 40){
        if(jet.PileupJetId() > 0.98) return true;
      }else if(jet.Pt() >= 30){
        if(jet.PileupJetId() > 0.96) return true;
      }else if(jet.Pt() >= 20){
        if(jet.PileupJetId() > 0.90) return true;
      }else{
        if(jet.PileupJetId() > 0.77) return true;
      }
    }
    else if(ID == "Medium"){
      if(jet.Pt() >= 40){
        if(jet.PileupJetId() > 0.96) return true;
      }else if(jet.Pt() >= 30){
        if(jet.PileupJetId() > 0.90) return true;
      }else if(jet.Pt() >= 20){
        if(jet.PileupJetId() > 0.68) return true;
      }else{
        if(jet.PileupJetId() > 0.26) return true;
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
        }else{
          if(jet.PileupJetId() > -0.95) return true;
        }
      }else if(fabs(jet.Eta()) < 2.75){
        if(jet.Pt() >= 40){
          if(jet.PileupJetId() > 0.22) return true;
        }else if(jet.Pt() >= 30){
          if(jet.PileupJetId() > -0.18) return true;
        }else if(jet.Pt() >= 20){
          if(jet.PileupJetId() > -0.55) return true;
        }else{
          if(jet.PileupJetId() > -0.72) return true;
        }
      }else if(fabs(jet.Eta()) < 3.0){
        if(jet.Pt() >= 40){
          if(jet.PileupJetId() > -0.13) return true;
        }else if(jet.Pt() >= 30){
          if(jet.PileupJetId() > -0.43) return true;
        }else if(jet.Pt() >= 20){
          if(jet.PileupJetId() > -0.60) return true;
        }else{
          if(jet.PileupJetId() > -0.68) return true;
        }
      }else{
        if(jet.Pt() >= 40){
          if(jet.PileupJetId() > -0.03) return true;
        }else if(jet.Pt() >= 30){
          if(jet.PileupJetId() > -0.24) return true;
        }else if(jet.Pt() >= 20){
          if(jet.PileupJetId() > -0.43) return true;
        }else{
          if(jet.PileupJetId() > -0.47) return true;
        }
      }
    }
  }
  else cout<<"[dybAnalyzer::PUJetIDPass] era is weird "<<endl;

  return false;
}

double dybAnalyzer::GetPUJetWeight(const vector<Jet>& jets, TString ID, int sys){
  if(IsDATA) return 1.;

  double weight(1.), weight_unc2(0.);
  for(unsigned int i=0; i<jets.size(); i++){
    double jetpt = jets.at(i).Pt();
    double jeteta = jets.at(i).Eta();
    if(jets.at(i).Pt() < 20) cout<<"jet pt < 20GeV, something wrong"<<endl;;
    if(jets.at(i).Pt() > 50) continue;
    if(abs(jets.at(i).Eta()) > 2.5) continue;

    double this_effSF = heff_sf->GetBinContent(heff_sf->FindBin(jetpt, jeteta));
    double this_effSF_unc = heff_sf_unc->GetBinContent(heff_sf_unc->FindBin(jetpt, jeteta));

    bool isRealJet = false;
    //isRealJet = (jets.at(i).GenHFHadronMatcherFlavour() >= 0.);
    vector<Gen> gens = GetGens();
    double deltaR = 0.4;
    for(unsigned int j=0; j<gens.size(); j++){
      if(!gens.at(j).isHardProcess()) continue;
      if(gens.at(j).DeltaR(jets.at(i)) > deltaR) continue;
      isRealJet = true;
      break;
    }
    bool isPassID = PUJetIDPass(jets.at(i), ID);

    if(isRealJet){
      if(isPassID){
        weight *= this_effSF;
        weight_unc2 += this_effSF_unc * this_effSF_unc;
      }else{
        //Prob_DATA *= 1.-this_DATA_eff;
        //Prob_MC *= 1.-this_MC_eff;
      }
    }
  }

  return weight + sys * sqrt(weight_unc2);
}

double dybAnalyzer::GetbChargeSFWeight(const vector<Jet>& jets, unsigned int mode, int sys, TString bChargeBins){
  double weight = 1.;
  if(IsDATA) return weight;

  vector<Gen> gens = GetGens();
  for(const auto& jet:jets){
    int genpid = 0;
    double dR = 99.;

    for(unsigned int i=0; i<gens.size(); i++){
      if(!gens.at(i).isPrompt()) continue;
      if(!gens.at(i).isHardProcess()) continue;
      if(fabs(gens.at(i).PID()) != 5) continue;
      if(gens.at(i).DeltaR(jet) > dR) continue;
      dR = gens.at(i).DeltaR(jet);
      if(dR < 0.4) genpid = gens.at(i).PID();
    }
    if(genpid == 0) continue; // No matched b-partons

    double Charge = jetCharge(jet); // No Normalization & TTJJ
    double alpha_plus_DATA_eff = 0.63181282; //0.63079885;
    double alpha_minus_DATA_eff = 0.61395198; //0.61295115;
    double alpha_plus_MC_eff = 0.65415066;
    double alpha_minus_MC_eff = 0.63769426;
    if(sys > 0){ // asym gets smaller
      alpha_plus_DATA_eff += -0.00088663; //-0.00086244;
      alpha_minus_DATA_eff += 0.00074632; //0.00072373;
    }else if(sys < 0){ // SF gets smaller
      alpha_plus_DATA_eff += 0.00250343; //0.00237556;
      alpha_minus_DATA_eff += 0.00297406; //0.00283084;
    }

    // Hyonsans's SF
    if(mode == 2){
      alpha_plus_DATA_eff = 0.63597127;
      alpha_minus_DATA_eff = 0.62062727;
      alpha_plus_MC_eff = 0.65465339;
      alpha_minus_MC_eff = 0.63993788;
      if(sys > 0){
        alpha_plus_DATA_eff += 0.00107148;
        alpha_minus_DATA_eff += -0.00107671;
      }else if(sys < 0){
        alpha_plus_DATA_eff += -0.00303157;
        alpha_minus_DATA_eff += -0.00304636;
      }
    }

    // 1D bChargeSF
    if(mode == 1){
      if(fabs(Charge) < afb_chbin[1]){// 0.1
        alpha_plus_DATA_eff = 0.52466892; //0.52455235;
        alpha_minus_DATA_eff = 0.51791587; //0.51779908;
        alpha_plus_MC_eff = 0.52968978;
        alpha_minus_MC_eff = 0.52540184;
        if(sys > 0 && bChargeBins.Contains("0")){
          alpha_plus_DATA_eff += -0.00159007; //-0.00157464;
          alpha_minus_DATA_eff += 0.0021877; //0.00209931;
        }else if(sys < 0 && bChargeBins.Contains("0")){
          alpha_plus_DATA_eff += 0.0029865; //0.00287364;
          alpha_minus_DATA_eff += 0.00217065; //0.00215546;
        }
      }else if(fabs(Charge) < afb_chbin[2]){// 0.2
        alpha_plus_DATA_eff = 0.57549775; //0.57503119;
        alpha_minus_DATA_eff = 0.56120889; //0.56072413;
        alpha_plus_MC_eff = 0.58939608;
        alpha_minus_MC_eff = 0.57641132;
        if(sys > 0 && bChargeBins.Contains("1")){
          alpha_plus_DATA_eff += -0.00163046; //-0.00176428;
          alpha_minus_DATA_eff += 0.00239976; //0.00222184;
        }else if(sys < 0 && bChargeBins.Contains("1")){
          alpha_plus_DATA_eff += 0.00365694; //0.00332861;
          alpha_minus_DATA_eff += 0.00248462; //0.00264313;
        }
      }else if(fabs(Charge) < afb_chbin[3]){// 0.6
        alpha_plus_DATA_eff = 0.66795153; //0.66685555;
        alpha_minus_DATA_eff = 0.64183908; //0.6406924;
        alpha_plus_MC_eff = 0.7016705;
        alpha_minus_MC_eff = 0.67580431;
        if(sys > 0 && bChargeBins.Contains("2")){
          alpha_plus_DATA_eff += -0.00198435; //-0.00192819;
          alpha_minus_DATA_eff += 0.00098174; //0.00094884;
        }else if(sys < 0 && bChargeBins.Contains("2")){
          alpha_plus_DATA_eff += 0.00278728; //0.00262054;
          alpha_minus_DATA_eff += 0.00563382; //0.00532535;
        }
      }else if(fabs(Charge) < afb_chbin[4]){// 1.0
        alpha_plus_DATA_eff = 0.78691713; //0.78271638;
        alpha_minus_DATA_eff = 0.74175893; //0.73778292;
        alpha_plus_MC_eff = 0.82685472;
        alpha_minus_MC_eff = 0.79243187;
        if(sys > 0 && bChargeBins.Contains("3")){
          alpha_plus_DATA_eff += -0.00277093; //-0.00284074;
          alpha_minus_DATA_eff += 0.00567753; //0.00513936;
        }else if(sys < 0 && bChargeBins.Contains("3")){
          alpha_plus_DATA_eff += 0.01097872; //0.01040963;
          alpha_minus_DATA_eff += 0.0053582; //0.00575384;
        }
      }else if(fabs(Charge) < afb_chbin[5]){// 3.0, soft muons
        alpha_plus_DATA_eff = 0.77490976; //0.77182254;
        alpha_minus_DATA_eff = 0.77647107; //0.77357429;
        alpha_plus_MC_eff = 0.78518778;
        alpha_minus_MC_eff = 0.78552888;
        if(sys > 0 && bChargeBins.Contains("4")){
          alpha_plus_DATA_eff += -0.00930939; //-0.00300315;
          alpha_minus_DATA_eff += 0.00629106; //0.01019522;
        }else if(sys < 0 && bChargeBins.Contains("4")){
          alpha_plus_DATA_eff += 0.00671166; //0.01005651;
          alpha_minus_DATA_eff += 0.00993179; //0.00296229;
        }
      }else{// 5.0, soft electrons
        alpha_plus_DATA_eff = 0.75926086; //0.75433632;
        alpha_minus_DATA_eff = 0.75920095; //0.75499201;
        alpha_plus_MC_eff = 0.76063023;
        alpha_minus_MC_eff = 0.76337732;
        if(sys > 0 && bChargeBins.Contains("5")){
          alpha_plus_DATA_eff += -0.00459789; //-0.00443982;
          alpha_minus_DATA_eff += 0.00433558; //0.00398256;
        }else if(sys < 0 && bChargeBins.Contains("5")){
          alpha_plus_DATA_eff += 0.00822754; //0.00811232;
          alpha_minus_DATA_eff += 0.00872532; //0.00904374;
        }
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

    /*
    int origin=jet.GenHFHadronMatcherOrigin();
    if(origin==-999) continue;
    if(origin*Charge<0){
      if(origin>0){
        weight *= alpha_minus_DATA_eff / alpha_minus_MC_eff;
      }else{
        weight *= alpha_plus_DATA_eff / alpha_plus_MC_eff;
      }
    }else{
      if(origin>0){
        weight *= (1. - alpha_minus_DATA_eff) / (1. - alpha_minus_MC_eff);
      }else{
        weight *= (1. - alpha_plus_DATA_eff) / (1. - alpha_plus_MC_eff);
      }
    }
    */
  }

  return weight;
}
