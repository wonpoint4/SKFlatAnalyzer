#include "ttljAnalyzer.h"

void ttljAnalyzer::initializeAnalyzer(){
  SMPAnalyzerCore::initializeAnalyzer(); //setup eff zpt roc z0
  IsNominalRun = !HasFlag("SYS") && !HasFlag("PDFSYS") && !HasFlag("LEPSYS");

  PDFbase = LHAPDF::mkPDF(306000);
  PDFnf4 = LHAPDF::mkPDF(325500);
  mcCorr->SetJetTaggingParameters({
    JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Tight, JetTagging::incl, JetTagging::comb),
  });
  SetupPUJetWeight();
  SetupJetVetoMap();
  SetupLikelihoods(0, 1);
}

void ttljAnalyzer::executeEvent(){
  ///////////////// GEN level /////////////////////
  if(!IsDATA) executeEventGen();

  ///////////////// RECO level /////////////////////
  jets_raw = GetAllJets(); // This can make event loops much slower in case of running over Unskimmed samples
  if(!IsDATA || DataStream.Contains("SingleMuon")){
    muons = MuonMomentumCorrection(SMPGetMuons("POGMediumWithLooseTrkIso", 8.0, 2.4), 0,0, true);
    executeEventWithParameter("m"+GetEraShort());
    if(HasFlag("SYS")){
      for(TString syst:{"_jet_scale_up", "_jet_scale_down", "_jet_smear_up", "_jet_smear_down"}){
        if(syst.Contains("scale") || !IsDATA) executeEventWithParameter("m"+GetEraShort(), syst);
      }
    }else if(HasFlag("LEPSYS")){
      for(unsigned int s=0; s<nmem_muon.size(); s++){
        for(unsigned int m=0; m<nmem_muon.at(s); m++){
          executeEventWithParameter("m"+GetEraShort(), Form("_MuonMomentum_s%dm%d", s, m), s,m);
        }
      }
    }else{
      for(TString syst:{"_jetpt25", "_jetpt30", "_jetpt35", "_jetpt45", "_jetpt50", "_jetpt55", "_jeteta5", "_jeteta1p5", "_HS", "_noSelQ"}){
        executeEventWithParameter("m"+GetEraShort(), syst);
      }
    }
  }
  if(!IsDATA || DataStream.Contains("SingleElectron") || DataStream.Contains("EGamma")){
    electrons = ElectronEnergyCorrection(SMPGetElectrons("passMediumID_SelQ", 8.0, 2.5), 0,0, true);
    executeEventWithParameter("e"+GetEraShort());
    if(HasFlag("SYS")){
      for(TString syst:{"_jet_scale_up", "_jet_scale_down", "_jet_smear_up", "_jet_smear_down"}){
        if(syst.Contains("scale") || !IsDATA) executeEventWithParameter("e"+GetEraShort(), syst);
      }
    }else if(HasFlag("LEPSYS")){
      for(unsigned int s=0; s<nmem_electron.size(); s++){
        for(unsigned int m=0; m<nmem_electron.at(s); m++){
          executeEventWithParameter("e"+GetEraShort(), Form("_ElectronEnergy_s%dm%d", s, m), s,m);
        }
      }
    }else{
      for(TString syst:{"_jetpt25", "_jetpt30", "_jetpt35", "_jetpt45", "_jetpt50", "_jetpt55", "_jeteta5", "_jeteta1p5"}){
        executeEventWithParameter("e"+GetEraShort(), syst);
      }
      electrons = ElectronEnergyCorrection(SMPGetElectrons("passMediumID", 8.0, 2.5), 0,0, true);
      for(TString syst:{"_HS", "_noSelQ"}){
        executeEventWithParameter("e"+GetEraShort(), syst);
      }
    }
  }
}

void ttljAnalyzer::executeEventWithParameter(TString channel, TString option, unsigned int set, unsigned int mem){

  lepton0 = NULL;
  prefix = channel+"/"+gprefix, hprefix = "", suffix = "";
  suffix += option;
  if((IsNominalRun || option != "") && set != 1) IsNominalLike = true;
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
  if(!option.Contains("HS")){
    if(!HasLeptons(channel, false, set, mem)) return;
  }else{
    if(!HasLeptons(channel, true, set, mem)) return;
  }

  // Jets
  vector<Jet> alljets = {};
  if(option.Contains("jet_scale_up")) alljets = SelectJets(ScaleJets(jets_raw, 1), "tightLepVeto", 25, (DataYear == 2016? 2.4: 2.5));
  else if(option.Contains("jet_scale_down")) alljets = SelectJets(ScaleJets(jets_raw, -1), "tightLepVeto", 25, (DataYear == 2016? 2.4: 2.5));
  else if(option.Contains("jet_smear_up")) alljets = SelectJets(SmearJets(jets_raw, 1), "tightLepVeto", 25, (DataYear == 2016? 2.4: 2.5));
  else if(option.Contains("jet_smear_down")) alljets = SelectJets(SmearJets(jets_raw, -1), "tightLepVeto", 25, (DataYear == 2016? 2.4: 2.5));
  else if(option.Contains("jeteta5")) alljets = SelectJets(SmearJets(jets_raw, 0), "tightLepVeto", 25, 5.0);
  else if(option.Contains("HS")) alljets = SelectJets(jets_raw, "tightLepVeto", 40, 2.4);
  else alljets = SelectJets(jets_raw, "tightLepVeto", 25, (DataYear == 2016? 2.4: 2.5));
  std::sort(alljets.begin(), alljets.end(), PtComparing);

  vector<Jet> lepvetojets = {}, realjets_before_vetomap = {}, realjets = {}, bjets = {}, ajets = {};
  for(const auto& jet:alljets){
    if(lepton0 && jet.DeltaR(*lepton0) < 0.4) continue;
    lepvetojets.push_back(jet);
  }
  for(const auto& jet:lepvetojets){
    if(!option.Contains("HS")){
      if(!PUJetIDPass(jet, "Loose")) continue;
    }
    realjets_before_vetomap.push_back(jet);
    if(IsNominalRun){
      FillHist(prefix+hprefix+"etaphi_realjets_before_vetomap"+suffix, jet.Eta(), jet.Phi(), 1., 200,-5.,5., 200,-3.2,3.2);
      FillHist("etaphi_realjets_before_vetomap", jet.Eta(), jet.Phi(), 1., 200,-5.,5., 200,-3.2,3.2);
    }
  }
  for(const auto& jet:realjets_before_vetomap){
    if(hvetomap->GetBinContent(hvetomap->FindBin(jet.Eta(), jet.Phi()))) continue;
    realjets.push_back(jet);
    if(IsNominalRun){
      FillHist(prefix+hprefix+"etaphi_realjets"+suffix, jet.Eta(), jet.Phi(), 1., 200,-5.,5., 200,-3.2,3.2);
      FillHist("etaphi_realjets", jet.Eta(), jet.Phi(), 1., 200,-5.,5., 200,-3.2,3.2);
    }
  }

  // b-tagging
  JetTagging::Parameters DeepJet_Tight = JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Tight, JetTagging::incl, JetTagging::comb);
  for(const auto& jet:realjets){
    //jet *= jet.BJetNNCorrection(); // full bJetEnergyCorrection (BBjetRegression)?
    if(jet.GetTaggerResult(DeepJet_Tight.j_Tagger) > mcCorr->GetJetTaggingCutValue(DeepJet_Tight.j_Tagger, DeepJet_Tight.j_WP) && fabs(jet.Eta()) < (DataYear == 2016? 2.4: 2.5)) bjets.push_back(jet);
    else{
      if(option.Contains("jetpt25")) ajets.push_back(jet);
      else if(option.Contains("jetpt30")){
        if(jet.Pt() > 30) ajets.push_back(jet);
      }else if(option.Contains("jetpt35")){
        if(jet.Pt() > 35) ajets.push_back(jet);
      }else if(option.Contains("jetpt45")){
        if(jet.Pt() > 45) ajets.push_back(jet);
      }else if(option.Contains("jetpt50")){
        if(jet.Pt() > 50) ajets.push_back(jet);
      }else if(option.Contains("jetpt55")){
        if(jet.Pt() > 55) ajets.push_back(jet);
      }else if(option.Contains("jeteta1p5")){
        if(jet.Pt() > 40 && fabs(jet.Eta()) < 1.5) ajets.push_back(jet);
      }else{ // default or "jeteta5"
        if(jet.Pt() > 40) ajets.push_back(jet);
      }
    }
  }

  // Jet related weights
  pujetSF = 1., btagSF = 1.;
  if(!IsDATA){
    pujetSF = GetPUJetWeight(lepvetojets, "Loose", 0);
    btagSF = mcCorr->GetBTaggingReweight_1a(realjets, DeepJet_Tight);
  }

  if(IsNominalLike){
    FillHist(prefix+hprefix+"njets_first"+suffix, realjets.size(), map_weight[""], 15,0,15);
    FillHist(prefix+hprefix+"nbjets_first"+suffix, bjets.size(), map_weight[""], 10,0,10);
    FillHist(prefix+hprefix+"najets_first"+suffix, ajets.size(), map_weight[""], 10,0,10);
  }

  if(!option.Contains("HS")){
    if(bjets.size() != 2) return;
    if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "Tight2b", map_weight[""]);
  }else{
    if(bjets.size() < 2) return;
    if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, ">=2b", map_weight[""]);
  }

  int bbCharges = 0;
  double b0_charge = jetCharge(bjets.at(0));
  double b1_charge = jetCharge(bjets.at(1));
  if(b0_charge > 0) bbCharges += 1<<0;
  if(b1_charge > 0) bbCharges += 1<<1;

  //bjets.at(0) *= bjets.at(0).BJetNNCorrection(); // bJetEnergyCorrection (BjetRegression)?
  //bjets.at(1) *= bjets.at(1).BJetNNCorrection(); // bJetEnergyCorrection (BjetRegression)?

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
    TString SelQ = "_SelQ";
    if(option.Contains("HS") || option.Contains("noSelQ")) SelQ = "";
    TString muonTrackingSF_key = "Muon_Tracking";
    TString muonRECOSF_key = "Muon_RECO";
    TString muonIDSF_key = "Muon_MediumID_trkIsoLoose";
    TString muonTriggerSF_key = "IsoMu24_MediumID_trkIsoLoose";
    TString electronRECOSF_key = "Electron_RECO";
    TString electronIDSF_key = "Electron"+SelQ+"_MediumID";
    TString electronTriggerSF_key = DataYear == 2018? "Ele28"+SelQ+"_MediumID": "Ele27"+SelQ+"_MediumID";

    if(HasFlag("LEPSYS") && option == ""){
      muonTrackingSF_sys = Make2DWeights(fEff->GetStructure(muonTrackingSF_key));
      muonRECOSF_sys = Make2DWeights(fEff->GetStructure(muonRECOSF_key));
      muonIDSF_sys = Make2DWeights(fEff->GetStructure(muonIDSF_key));
      muonTriggerSF_sys = Make2DWeights(fEff->GetStructure(muonTriggerSF_key));
      electronRECOSF_sys = Make2DWeights(fEff->GetStructure(electronRECOSF_key));
      electronIDSF_sys = Make2DWeights(fEff->GetStructure(electronIDSF_key));
      electronTriggerSF_sys = Make2DWeights(fEff->GetStructure(electronTriggerSF_key));
    }

    // Tracking, RECO, ID SF per lepton
    for(const Lepton* lepton:leptons){
      if(lepton->LeptonFlavour() == Lepton::MUON){
        muonTrackingSF *= fEff->GetEfficiencySF(muonTrackingSF_key, lepton, 0,0);
        muonRECOSF *= fEff->GetEfficiencySF(muonRECOSF_key, lepton, 0,0);
        muonIDSF *= fEff->GetEfficiencySF(muonIDSF_key, lepton, 0,0);

        if(HasFlag("LEPSYS") && option == ""){
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

        if(HasFlag("LEPSYS") && option == ""){
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
    if(channel.Contains("m"+GetEraShort())){
      if(DataYear != 2017){
        muonTriggerSF *= GetLeptonTriggerSF(muonTriggerSF_key, leptons, 0,0);
        if(HasFlag("LEPSYS") && option == ""){
          for(unsigned int s=0; s<muonTriggerSF_sys.size(); s++){
            for(unsigned int m=0; m<muonTriggerSF_sys[s].size(); m++){
              muonTriggerSF_sys[s][m] *= GetLeptonTriggerSF(muonTriggerSF_key, leptons, s,m);
            }
          }
        }
      }else{
        muonTriggerSF *= GetLeptonTriggerORSF({"IsoMu24_MediumID_trkIsoLoose", "IsoMu27_MediumID_trkIsoLoose"}, leptons, 0,0);
        if(HasFlag("LEPSYS") && option == ""){
          for(unsigned int s=0; s<muonTriggerSF_sys.size(); s++){
            for(unsigned int m=0; m<muonTriggerSF_sys[s].size(); m++){
              muonTriggerSF_sys[s][m] *= GetLeptonTriggerORSF({"IsoMu24_MediumID_trkIsoLoose", "IsoMu27_MediumID_trkIsoLoose"}, leptons, s,m);
            }
          }
        }
      }
    }
    if(channel.Contains("e"+GetEraShort())){
      if(DataYear == 2016){
        electronTriggerSF *= GetLeptonTriggerSF(electronTriggerSF_key, leptons, 0,0);
        if(HasFlag("LEPSYS") && option == ""){
          for(unsigned int s=0; s<electronTriggerSF_sys.size(); s++){
            for(unsigned int m=0; m<electronTriggerSF_sys[s].size(); m++){
              electronTriggerSF_sys[s][m] *= GetLeptonTriggerSF(electronTriggerSF_key, leptons, s,m);
            }
          }
        }
      }else if(DataYear == 2017){
        electronTriggerSF *= GetLeptonTriggerORSF({"Ele27"+SelQ+"_MediumID", "Ele32"+SelQ+"_MediumID"}, leptons, 0,0);
        if(HasFlag("LEPSYS") && option == ""){
          for(unsigned int s=0; s<electronTriggerSF_sys.size(); s++){
            for(unsigned int m=0; m<electronTriggerSF_sys[s].size(); m++){
              electronTriggerSF_sys[s][m] *= GetLeptonTriggerORSF({"Ele27"+SelQ+"_MediumID", "Ele32"+SelQ+"_MediumID"}, leptons, s,m);
            }
          }
        }
      }else{
        electronTriggerSF *= GetLeptonTriggerORSF({"Ele28"+SelQ+"_MediumID", "Ele32"+SelQ+"_MediumID"}, leptons, 0,0);
        if(HasFlag("LEPSYS") && option == ""){
          for(unsigned int s=0; s<electronTriggerSF_sys.size(); s++){
            for(unsigned int m=0; m<electronTriggerSF_sys[s].size(); m++){
              electronTriggerSF_sys[s][m] *= GetLeptonTriggerORSF({"Ele28"+SelQ+"_MediumID", "Ele32"+SelQ+"_MediumID"}, leptons, s,m);
            }
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

  double bChargeSF0 = GetbChargeSFWeight(bjets, 0, 0);
  double bChargeSFHS = GetbChargeSFWeight(bjets, 2, 0);
  double bChargeSF1 = GetbChargeSFWeight(bjets, 1, 0);

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
    map_weight["_btagSF_huncorr"+GetEraShort()] = map_weight[""] / btagSF * mcCorr->GetBTaggingReweight_1a(realjets, DeepJet_Tight, "SystUpHTagUnCorr");
    map_weight["_btagSF_lup"] = map_weight[""] / btagSF * mcCorr->GetBTaggingReweight_1a(realjets, DeepJet_Tight, "SystUpLTag");
    map_weight["_btagSF_ldown"] = map_weight[""] / btagSF * mcCorr->GetBTaggingReweight_1a(realjets, DeepJet_Tight, "SystDownLTag");
    map_weight["_btagSF_lcorr"] = map_weight[""] / btagSF * mcCorr->GetBTaggingReweight_1a(realjets, DeepJet_Tight, "SystUpLTagCorr");;
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
    for(TString bCh:{"0", "1", "2", "3", "4", "5"}){
      map_weight["_bChargeSF1_up"+bCh] = map_weight[""] * GetbChargeSFWeight(bjets, 1, 1, bCh);
      map_weight["_bChargeSF1_down"+bCh] = map_weight[""] * GetbChargeSFWeight(bjets, 1, -1, bCh);
    }
    map_weight["_bChargeSFHS"] = map_weight[""] * GetbChargeSFWeight(bjets, 2, 0);
    map_weight["_bChargeSFHS_up"] = map_weight[""] * GetbChargeSFWeight(bjets, 2, 1);
    map_weight["_bChargeSFHS_down"] = map_weight[""] * GetbChargeSFWeight(bjets, 2, -1);

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
  }else if(!IsDATA && HasFlag("LEPSYS") && option == ""){
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
    if(channel.Contains("m"+GetEraShort())){
      for(unsigned int i=0; i<nmem_electron.size(); i++){
        for(unsigned int j=0; j<nmem_electron.at(i); j++){
          map_weight[Form("_ElectronEnergy_s%dm%d", i, j)] = map_weight[""];
        }
      }
    }
    if(channel.Contains("e"+GetEraShort())){
      for(unsigned int i=0; i<nmem_muon.size(); i++){
        for(unsigned int j=0; j<nmem_muon.at(i); j++){
          map_weight[Form("_MuonMomentum_s%dm%d", i, j)] = map_weight[""];
        }
      }
    }
  }else if(!IsDATA && IsNominalRun) map_weight["_bChargeSF1"] = map_weight[""] * bChargeSF1;

  if((HasFlag("SYS") || HasFlag("PDFSYS") || HasFlag("LEPSYS")) && option == "") map_weight.erase("");

  if(IsNominalLike){
    FillHist(prefix+hprefix+"njets_incTT"+suffix, realjets.size(), map_weight, 15,0,15);
    FillHist(prefix+hprefix+"nbjets_incTT"+suffix, bjets.size(), map_weight, 10,0,10);
    FillHist(prefix+hprefix+"najets_incTT"+suffix, ajets.size(), map_weight, 10,0,10);

    FillHist(prefix+hprefix+"lpt_incTT"+suffix, lepton0->Pt(), map_weight, 100,0,200);
    FillHist(prefix+hprefix+"leta_incTT"+suffix, lepton0->Eta(), map_weight, 100,-5,5);
    FillHist(prefix+hprefix+"bpt_incTT"+suffix, bjets.at(0).Pt(), map_weight, 100,0,200);
    FillHist(prefix+hprefix+"bpt_incTT"+suffix, bjets.at(1).Pt(), map_weight, 100,0,200);
    FillHist(prefix+hprefix+"beta_incTT"+suffix, bjets.at(0).Eta(), map_weight, 100,-5,5);
    FillHist(prefix+hprefix+"beta_incTT"+suffix, bjets.at(1).Eta(), map_weight, 100,-5,5);
    FillHist(prefix+hprefix+"bphi_incTT"+suffix, bjets.at(0).Phi(), map_weight, 64,-3.2,3.2);
    FillHist(prefix+hprefix+"bphi_incTT"+suffix, bjets.at(1).Phi(), map_weight, 64,-3.2,3.2);
    FillHist(prefix+hprefix+"met_incTT"+suffix, PuppiMET_Type1_pt, map_weight, 200,0,200);
    FillHist(prefix+hprefix+"metphi_incTT"+suffix, PuppiMET_Type1_phi, map_weight, 64,-3.2,3.2);
    FillHist(prefix+hprefix+"met_metphi_incTT"+suffix, PuppiMET_Type1_pt, PuppiMET_Type1_phi, map_weight[""], 40,0,200, 32,-3.2,3.2);
    FillHist(prefix+hprefix+"bbCharges_incTT"+suffix, bbCharges, map_weight, 4,0,4);
  }

  if(!option.Contains("HS")){
    if(ajets.size() < 2) return;
    if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, ">=2j", map_weight[""]);
  }else return;

  //==== Making Likelihood
  if(IsNominalLike){
    if(!IsDATA && MCSample.Contains("TTLJ")) FillingLikelihood(bjets, ajets, map_weight[""], suffix, 0, 0); // ByungHun Oh's method - drop events with ambiguity
    if(!IsDATA && MCSample.Contains("TTLJ")) FillingLikelihood(bjets, ajets, map_weight[""], suffix, 1, 0); // Charmonium guy's method - match smaller dR < 0.3
  }

  if(IsNominalLike){
    FillHist(prefix+hprefix+"lpt_incTTLJ"+suffix, lepton0->Pt(), map_weight, 100,0,200);
    FillHist(prefix+hprefix+"leta_incTTLJ"+suffix, lepton0->Eta(), map_weight, 100,-5,5);
    FillHist(prefix+hprefix+"bpt_incTTLJ"+suffix, bjets.at(0).Pt(), map_weight, 100,0,200);
    FillHist(prefix+hprefix+"bpt_incTTLJ"+suffix, bjets.at(1).Pt(), map_weight, 100,0,200);
    FillHist(prefix+hprefix+"beta_incTTLJ"+suffix, bjets.at(0).Eta(), map_weight, 100,-5,5);
    FillHist(prefix+hprefix+"beta_incTTLJ"+suffix, bjets.at(1).Eta(), map_weight, 100,-5,5);
    FillHist(prefix+hprefix+"bphi_incTTLJ"+suffix, bjets.at(0).Phi(), map_weight, 64,-3.2,3.2);
    FillHist(prefix+hprefix+"bphi_incTTLJ"+suffix, bjets.at(1).Phi(), map_weight, 64,-3.2,3.2);
    FillHist(prefix+hprefix+"jpt_incTTLJ"+suffix, ajets.at(0).Pt(), map_weight, 100,0,200);
    FillHist(prefix+hprefix+"jpt_incTTLJ"+suffix, ajets.at(1).Pt(), map_weight, 100,0,200);
    FillHist(prefix+hprefix+"jeta_incTTLJ"+suffix, ajets.at(0).Eta(), map_weight, 100,-5,5);
    FillHist(prefix+hprefix+"jeta_incTTLJ"+suffix, ajets.at(1).Eta(), map_weight, 100,-5,5);
    FillHist(prefix+hprefix+"jphi_incTTLJ"+suffix, ajets.at(0).Phi(), map_weight, 64,-3.2,3.2);
    FillHist(prefix+hprefix+"jphi_incTTLJ"+suffix, ajets.at(1).Phi(), map_weight, 64,-3.2,3.2);
    FillHist(prefix+hprefix+"met_incTTLJ"+suffix, PuppiMET_Type1_pt, map_weight, 200,0,200);
    FillHist(prefix+hprefix+"metphi_incTTLJ"+suffix, PuppiMET_Type1_phi, map_weight, 64,-3.2,3.2);
    FillHist(prefix+hprefix+"met_metphi_incTTLJ"+suffix, PuppiMET_Type1_pt, PuppiMET_Type1_phi, map_weight[""], 40,0,200, 32,-3.2,3.2);
    FillHist(prefix+hprefix+"bbCharges_incTTLJ"+suffix, bbCharges, map_weight, 4,0,4);
  }

  //==== Finding the correct bbjj combination
  vector<unsigned int> idx_bbjj = {0, 0, 0, 0};
  vector<double> LRs = {0, 0, 0, 0, 0, 0};
  idx_bbjj = Finding_bbjj_byLikelihood(channel, bjets, ajets, LRs);

  bool goodKinematic = true;
  if((idx_bbjj.at(0) == idx_bbjj.at(1)) || (idx_bbjj.at(2) == idx_bbjj.at(3))) goodKinematic = false;
  if(IsNominalLike) FillHist(prefix+hprefix+"LR_efficiency"+suffix, (goodKinematic? 1: 0), map_weight[""], 2,0,2);
  if(!goodKinematic) return;
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "Kin. cuts", map_weight[""]);

  if(IsNominalLike){
    FillHist(prefix+hprefix+"likelihood_ratio_beforeLRcut"+suffix, LRs.at(0) * LRs.at(1) * LRs.at(2) * LRs.at(3), map_weight, 1000,0,50);
    FillHist(prefix+hprefix+"likelihood_ratio_Mbl_beforeLRcut"+suffix, LRs.at(0), map_weight, 100,0,5);
    FillHist(prefix+hprefix+"likelihood_ratio_MblMET_beforeLRcut"+suffix, LRs.at(1), map_weight, 100,0,5);
    FillHist(prefix+hprefix+"likelihood_ratio_Mbjj_beforeLRcut"+suffix, LRs.at(2), map_weight, 100,0,5);
    FillHist(prefix+hprefix+"likelihood_ratio_Mjj_beforeLRcut"+suffix, LRs.at(3), map_weight, 100,0,5);
    FillHist(prefix+hprefix+"likelihood_ratio_dRtt_beforeLRcut"+suffix, LRs.at(4), map_weight, 100,0,5);
    FillHist(prefix+hprefix+"likelihood_ratio_dPhitt_beforeLRcut"+suffix, LRs.at(5), map_weight, 100,0,5);

    FillHist(prefix+hprefix+"lpt_beforeLRcut"+suffix, lepton0->Pt(), map_weight, 100,0,200);
    FillHist(prefix+hprefix+"leta_beforeLRcut"+suffix, lepton0->Eta(), map_weight, 100,-5,5);
    FillHist(prefix+hprefix+"bpt_beforeLRcut"+suffix, bjets.at(0).Pt(), map_weight, 100,0,200);
    FillHist(prefix+hprefix+"bpt_beforeLRcut"+suffix, bjets.at(1).Pt(), map_weight, 100,0,200);
    FillHist(prefix+hprefix+"beta_beforeLRcut"+suffix, bjets.at(0).Eta(), map_weight, 100,-5,5);
    FillHist(prefix+hprefix+"beta_beforeLRcut"+suffix, bjets.at(1).Eta(), map_weight, 100,-5,5);
    FillHist(prefix+hprefix+"bphi_beforeLRcout"+suffix, bjets.at(0).Phi(), map_weight, 64,-3.2,3.2);
    FillHist(prefix+hprefix+"bphi_beforeLRcout"+suffix, bjets.at(1).Phi(), map_weight, 64,-3.2,3.2);
    FillHist(prefix+hprefix+"jpt_beforeLRcut"+suffix, ajets.at(0).Pt(), map_weight, 100,0,200);
    FillHist(prefix+hprefix+"jpt_beforeLRcut"+suffix, ajets.at(1).Pt(), map_weight, 100,0,200);
    FillHist(prefix+hprefix+"jeta_beforeLRcut"+suffix, ajets.at(0).Eta(), map_weight, 100,-5,5);
    FillHist(prefix+hprefix+"jeta_beforeLRcut"+suffix, ajets.at(1).Eta(), map_weight, 100,-5,5);
    FillHist(prefix+hprefix+"jphi_beforeLRcout"+suffix, ajets.at(0).Phi(), map_weight, 64,-3.2,3.2);
    FillHist(prefix+hprefix+"jphi_beforeLRcout"+suffix, ajets.at(1).Phi(), map_weight, 64,-3.2,3.2);
    FillHist(prefix+hprefix+"met_beforeLRcout"+suffix, PuppiMET_Type1_pt, map_weight, 200,0,200);
    FillHist(prefix+hprefix+"metphi_beforeLRcout"+suffix, PuppiMET_Type1_phi, map_weight, 64,-3.2,3.2);
    FillHist(prefix+hprefix+"met_metphi_beforeLRcout"+suffix, PuppiMET_Type1_pt, PuppiMET_Type1_phi, map_weight[""], 40,0,200, 32,-3.2,3.2);
    FillHist(prefix+hprefix+"bbCharges_beforeLRcut"+suffix, bbCharges, map_weight, 4,0,4);
  }

  if(LRs.at(0) * LRs.at(1) * LRs.at(2) * LRs.at(3) < 0.5) return;
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "LR0p5 cuts", map_weight[""]);

  if(IsNominalLike){
    FillHist(prefix+hprefix+"lpt_afterLRcut"+suffix, lepton0->Pt(), map_weight, 100,0,200);
    FillHist(prefix+hprefix+"leta_afterLRcut"+suffix, lepton0->Eta(), map_weight, 100,-5,5);
    FillHist(prefix+hprefix+"bpt_afterLRcut"+suffix, bjets.at(0).Pt(), map_weight, 100,0,200);
    FillHist(prefix+hprefix+"bpt_afterLRcut"+suffix, bjets.at(1).Pt(), map_weight, 100,0,200);
    FillHist(prefix+hprefix+"beta_afterLRcut"+suffix, bjets.at(0).Eta(), map_weight, 100,-5,5);
    FillHist(prefix+hprefix+"beta_afterLRcut"+suffix, bjets.at(1).Eta(), map_weight, 100,-5,5);
    FillHist(prefix+hprefix+"bphi_afterLRcout"+suffix, bjets.at(0).Phi(), map_weight, 64,-3.2,3.2);
    FillHist(prefix+hprefix+"bphi_afterLRcout"+suffix, bjets.at(1).Phi(), map_weight, 64,-3.2,3.2);
    FillHist(prefix+hprefix+"jpt_afterLRcut"+suffix, ajets.at(0).Pt(), map_weight, 100,0,200);
    FillHist(prefix+hprefix+"jpt_afterLRcut"+suffix, ajets.at(1).Pt(), map_weight, 100,0,200);
    FillHist(prefix+hprefix+"jeta_afterLRcut"+suffix, ajets.at(0).Eta(), map_weight, 100,-5,5);
    FillHist(prefix+hprefix+"jeta_afterLRcut"+suffix, ajets.at(1).Eta(), map_weight, 100,-5,5);
    FillHist(prefix+hprefix+"jphi_afterLRcout"+suffix, ajets.at(0).Phi(), map_weight, 64,-3.2,3.2);
    FillHist(prefix+hprefix+"jphi_afterLRcout"+suffix, ajets.at(1).Phi(), map_weight, 64,-3.2,3.2);
    FillHist(prefix+hprefix+"met_afterLRcout"+suffix, PuppiMET_Type1_pt, map_weight, 200,0,200);
    FillHist(prefix+hprefix+"metphi_afterLRcout"+suffix, PuppiMET_Type1_phi, map_weight, 64,-3.2,3.2);
    FillHist(prefix+hprefix+"met_metphi_afterLRcout"+suffix, PuppiMET_Type1_pt, PuppiMET_Type1_phi, map_weight[""], 40,0,200, 32,-3.2,3.2);
    FillHist(prefix+hprefix+"bbCharges_afterLRcut"+suffix, bbCharges, map_weight, 4,0,4);
  }

  TString prefix_nPV = prefix+hprefix;
  if(nPV <= 10) prefix_nPV.ReplaceAll(GetEraShort(), GetEraShort()+"F");      // few
  else if(nPV <= 20) prefix_nPV.ReplaceAll(GetEraShort(), GetEraShort()+"S"); // some
  else if(nPV <= 30) prefix_nPV.ReplaceAll(GetEraShort(), GetEraShort()+"L"); // little
  else if(nPV <= 40) prefix_nPV.ReplaceAll(GetEraShort(), GetEraShort()+"M"); // middle
  else if(nPV <= 50) prefix_nPV.ReplaceAll(GetEraShort(), GetEraShort()+"H"); // high
  else prefix_nPV.ReplaceAll(GetEraShort(), GetEraShort()+"V");               // very high

  Jet lepb = bjets.at(idx_bbjj.at(0));
  Jet hadb = bjets.at(idx_bbjj.at(1));
  Jet Wj0 = ajets.at(idx_bbjj.at(2));
  Jet Wj1 = ajets.at(idx_bbjj.at(3));
  double lepb_charge = idx_bbjj.at(0) == 0? b0_charge: b1_charge;
  double hadb_charge = idx_bbjj.at(1) == 0? b0_charge: b1_charge;
  double a0charge = jetCharge(Wj0);
  double a1charge = jetCharge(Wj1);

  if(!IsDATA && MCSample.Contains("TTLJ")){
    bool isTTLJinAcceptance = true;
    for(Gen gen:{gen_b0, gen_b1, gen_j0, gen_j1}){
      if(gen.Pt() < 30) isTTLJinAcceptance = false;
      if(fabs(gen.Eta()) > (DataYear == 2016? 2.4: 2.5)) isTTLJinAcceptance = false;
    }
    FillHist(prefix+hprefix+"LR_acceptance"+suffix, (isTTLJinAcceptance? 1: 0), map_weight, 2,0,2);

    double match_dR = 0.4;
    bool Gen_LR_match_onlyb = false;
    bool Gen_LR_match_Whad = false;
    bool Gen_LR_match_full = false;
    if((gen_j0.DeltaR(Wj0) < match_dR && gen_j1.DeltaR(Wj1) < match_dR) || (gen_j0.DeltaR(Wj1) < match_dR && gen_j1.DeltaR(Wj0) < match_dR)) Gen_LR_match_Whad = true;

    // ChargeEasy => 0:negative, 1:positive
    // purity => -1:UnMatched, 0:Wrong, 1:Correct in TTLJ
    // correct => 0:Wrong, 1:Correct in TTLJ

    //When lepb = b, hadb = B (bbar) -> lep+
    if(gen_b0.DeltaR(lepb) < match_dR && gen_b1.DeltaR(hadb) < match_dR){
      Gen_LR_match_onlyb = true;
      if(Gen_LR_match_Whad) Gen_LR_match_full = true;
      if(lepton0->Charge() > 0){
        FillHist(prefix+hprefix+"LR_purity"+suffix, 1, map_weight, 4,-2,2);
        FillHist(prefix+hprefix+"LR_correct"+suffix, 1, map_weight, 2,0,2);
        FillHist(prefix_nPV+"LR_purity"+suffix, 1, map_weight, 4,-2,2);
        FillHist(prefix_nPV+"LR_correct"+suffix, 1, map_weight, 2,0,2);
        if(isTTLJinAcceptance) FillHist(prefix+hprefix+"LR_purity_inAcceptance"+suffix, 1, map_weight, 2,0,2);
        if(isTTLJinAcceptance) FillHist(prefix+hprefix+"LR_correct_inAcceptance"+suffix, 1, map_weight, 2,0,2);
        if(isTTLJinAcceptance) FillHist(prefix_nPV+"LR_purity_inAcceptance"+suffix, 1, map_weight, 2,0,2);
        if(isTTLJinAcceptance) FillHist(prefix_nPV+"LR_correct_inAcceptance"+suffix, 1, map_weight, 2,0,2);
        FillHist(prefix+hprefix+"LR_purity_Lp"+suffix, 1, map_weight, 4,-2,2);
        FillHist(prefix+hprefix+"LR_correct_Lp"+suffix, 1, map_weight, 2,0,2);

        // For bCharge Accuracy (gen-info)
        FillHist(prefix+hprefix+"genbCharge_Lp"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        FillHist(prefix+hprefix+"genBCharge_Lp"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        FillHist(prefix_nPV+"genbCharge_Lp"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        FillHist(prefix_nPV+"genBCharge_Lp"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        if(lepb.Pt() < 35) FillHist(prefix+hprefix+"genbCharge_pt0_Lp"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        else if(lepb.Pt() < 50) FillHist(prefix+hprefix+"genbCharge_pt1_Lp"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        else if(lepb.Pt() < 80) FillHist(prefix+hprefix+"genbCharge_pt2_Lp"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        else if(lepb.Pt() < 120) FillHist(prefix+hprefix+"genbCharge_pt3_Lp"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        else FillHist(prefix+hprefix+"genbCharge_pt4_Lp"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        if(hadb.Pt() < 35) FillHist(prefix+hprefix+"genBCharge_pt0_Lp"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        else if(hadb.Pt() < 50) FillHist(prefix+hprefix+"genBCharge_pt1_Lp"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        else if(hadb.Pt() < 80) FillHist(prefix+hprefix+"genBCharge_pt2_Lp"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        else if(hadb.Pt() < 120) FillHist(prefix+hprefix+"genBCharge_pt3_Lp"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        else FillHist(prefix+hprefix+"genBCharge_pt4_Lp"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);

        for(unsigned int i=1; i<afb_chbinnum+1; i++){
          if(afb_chbin[i-1] < abs(lepb_charge) && abs(lepb_charge) < afb_chbin[i]){
            FillHist(Form(prefix+hprefix+"genbCharge%d_Lp"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            FillHist(Form(prefix_nPV+"genbCharge%d_Lp"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            if(lepb.Pt() < 35) FillHist(Form(prefix+hprefix+"genbCharge%d_pt0_Lp"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            else if(lepb.Pt() < 50) FillHist(Form(prefix+hprefix+"genbCharge%d_pt1_Lp"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            else if(lepb.Pt() < 80) FillHist(Form(prefix+hprefix+"genbCharge%d_pt2_Lp"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            else if(lepb.Pt() < 120) FillHist(Form(prefix+hprefix+"genbCharge%d_pt3_Lp"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            else FillHist(Form(prefix+hprefix+"genbCharge%d_pt4_Lp"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
          }
          if(afb_chbin[i-1] < abs(hadb_charge) && abs(hadb_charge) < afb_chbin[i]){
            FillHist(Form(prefix+hprefix+"genBCharge%d_Lp"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            FillHist(Form(prefix_nPV+"genBCharge%d_Lp"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            if(hadb.Pt() < 35) FillHist(Form(prefix+hprefix+"genBCharge%d_pt0_Lp"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            else if(hadb.Pt() < 50) FillHist(Form(prefix+hprefix+"genBCharge%d_pt1_Lp"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            else if(hadb.Pt() < 80) FillHist(Form(prefix+hprefix+"genBCharge%d_pt2_Lp"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            else if(hadb.Pt() < 120) FillHist(Form(prefix+hprefix+"genBCharge%d_pt3_Lp"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            else FillHist(Form(prefix+hprefix+"genBCharge%d_pt4_Lp"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
          }
        }

        prefix += "Correct"; // lep+
        prefix_nPV += "Correct";
      }else{
        FillHist(prefix+hprefix+"LR_purity"+suffix, 0, map_weight, 4,-2,2);
        FillHist(prefix+hprefix+"LR_correct"+suffix, 0, map_weight, 2,0,2);
        FillHist(prefix_nPV+"LR_purity"+suffix, 0, map_weight, 4,-2,2);
        FillHist(prefix_nPV+"LR_correct"+suffix, 0, map_weight, 2,0,2);
        if(isTTLJinAcceptance) FillHist(prefix+hprefix+"LR_purity_inAcceptance"+suffix, 0, map_weight, 2,0,2);
        if(isTTLJinAcceptance) FillHist(prefix+hprefix+"LR_correct_inAcceptance"+suffix, 0, map_weight, 2,0,2);
        if(isTTLJinAcceptance) FillHist(prefix_nPV+"LR_purity_inAcceptance"+suffix, 0, map_weight, 2,0,2);
        if(isTTLJinAcceptance) FillHist(prefix_nPV+"LR_correct_inAcceptance"+suffix, 0, map_weight, 2,0,2);
        FillHist(prefix+hprefix+"LR_purity_Lm"+suffix, 0, map_weight, 4,-2,2);
        FillHist(prefix+hprefix+"LR_correct_Lm"+suffix, 0, map_weight, 2,0,2);

        // For bCharge Accuracy (gen-info)
        FillHist(prefix+hprefix+"genbCharge_Lm"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        FillHist(prefix+hprefix+"genBCharge_Lm"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        FillHist(prefix_nPV+"genbCharge_Lm"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        FillHist(prefix_nPV+"genBCharge_Lm"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        if(lepb.Pt() < 35) FillHist(prefix+hprefix+"genbCharge_pt0_Lm"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        else if(lepb.Pt() < 50) FillHist(prefix+hprefix+"genbCharge_pt1_Lm"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        else if(lepb.Pt() < 80) FillHist(prefix+hprefix+"genbCharge_pt2_Lm"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        else if(lepb.Pt() < 120) FillHist(prefix+hprefix+"genbCharge_pt3_Lm"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        else FillHist(prefix+hprefix+"genbCharge_pt4_Lm"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        if(hadb.Pt() < 35) FillHist(prefix+hprefix+"genBCharge_pt0_Lm"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        else if(hadb.Pt() < 50) FillHist(prefix+hprefix+"genBCharge_pt1_Lm"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        else if(hadb.Pt() < 80) FillHist(prefix+hprefix+"genBCharge_pt2_Lm"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        else if(hadb.Pt() < 120) FillHist(prefix+hprefix+"genBCharge_pt3_Lm"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        else FillHist(prefix+hprefix+"genBCharge_pt4_Lm"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);

        for(unsigned int i=1; i<afb_chbinnum+1; i++){
          if(afb_chbin[i-1] < abs(lepb_charge) && abs(lepb_charge) < afb_chbin[i]){
            FillHist(Form(prefix+hprefix+"genbCharge%d_Lm"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            FillHist(Form(prefix_nPV+"genbCharge%d_Lm"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            if(lepb.Pt() < 35) FillHist(Form(prefix+hprefix+"genbCharge%d_pt0_Lm"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            else if(lepb.Pt() < 50) FillHist(Form(prefix+hprefix+"genbCharge%d_pt1_Lm"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            else if(lepb.Pt() < 80) FillHist(Form(prefix+hprefix+"genbCharge%d_pt2_Lm"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            else if(lepb.Pt() < 120) FillHist(Form(prefix+hprefix+"genbCharge%d_pt3_Lm"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            else FillHist(Form(prefix+hprefix+"genbCharge%d_pt4_Lm"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
          }
          if(afb_chbin[i-1] < abs(hadb_charge) && abs(hadb_charge) < afb_chbin[i]){
            FillHist(Form(prefix+hprefix+"genBCharge%d_Lm"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            FillHist(Form(prefix_nPV+"genBCharge%d_Lm"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            if(hadb.Pt() < 35) FillHist(Form(prefix+hprefix+"genBCharge%d_pt0_Lm"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            else if(hadb.Pt() < 50) FillHist(Form(prefix+hprefix+"genBCharge%d_pt1_Lm"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            else if(hadb.Pt() < 80) FillHist(Form(prefix+hprefix+"genBCharge%d_pt2_Lm"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            else if(hadb.Pt() < 120) FillHist(Form(prefix+hprefix+"genBCharge%d_pt3_Lm"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            else FillHist(Form(prefix+hprefix+"genBCharge%d_pt4_Lm"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
          }
        }

        prefix += "Wrong"; // lep-
        prefix_nPV += "Wrong";
      }
    }
    //When lepb = B (bbar), hadb = b => lep-
    else if(gen_b0.DeltaR(hadb) < match_dR && gen_b1.DeltaR(lepb) < match_dR){
      Gen_LR_match_onlyb = true;
      if(Gen_LR_match_Whad) Gen_LR_match_full = true;
      if(lepton0->Charge() < 0){
        FillHist(prefix+hprefix+"LR_purity"+suffix, 1, map_weight, 4,-2,2);
        FillHist(prefix+hprefix+"LR_correct"+suffix, 1, map_weight, 2,0,2);
        FillHist(prefix_nPV+"LR_purity"+suffix, 1, map_weight, 4,-2,2);
        FillHist(prefix_nPV+"LR_correct"+suffix, 1, map_weight, 2,0,2);
        if(isTTLJinAcceptance) FillHist(prefix+hprefix+"LR_purity_inAcceptance"+suffix, 1, map_weight, 2,0,2);
        if(isTTLJinAcceptance) FillHist(prefix+hprefix+"LR_correct_inAcceptance"+suffix, 1, map_weight, 2,0,2);
        if(isTTLJinAcceptance) FillHist(prefix_nPV+"LR_purity_inAcceptance"+suffix, 1, map_weight, 2,0,2);
        if(isTTLJinAcceptance) FillHist(prefix_nPV+"LR_correct_inAcceptance"+suffix, 1, map_weight, 2,0,2);
        FillHist(prefix+hprefix+"LR_purity_Lm"+suffix, 1, map_weight, 4,-2,2);
        FillHist(prefix+hprefix+"LR_correct_Lm"+suffix, 1, map_weight, 2,0,2);

        // For bCharge Accuracy (gen-info)
        FillHist(prefix+hprefix+"genbCharge_Lm"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        FillHist(prefix+hprefix+"genBCharge_Lm"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        FillHist(prefix_nPV+"genbCharge_Lm"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        FillHist(prefix_nPV+"genBCharge_Lm"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        if(hadb.Pt() < 35) FillHist(prefix+hprefix+"genbCharge_pt0_Lm"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        else if(hadb.Pt() < 50) FillHist(prefix+hprefix+"genbCharge_pt1_Lm"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        else if(hadb.Pt() < 80) FillHist(prefix+hprefix+"genbCharge_pt2_Lm"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        else if(hadb.Pt() < 120) FillHist(prefix+hprefix+"genbCharge_pt3_Lm"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        else FillHist(prefix+hprefix+"genbCharge_pt4_Lm"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        if(lepb.Pt() < 35) FillHist(prefix+hprefix+"genBCharge_pt0_Lm"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        else if(lepb.Pt() < 50) FillHist(prefix+hprefix+"genBCharge_pt1_Lm"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        else if(lepb.Pt() < 80) FillHist(prefix+hprefix+"genBCharge_pt2_Lm"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        else if(lepb.Pt() < 120) FillHist(prefix+hprefix+"genBCharge_pt3_Lm"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        else FillHist(prefix+hprefix+"genBCharge_pt4_Lm"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);

        for(unsigned int i=1; i<afb_chbinnum+1; i++){
          if(afb_chbin[i-1] < abs(hadb_charge) && abs(hadb_charge) < afb_chbin[i]){
            FillHist(Form(prefix+hprefix+"genbCharge%d_Lm"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            FillHist(Form(prefix_nPV+"genbCharge%d_Lm"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            if(hadb.Pt() < 35) FillHist(Form(prefix+hprefix+"genbCharge%d_pt0_Lm"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            else if(hadb.Pt() < 50) FillHist(Form(prefix+hprefix+"genbCharge%d_pt1_Lm"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            else if(hadb.Pt() < 80) FillHist(Form(prefix+hprefix+"genbCharge%d_pt2_Lm"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            else if(hadb.Pt() < 120) FillHist(Form(prefix+hprefix+"genbCharge%d_pt3_Lm"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            else FillHist(Form(prefix+hprefix+"genbCharge%d_pt4_Lm"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
          }
          if(afb_chbin[i-1] < abs(lepb_charge) && abs(lepb_charge) < afb_chbin[i]){
            FillHist(Form(prefix+hprefix+"genBCharge%d_Lm"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            FillHist(Form(prefix_nPV+"genBCharge%d_Lm"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            if(lepb.Pt() < 35) FillHist(Form(prefix+hprefix+"genBCharge%d_pt0_Lm"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            else if(lepb.Pt() < 50) FillHist(Form(prefix+hprefix+"genBCharge%d_pt1_Lm"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            else if(lepb.Pt() < 80) FillHist(Form(prefix+hprefix+"genBCharge%d_pt2_Lm"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            else if(lepb.Pt() < 120) FillHist(Form(prefix+hprefix+"genBCharge%d_pt3_Lm"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            else FillHist(Form(prefix+hprefix+"genBCharge%d_pt4_Lm"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
          }
        }

        prefix += "Correct"; // lep-
        prefix_nPV += "Correct";
      }else{
        FillHist(prefix+hprefix+"LR_purity"+suffix, 0, map_weight, 4,-2,2);
        FillHist(prefix+hprefix+"LR_correct"+suffix, 0, map_weight, 2,0,2);
        FillHist(prefix_nPV+"LR_purity"+suffix, 0, map_weight, 4,-2,2);
        FillHist(prefix_nPV+"LR_correct"+suffix, 0, map_weight, 2,0,2);
        if(isTTLJinAcceptance) FillHist(prefix+hprefix+"LR_purity_inAcceptance"+suffix, 0, map_weight, 2,0,2);
        if(isTTLJinAcceptance) FillHist(prefix+hprefix+"LR_correct_inAcceptance"+suffix, 0, map_weight, 2,0,2);
        if(isTTLJinAcceptance) FillHist(prefix_nPV+"LR_purity_inAcceptance"+suffix, 0, map_weight, 2,0,2);
        if(isTTLJinAcceptance) FillHist(prefix_nPV+"LR_correct_inAcceptance"+suffix, 0, map_weight, 2,0,2);
        FillHist(prefix+hprefix+"LR_purity_Lp"+suffix, 0, map_weight, 4,-2,2);
        FillHist(prefix+hprefix+"LR_correct_Lp"+suffix, 0, map_weight, 2,0,2);

        // For bCharge Accuracy (gen-info)
        FillHist(prefix+hprefix+"genbCharge_Lp"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        FillHist(prefix+hprefix+"genBCharge_Lp"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        FillHist(prefix_nPV+"genbCharge_Lp"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        FillHist(prefix_nPV+"genBCharge_Lp"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        if(hadb.Pt() < 35) FillHist(prefix+hprefix+"genbCharge_pt0_Lp"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        else if(hadb.Pt() < 50) FillHist(prefix+hprefix+"genbCharge_pt1_Lp"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        else if(hadb.Pt() < 80) FillHist(prefix+hprefix+"genbCharge_pt2_Lp"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        else if(hadb.Pt() < 120) FillHist(prefix+hprefix+"genbCharge_pt3_Lp"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        else FillHist(prefix+hprefix+"genbCharge_pt4_Lp"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        if(lepb.Pt() < 35) FillHist(prefix+hprefix+"genBCharge_pt0_Lp"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        else if(lepb.Pt() < 50) FillHist(prefix+hprefix+"genBCharge_pt1_Lp"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        else if(lepb.Pt() < 80) FillHist(prefix+hprefix+"genBCharge_pt2_Lp"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        else if(lepb.Pt() < 120) FillHist(prefix+hprefix+"genBCharge_pt3_Lp"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        else FillHist(prefix+hprefix+"genBCharge_pt4_Lp"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);

        for(unsigned int i=1; i<afb_chbinnum+1; i++){
          if(afb_chbin[i-1] < abs(hadb_charge) && abs(hadb_charge) < afb_chbin[i]){
            FillHist(Form(prefix+hprefix+"genbCharge%d_Lp"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            FillHist(Form(prefix_nPV+"genbCharge%d_Lp"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            if(hadb.Pt() < 35) FillHist(Form(prefix+hprefix+"genbCharge%d_pt0_Lp"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            else if(hadb.Pt() < 50) FillHist(Form(prefix+hprefix+"genbCharge%d_pt1_Lp"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            else if(hadb.Pt() < 80) FillHist(Form(prefix+hprefix+"genbCharge%d_pt2_Lp"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            else if(hadb.Pt() < 120) FillHist(Form(prefix+hprefix+"genbCharge%d_pt3_Lp"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            else FillHist(Form(prefix+hprefix+"genbCharge%d_pt4_Lp"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
          }
          if(afb_chbin[i-1] < abs(lepb_charge) && abs(lepb_charge) < afb_chbin[i]){
            FillHist(Form(prefix+hprefix+"genBCharge%d_Lp"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            FillHist(Form(prefix_nPV+"genBCharge%d_Lp"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            if(lepb.Pt() < 35) FillHist(Form(prefix+hprefix+"genBCharge%d_pt0_Lp"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            else if(lepb.Pt() < 50) FillHist(Form(prefix+hprefix+"genBCharge%d_pt1_Lp"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            else if(lepb.Pt() < 80) FillHist(Form(prefix+hprefix+"genBCharge%d_pt2_Lp"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            else if(lepb.Pt() < 120) FillHist(Form(prefix+hprefix+"genBCharge%d_pt3_Lp"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
            else FillHist(Form(prefix+hprefix+"genBCharge%d_pt4_Lp"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
          }
        }

        prefix += "Wrong"; // lep+
        prefix_nPV += "Wrong";
      }
    }
    else{
      FillHist(prefix+hprefix+"LR_purity"+suffix, -1, map_weight, 4,-2,2);
      if(isTTLJinAcceptance) FillHist(prefix+hprefix+"LR_purity_inAcceptance"+suffix, -1, map_weight, 2,0,2);
      if(isTTLJinAcceptance) FillHist(prefix_nPV+"LR_purity_inAcceptance"+suffix, -1, map_weight, 2,0,2);
      if(lepton0->Charge() > 0) FillHist(prefix+hprefix+"LR_purity_Lp"+suffix, -1, map_weight, 4,-2,2);
      if(lepton0->Charge() < 0) FillHist(prefix+hprefix+"LR_purity_Lm"+suffix, -1, map_weight, 4,-2,2);
      prefix += "UnMatched";
      prefix_nPV += "UnMatched";
    }
    // Too many histograms
    //prefix += (isTTLJinAcceptance? "0_": "1_");
    //prefix_nPV += (isTTLJinAcceptance? "0_": "1_");
    prefix += "_";
    prefix_nPV += "_";

    FillHist("Gen_LR_match_onlyb"+suffix, Gen_LR_match_onlyb, map_weight, 2,0,2);
    FillHist("Gen_LR_match_Whad"+suffix, Gen_LR_match_Whad, map_weight, 2,0,2);
    FillHist("Gen_LR_match_full"+suffix, Gen_LR_match_full, map_weight, 2,0,2);
    FillHist(prefix+hprefix+"Gen_LR_match_onlyb"+suffix, Gen_LR_match_onlyb, map_weight, 2,0,2);
    FillHist(prefix+hprefix+"Gen_LR_match_Whad"+suffix, Gen_LR_match_Whad, map_weight, 2,0,2);
    FillHist(prefix+hprefix+"Gen_LR_match_full"+suffix, Gen_LR_match_full, map_weight, 2,0,2);
    FillHist(channel+"/Gen_LR_match_onlyb"+suffix, Gen_LR_match_onlyb, map_weight, 2,0,2);
    FillHist(channel+"/Gen_LR_match_Whad"+suffix, Gen_LR_match_Whad, map_weight, 2,0,2);
    FillHist(channel+"/Gen_LR_match_full"+suffix, Gen_LR_match_full, map_weight, 2,0,2);
    if(isTTLJinAcceptance){
      FillHist("Gen_LR_match_onlyb_inAcceptance"+suffix, Gen_LR_match_onlyb, map_weight, 2,0,2);
      FillHist("Gen_LR_match_Whad_inAcceptance"+suffix, Gen_LR_match_Whad, map_weight, 2,0,2);
      FillHist("Gen_LR_match_full_inAcceptance"+suffix, Gen_LR_match_full, map_weight, 2,0,2);
      FillHist(prefix+hprefix+"Gen_LR_match_onlyb_inAcceptance"+suffix, Gen_LR_match_onlyb, map_weight, 2,0,2);
      FillHist(prefix+hprefix+"Gen_LR_match_Whad_inAcceptance"+suffix, Gen_LR_match_Whad, map_weight, 2,0,2);
      FillHist(prefix+hprefix+"Gen_LR_match_full_inAcceptance"+suffix, Gen_LR_match_full, map_weight, 2,0,2);
      FillHist(channel+"/Gen_LR_match_onlyb_inAcceptance"+suffix, Gen_LR_match_onlyb, map_weight, 2,0,2);
      FillHist(channel+"/Gen_LR_match_Whad_inAcceptance"+suffix, Gen_LR_match_Whad, map_weight, 2,0,2);
      FillHist(channel+"/Gen_LR_match_full_inAcceptance"+suffix, Gen_LR_match_full, map_weight, 2,0,2);
    }
  }

  //==========================
  //==== Now reco fill histograms
  //==========================
  FillHist(prefix+hprefix+"mass_jj"+suffix, (Wj0 + Wj1).M(), map_weight, 40,0,200);
  FillHist(prefix+hprefix+"mass_bjj"+suffix, (hadb + Wj0 + Wj1).M(), map_weight, 40,100,300);
  //FillHist(prefix+hprefix+"mass_Toplep"+suffix, (lepb + *lepton0 + neutrino).M(), map_weight, 40,100,300);
  FillHist(prefix+hprefix+"mass_blMET"+suffix, (lepb + *lepton0 + met).M(), map_weight, 40,100,300);
  FillHist(prefix+hprefix+"mass_bl"+suffix, (lepb + *lepton0).M(), map_weight, 40,0,200);

  FillHist(prefix+hprefix+"dR_bjj_blMET"+suffix, (hadb + Wj0 + Wj1).DeltaR(lepb + *lepton0 + met), map_weight, 200,0,10);
  FillHist(prefix+hprefix+"dR_bjj_lMET"+suffix, (hadb + Wj0 + Wj1).DeltaR(*lepton0 + met), map_weight, 200,0,10);
  FillHist(prefix+hprefix+"dR_bjj_jj"+suffix, (hadb + Wj0 + Wj1).DeltaR(Wj0 + Wj1), map_weight, 200,0,10);
  FillHist(prefix+hprefix+"dR_blMET_lMET"+suffix, (lepb + *lepton0 + met).DeltaR(*lepton0 + met), map_weight, 200,0,10);
  FillHist(prefix+hprefix+"dR_blMET_jj"+suffix, (lepb + *lepton0 + met).DeltaR(Wj0 + Wj1), map_weight, 200,0,10);
  FillHist(prefix+hprefix+"dR_lMET_jj"+suffix, (*lepton0 + met).DeltaR(Wj0 + Wj1), map_weight, 200,0,10);

  FillHist(prefix+hprefix+"dPhi_bjj_blMET"+suffix, fabs((hadb + Wj0 + Wj1).DeltaPhi(lepb + *lepton0 + met)), map_weight, 200,0,10);
  FillHist(prefix+hprefix+"dPhi_bjj_lMET"+suffix, fabs((hadb + Wj0 + Wj1).DeltaPhi(*lepton0 + met)), map_weight, 200,0,10);
  FillHist(prefix+hprefix+"dPhi_bjj_jj"+suffix, fabs((hadb + Wj0 + Wj1).DeltaPhi(Wj0 + Wj1)), map_weight, 200,0,10);
  FillHist(prefix+hprefix+"dPhi_blMET_lMET"+suffix, fabs((lepb + *lepton0 + met).DeltaPhi(*lepton0 + met)), map_weight, 200,0,10);
  FillHist(prefix+hprefix+"dPhi_blMET_jj"+suffix, fabs((lepb + *lepton0 + met).DeltaPhi(Wj0 + Wj1)), map_weight, 200,0,10);
  FillHist(prefix+hprefix+"dPhi_lMET_jj"+suffix, fabs((*lepton0 + met).DeltaPhi(Wj0 + Wj1)), map_weight, 200,0,10);

  FillHist(prefix+hprefix+"likelihood_ratio"+suffix, LRs.at(0) * LRs.at(1) * LRs.at(2) * LRs.at(3), map_weight, 1000,0,50);
  FillHist(prefix+hprefix+"likelihood_ratio_Mbl"+suffix, LRs.at(0), map_weight, 100,0,5);
  FillHist(prefix+hprefix+"likelihood_ratio_MblMET"+suffix, LRs.at(1), map_weight, 100,0,5);
  FillHist(prefix+hprefix+"likelihood_ratio_Mbjj"+suffix, LRs.at(2), map_weight, 100,0,5);
  FillHist(prefix+hprefix+"likelihood_ratio_Mjj"+suffix, LRs.at(3), map_weight, 100,0,5);
  FillHist(prefix+hprefix+"likelihood_ratio_dRtt"+suffix, LRs.at(4), map_weight, 100,0,5);
  FillHist(prefix+hprefix+"likelihood_ratio_dPhitt"+suffix, LRs.at(5), map_weight, 100,0,5);

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

  FillHist(prefix+hprefix+"leta"+suffix, lepton0->Eta(), map_weight, 100,-5,5);
  FillHist(prefix+hprefix+"beta"+suffix, lepb.Eta(), map_weight, 100,-5,5);
  FillHist(prefix+hprefix+"beta"+suffix, hadb.Eta(), map_weight, 100,-5,5);
  FillHist(prefix+hprefix+"lepbeta"+suffix, lepb.Eta(), map_weight, 100,-5,5);
  FillHist(prefix+hprefix+"hadbeta"+suffix, hadb.Eta(), map_weight, 100,-5,5);
  FillHist(prefix+hprefix+"bphi"+suffix, lepb.Phi(), map_weight, 64,-3.2,3.2);
  FillHist(prefix+hprefix+"bphi"+suffix, hadb.Phi(), map_weight, 64,-3.2,3.2);
  FillHist(prefix+hprefix+"lepbphi"+suffix, lepb.Phi(), map_weight, 64,-3.2,3.2);
  FillHist(prefix+hprefix+"hadbphi"+suffix, hadb.Phi(), map_weight, 64,-3.2,3.2);
  FillHist(prefix+hprefix+"jeta"+suffix, Wj0.Eta(), map_weight, 100,-5,5);
  FillHist(prefix+hprefix+"jeta"+suffix, Wj1.Eta(), map_weight, 100,-5,5);
  FillHist(prefix+hprefix+"j0eta"+suffix, Wj0.Eta(), map_weight, 100,-5,5);
  FillHist(prefix+hprefix+"j1eta"+suffix, Wj1.Eta(), map_weight, 100,-5,5);
  FillHist(prefix+hprefix+"jphi"+suffix, Wj0.Phi(), map_weight, 64,-3.2,3.2);
  FillHist(prefix+hprefix+"jphi"+suffix, Wj1.Phi(), map_weight, 64,-3.2,3.2);
  FillHist(prefix+hprefix+"j0phi"+suffix, Wj0.Phi(), map_weight, 64,-3.2,3.2);
  FillHist(prefix+hprefix+"j1phi"+suffix, Wj1.Phi(), map_weight, 64,-3.2,3.2);
  FillHist(prefix+hprefix+"met"+suffix, PuppiMET_Type1_pt, map_weight, 200,0,200);
  FillHist(prefix+hprefix+"metphi"+suffix, PuppiMET_Type1_phi, map_weight, 64,-3.2,3.2);
  FillHist(prefix+hprefix+"met_metphi"+suffix, PuppiMET_Type1_pt, PuppiMET_Type1_phi, map_weight[""], 40,0,200, 32,-3.2,3.2);

  // bCharges
  FillHist(prefix+hprefix+"bbCharges"+suffix, bbCharges, map_weight, 4,0,4);
  FillHist(prefix+hprefix+"lepbChargeRaw"+suffix, lepb_charge, map_weight, 200,-5,5);
  FillHist(prefix+hprefix+"lepbChargeRaw2"+suffix, lepb_charge, map_weight, 100,-5,5);
  FillHist(prefix+hprefix+"lepbCharge"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
  FillHist(prefix+hprefix+"hadbChargeRaw"+suffix, hadb_charge, map_weight, 200,-5,5);
  FillHist(prefix+hprefix+"hadbChargeRaw2"+suffix, hadb_charge, map_weight, 100,-5,5);
  FillHist(prefix+hprefix+"hadbCharge"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
  FillHist(prefix+hprefix+"a0ChargeRaw"+suffix, a0charge, map_weight, 200,-5,5);
  FillHist(prefix+hprefix+"a0Charge"+suffix, (a0charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
  FillHist(prefix+hprefix+"a1ChargeRaw"+suffix, a1charge, map_weight, 200,-5,5);
  FillHist(prefix+hprefix+"a1Charge"+suffix, (a1charge < 0? -0.5: 0.5), map_weight, 2,-1,1);

  if(lepton0->Charge() < 0){
    FillHist(prefix+hprefix+"lepbChargeRaw_Lm"+suffix, lepb_charge, map_weight, 200,-5,5);
    FillHist(prefix+hprefix+"lepbCharge_Lm"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
    FillHist(prefix+hprefix+"hadbChargeRaw_Lm"+suffix, hadb_charge, map_weight, 200,-5,5);
    FillHist(prefix+hprefix+"hadbCharge_Lm"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
    FillHist(prefix+hprefix+"recoBChargeRaw"+suffix, lepb_charge, map_weight, 200,-5,5);
    FillHist(prefix+hprefix+"recoBChargeRaw2"+suffix, lepb_charge, map_weight, 100,-5,5);
    FillHist(prefix+hprefix+"recoBCharge"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
    FillHist(prefix+hprefix+"recobChargeRaw"+suffix, hadb_charge, map_weight, 200,-5,5);
    FillHist(prefix+hprefix+"recobChargeRaw2"+suffix, hadb_charge, map_weight, 100,-5,5);
    FillHist(prefix+hprefix+"recobCharge"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
    if(lepb.Pt() < 35){
      FillHist(prefix+hprefix+"recoBChargeRaw_pt0"+suffix, lepb_charge, map_weight, 200,-5,5);
      FillHist(prefix+hprefix+"recoBCharge_pt0"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
    }else if(lepb.Pt() < 50){
      FillHist(prefix+hprefix+"recoBChargeRaw_pt1"+suffix, lepb_charge, map_weight, 200,-5,5);
      FillHist(prefix+hprefix+"recoBCharge_pt1"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
    }else if(lepb.Pt() < 80){
      FillHist(prefix+hprefix+"recoBChargeRaw_pt2"+suffix, lepb_charge, map_weight, 200,-5,5);
      FillHist(prefix+hprefix+"recoBCharge_pt2"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
    }else if(lepb.Pt() < 120){
      FillHist(prefix+hprefix+"recoBChargeRaw_pt3"+suffix, lepb_charge, map_weight, 200,-5,5);
      FillHist(prefix+hprefix+"recoBCharge_pt3"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
    }else{
      FillHist(prefix+hprefix+"recoBChargeRaw_pt4"+suffix, lepb_charge, map_weight, 200,-5,5);
      FillHist(prefix+hprefix+"recoBCharge_pt4"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
    }
    if(hadb.Pt() < 35){
      FillHist(prefix+hprefix+"recobChargeRaw_pt0"+suffix, hadb_charge, map_weight, 200,-5,5);
      FillHist(prefix+hprefix+"recobCharge_pt0"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
    }else if(hadb.Pt() < 50){
      FillHist(prefix+hprefix+"recobChargeRaw_pt1"+suffix, hadb_charge, map_weight, 200,-5,5);
      FillHist(prefix+hprefix+"recobCharge_pt1"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
    }else if(hadb.Pt() < 80){
      FillHist(prefix+hprefix+"recobChargeRaw_pt2"+suffix, hadb_charge, map_weight, 200,-5,5);
      FillHist(prefix+hprefix+"recobCharge_pt2"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
    }else if(hadb.Pt() < 120){
      FillHist(prefix+hprefix+"recobChargeRaw_pt3"+suffix, hadb_charge, map_weight, 200,-5,5);
      FillHist(prefix+hprefix+"recobCharge_pt3"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
    }else{
      FillHist(prefix+hprefix+"recobChargeRaw_pt4"+suffix, hadb_charge, map_weight, 200,-5,5);
      FillHist(prefix+hprefix+"recobCharge_pt4"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
    }
  }else{
    FillHist(prefix+hprefix+"lepbChargeRaw_Lp"+suffix, lepb_charge, map_weight, 200,-5,5);
    FillHist(prefix+hprefix+"lepbCharge_Lp"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
    FillHist(prefix+hprefix+"hadbChargeRaw_Lp"+suffix, hadb_charge, map_weight, 200,-5,5);
    FillHist(prefix+hprefix+"hadbCharge_Lp"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
    FillHist(prefix+hprefix+"recobChargeRaw"+suffix, lepb_charge, map_weight, 200,-5,5);
    FillHist(prefix+hprefix+"recobChargeRaw2"+suffix, lepb_charge, map_weight, 100,-5,5);
    FillHist(prefix+hprefix+"recobCharge"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
    FillHist(prefix+hprefix+"recoBChargeRaw"+suffix, hadb_charge, map_weight, 200,-5,5);
    FillHist(prefix+hprefix+"recoBChargeRaw2"+suffix, hadb_charge, map_weight, 100,-5,5);
    FillHist(prefix+hprefix+"recoBCharge"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
    if(lepb.Pt() < 35){
      FillHist(prefix+hprefix+"recobChargeRaw_pt0"+suffix, lepb_charge, map_weight, 200,-5,5);
      FillHist(prefix+hprefix+"recobCharge_pt0"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
    }else if(lepb.Pt() < 50){
      FillHist(prefix+hprefix+"recobChargeRaw_pt1"+suffix, lepb_charge, map_weight, 200,-5,5);
      FillHist(prefix+hprefix+"recobCharge_pt1"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
    }else if(lepb.Pt() < 80){
      FillHist(prefix+hprefix+"recobChargeRaw_pt2"+suffix, lepb_charge, map_weight, 200,-5,5);
      FillHist(prefix+hprefix+"recobCharge_pt2"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
    }else if(lepb.Pt() < 120){
      FillHist(prefix+hprefix+"recobChargeRaw_pt3"+suffix, lepb_charge, map_weight, 200,-5,5);
      FillHist(prefix+hprefix+"recobCharge_pt3"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
    }else{
      FillHist(prefix+hprefix+"recobChargeRaw_pt4"+suffix, lepb_charge, map_weight, 200,-5,5);
      FillHist(prefix+hprefix+"recobCharge_pt4"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
    }
    if(hadb.Pt() < 35){
      FillHist(prefix+hprefix+"recoBChargeRaw_pt0"+suffix, hadb_charge, map_weight, 200,-5,5);
      FillHist(prefix+hprefix+"recoBCharge_pt0"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
    }else if(hadb.Pt() < 50){
      FillHist(prefix+hprefix+"recoBChargeRaw_pt1"+suffix, hadb_charge, map_weight, 200,-5,5);
      FillHist(prefix+hprefix+"recoBCharge_pt1"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
    }else if(hadb.Pt() < 80){
      FillHist(prefix+hprefix+"recoBChargeRaw_pt2"+suffix, hadb_charge, map_weight, 200,-5,5);
      FillHist(prefix+hprefix+"recoBCharge_pt2"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
    }else if(hadb.Pt() < 120){
      FillHist(prefix+hprefix+"recoBChargeRaw_pt3"+suffix, hadb_charge, map_weight, 200,-5,5);
      FillHist(prefix+hprefix+"recoBCharge_pt3"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
    }else{
      FillHist(prefix+hprefix+"recoBChargeRaw_pt4"+suffix, hadb_charge, map_weight, 200,-5,5);
      FillHist(prefix+hprefix+"recoBCharge_pt4"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
    }
  }

  for(int i=1; i<afb_chbinnum+1; i++){
    if(lepton0->Charge() < 0){
      if(afb_chbin[i-1] < abs(lepb_charge) && abs(lepb_charge) < afb_chbin[i]){
        FillHist(Form(prefix+hprefix+"lepbCharge%dRaw_Lm"+suffix, i-1), lepb_charge, map_weight, 200,-5,5);
        FillHist(Form(prefix+hprefix+"lepbCharge%d_Lm"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        FillHist(prefix+hprefix+"lepbChargebin_Lm"+suffix, LHAPDF::sgn(lepb_charge) * (i - 0.5), map_weight, 12,-6,6);
        FillHist(Form(prefix+hprefix+"recoBCharge%dRaw"+suffix, i-1), lepb_charge, map_weight, 200,-5,5);
        FillHist(Form(prefix+hprefix+"recoBCharge%dRaw2"+suffix, i-1), lepb_charge, map_weight, 100,-5,5);
        FillHist(Form(prefix+hprefix+"recoBCharge%d"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        FillHist(prefix+hprefix+"recoBCharges"+suffix, (lepb_charge < 0? (2 * i) - 1.5: (2 * i) - 0.5), map_weight, 12,0,12);
        FillHist(prefix+hprefix+"recoBCharges2"+suffix, (lepb_charge < 0? -i + 0.5: i - 0.5), map_weight, 12,-6,6);
        if(lepb.Pt() < 35){
          FillHist(Form(prefix+hprefix+"recoBCharge%dRaw_pt0"+suffix, i-1), lepb_charge, map_weight, 200,-5,5);
          FillHist(Form(prefix+hprefix+"recoBCharge%d_pt0"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        }else if(lepb.Pt() < 50){
          FillHist(Form(prefix+hprefix+"recoBCharge%dRaw_pt1"+suffix, i-1), lepb_charge, map_weight, 200,-5,5);
          FillHist(Form(prefix+hprefix+"recoBCharge%d_pt1"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        }else if(lepb.Pt() < 80){
          FillHist(Form(prefix+hprefix+"recoBCharge%dRaw_pt2"+suffix, i-1), lepb_charge, map_weight, 200,-5,5);
          FillHist(Form(prefix+hprefix+"recoBCharge%d_pt2"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        }else if(lepb.Pt() < 120){
          FillHist(Form(prefix+hprefix+"recoBCharge%dRaw_pt3"+suffix, i-1), lepb_charge, map_weight, 200,-5,5);
          FillHist(Form(prefix+hprefix+"recoBCharge%d_pt3"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        }else{
          FillHist(Form(prefix+hprefix+"recoBCharge%dRaw_pt4"+suffix, i-1), lepb_charge, map_weight, 200,-5,5);
          FillHist(Form(prefix+hprefix+"recoBCharge%d_pt4"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        }
        FillHist(prefix+hprefix+"recoBChargebin"+suffix, LHAPDF::sgn(lepb_charge) * (i - 0.5), map_weight, 12,-6,6);
        FillHist(prefix+hprefix+"recoBChargeAbsbin"+suffix, i - 0.5, map_weight, 6,0,6);
        FillHist(prefix+hprefix+"recoChargebin"+suffix, LHAPDF::sgn(lepb_charge) * (i - 0.5), map_weight, 12,-6,6);
        FillHist(prefix+hprefix+"recoChargeAbsbin"+suffix, i - 0.5, map_weight, 6,0,6);
      }
      if(afb_chbin[i-1] < abs(hadb_charge) && abs(hadb_charge) < afb_chbin[i]){
        FillHist(Form(prefix+hprefix+"hadbCharge%dRaw_Lm"+suffix, i-1), hadb_charge, map_weight, 200,-5,5);
        FillHist(Form(prefix+hprefix+"hadbCharge%d_Lm"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        FillHist(prefix+hprefix+"hadbChargebin_Lm"+suffix, LHAPDF::sgn(hadb_charge) * (i - 0.5), map_weight, 12,-6,6);
        FillHist(Form(prefix+hprefix+"recobCharge%dRaw"+suffix, i-1), hadb_charge, map_weight, 200,-5,5);
        FillHist(Form(prefix+hprefix+"recobCharge%dRaw2"+suffix, i-1), hadb_charge, map_weight, 100,-5,5);
        FillHist(Form(prefix+hprefix+"recobCharge%d"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        FillHist(prefix+hprefix+"recobCharges"+suffix, (hadb_charge < 0? (2 * i) - 1.5: (2 * i) - 0.5), map_weight, 12,0,12);
        FillHist(prefix+hprefix+"recobCharges2"+suffix, (hadb_charge < 0? -i + 0.5: i - 0.5), map_weight, 12,-6,6);
        if(hadb.Pt() < 35){
          FillHist(Form(prefix+hprefix+"recobCharge%dRaw_pt0"+suffix, i-1), hadb_charge, map_weight, 200,-5,5);
          FillHist(Form(prefix+hprefix+"recobCharge%d_pt0"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        }else if(hadb.Pt() < 50){
          FillHist(Form(prefix+hprefix+"recobCharge%dRaw_pt1"+suffix, i-1), hadb_charge, map_weight, 200,-5,5);
          FillHist(Form(prefix+hprefix+"recobCharge%d_pt1"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        }else if(hadb.Pt() < 80){
          FillHist(Form(prefix+hprefix+"recobCharge%dRaw_pt2"+suffix, i-1), hadb_charge, map_weight, 200,-5,5);
          FillHist(Form(prefix+hprefix+"recobCharge%d_pt2"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        }else if(hadb.Pt() < 120){
          FillHist(Form(prefix+hprefix+"recobCharge%dRaw_pt3"+suffix, i-1), hadb_charge, map_weight, 200,-5,5);
          FillHist(Form(prefix+hprefix+"recobCharge%d_pt3"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        }else{
          FillHist(Form(prefix+hprefix+"recobCharge%dRaw_pt4"+suffix, i-1), hadb_charge, map_weight, 200,-5,5);
          FillHist(Form(prefix+hprefix+"recobCharge%d_pt4"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        }
        FillHist(prefix+hprefix+"recobChargebin"+suffix, LHAPDF::sgn(hadb_charge) * (i - 0.5), map_weight, 12,-6,6);
        FillHist(prefix+hprefix+"recobChargeAbsbin"+suffix, i - 0.5, map_weight, 6,0,6);
        FillHist(prefix+hprefix+"recoChargebin"+suffix, LHAPDF::sgn(hadb_charge) * (i - 0.5), map_weight, 12,-6,6);
        FillHist(prefix+hprefix+"recoChargeAbsbin"+suffix, i - 0.5, map_weight, 6,0,6);
      }
    }else{
      if(afb_chbin[i-1] < abs(lepb_charge) && abs(lepb_charge) < afb_chbin[i]){
        FillHist(Form(prefix+hprefix+"lepbCharge%dRaw_Lp"+suffix, i-1), lepb_charge, map_weight, 200,-5,5);
        FillHist(Form(prefix+hprefix+"lepbCharge%d_Lp"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        FillHist(prefix+hprefix+"lepbChargebin_Lp"+suffix, LHAPDF::sgn(lepb_charge) * (i - 0.5), map_weight, 12,-6,6);
        FillHist(Form(prefix+hprefix+"recobCharge%dRaw"+suffix, i-1), lepb_charge, map_weight, 200,-5,5);
        FillHist(Form(prefix+hprefix+"recobCharge%dRaw2"+suffix, i-1), lepb_charge, map_weight, 100,-5,5);
        FillHist(Form(prefix+hprefix+"recobCharge%d"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        FillHist(prefix+hprefix+"recobCharges"+suffix, (lepb_charge < 0? (2 * i) - 1.5: (2 * i) - 0.5), map_weight, 12,0,12);
        FillHist(prefix+hprefix+"recobCharges2"+suffix, (lepb_charge < 0? -i + 0.5: i - 0.5), map_weight, 12,-6,6);
        if(lepb.Pt() < 35){
          FillHist(Form(prefix+hprefix+"recobCharge%dRaw_pt0"+suffix, i-1), lepb_charge, map_weight, 200,-5,5);
          FillHist(Form(prefix+hprefix+"recobCharge%d_pt0"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        }else if(lepb.Pt() < 50){
          FillHist(Form(prefix+hprefix+"recobCharge%dRaw_pt1"+suffix, i-1), lepb_charge, map_weight, 200,-5,5);
          FillHist(Form(prefix+hprefix+"recobCharge%d_pt1"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        }else if(lepb.Pt() < 80){
          FillHist(Form(prefix+hprefix+"recobCharge%dRaw_pt2"+suffix, i-1), lepb_charge, map_weight, 200,-5,5);
          FillHist(Form(prefix+hprefix+"recobCharge%d_pt2"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        }else if(lepb.Pt() < 120){
          FillHist(Form(prefix+hprefix+"recobCharge%dRaw_pt3"+suffix, i-1), lepb_charge, map_weight, 200,-5,5);
          FillHist(Form(prefix+hprefix+"recobCharge%d_pt3"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        }else{
          FillHist(Form(prefix+hprefix+"recobCharge%dRaw_pt4"+suffix, i-1), lepb_charge, map_weight, 200,-5,5);
          FillHist(Form(prefix+hprefix+"recobCharge%d_pt4"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        }
        FillHist(prefix+hprefix+"recobChargebin"+suffix, LHAPDF::sgn(lepb_charge) * (i - 0.5), map_weight, 12,-6,6);
        FillHist(prefix+hprefix+"recobChargeAbsbin"+suffix, i - 0.5, map_weight, 6,0,6);
        FillHist(prefix+hprefix+"recoChargebin"+suffix, LHAPDF::sgn(lepb_charge) * (i - 0.5), map_weight, 12,-6,6);
        FillHist(prefix+hprefix+"recoChargeAbsbin"+suffix, i - 0.5, map_weight, 6,0,6);
      }
      if(afb_chbin[i-1] < abs(hadb_charge) && abs(hadb_charge) < afb_chbin[i]){
        FillHist(Form(prefix+hprefix+"hadbCharge%dRaw_Lp"+suffix, i-1), hadb_charge, map_weight, 200,-5,5);
        FillHist(Form(prefix+hprefix+"hadbCharge%d_Lp"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        FillHist(prefix+hprefix+"hadbChargebin_Lp"+suffix, LHAPDF::sgn(hadb_charge) * (i - 0.5), map_weight, 12,-6,6);
        FillHist(Form(prefix+hprefix+"recoBCharge%dRaw"+suffix, i-1), hadb_charge, map_weight, 200,-5,5);
        FillHist(Form(prefix+hprefix+"recoBCharge%dRaw2"+suffix, i-1), hadb_charge, map_weight, 100,-5,5);
        FillHist(Form(prefix+hprefix+"recoBCharge%d"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        FillHist(prefix+hprefix+"recoBCharges"+suffix, (hadb_charge < 0? (2 * i) - 1.5: (2 * i) - 0.5), map_weight, 12,0,12);
        FillHist(prefix+hprefix+"recoBCharges2"+suffix, (hadb_charge < 0? -i + 0.5: i - 0.5), map_weight, 12,-6,6);
        if(hadb.Pt() < 35){
          FillHist(Form(prefix+hprefix+"recoBCharge%dRaw_pt0"+suffix, i-1), hadb_charge, map_weight, 200,-5,5);
          FillHist(Form(prefix+hprefix+"recoBCharge%d_pt0"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        }else if(hadb.Pt() < 50){
          FillHist(Form(prefix+hprefix+"recoBCharge%dRaw_pt1"+suffix, i-1), hadb_charge, map_weight, 200,-5,5);
          FillHist(Form(prefix+hprefix+"recoBCharge%d_pt1"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        }else if(hadb.Pt() < 80){
          FillHist(Form(prefix+hprefix+"recoBCharge%dRaw_pt2"+suffix, i-1), hadb_charge, map_weight, 200,-5,5);
          FillHist(Form(prefix+hprefix+"recoBCharge%d_pt2"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        }else if(hadb.Pt() < 120){
          FillHist(Form(prefix+hprefix+"recoBCharge%dRaw_pt3"+suffix, i-1), hadb_charge, map_weight, 200,-5,5);
          FillHist(Form(prefix+hprefix+"recoBCharge%d_pt3"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        }else{
          FillHist(Form(prefix+hprefix+"recoBCharge%dRaw_pt4"+suffix, i-1), hadb_charge, map_weight, 200,-5,5);
          FillHist(Form(prefix+hprefix+"recoBCharge%d_pt4"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        }
        FillHist(prefix+hprefix+"recoBChargebin"+suffix, LHAPDF::sgn(hadb_charge) * (i - 0.5), map_weight, 12,-6,6);
        FillHist(prefix+hprefix+"recoBChargeAbsbin"+suffix, i - 0.5, map_weight, 6,0,6);
        FillHist(prefix+hprefix+"recoChargebin"+suffix, LHAPDF::sgn(hadb_charge) * (i - 0.5), map_weight, 12,-6,6);
        FillHist(prefix+hprefix+"recoChargeAbsbin"+suffix, i - 0.5, map_weight, 6,0,6);
      }
    }
  }

  // bCharges vs. nPV
  FillHist(prefix+hprefix+"nPV"+suffix, nPV, map_weight, 100,0,100);
  FillHist(prefix_nPV+"bbCharges"+suffix, bbCharges, map_weight, 4,0,4);
  if(lepton0->Charge() < 0){
    FillHist(prefix_nPV+"lepbChargeRaw_Lm"+suffix, lepb_charge, map_weight, 200,-5,5);
    FillHist(prefix_nPV+"lepbCharge_Lm"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
    FillHist(prefix_nPV+"hadbChargeRaw_Lm"+suffix, hadb_charge, map_weight, 200,-5,5);
    FillHist(prefix_nPV+"hadbCharge_Lm"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
    FillHist(prefix_nPV+"recoBChargeRaw"+suffix, lepb_charge, map_weight, 200,-5,5);
    FillHist(prefix_nPV+"recoBChargeRaw2"+suffix, lepb_charge, map_weight, 100,-5,5);
    FillHist(prefix_nPV+"recoBCharge"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
    FillHist(prefix_nPV+"recobChargeRaw"+suffix, hadb_charge, map_weight, 200,-5,5);
    FillHist(prefix_nPV+"recobChargeRaw2"+suffix, hadb_charge, map_weight, 100,-5,5);
    FillHist(prefix_nPV+"recobCharge"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
  }else{
    FillHist(prefix_nPV+"lepbChargeRaw_Lp"+suffix, lepb_charge, map_weight, 200,-5,5);
    FillHist(prefix_nPV+"lepbCharge_Lp"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
    FillHist(prefix_nPV+"hadbChargeRaw_Lp"+suffix, hadb_charge, map_weight, 200,-5,5);
    FillHist(prefix_nPV+"hadbCharge_Lp"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
    FillHist(prefix_nPV+"recobChargeRaw"+suffix, lepb_charge, map_weight, 200,-5,5);
    FillHist(prefix_nPV+"recobChargeRaw2"+suffix, lepb_charge, map_weight, 100,-5,5);
    FillHist(prefix_nPV+"recobCharge"+suffix, (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
    FillHist(prefix_nPV+"recoBChargeRaw"+suffix, hadb_charge, map_weight, 200,-5,5);
    FillHist(prefix_nPV+"recoBChargeRaw2"+suffix, hadb_charge, map_weight, 100,-5,5);
    FillHist(prefix_nPV+"recoBCharge"+suffix, (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
  }
  for(int i=1; i<afb_chbinnum+1; i++){
    if(lepton0->Charge() < 0){
      if(afb_chbin[i-1] < abs(lepb_charge) && abs(lepb_charge) < afb_chbin[i]){
        FillHist(Form(prefix_nPV+"lepbCharge%dRaw_Lm"+suffix, i-1), lepb_charge, map_weight, 200,-5,5);
        FillHist(Form(prefix_nPV+"lepbCharge%d_Lm"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        FillHist(prefix_nPV+"lepbChargebin_Lm"+suffix, LHAPDF::sgn(lepb_charge) * (i - 0.5), map_weight, 12,-6,6);
        FillHist(Form(prefix_nPV+"recoBCharge%dRaw"+suffix, i-1), lepb_charge, map_weight, 200,-5,5);
        FillHist(Form(prefix_nPV+"recoBCharge%d"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        FillHist(prefix_nPV+"recoBChargebin"+suffix, LHAPDF::sgn(lepb_charge) * (i - 0.5), map_weight, 12,-6,6);
        FillHist(prefix_nPV+"recoBChargeAbsbin"+suffix, i - 0.5, map_weight, 6,0,6);
        FillHist(prefix_nPV+"recoChargebin"+suffix, LHAPDF::sgn(lepb_charge) * (i - 0.5), map_weight, 12,-6,6);
        FillHist(prefix_nPV+"recoChargeAbsbin"+suffix, i - 0.5, map_weight, 6,0,6);
      }
      if(afb_chbin[i-1] < abs(hadb_charge) && abs(hadb_charge) < afb_chbin[i]){
        FillHist(Form(prefix_nPV+"hadbCharge%dRaw_Lm"+suffix, i-1), hadb_charge, map_weight, 200,-5,5);
        FillHist(Form(prefix_nPV+"hadbCharge%d_Lm"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        FillHist(prefix_nPV+"hadbChargebin_Lm"+suffix, LHAPDF::sgn(hadb_charge) * (i - 0.5), map_weight, 12,-6,6);
        FillHist(Form(prefix_nPV+"recobCharge%dRaw"+suffix, i-1), hadb_charge, map_weight, 200,-5,5);
        FillHist(Form(prefix_nPV+"recobCharge%d"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        FillHist(prefix_nPV+"recobChargebin"+suffix, LHAPDF::sgn(hadb_charge) * (i - 0.5), map_weight, 12,-6,6);
        FillHist(prefix_nPV+"recobChargeAbsbin"+suffix, i - 0.5, map_weight, 6,0,6);
        FillHist(prefix_nPV+"recoChargebin"+suffix, LHAPDF::sgn(hadb_charge) * (i - 0.5), map_weight, 12,-6,6);
        FillHist(prefix_nPV+"recoChargeAbsbin"+suffix, i - 0.5, map_weight, 6,0,6);
      }
    }else{
      if(afb_chbin[i-1] < abs(lepb_charge) && abs(lepb_charge) < afb_chbin[i]){
        FillHist(Form(prefix_nPV+"lepbCharge%dRaw_Lp"+suffix, i-1), lepb_charge, map_weight, 200,-5,5);
        FillHist(Form(prefix_nPV+"lepbCharge%d_Lp"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        FillHist(prefix_nPV+"lepbChargebin_Lp"+suffix, LHAPDF::sgn(lepb_charge) * (i - 0.5), map_weight, 12,-6,6);
        FillHist(Form(prefix_nPV+"recobCharge%dRaw"+suffix, i-1), lepb_charge, map_weight, 200,-5,5);
        FillHist(Form(prefix_nPV+"recobCharge%d"+suffix, i-1), (lepb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        FillHist(prefix_nPV+"recobChargebin"+suffix, LHAPDF::sgn(lepb_charge) * (i - 0.5), map_weight, 12,-6,6);
        FillHist(prefix_nPV+"recobChargeAbsbin"+suffix, i - 0.5, map_weight, 6,0,6);
        FillHist(prefix_nPV+"recoChargebin"+suffix, LHAPDF::sgn(lepb_charge) * (i - 0.5), map_weight, 12,-6,6);
        FillHist(prefix_nPV+"recoChargeAbsbin"+suffix, i - 0.5, map_weight, 6,0,6);
      }
      if(afb_chbin[i-1] < abs(hadb_charge) && abs(hadb_charge) < afb_chbin[i]){
        FillHist(Form(prefix_nPV+"hadbCharge%dRaw_Lp"+suffix, i-1), hadb_charge, map_weight, 200,-5,5);
        FillHist(Form(prefix_nPV+"hadbCharge%d_Lp"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        FillHist(prefix_nPV+"hadbChargebin_Lp"+suffix, LHAPDF::sgn(hadb_charge) * (i - 0.5), map_weight, 12,-6,6);
        FillHist(Form(prefix_nPV+"recoBCharge%dRaw"+suffix, i-1), hadb_charge, map_weight, 200,-5,5);
        FillHist(Form(prefix_nPV+"recoBCharge%d"+suffix, i-1), (hadb_charge < 0? -0.5: 0.5), map_weight, 2,-1,1);
        FillHist(prefix_nPV+"recoBChargebin"+suffix, LHAPDF::sgn(hadb_charge) * (i - 0.5), map_weight, 12,-6,6);
        FillHist(prefix_nPV+"recoBChargeAbsbin"+suffix, i - 0.5, map_weight, 6,0,6);
        FillHist(prefix_nPV+"recoChargebin"+suffix, LHAPDF::sgn(hadb_charge) * (i - 0.5), map_weight, 12,-6,6);
        FillHist(prefix_nPV+"recoChargeAbsbin"+suffix, i - 0.5, map_weight, 6,0,6);
      }
    }
  }

  FillHist(prefix+hprefix+"asChargeSum"+suffix, a0charge + a1charge, map_weight, 400,-10,10);
  FillHist(prefix+hprefix+"asChargeAbsSum"+suffix, (a0charge < 0? -1: 1) + (a1charge < 0? -1: 1), map_weight, 8,-4,4);
  FillHist(prefix+hprefix+"bsChargeSum"+suffix, lepb_charge + hadb_charge, map_weight, 400,-10,10);
  FillHist(prefix+hprefix+"bsChargeAbsSum"+suffix, (lepb_charge < 0? -1: 1) + (hadb_charge < 0? -1: 1), map_weight, 8,-4,4);
}

bool ttljAnalyzer::HasLeptons(TString channel, bool leps, unsigned int s, unsigned int m){
  bool moreleptons = false;
  double l0pt = 26.;
  if(channel.Contains("m"+GetEraShort())){
    muons = MuonMomentumCorrection(muons, s,m, false);
    if(muons.size() > 0) lepton0 = &muons.at(0);
    if(muons.size() > 1) moreleptons = true;
  }else if(channel.Contains("e"+GetEraShort())){
    l0pt = 30.;
    electrons = ElectronEnergyCorrection(electrons, s,m, false);
    if(electrons.size() > 0) lepton0 = &electrons.at(0);
    if(electrons.size() > 1) moreleptons = true;
  }else{
    cout<<"[ttljAnalyzer::Hasleptons] channel="<<channel<<" is weird"<<endl;
  }

  if(!lepton0) return false;
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "1Leptons", map_weight[""]);
  if(lepton0->Pt() < l0pt) return false;
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "LepPt", map_weight[""]);
  if(!leps){
    if(moreleptons) return false;
    if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "1Lepton", map_weight[""]);

    met = GetEvent().GetMETVector();
    if(met.Pt() < 20.) return false;
    if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "MET20", map_weight[""]);
  }

  leptons ={};
  leptons.push_back(lepton0);

  return true;
}

void ttljAnalyzer::executeEventGen(){
  FillHist("gen/executeEventGen", 1, 1, 2,0,2);
  gens=GetGens();
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
