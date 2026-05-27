#include "dybAnalyzer.h"

void dybAnalyzer::initializeAnalyzer(){
  SMPAnalyzerCore::initializeAnalyzer(); //setup eff zpt roc z0 cf
  IsNominalRun = !HasFlag("SYS") && !HasFlag("PDFSYS");

  PDFbase = LHAPDF::mkPDF(306000);
  PDFnf4 = LHAPDF::mkPDF(325500);
  mcCorr->SetJetTaggingParameters({
    JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Tight, JetTagging::incl, JetTagging::comb),
    JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Loose, JetTagging::incl, JetTagging::comb),
  });
  SetupPUJetWeight();
  SetupJetVetoMap();
}

void dybAnalyzer::executeEvent(){
  //// FIXME some events of DYJets has nan PDF weights. I don't know why...
  if(MCSample == "DYJets" && !isnormal(weight_Scale->at(0))) return;

  ///////////////// GEN level /////////////////////
  zptweight = 1., topptweight = 1., weakweight = 1.;
  if(IsDYSample || IsTTSample) executeEventGen();

  ///////////////// RECO level /////////////////////
  muons_raw = GetAllMuons(); // This can make event loops much slower in case of running over Unskimmed samples
  electrons_raw = GetAllElectrons();
  jets_raw = GetAllJets();
  if(!IsDATA || DataStream.Contains("DoubleMuon")){
    vector<Muon> muons_2p4 = MuonMomentumCorrection(SMPGetMuons("POGMediumWithLooseTrkIso", 8.0, 2.4), 0,0, true);
    muons = muons_2p4;
    executeEventWithParameter("mm"+GetEraShort());
    /*
    executeEventWithParameter("mm"+GetEraShort(), "_jet_dR_0p01");
    executeEventWithParameter("mm"+GetEraShort(), "_jet_dR_0p05");
    executeEventWithParameter("mm"+GetEraShort(), "_jet_dR_0p10");
    executeEventWithParameter("mm"+GetEraShort(), "_jet_dR_0p15");
    executeEventWithParameter("mm"+GetEraShort(), "_jet_dR_0p20");
    executeEventWithParameter("mm"+GetEraShort(), "_jet_dR_0p25");
    executeEventWithParameter("mm"+GetEraShort(), "_jet_dR_0p30");
    executeEventWithParameter("mm"+GetEraShort(), "_jet_dR_0p35");
    executeEventWithParameter("mm"+GetEraShort(), "_jet_dR_0p50");
    executeEventWithParameter("mm"+GetEraShort(), "_jet_dR_0p60");
    executeEventWithParameter("mm"+GetEraShort(), "_jet_dR_0p70");
    executeEventWithParameter("mm"+GetEraShort(), "_jet_dR_0p80");
    executeEventWithParameter("mm"+GetEraShort(), "_jet_dR_0p90");
    executeEventWithParameter("mm"+GetEraShort(), "_jet_dR_1p00");
    */
    if(HasFlag("SYS")){
      vector<Muon> muons_2p1 = MuonMomentumCorrection(SMPGetMuons("POGMediumWithLooseTrkIso", 8.0, 2.1), 0,0, true);
      for(TString syst:{"_jet_scale_up", "_jet_scale_down", "_jet_smear_up", "_jet_smear_down", "_jet_dR_0p20", "_jet_dR_0p30", "_jet_dR_0p60", "_jet_dR_0p80"}){//#### , "_jet_pt_up", "_jet_pt_down", "_jet_eta_down", "_lep_pt_up", "_lep_pt_down", "_lep_eta_down"}){
        muons = syst.Contains("lep_eta_down")? muons_2p1: muons_2p4;
        if(!(syst.Contains("smear") && IsDATA)) executeEventWithParameter("mm"+GetEraShort(), syst);
      }
      muons = muons_2p4;
      for(unsigned int s=0; s<nmem_muon.size(); s++){
        for(unsigned int m=0; m<nmem_muon.at(s); m++){
          executeEventWithParameter("mm"+GetEraShort(), Form("_MuonMomentum_s%dm%d", s, m), s,m);
        }
      }
    }
  }
  if(!IsDATA || DataStream.Contains("DoubleEG") || DataStream.Contains("EGamma")){
    vector<Electron> electrons_2p5 = ElectronEnergyCorrection(SMPGetElectrons("passMediumID", 8.0, 2.5), 0,0, true);
    electrons = electrons_2p5;
    executeEventWithParameter("ee"+GetEraShort());
    /*
    executeEventWithParameter("ee"+GetEraShort(), "_jet_dR_0p01");
    executeEventWithParameter("ee"+GetEraShort(), "_jet_dR_0p05");
    executeEventWithParameter("ee"+GetEraShort(), "_jet_dR_0p10");
    executeEventWithParameter("ee"+GetEraShort(), "_jet_dR_0p15");
    executeEventWithParameter("ee"+GetEraShort(), "_jet_dR_0p20");
    executeEventWithParameter("ee"+GetEraShort(), "_jet_dR_0p25");
    executeEventWithParameter("ee"+GetEraShort(), "_jet_dR_0p30");
    executeEventWithParameter("ee"+GetEraShort(), "_jet_dR_0p35");
    executeEventWithParameter("ee"+GetEraShort(), "_jet_dR_0p50");
    executeEventWithParameter("ee"+GetEraShort(), "_jet_dR_0p60");
    executeEventWithParameter("ee"+GetEraShort(), "_jet_dR_0p70");
    executeEventWithParameter("ee"+GetEraShort(), "_jet_dR_0p80");
    executeEventWithParameter("ee"+GetEraShort(), "_jet_dR_0p90");
    executeEventWithParameter("ee"+GetEraShort(), "_jet_dR_1p00");
    */
    if(HasFlag("SYS")){
      vector<Electron> electrons_2p1 = ElectronEnergyCorrection(SMPGetElectrons("passMediumID", 8.0, 2.1), 0,0, true);
      for(TString syst:{"_jet_scale_up", "_jet_scale_down", "_jet_smear_up", "_jet_smear_down", "_jet_dR_0p20", "_jet_dR_0p30", "_jet_dR_0p60", "_jet_dR_0p80"}){//#### , "_jet_pt_up", "_jet_pt_down", "_jet_eta_down", "_lep_pt_up", "_lep_pt_down", "_lep_eta_down"}){
        electrons = syst.Contains("lep_eta_down")? electrons_2p1: electrons_2p5;
        if(!(syst.Contains("smear") && IsDATA)) executeEventWithParameter("ee"+GetEraShort(), syst);
      }
      electrons = electrons_2p5;
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
  prefix = channel+"/"+gprefix, hprefix = "", suffix = "";
  suffix += option;
  if((IsNominalRun || option != "") && set != 1) IsNominalLike = true;
  else IsNominalLike = false;

  // Weights Setup
  lumiweight = 1., PUweight = 1., prefireweight = 1.;
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

  //if(!IsDATA && weight_Scale->size() > 0){
  //  map_weight[""] *= TMath::Range(-10., 10., weight_Scale->at(8));
  //  if(IsNominalLike){
  //    FillHist(prefix+hprefix+"weight_Scale8"+suffix, lumiweight, map_weight[""], 200,-5,5);
  //    FillCutflow(prefix+hprefix+"cutflow"+suffix, "Scale8", map_weight[""]);
  //  }
  //}

  // Trigger
  if(!IsFiredTriggers(channel)) return;
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "PassTrig", map_weight[""]);

  // MET Filter
  if(!PassMETFilter()) return;
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "METfilter", map_weight[""]);

  // Dilepton + pT + OS + Mass
  int lep_pt_sys = 0;
  if(option.Contains("lep_pt_up")) lep_pt_sys = 1;
  else if(option.Contains("lep_pt_down")) lep_pt_sys = -1;
  if(!HasDileptons(channel, set, mem, lep_pt_sys)) return;

  // Jets
  vector<Jet> alljets = {};
  if(option.Contains("jet_scale_up")) alljets = SelectJets(ScaleJets(jets_raw, 1), "tightLepVeto", 30, (DataYear == 2016? 2.4: 2.5));
  else if(option.Contains("jet_scale_down")) alljets = SelectJets(ScaleJets(jets_raw, -1), "tightLepVeto", 30, (DataYear == 2016? 2.4: 2.5));
  else if(option.Contains("jet_smear_up")) alljets = SelectJets(SmearJets(jets_raw, 1), "tightLepVeto", 30, (DataYear == 2016? 2.4: 2.5));
  else if(option.Contains("jet_smear_down")) alljets = SelectJets(SmearJets(jets_raw, -1), "tightLepVeto", 30, (DataYear == 2016? 2.4: 2.5));
  else if(option.Contains("jet_pt_up")) alljets = SelectJets(jets_raw, "tightLepVeto", 30, (DataYear == 2016? 2.4: 2.5));
  else if(option.Contains("jet_pt_down")) alljets = SelectJets(jets_raw, "tightLepVeto", 20, (DataYear == 2016? 2.4: 2.5));
  else if(option.Contains("jet_eta_down")) alljets = SelectJets(jets_raw, "tightLepVeto", 30, 2.1);
  else alljets = SelectJets(jets_raw, "tightLepVeto", 30, (DataYear == 2016? 2.4: 2.5));
  std::sort(alljets.begin(), alljets.end(), PtComparing);

  vector<Jet> lepvetojets = {}, realjets_before_vetomap = {}, realjets = {}, bjets = {}, ajets = {};
  double dR_lj = 0.4;
  if(option.Contains("jet_dR_0p01")) dR_lj = 0.01;
  else if(option.Contains("jet_dR_0p05")) dR_lj = 0.05;
  else if(option.Contains("jet_dR_0p10")) dR_lj = 0.10;
  else if(option.Contains("jet_dR_0p15")) dR_lj = 0.15;
  else if(option.Contains("jet_dR_0p20")) dR_lj = 0.2;
  else if(option.Contains("jet_dR_0p25")) dR_lj = 0.25;
  else if(option.Contains("jet_dR_0p30")) dR_lj = 0.3;
  else if(option.Contains("jet_dR_0p35")) dR_lj = 0.35;
  else if(option.Contains("jet_dR_0p50")) dR_lj = 0.5;
  else if(option.Contains("jet_dR_0p60")) dR_lj = 0.6;
  else if(option.Contains("jet_dR_0p70")) dR_lj = 0.7;
  else if(option.Contains("jet_dR_0p80")) dR_lj = 0.8;
  else if(option.Contains("jet_dR_0p90")) dR_lj = 0.9;
  else if(option.Contains("jet_dR_1p00")) dR_lj = 1.0;
  for(const auto& jet:alljets){
    if(lepton0 && jet.DeltaR(*lepton0) < dR_lj) continue;
    if(lepton1 && jet.DeltaR(*lepton1) < dR_lj) continue;
    lepvetojets.push_back(jet);
  }
  for(const auto& jet:lepvetojets){
    if(!PUJetIDPass(jet, "Loose")) continue;
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
  JetTagging::Parameters DeepJet_Loose = JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Loose, JetTagging::incl, JetTagging::comb);

  // Check jetcharge vs. b-tagging score, pt, eta, nPV correlations in MC
  if(!IsDATA && IsNominalRun) Checks_bjet_information(prefix+hprefix, realjets, map_weight[""]);

  for(auto& jet:realjets){
    //jet *= jet.BJetNNCorrection(); // full bJetEnergyCorrection (BBjetRegression)?
    if(jet.GetTaggerResult(DeepJet_Tight.j_Tagger) > mcCorr->GetJetTaggingCutValue(DeepJet_Tight.j_Tagger, DeepJet_Tight.j_WP)) bjets.push_back(jet);
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
  //map_weight[""] *= zptweight;
  //if(IsNominalLike){
  //  FillHist(prefix+hprefix+"weight_Zpt"+suffix, zptweight, map_weight[""], 200,-5,5);
  //  FillCutflow(prefix+hprefix+"cutflow"+suffix, "Zpt", map_weight[""]);
  //}
  map_weight[""] *= topptweight;
  if(IsNominalLike){
    FillHist(prefix+hprefix+"weight_Toppt"+suffix, topptweight, map_weight[""], 200,-5,5);
    FillCutflow(prefix+hprefix+"cutflow"+suffix, "Toppt", map_weight[""]);
  }

  // New Weak corrections with sin2w variations
  double lhe_mass = -2.;
  double lhe_costheta_CS = -2.;
  double lhe_costheta_Recoil = -2.;
  if(IsDYSample){
    lhe_mass = ((Particle)lhe_l0 + (Particle)lhe_l1).M();
    if(lhe_p0.ID() == 21 || lhe_p0.ID() == 22){
      if(lhe_p1.ID() == 21 || lhe_p1.ID() == 22) lhe_costheta_CS = GetCosThetaCS(&lhe_l0, &lhe_l1, 0);
      else if(lhe_p1.ID() > 0) lhe_costheta_CS = GetCosThetaCS(&lhe_l0, &lhe_l1, -1);
      else if(lhe_p1.ID() < 0) lhe_costheta_CS = GetCosThetaCS(&lhe_l0, &lhe_l1, 1);
    }else if(lhe_p0.ID() > 0){
      if(lhe_p1.ID() == 21 || lhe_p1.ID() == 22) lhe_costheta_CS = GetCosThetaCS(&lhe_l0, &lhe_l1, 1);
      else if(lhe_p1.ID() > 0) lhe_costheta_CS = GetCosThetaCS(&lhe_l0, &lhe_l1, 0);
      else if(lhe_p1.ID() < 0) lhe_costheta_CS = GetCosThetaCS(&lhe_l0, &lhe_l1, 1);
    }else if(lhe_p0.ID() < 0){
      if(lhe_p1.ID() == 21 || lhe_p1.ID() == 22) lhe_costheta_CS = GetCosThetaCS(&lhe_l0, &lhe_l1, -1);
      else if(lhe_p1.ID() > 0) lhe_costheta_CS = GetCosThetaCS(&lhe_l0, &lhe_l1, -1);
      else if(lhe_p1.ID() < 0) lhe_costheta_CS = GetCosThetaCS(&lhe_l0, &lhe_l1, 0);
    }
    if(lhe_costheta_CS == -2.){
      cout<<"wrong pid for parton: "<<lhe_p0.ID()<<" "<<lhe_p1.ID()<<endl;
      exit(EXIT_FAILURE);
    }

    if(lhe_j0.ID() != 0) lhe_costheta_Recoil = GetCosThetaRecoil((Particle*)&lhe_l0, (Particle*)&lhe_l1, (Particle*)&lhe_j0, (lhe_j0.ID() < 0? 1: -1));

    if(IsNominalRun){
      FillHist(prefix+hprefix+"lhe_mass"+suffix, lhe_mass, map_weight[""], 205,-5,200);
      FillHist(prefix+hprefix+"lhe_costheta_CS"+suffix, lhe_costheta_CS, map_weight[""], 40,-2,2);
      FillHist(prefix+hprefix+"lhe_costheta_Recoil"+suffix, lhe_costheta_Recoil, map_weight[""], 40,-2,2);
      FillHist(prefix+hprefix+"lhe_mass_costheta_Recoil"+suffix, lhe_mass, lhe_costheta_Recoil, map_weight[""], 205,-5,200, 40,-2,2);
      FillHist(prefix+hprefix+"lhe_j0_ID_costheta_Recoil"+suffix, lhe_j0.ID(), lhe_costheta_Recoil, map_weight[""], 50,-25,25, 40,-2,2);
    }
    //weakweight = GetDYWeakWeight(lhe_mass, lhe_costheta_CS, 0, 5);
    weakweight = GetDYWeakWeight(lhe_mass, lhe_costheta_Recoil, 1, 5, lhe_j0.ID());
    //weakweight = GetDYWeakWeight(lhe_mass, lhe_costheta_CS, 2);

    if(IsNominalRun){
      FillHist(prefix+hprefix+"lhe_mass_weakweight"+suffix, lhe_mass, weakweight, map_weight[""], 205,-5,200, 200,0,5);
      FillHist(prefix+hprefix+"lhe_costheta_Recoil_weakweight"+suffix, lhe_costheta_Recoil, weakweight, map_weight[""], 40,-2,2, 200,0,5);
      FillHist(prefix+hprefix+"lhe_j0_ID_weakweight"+suffix, lhe_j0.ID(), weakweight, map_weight[""], 50,-25,25, 200,0,5);

      FillHist(prefix+hprefix+"lhe_mass_weakweight2"+suffix, lhe_mass, GetDYWeakWeight(lhe_mass, lhe_costheta_Recoil, 1, 5, lhe_j0.ID()), map_weight[""], 205,-5,200, 200,0,5);
      FillHist(prefix+hprefix+"lhe_costheta_Recoil_weakweight2"+suffix, lhe_costheta_Recoil, GetDYWeakWeight(lhe_mass, lhe_costheta_Recoil, 1, 5, lhe_j0.ID()), map_weight[""], 40,-2,2, 200,0,5);
      FillHist(prefix+hprefix+"lhe_j0_ID_weakweight2"+suffix, lhe_j0.ID(), GetDYWeakWeight(lhe_mass, lhe_costheta_Recoil, 1, 5, lhe_j0.ID()), map_weight[""], 50,-25,25, 200,0,5);
    }
  }
  map_weight[""] *= weakweight;
  if(IsNominalLike){
    FillHist(prefix+hprefix+"weight_Weak"+suffix, weakweight, map_weight[""], 200,-5,5);
    FillCutflow(prefix+hprefix+"cutflow"+suffix, "Weak", map_weight[""]);
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

    if(HasFlag("SYS") && option == ""){
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

        if(HasFlag("SYS") && option == ""){
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

        if(HasFlag("SYS") && option == ""){
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
      if(HasFlag("SYS") && option == ""){
        for(unsigned int s=0; s<muonTriggerSF_sys.size(); s++){
          for(unsigned int m=0; m<muonTriggerSF_sys[s].size(); m++){
            muonTriggerSF_sys[s][m] *= GetDileptonTriggerSF(muonTriggerLeg1SF_key, muonTriggerLeg2SF_key, muonTriggerDZSF_key, leptons, s,m);
          }
        }
      }
    }
    if(channel.Contains("ee"+GetEraShort())){
      electronTriggerSF *= GetDileptonTriggerSF(electronTriggerLeg1SF_key, electronTriggerLeg2SF_key, electronTriggerDZSF_key, leptons, 0,0);
      if(HasFlag("SYS") && option == ""){
        for(unsigned int s=0; s<electronTriggerSF_sys.size(); s++){
          for(unsigned int m=0; m<electronTriggerSF_sys[s].size(); m++){
            electronTriggerSF_sys[s][m] *= GetDileptonTriggerSF(electronTriggerLeg1SF_key, electronTriggerLeg2SF_key, electronTriggerDZSF_key, leptons, s,m);
          }
        }
      }
    }
  }

  // electron charge-flip SF
  chargeflipSF = 1.;
  if(!IsDATA){
    truth_lepton0 = SMPGetGenMatchedLepton(*lepton0, gens);
    truth_lepton1 = SMPGetGenMatchedLepton(*lepton1, gens);
    chargeflipSF = GetCFSF(0);
  }
  map_weight[""] *= chargeflipSF;
  if(IsNominalLike){
    FillHist(prefix+hprefix+"weight_CFSF"+suffix, chargeflipSF, map_weight[""], 200,-5,5);
    FillCutflow(prefix+hprefix+"cutflow"+suffix, "CFSF", map_weight[""]);
  }

  // lepton efficiency SF
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

  //==== Weights of Systematics
  if(!IsDATA && HasFlag("SYS") && option == ""){
    // Prefiring weight
    map_weight["_noprefireweight"] =  map_weight[""] / prefireweight;
    map_weight["_prefireweight_up"] =  map_weight[""] / prefireweight * L1PrefireReweight_Up;
    map_weight["_prefireweight_down"] = map_weight[""] / prefireweight * L1PrefireReweight_Down;

    // PU reweight
    map_weight["_noPUweight"] = map_weight[""] / PUweight;
    map_weight["_PUweight_up"] = map_weight[""] / PUweight * GetPileUpWeight(nPileUp, 1);
    map_weight["_PUweight_down"] = map_weight[""] / PUweight * GetPileUpWeight(nPileUp, -1);

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
  }else if(!IsDATA && HasFlag("PDFSYS") && option == ""){
    // Zpt Reweight
    //map_weight["_noZpt"] =  map_weight[""] / zptweight;
    //map_weight["_Zpt_gym"] =  map_weight[""] / zptweight * zptweight_gym;
    map_weight["_Zpt"] =  map_weight[""] * zptweight;
    map_weight["_Zpt_gym"] =  map_weight[""] * zptweight_gym;

    // Weak Reweight
    map_weight["_noWeak"] =  map_weight[""] / weakweight;
    if(IsDYSample) map_weight["_oldWeak"] =  map_weight[""] / weakweight * GetDYWeakWeight(lhe_mass, 0., 2);

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
  }else if(!IsDATA && MCSample.Contains("MiNNLO") && IsNominalRun && !hprefix.Contains("ss_")){
    for(unsigned int i=0; i<weight_sthw2->size(); i++) map_weight[Form("_sthw2_%d", i)] = map_weight[""] * weight_sthw2->at(i);

    // New Weak corrections with sin2w variations
    for(unsigned int mem=0; mem<NWEIGHTS; mem++){
      // map_weight["_CS_"+TString(WEIGHT_NAMES[mem])] = map_weight[""] / weakweight * GetDYWeakWeight(lhe_mass, lhe_costheta_CS, 0, mem);
      map_weight["_Recoil_"+TString(WEIGHT_NAMES[mem])] = map_weight[""] / weakweight * GetDYWeakWeight(lhe_mass, lhe_costheta_Recoil, 1, mem, lhe_j0.ID());
    }
  }

  // Electron charge flip SF
  if(!IsDATA && IsNominalRun){
    map_weight["_noCFSF"] = map_weight[""] / chargeflipSF;
    map_weight["_CFSF_up"] = map_weight[""] / chargeflipSF * GetCFSF(1);
    map_weight["_CFSF_down"] = map_weight[""] / chargeflipSF * GetCFSF(-1);
  }

  double weight_default = map_weight[""];
  if((HasFlag("SYS") || HasFlag("PDFSYS")) && option == "") map_weight.erase("");

  //==== Inclusive DY
  double dimass = (*lepton0 + *lepton1).M();
  double dirap = (*lepton0 + *lepton1).Rapidity();
  double dipt = (*lepton0 + *lepton1).Pt();
  double costhetaCS = GetCosThetaCS(lepton0, lepton1);

  if(bjets.size() == 0) return;
  //if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "Tightnb", map_weight[""]);
  //bjets.at(0) *= bjets.at(0).BJetNNCorrection(); // bJetEnergyCorrection (BjetRegression)?
  jet0 = &bjets.at(0);
  double bcharge = jetCharge(*jet0);
  bchargeSF = GetbChargeSFWeight(bjets, 1, 0);
  double costhetaRecoil = GetCosThetaRecoil(lepton0, lepton1, jet0, bcharge);
  double costhetaRecoil_nobch = GetCosThetaRecoil(lepton0, lepton1, jet0, 0, -1);

  if(IsDYSample){
    LHE* lhe_b0 = NULL;
    double mindR = 0.4;
    for(unsigned int i=0; i<HSb.size(); i++){
      if(HSb.at(i).DeltaR(*jet0) > mindR) continue;
      lhe_b0 = &HSb.at(i);
      mindR = lhe_b0->DeltaR(*jet0);
    }

    bool isSignal = lhe_prefix.Contains("sig") && (lhe_b0);
    if(isSignal){
      if(jet0->partonFlavour() == 5) prefix += "dyb_";
      else if(jet0->partonFlavour() == -5) prefix += "dyB_";
    }
  }

  weight_default *= pujetSF;
  map_weight = map_weight * pujetSF;
  if(IsNominalLike){
    FillHist(prefix+hprefix+"weight_PUjetSF"+suffix, pujetSF, map_weight[""], 200,-5,5);
    FillCutflow(prefix+hprefix+"cutflow"+suffix, "PUjetSF", map_weight[""]);
  }
  weight_default *= btagSF;
  map_weight = map_weight * btagSF;
  if(IsNominalLike){
    FillHist(prefix+hprefix+"weight_btagSF"+suffix, btagSF, map_weight[""], 200,-5,5);
    FillCutflow(prefix+hprefix+"cutflow"+suffix, "btagSF", map_weight[""]);
  }
  weight_default *= bchargeSF;
  map_weight = map_weight * bchargeSF;
  if(IsNominalLike){
    FillHist(prefix+hprefix+"weight_bChargeSF1"+suffix, bchargeSF, map_weight[""], 200,-5,5);
    FillCutflow(prefix+hprefix+"cutflow"+suffix, "bChargeSF1", map_weight[""]);
  }
  //map_weight = map_weight * GetAdhocbChargeSFWeight(bcharge);
  //if(IsNominalLike){
  //  FillHist(prefix+hprefix+"weight_bChargeSF_adhoc"+suffix, bchargeSF, map_weight[""], 200,-5,5);
  //  FillCutflow(prefix+hprefix+"cutflow"+suffix, "bChargeSF_adhoc", map_weight[""]);
  //}


  //==== Weights of Systematics (regarding jets)
  if(!IsDATA && HasFlag("SYS") && option == ""){
    // b-tagging SF
    map_weight["_nobtagSF"] =  weight_default / btagSF;
    map_weight["_btagSF_hup"] = weight_default / btagSF * GetBTaggingReweight_1a_2WP(realjets, DeepJet_Tight, DeepJet_Loose, "SystUpHTag");
    map_weight["_btagSF_hdown"] = weight_default / btagSF * GetBTaggingReweight_1a_2WP(realjets, DeepJet_Tight, DeepJet_Loose, "SystDownHTag");
    map_weight["_btagSF_hcorr"] = weight_default / btagSF * GetBTaggingReweight_1a_2WP(realjets, DeepJet_Tight, DeepJet_Loose, "SystUpHTagCorr");
    map_weight["_btagSF_huncorr"] = weight_default / btagSF * GetBTaggingReweight_1a_2WP(realjets, DeepJet_Tight, DeepJet_Loose, "SystUpHTagUnCorr");
    map_weight["_btagSF_lup"] = weight_default / btagSF * GetBTaggingReweight_1a_2WP(realjets, DeepJet_Tight, DeepJet_Loose, "SystUpLTag");
    map_weight["_btagSF_ldown"] = weight_default / btagSF * GetBTaggingReweight_1a_2WP(realjets, DeepJet_Tight, DeepJet_Loose, "SystDownLTag");
    map_weight["_btagSF_lcorr"] = weight_default / btagSF * GetBTaggingReweight_1a_2WP(realjets, DeepJet_Tight, DeepJet_Loose, "SystUpLTagCorr");;
    map_weight["_btagSF_luncorr"] = weight_default / btagSF * GetBTaggingReweight_1a_2WP(realjets, DeepJet_Tight, DeepJet_Loose, "SystUpLTagUnCorr");

    // PUjetID SF
    map_weight["_noPUjetSF"] =  weight_default / pujetSF;
    map_weight["_PUjetSF_up"] =  weight_default / pujetSF * GetPUJetWeight(lepvetojets, "Loose", 1);
    map_weight["_PUjetSF_down"] = weight_default / pujetSF * GetPUJetWeight(lepvetojets, "Loose", -1);

    // bChargeID SF
    map_weight["_bChargeSF0"] = weight_default / bchargeSF * GetbChargeSFWeight(bjets, 0, 0);
    map_weight["_bChargeSF0_up"] = weight_default / bchargeSF * GetbChargeSFWeight(bjets, 0, 1);
    map_weight["_bChargeSF0_down"] = weight_default / bchargeSF * GetbChargeSFWeight(bjets, 0, -1);
    map_weight["_nobChargeSF1"] = weight_default / bchargeSF;
    for(TString bCh:{"0", "1", "2", "3", "4", "5"}){
      map_weight["_bChargeSF1_up"+bCh] = weight_default / bchargeSF * GetbChargeSFWeight(bjets, 1, 1, bCh);
      map_weight["_bChargeSF1_down"+bCh] = weight_default / bchargeSF * GetbChargeSFWeight(bjets, 1, -1, bCh);
    }
    map_weight["_bChargeSFHS"] = weight_default / bchargeSF * GetbChargeSFWeight(bjets, 2, 0);
    map_weight["_bChargeSFHS_up"] = weight_default / bchargeSF * GetbChargeSFWeight(bjets, 2, 1);
    map_weight["_bChargeSFHS_down"] = weight_default / bchargeSF * GetbChargeSFWeight(bjets, 2, -1);
  }else if(!IsDATA && IsNominalRun && !hprefix.Contains("ss_")){
    map_weight["_bChargeSF_adhoc"] = weight_default * GetAdhocbChargeSFWeight(bcharge);
    //map_weight["_nobChargeSF_adhoc"] = weight_default / GetAdhocbChargeSFWeight(bcharge);

    map_weight["_Zpt_adhoc"] = weight_default * GetAdhocZptWeight(dipt);
    map_weight["_njets_adhoc"] = weight_default * GetAdhocNjetsWeight(realjets.size());
  }

  if(bjets.size() != 1) return;
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "Tight1b", map_weight[""]);

  if(ajets.size() > 0) return;
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "Veto2b", map_weight[""]);

  if(PuppiMET_Type1_pt > 60) return;
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "MET60", map_weight[""]);

  if((*lepton0 + *lepton1).Pt() < 15) return;
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "ZpT15", map_weight[""]);

  // Checks QCD uncertainty controls
  FillHist(prefix+hprefix+"nalljets_zpt15"+suffix, min((int)alljets.size(), 5), map_weight, 6,0,6);
  FillHist(prefix+hprefix+"nlepvetojets_zpt15"+suffix, min((int)lepvetojets.size(), 5), map_weight, 6,0,6);
  FillHist(prefix+hprefix+"nrealjets_before_vetomap_zpt15"+suffix, min((int)realjets_before_vetomap.size(), 5), map_weight, 6,0,6);
  FillHist(prefix+hprefix+"nrealjets_zpt15"+suffix, min((int)realjets.size(), 5), map_weight, 6,0,6);

  int nnearjets = 0;
  for(int i=0; i<(int)realjets.size(); i++){
    if(realjets.at(i) == *jet0){
      FillHist(prefix+hprefix+"bidx"+suffix, min(i, 6), map_weight, 6,0,6);
      continue;
    }
    if(realjets.at(i).DeltaR(*jet0) > 1.5) continue;

    FillHist(prefix+hprefix+"jidx"+suffix, min(i, 6), map_weight, 6,0,6);
    FillHist(prefix+hprefix+"jpt"+suffix, realjets.at(i).Pt(), map_weight, 200,0,200);
    FillHist(prefix+hprefix+"jeta"+suffix, realjets.at(i).Eta(), map_weight, 50,-2.5,2.5);
    FillHist(prefix+hprefix+"jphi"+suffix, realjets.at(i).Phi(), map_weight, 64,-3.2,3.2);
    FillHist(prefix+hprefix+"bjdR"+suffix, realjets.at(i).DeltaR(*jet0), map_weight, 30,0,1.5);
    FillHist(prefix+hprefix+"bjdPhi"+suffix, realjets.at(i).DeltaPhi(*jet0), map_weight, 60,-1.5,1.5);
    FillHist(prefix+hprefix+"bjdpt"+suffix, jet0->Pt() - realjets.at(i).Pt(), map_weight, 400,-200,200);
    if(nnearjets == 0){
      FillHist(prefix+hprefix+"j0idx"+suffix, min(i, 6), map_weight, 6,0,6);
      FillHist(prefix+hprefix+"j0pt"+suffix, realjets.at(i).Pt(), map_weight, 200,0,200);
      FillHist(prefix+hprefix+"j0eta"+suffix, realjets.at(i).Eta(), map_weight, 50,-2.5,2.5);
      FillHist(prefix+hprefix+"j0phi"+suffix, realjets.at(i).Phi(), map_weight, 64,-3.2,3.2);
      FillHist(prefix+hprefix+"bj0dR"+suffix, realjets.at(i).DeltaR(*jet0), map_weight, 30,0,1.5);
      FillHist(prefix+hprefix+"bj0dPhi"+suffix, realjets.at(i).DeltaPhi(*jet0), map_weight, 60,-1.5,1.5);
      FillHist(prefix+hprefix+"bj0dpt"+suffix, jet0->Pt() - realjets.at(i).Pt(), map_weight, 400,-200,200);
    }
    if(fabs(jet0->Eta()) < 0.3){
      FillHist(prefix+hprefix+"bjdR_beta0p3"+suffix, realjets.at(i).DeltaR(*jet0), map_weight, 30,0,1.5);
      FillHist(prefix+hprefix+"bjdPhi_beta0p3"+suffix, realjets.at(i).DeltaPhi(*jet0), map_weight, 60,-1.5,1.5);
    }else if(fabs(jet0->Eta()) < 0.6){
      FillHist(prefix+hprefix+"bjdR_beta0p6"+suffix, realjets.at(i).DeltaR(*jet0), map_weight, 30,0,1.5);
      FillHist(prefix+hprefix+"bjdPhi_beta0p6"+suffix, realjets.at(i).DeltaPhi(*jet0), map_weight, 60,-1.5,1.5);
    }else if(fabs(jet0->Eta()) < 0.9){
      FillHist(prefix+hprefix+"bjdR_beta0p9"+suffix, realjets.at(i).DeltaR(*jet0), map_weight, 30,0,1.5);
      FillHist(prefix+hprefix+"bjdPhi_beta0p9"+suffix, realjets.at(i).DeltaPhi(*jet0), map_weight, 60,-1.5,1.5);
    }else if(fabs(jet0->Eta()) < 1.2){
      FillHist(prefix+hprefix+"bjdR_beta1p2"+suffix, realjets.at(i).DeltaR(*jet0), map_weight, 30,0,1.5);
      FillHist(prefix+hprefix+"bjdPhi_beta1p2"+suffix, realjets.at(i).DeltaPhi(*jet0), map_weight, 60,-1.5,1.5);
    }else if(fabs(jet0->Eta()) < 1.5){
      FillHist(prefix+hprefix+"bjdR_beta1p5"+suffix, realjets.at(i).DeltaR(*jet0), map_weight, 30,0,1.5);
      FillHist(prefix+hprefix+"bjdPhi_beta1p5"+suffix, realjets.at(i).DeltaPhi(*jet0), map_weight, 60,-1.5,1.5);
    }else if(fabs(jet0->Eta()) < 1.8){
      FillHist(prefix+hprefix+"bjdR_beta1p8"+suffix, realjets.at(i).DeltaR(*jet0), map_weight, 30,0,1.5);
      FillHist(prefix+hprefix+"bjdPhi_beta1p8"+suffix, realjets.at(i).DeltaPhi(*jet0), map_weight, 60,-1.5,1.5);
    }else if(fabs(jet0->Eta()) < 2.1){
      FillHist(prefix+hprefix+"bjdR_beta2p1"+suffix, realjets.at(i).DeltaR(*jet0), map_weight, 30,0,1.5);
      FillHist(prefix+hprefix+"bjdPhi_beta2p1"+suffix, realjets.at(i).DeltaPhi(*jet0), map_weight, 60,-1.5,1.5);
    }else{
      FillHist(prefix+hprefix+"bjdR_beta2p5"+suffix, realjets.at(i).DeltaR(*jet0), map_weight, 30,0,1.5);
      FillHist(prefix+hprefix+"bjdPhi_beta2p5"+suffix, realjets.at(i).DeltaPhi(*jet0), map_weight, 60,-1.5,1.5);
    }
    nnearjets++;
  }
  FillHist(prefix+hprefix+"nj"+suffix, min(nnearjets, 6), map_weight, 6,0,6);

  if(realjets.size() > 1) return;
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "Veto2j", map_weight[""]);

  // w : 0, w^2 : 1, w*v : 2, w*(w*v) : 3, (w*v)^2 : 4
  FillHist(prefix+hprefix+"statCorr"+suffix, 0, map_weight / map_weight * weight_default, 5,0,5);
  FillHist(prefix+hprefix+"statCorr"+suffix, 1, map_weight / map_weight * (weight_default * weight_default), 5,0,5);
  FillHist(prefix+hprefix+"statCorr"+suffix, 2, map_weight, 5,0,5);
  FillHist(prefix+hprefix+"statCorr"+suffix, 3, map_weight * weight_default, 5,0,5);
  FillHist(prefix+hprefix+"statCorr"+suffix, 4, map_weight * map_weight, 5,0,5);
  FillHist(prefix+hprefix+"mll"+suffix, dimass, map_weight, 80,70,110);
  FillHist(prefix+hprefix+"yll"+suffix, dirap, map_weight, 50,-2.5,2.5);
  FillHist(prefix+hprefix+"ptll"+suffix, dipt, map_weight, 200,0,200);
  FillHist(prefix+hprefix+"lpt"+suffix, lepton0->Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"leta"+suffix, lepton0->Eta(), map_weight, 50,-2.5,2.5);
  FillHist(prefix+hprefix+"l0pt"+suffix, lepton0->Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"l0eta"+suffix, lepton0->Eta(), map_weight, 50,-2.5,2.5);
  FillHist(prefix+hprefix+"lpt"+suffix, lepton1->Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"leta"+suffix, lepton1->Eta(), map_weight, 50,-2.5,2.5);
  FillHist(prefix+hprefix+"l1pt"+suffix, lepton1->Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"l1eta"+suffix, lepton1->Eta(), map_weight, 50,-2.5,2.5);
  FillHist(prefix+hprefix+"bpt"+suffix, jet0->Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"beta"+suffix, jet0->Eta(), map_weight, 50,-2.5,2.5);
  FillHist(prefix+hprefix+"bphi"+suffix, jet0->Phi(), map_weight, 64,-3.2,3.2);
  if(IsNominalRun){
    FillHist(prefix+hprefix+"betaphi"+suffix, jet0->Eta(), jet0->Phi(), map_weight[""], 200,-2.5,2.5, 200,-3.2,3.2);
    if(!hprefix.Contains("ss_"))FillHist("betaphi", jet0->Eta(), jet0->Phi(), map_weight[""], 200,-2.5,2.5, 200,-3.2,3.2);
  }
  FillHist(prefix+hprefix+"bChargeRaw"+suffix, bcharge, map_weight, 200,-5,5);
  FillHist(prefix+hprefix+"bChargeRaw2"+suffix, bcharge, map_weight, 100,-5,5);
  FillHist(prefix+hprefix+"bCharge"+suffix, (bcharge < 0? -0.5: 0.5), map_weight, 2,-1,1);
  FillHist(prefix+hprefix+"met"+suffix, PuppiMET_Type1_pt, map_weight, 200,0,200);
  FillHist(prefix+hprefix+"metphi"+suffix, PuppiMET_Type1_phi, map_weight, 64,-3.2,3.2);
  if(IsNominalRun) FillHist(prefix+hprefix+"met_metphi"+suffix, PuppiMET_Type1_pt, PuppiMET_Type1_phi, map_weight[""], 40,0,200, 32,-3.2,3.2);
  FillHist(prefix+hprefix+"met2"+suffix, PuppiMET_Type1_PhiCor_pt, map_weight, 200,0,200);
  FillHist(prefix+hprefix+"met2phi"+suffix, PuppiMET_Type1_PhiCor_phi, map_weight, 64,-3.2,3.2);
  FillHist(prefix+hprefix+"met_met2"+suffix, PuppiMET_Type1_pt - PuppiMET_Type1_PhiCor_pt, map_weight, 200,-100,100);
  if(IsNominalRun) FillHist(prefix+hprefix+"met2_metphi2"+suffix, PuppiMET_Type1_PhiCor_pt, PuppiMET_Type1_PhiCor_phi, map_weight[""], 40,0,200, 32,-3.2,3.2);
  FillHist(prefix+hprefix+"ZbdR"+suffix, (*lepton0 + *lepton1).DeltaR(*jet0), map_weight, 200,0,10);
  FillHist(prefix+hprefix+"ZbdPhi"+suffix, (*lepton0 + *lepton1).DeltaPhi(*jet0), map_weight, 128,-3.2,3.2);
  FillHist(prefix+hprefix+"Zbpt"+suffix, (*lepton0 + *lepton1 + *jet0).Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"costhetaRecoil"+suffix, dimass, fabs(bcharge), costhetaRecoil, map_weight, afb_mbinnum,(double*)afb_mbin, afb_chbinnum,(double*)afb_chbin, 20,-1,1);
  FillHist(prefix+hprefix+"costhetaRecoil2"+suffix, dimass, fabs(bcharge), costhetaRecoil, map_weight, afb_mbinnum_HS,(double*)afb_mbin_HS, afb_chbinnum,(double*)afb_chbin, 20,-1,1);
  FillHist(prefix+hprefix+"costhetaRecoil3"+suffix, dimass, fabs(bcharge), costhetaRecoil, map_weight, afb_mbinnum_original,(double*)afb_mbin_original, afb_chbinnum,(double*)afb_chbin, 20,-1,1);
  FillHist(prefix+hprefix+"AbscosthetaRecoil"+suffix, dimass, fabs(bcharge), fabs(costhetaRecoil), map_weight, afb_mbinnum,(double*)afb_mbin, afb_chbinnum,(double*)afb_chbin, 10,0,1);
  FillHist(prefix+hprefix+"costhetaRecoil_nobch"+suffix, dimass, fabs(bcharge), costhetaRecoil_nobch, map_weight, afb_mbinnum,(double*)afb_mbin, afb_chbinnum,(double*)afb_chbin, 10,0,1);
  FillHist(prefix+hprefix+"costhetaCS"+suffix, dimass, dirap, costhetaCS, map_weight, afb_mbinnum,(double*)afb_mbin, afb_ybinnum,(double*)afb_ybin, 20,-1,1);
  if(fabs(jet0->Eta()) < 0.3){
    FillHist(prefix+hprefix+"ZbdR_beta0p3"+suffix, (*lepton0 + *lepton1).DeltaR(*jet0), map_weight, 200,0,10);
    FillHist(prefix+hprefix+"ZbdPhi_beta0p3"+suffix, (*lepton0 + *lepton1).DeltaPhi(*jet0), map_weight, 128,-3.2,3.2);
  }else if(fabs(jet0->Eta()) < 0.6){
    FillHist(prefix+hprefix+"ZbdR_beta0p6"+suffix, (*lepton0 + *lepton1).DeltaR(*jet0), map_weight, 200,0,10);
    FillHist(prefix+hprefix+"ZbdPhi_beta0p6"+suffix, (*lepton0 + *lepton1).DeltaPhi(*jet0), map_weight, 128,-3.2,3.2);
  }else if(fabs(jet0->Eta()) < 0.9){
    FillHist(prefix+hprefix+"ZbdR_beta0p9"+suffix, (*lepton0 + *lepton1).DeltaR(*jet0), map_weight, 200,0,10);
    FillHist(prefix+hprefix+"ZbdPhi_beta0p9"+suffix, (*lepton0 + *lepton1).DeltaPhi(*jet0), map_weight, 128,-3.2,3.2);
  }else if(fabs(jet0->Eta()) < 1.2){
    FillHist(prefix+hprefix+"ZbdR_beta1p2"+suffix, (*lepton0 + *lepton1).DeltaR(*jet0), map_weight, 200,0,10);
    FillHist(prefix+hprefix+"ZbdPhi_beta1p2"+suffix, (*lepton0 + *lepton1).DeltaPhi(*jet0), map_weight, 128,-3.2,3.2);
  }else if(fabs(jet0->Eta()) < 1.5){
    FillHist(prefix+hprefix+"ZbdR_beta1p5"+suffix, (*lepton0 + *lepton1).DeltaR(*jet0), map_weight, 200,0,10);
    FillHist(prefix+hprefix+"ZbdPhi_beta1p5"+suffix, (*lepton0 + *lepton1).DeltaPhi(*jet0), map_weight, 128,-3.2,3.2);
  }else if(fabs(jet0->Eta()) < 1.8){
    FillHist(prefix+hprefix+"ZbdR_beta1p8"+suffix, (*lepton0 + *lepton1).DeltaR(*jet0), map_weight, 200,0,10);
    FillHist(prefix+hprefix+"ZbdPhi_beta1p8"+suffix, (*lepton0 + *lepton1).DeltaPhi(*jet0), map_weight, 128,-3.2,3.2);
  }else if(fabs(jet0->Eta()) < 2.1){
    FillHist(prefix+hprefix+"ZbdR_beta2p1"+suffix, (*lepton0 + *lepton1).DeltaR(*jet0), map_weight, 200,0,10);
    FillHist(prefix+hprefix+"ZbdPhi_beta2p1"+suffix, (*lepton0 + *lepton1).DeltaPhi(*jet0), map_weight, 128,-3.2,3.2);
  }else{
    FillHist(prefix+hprefix+"ZbdR_beta2p5"+suffix, (*lepton0 + *lepton1).DeltaR(*jet0), map_weight, 200,0,10);
    FillHist(prefix+hprefix+"ZbdPhi_beta2p5"+suffix, (*lepton0 + *lepton1).DeltaPhi(*jet0), map_weight, 128,-3.2,3.2);
  }

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
      FillHist(Form(prefix+hprefix+"bCharge%dRaw2"+suffix, i-1), bcharge, map_weight, 100,-5,5);
      FillHist(Form(prefix+hprefix+"bCharge%d"+suffix, i-1), (bcharge < 0? -0.5: 0.5), map_weight, 2,-1,1);
      FillHist(prefix+hprefix+"bChargebin"+suffix, LHAPDF::sgn(bcharge) * (i - 0.5), map_weight, 12,-6,6);
      FillHist(prefix+hprefix+"bChargeAbsbin"+suffix, i - 0.5, map_weight, 6,0,6);
      FillHist(Form(prefix_nPV+"bCharge%dRaw"+suffix, i-1), bcharge, map_weight, 200,-5,5);
      FillHist(Form(prefix_nPV+"bCharge%dRaw2"+suffix, i-1), bcharge, map_weight, 100,-5,5);
      FillHist(Form(prefix_nPV+"bCharge%d"+suffix, i-1), (bcharge < 0? -0.5: 0.5), map_weight, 2,-1,1);
      FillHist(prefix_nPV+"bChargebin"+suffix, LHAPDF::sgn(bcharge) * (i - 0.5), map_weight, 12,-6,6);
      FillHist(prefix_nPV+"bChargeAbsbin"+suffix, i - 0.5, map_weight, 6,0,6);
    }
  }
}

bool dybAnalyzer::HasDileptons(TString channel, unsigned int s, unsigned int m, int sys){
  double l0pt = 0., l1pt = 0.;
  if(channel.Contains("mm"+GetEraShort())){
    l0pt = 20. + 3 * sys;
    l1pt = 10. + 2 * sys;
    muons = MuonMomentumCorrection(muons, s,m, false);
    if(muons.size() > 0) lepton0 = &muons.at(0);
    if(muons.size() > 1) lepton1 = &muons.at(1);
  }else if(channel.Contains("ee"+GetEraShort())){
    l0pt = 25. + 2 * sys;
    l1pt = 15. + 3 * sys;
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

  if(hprefix.Contains("ss_") && !IsNominalRun) return false; // To reduce SS histograms

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

  // Selections from Hyonsan's AN (AN-20-216 v3)
  vector<Muon> bmuon;
  if(muons_raw.size() == 0) muons_raw = GetAllMuons();
  for(const auto& mu:muons_raw){
    if(jet.DeltaR(mu) > 0.3) continue;
    if(mu.P() * sin(mu.Angle(jet.Vect())) < 1.0) continue;
    if(mu.IsType(Muon::Type::PFMuon)) bmuon.push_back(mu);
  }
  vector<Electron> belectron;
  if(electrons_raw.size() == 0) electrons_raw = GetAllElectrons();
  for(const auto& el:electrons_raw){
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

void dybAnalyzer::Checks_bjet_information(TString prefix_hist, const vector<Jet> realjets, double weight){
  if(IsDATA) return;

  for(const auto& jet:realjets){
    if(fabs(jet.partonFlavour()) != 5) continue;

    vector<Muon> bmuon;
    if(muons_raw.size() == 0) muons_raw = GetAllMuons();
    for(const auto& mu:muons_raw){
      if(jet.DeltaR(mu) > 0.3) continue;
      if(mu.P() * sin(mu.Angle(jet.Vect())) < 1.0) continue;
      if(mu.IsType(Muon::Type::PFMuon)) bmuon.push_back(mu);
    }
    vector<Electron> belectron;
    if(electrons_raw.size() == 0) electrons_raw = GetAllElectrons();
    for(const auto& el:electrons_raw){
      if(jet.DeltaR(el) > 0.3) continue;
      if(el.P() * sin(el.Angle(jet.Vect())) < 1.0) continue;
      if(!PassID(&el, "SoftElectronID")) continue;
      if(!el.PassConversionVeto()) continue;
      if(el.IsGsfCtfScPixChargeConsistent()) belectron.push_back(el);
    }
    double jetCharge = jet.Charge();
    if(bmuon.size() > 0) jetCharge += 2 * bmuon.at(0).Charge();
    else if(belectron.size() > 0) jetCharge += 4 * belectron.at(0).Charge();

    TString samples = "etc/";
    if(IsDYSample) samples = "DY/";
    else if(IsTTSample) samples = "ttbar/";
    TString suf = "";
    for(TString pre:{(TString)"PID5/", GetEraShort()+"/", prefix_hist, samples}){
      FillHists_bjet_information(pre, suf, jet, bmuon, belectron, weight);

      if(jet.Pt() < 35) suf = "_pt0";
      else if(jet.Pt() < 50) suf = "_pt1";
      else if(jet.Pt() < 80) suf = "_pt2";
      else if(jet.Pt() < 120) suf = "_pt3";
      else suf = "_pt4";
      FillHists_bjet_information(pre, suf, jet, bmuon, belectron, weight);

      if(fabs(jet.Eta()) < 0.6) suf = "_eta0";
      else if(fabs(jet.Eta()) < 1.2) suf = "_eta1";
      else if(fabs(jet.Eta()) < 1.8) suf = "_eta2";
      else suf = "_eta3";
      FillHists_bjet_information(pre, suf, jet, bmuon, belectron, weight);

      for(double ch:{0.7, 0.8, 0.9, 0.95}){
        if(fabs(jet.Charge()) < ch) break;
        suf = Form("_charge0p%d", int(ch * 100));
        FillHists_bjet_information(pre, suf, jet, bmuon, belectron, weight);
      }

      TString PID5 = (jet.partonFlavour() > 0? "b": "B");
      if(jetCharge * jet.partonFlavour() < 0) suf = "_Correct";
      else suf = "_Incorrect";
      FillHist(pre+"PID5_jetcharge"+suf, jet.Charge(), weight, 200,-1,1);
      FillHist(pre+"PID5_jetCharge"+suf, jetCharge, weight, 200,-5,5);
      FillHist(pre+PID5+"_jetcharge"+suf, jet.Charge(), weight, 200,-1,1);
      FillHist(pre+PID5+"_jetCharge"+suf, jetCharge, weight, 200,-5,5);

      double bScore = jet.GetTaggerResult(JetTagging::DeepJet);
      FillHist(pre+"PID5_accuracy", (jet.Charge() * jet.partonFlavour() < 0? 1: 0), weight, 2,0,2);
      FillHist(pre+"PID5_score_jetcharge", bScore, jet.Charge(), weight, 100,0,1, 200,-1,1);
      FillHist(pre+"PID5_score_accuracy", bScore, (jet.Charge() * jet.partonFlavour() < 0? 1: 0), weight, 100,0,1, 2,0,2);
      FillHist(pre+"PID5_Accuracy", (jetCharge * jet.partonFlavour() < 0? 1: 0), weight, 2,0,2);
      FillHist(pre+"PID5_score_jetCharge", bScore, jetCharge, weight, 100,0,1, 200,-5,5);
      FillHist(pre+"PID5_score_Accuracy", bScore, (jetCharge * jet.partonFlavour() < 0? 1: 0), weight, 100,0,1, 2,0,2);
      FillHist(pre+PID5+"_accuracy", (jet.Charge() * jet.partonFlavour() < 0? 1: 0), weight, 2,0,2);
      FillHist(pre+PID5+"_score_jetcharge", bScore, jet.Charge(), weight, 100,0,1, 200,-1,1);
      FillHist(pre+PID5+"_score_accuracy", bScore, (jet.Charge() * jet.partonFlavour() < 0? 1: 0), weight, 100,0,1, 2,0,2);
      FillHist(pre+PID5+"_Accuracy", (jetCharge * jet.partonFlavour() < 0? 1: 0), weight, 2,0,2);
      FillHist(pre+PID5+"_score_jetCharge", bScore, jetCharge, weight, 100,0,1, 200,-5,5);
      FillHist(pre+PID5+"_score_Accuracy", bScore, (jetCharge * jet.partonFlavour() < 0? 1: 0), weight, 100,0,1, 2,0,2);
    }
  }
}

void dybAnalyzer::FillHists_bjet_information(TString prefix_hist, TString suffix_hist, const Jet jet, const vector<Muon> bmuon, const vector<Electron> belectron, double weight){
  TString PID5 = (jet.partonFlavour() > 0? "b": "B");
  double bScore = jet.GetTaggerResult(JetTagging::DeepJet);
  double jetCharge = jet.Charge();
  if(bmuon.size() > 0) jetCharge += 2 * bmuon.at(0).Charge();
  else if(belectron.size() > 0) jetCharge += 4 * belectron.at(0).Charge();
  double chargedAll = jet.chargedHadronEnergyFraction() + jet.chargedEmEnergyFraction() + jet.muonEnergyFraction();
  double neutralAll = jet.neutralHadronEnergyFraction() + jet.neutralEmEnergyFraction();

  for(TString pre:{prefix_hist+"PID5_", prefix_hist+PID5+"_"}){
    for(TString suf:{suffix_hist, suffix_hist+"_correct", suffix_hist+"_incorrect"}){
      if(suf.Contains("_correct") && jet.Charge() * jet.partonFlavour() > 0) continue;
      else if(suf.Contains("_incorrect") && jet.Charge() * jet.partonFlavour() < 0) continue;

      FillHist(pre+"pt"+suf, jet.Pt(), weight, 200,0,200);
      FillHist(pre+"eta"+suf, jet.Eta(), weight, 50,-2.5,2.5);
      FillHist(pre+"phi"+suf, jet.Phi(), weight, 64,-3.2,3.2);
      FillHist(pre+"score"+suf, bScore, weight, 100,0,1);
      FillHist(pre+"nPV"+suf, nPV, weight, 100,0,100);
      FillHist(pre+"jetcharge"+suf, jet.Charge(), weight, 200,-1,1);
      FillHist(pre+"jetCharge"+suf, jetCharge, weight, 200,-5,5);
      FillHist(pre+"nmuons"+suf, bmuon.size(), weight, 10,0,10);
      FillHist(pre+"nelectrons"+suf, belectron.size(), weight, 10,0,10);
      FillHist(pre+"fbmuon"+suf, (bmuon.size() > 0? 1: 0), weight, 2,0,2);
      FillHist(pre+"fbelectron"+suf, (belectron.size() > 0? 1: 0), weight, 2,0,2);
      FillHist(pre+"fblepton"+suf, (bmuon.size() * belectron.size() > 0? 1: 0), weight, 2,0,2);
      FillHist(pre+"chargedHadron"+suf, jet.chargedHadronEnergyFraction(), weight, 100,0,1);
      FillHist(pre+"neutralHadron"+suf, jet.neutralHadronEnergyFraction(), weight, 100,0,1);
      FillHist(pre+"chargedEm"+suf, jet.chargedEmEnergyFraction(), weight, 100,0,1);
      FillHist(pre+"neutralEm"+suf, jet.neutralEmEnergyFraction(), weight, 100,0,1);
      FillHist(pre+"muon"+suf, jet.muonEnergyFraction(), weight, 100,0,1);
      FillHist(pre+"chargedAll"+suf, chargedAll, weight, 100,0,1);
      FillHist(pre+"neutralAll"+suf, neutralAll, weight, 100,0,1);
      FillHist(pre+"chargedMultiplicity"+suf, jet.chargedMultiplicity(), weight, 100,0,100);
      FillHist(pre+"neutralMultiplicity"+suf, jet.neutralMultiplicity(), weight, 100,0,100);
      FillHist(pre+"pt_chargedHadron"+suf, jet.Pt(), jet.chargedHadronEnergyFraction(), weight, 40,0,200, 20,0,1);
      FillHist(pre+"pt_neutralHadron"+suf, jet.Pt(), jet.neutralHadronEnergyFraction(), weight, 40,0,200, 20,0,1);
      FillHist(pre+"pt_chargedEm"+suf, jet.Pt(), jet.chargedEmEnergyFraction(), weight, 40,0,200, 20,0,1);
      FillHist(pre+"pt_neutralEm"+suf, jet.Pt(), jet.neutralEmEnergyFraction(), weight, 40,0,200, 20,0,1);
      FillHist(pre+"pt_muon"+suf, jet.Pt(), jet.muonEnergyFraction(), weight, 40,0,200, 20,0,1);
      FillHist(pre+"pt_chargedAll"+suf, jet.Pt(), chargedAll, weight, 40,0,200, 20,0,1);
      FillHist(pre+"pt_neutralAll"+suf, jet.Pt(), neutralAll, weight, 40,0,200, 20,0,1);
      FillHist(pre+"eta_pt"+suf, jet.Eta(), jet.Pt(), weight, 50,-2.5,2.5, 40,0,200);
      FillHist(pre+"eta_phi"+suf, jet.Eta(), jet.Phi(), weight, 50,-2.5,2.5, 64,-3.2,3.2);
      FillHist(pre+"eta_chargedHadron"+suf, jet.Eta(), jet.chargedHadronEnergyFraction(), weight, 50,-2.5,2.5, 20,0,1);
      FillHist(pre+"eta_chargedAll"+suf, jet.Eta(), chargedAll, weight, 50,-2.5,2.5, 20,0,1);
      FillHist(pre+"eta_neutralAll"+suf, jet.Eta(), neutralAll, weight, 50,-2.5,2.5, 20,0,1);
      FillHist(pre+"eta_chargedMultiplicity"+suf, jet.Eta(), jet.chargedMultiplicity(), weight, 50,-2.5,2.5, 30,0,30);
      FillHist(pre+"chargedHadron_chargedMultiplicity"+suf, jet.chargedHadronEnergyFraction(), jet.chargedMultiplicity(), weight, 20,0,1, 30,0,30);
      FillHist(pre+"muon_chargedMultiplicity"+suf, jet.muonEnergyFraction(), jet.chargedMultiplicity(), weight, 20,0,1, 30,0,30);
      FillHist(pre+"chargedAll_chargedMultiplicity"+suf, chargedAll, jet.chargedMultiplicity(), weight, 20,0,1, 30,0,30);
      FillHist(pre+"neutralHadron_neutralMultiplicity"+suf, jet.neutralHadronEnergyFraction(), jet.neutralMultiplicity(), weight, 20,0,1, 30,0,30);
      FillHist(pre+"neutralEm_neutralMultiplicity"+suf, jet.neutralEmEnergyFraction(), jet.neutralMultiplicity(), weight, 20,0,1, 30,0,30);
      FillHist(pre+"neutralAll_neutralMultiplicity"+suf, neutralAll, jet.neutralMultiplicity(), weight, 20,0,1, 30,0,30);
      FillHist(pre+"chargedHadron_muon"+suf, jet.chargedHadronEnergyFraction(),jet.muonEnergyFraction(), weight, 20,0,1, 20,0,1);
      FillHist(pre+"chargedHadron_neutralHadron"+suf, jet.chargedHadronEnergyFraction(), jet.neutralHadronEnergyFraction(), weight, 20,0,1, 20,0,1);
      FillHist(pre+"chargedHadron_neutralEm"+suf, jet.chargedHadronEnergyFraction(), jet.neutralEmEnergyFraction(), weight, 20,0,1, 20,0,1);
      FillHist(pre+"chargedHadron_neutralAll"+suf, jet.chargedHadronEnergyFraction(), neutralAll, weight, 20,0,1, 20,0,1);
      FillHist(pre+"chargedHadron_neutralMultiplicity"+suf, jet.chargedHadronEnergyFraction(), jet.neutralMultiplicity(), weight, 20,0,1, 30,0,30);
      FillHist(pre+"chargedMultiplicity_neutralMultiplicity"+suf, jet.chargedMultiplicity(), jet.neutralMultiplicity(), weight, 30,0,30, 30,0,30);
      FillHist(pre+"muon_nmuons"+suf, jet.muonEnergyFraction(), bmuon.size(), weight, 20,0,1, 10,0,10);
      FillHist(pre+"muon_fbmuon"+suf, jet.muonEnergyFraction(), (bmuon.size() > 0? 1: 0), weight, 20,0,1, 2,0,2);

      FillHist(pre+"pt_jetcharge"+suf, jet.Pt(), fabs(jet.Charge()), weight, 200,0,200, 20,0,1);
      FillHist(pre+"eta_jetcharge"+suf, jet.Eta(), fabs(jet.Charge()), weight, 50,-2.5,2.5, 20,0,1);
      FillHist(pre+"phi_jetcharge"+suf, jet.Phi(), fabs(jet.Charge()), weight, 64,-3.2,3.2, 20,0,1);
      FillHist(pre+"score_jetcharge"+suf, bScore, fabs(jet.Charge()), weight, 20,0,1, 20,0,1);
      FillHist(pre+"nPV_jetcharge"+suf, nPV, fabs(jet.Charge()), weight, 20,0,100, 20,0,1);
      FillHist(pre+"nmuons_jetcharge"+suf, bmuon.size(), fabs(jet.Charge()), weight, 10,0,10, 20,0,1);
      FillHist(pre+"nelectrons_jetcharge"+suf, belectron.size(), fabs(jet.Charge()), weight, 10,0,10, 20,0,1);
      FillHist(pre+"fbmuon_jetcharge"+suf, (bmuon.size() > 0? 1: 0), fabs(jet.Charge()), weight, 2,0,2, 20,0,1);
      FillHist(pre+"fbelectron_jetcharge"+suf, (belectron.size() > 0? 1: 0), fabs(jet.Charge()), weight, 2,0,2, 20,0,1);
      FillHist(pre+"fblepton_jetcharge"+suf, (bmuon.size() * belectron.size() > 0? 1: 0), fabs(jet.Charge()), weight, 2,0,2, 20,0,1);
      FillHist(pre+"chargedHadron_jetcharge"+suf, jet.chargedHadronEnergyFraction(), fabs(jet.Charge()), weight, 20,0,1, 20,0,1);
      FillHist(pre+"neutralHadron_jetcharge"+suf, jet.neutralHadronEnergyFraction(), fabs(jet.Charge()), weight, 20,0,1, 20,0,1);
      FillHist(pre+"chargedEm_jetcharge"+suf, jet.chargedEmEnergyFraction(), fabs(jet.Charge()), weight, 20,0,1, 20,0,1);
      FillHist(pre+"neutralEm_jetcharge"+suf, jet.neutralEmEnergyFraction(), fabs(jet.Charge()), weight, 20,0,1, 20,0,1);
      FillHist(pre+"muon_jetcharge"+suf, jet.muonEnergyFraction(), fabs(jet.Charge()), weight, 20,0,1, 20,0,1);
      FillHist(pre+"chargedAll_jetcharge"+suf, chargedAll, fabs(jet.Charge()), weight, 20,0,1, 20,0,1);
      FillHist(pre+"neutralAll_jetcharge"+suf, neutralAll, fabs(jet.Charge()), weight, 20,0,1, 20,0,1);
      FillHist(pre+"chargedMultiplicity_jetcharge"+suf, jet.chargedMultiplicity(), fabs(jet.Charge()), weight, 30,0,30, 20,0,1);
      FillHist(pre+"neutralMultiplicity_jetcharge"+suf, jet.neutralMultiplicity(), fabs(jet.Charge()), weight, 30,0,30, 20,0,1);
    }
  }
}

dybAnalyzer::dybAnalyzer(){}
dybAnalyzer::~dybAnalyzer(){}

// From Hyonsan's functions in SMPAnalyzerCore
void dybAnalyzer::executeEventGen(){
  lhe_prefix = "";
  gprefix = "";
  if(IsData) return;

  gens = GetGens();
  if(IsDYSample){
    // LHE Setting
    lhes = GetLHEs();
    LHE lhe_j0_qg = LHE();
    double lhe_j0_qg_pt = 0.1;
    lhe_l0 = LHE(), lhe_l1 = LHE(), lhe_p0 = LHE(), lhe_p1 = LHE(), lhe_j0 = LHE();
    HSb = {}, HSc = {};
    for(int i=0; i<(int)lhes.size(); i++){
      if(lhe_l0.ID() == 0 && (abs(lhes[i].ID()) == 11 || abs(lhes[i].ID()) == 13 || abs(lhes[i].ID()) == 15)) lhe_l0 = lhes[i];
      if(lhe_l0.ID()      && (abs(lhes[i].ID()) == 11 || abs(lhes[i].ID()) == 13 || abs(lhes[i].ID()) == 15)) lhe_l1 = lhes[i];

      if(lhe_p0.ID() == 0 && lhes[i].Status() == -1 && lhes[i].Pz() > 0) lhe_p0 = lhes[i];
      if(lhe_p1.ID() == 0 && lhes[i].Status() == -1 && lhes[i].Pz() < 0) lhe_p1 = lhes[i];

      if(fabs(lhes[i].ID()) == 5 && lhes[i].Status() == 1) HSb.push_back(lhes[i]);
      if(fabs(lhes[i].ID()) == 4 && lhes[i].Status() == 1) HSc.push_back(lhes[i]);
      if((lhes[i].ID() == 21 || abs(lhes[i].ID()) < 7) && lhes[i].Status() == 1){
        if(lhe_j0.ID() == 0) lhe_j0 = lhes[i];
        else if(lhes[i].Pt() > lhe_j0.Pt()) lhe_j0 = lhes[i];
      }
      // for qg collisions
      if(lhe_p0.ID() && lhe_p1.ID()){
        if(lhe_p0.ID() == 21 && abs(lhe_p1.ID()) < 7 && lhes[i].ID() == lhe_p1.ID() && lhes[i].Pt() > lhe_j0_qg_pt){
          lhe_j0_qg = lhes[i];
          lhe_j0_qg_pt = lhes[i].Pt();
        }else if(lhe_p1.ID() == 21 && abs(lhe_p0.ID()) < 7 && lhes[i].ID() == lhe_p0.ID() && lhes[i].Pt() > lhe_j0_qg_pt){
          lhe_j0_qg = lhes[i];
          lhe_j0_qg_pt = lhes[i].Pt();
        }
      }
    }

    if(lhe_p0.ID() == 0 || lhe_p1.ID() == 0){
      cout <<"[dybAnalyzer::executeEventGen] something lhe is wrong in partons"<<endl;
      exit(EXIT_FAILURE);
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

    lhe_prefix = "";
    double dimass_lhe = ((Particle)lhe_l0 + (Particle)lhe_l1).M();
    double dirap_lhe = ((Particle)lhe_l0 + (Particle)lhe_l1).Rapidity();
    double dipt_lhe = ((Particle)lhe_l0 + (Particle)lhe_l1).Pt();
    double genweight = reductionweight * MCweight() * _event.GetTriggerLumi("Full");
    if(IsNominalRun){
      FillHist("lhe/"+lhe_prefix+"pZ_p0", lhe_p0.Pz(), genweight, 1000,-14000,14000);
      FillHist("lhe/"+lhe_prefix+"ID_p0", lhe_p0.ID(), genweight, 60,-30,30);
      FillHist("lhe/"+lhe_prefix+"pZ_p1", lhe_p1.Pz(), genweight, 1000,-14000,14000);
      FillHist("lhe/"+lhe_prefix+"ID_p1", lhe_p1.ID(), genweight, 60,-30,30);
      FillHist("lhe/"+lhe_prefix+"p0p1_ID", lhe_p0.ID() * lhe_p1.ID(), genweight, 2000,-1000,1000);
      FillHist("lhe/"+lhe_prefix+"nHSb", HSb.size(), genweight, 5,0,5);
      FillHist("lhe/"+lhe_prefix+"nHSc", HSc.size(), genweight, 5,0,5);
      FillHist("lhe/"+lhe_prefix+"nHSb_nHSc", HSb.size(), HSc.size(), genweight, 5,0,5, 5,0,5);
      FillHist("lhe/"+lhe_prefix+"p0p1_ID_nHSb", lhe_p0.ID() * lhe_p1.ID(), HSb.size(), genweight, 600,-150,450, 5,0,5);
      FillHist("lhe/"+lhe_prefix+"p0p1_ID_nHSc", lhe_p0.ID() * lhe_p1.ID(), HSc.size(), genweight, 600,-150,450, 5,0,5);
      FillHist("lhe/"+lhe_prefix+"p0p1_ID_nHSbc", lhe_p0.ID() * lhe_p1.ID(), HSb.size() + HSc.size() * 5, genweight, 600,-150,450, 20,0,20);
    }

    if(HSb.size() == 1){ // bg, bq collisions
      LHE b = HSb.at(0);
      if(fabs(lhe_p0.ID() * lhe_p1.ID()) == 105) lhe_prefix = "bg_";
      else lhe_prefix = "bq_";

      double costhetaRecoil_lhe_HSb = GetCosThetaRecoil((Particle*)&lhe_l0, (Particle*)&lhe_l1, (Particle*)&b, (b.ID() < 0? 1: -1));
      if(IsNominalRun){
        FillHist("lhe/"+lhe_prefix+"HSb_pt", b.Pt(), genweight, 200,0,200);
        FillHist("lhe/"+lhe_prefix+"HSb_eta", b.Eta(), genweight, 200,-10,10);
        FillHist("lhe/"+lhe_prefix+"HSbZ_pt", ((Particle)b + (Particle)lhe_l0 + (Particle)lhe_l1).Pt(), genweight, 200,0,200);
        FillHist("lhe/"+lhe_prefix+"HSbZ_dR", b.DeltaR(((Particle)lhe_l0 + (Particle)lhe_l1)), genweight, 200,0,10);
        FillHist("lhe/"+lhe_prefix+"HSbZ_dPhi", fabs(b.DeltaPhi(((Particle)lhe_l0 + (Particle)lhe_l1))), genweight, 200,0,10);
        FillHist("lhe/"+lhe_prefix+"HSbZ_dEta", fabs(b.Eta() - ((Particle)lhe_l0 + (Particle)lhe_l1).Eta()), genweight, 200,0,10);
        FillHist("lhe/"+lhe_prefix+"costhetaRecoil_HSb", dimass_lhe, (b.ID() < 0? 1: -1), costhetaRecoil_lhe_HSb, genweight, afb_mbinnum,(double*)afb_mbin, afb_chbinnum,(double*)afb_chbin, 20,-1,1);
      }
    }else if(HSb.size() > 1){ // bb, bB, gg, qQ collisions
      LHE b0 = HSb.at(0), b1 = HSb.at(1);
      if(b0.Pt() < b1.Pt()){
        b0 = HSb.at(1);
        b1 = HSb.at(0);
      }
      if(lhe_p0.ID() * lhe_p1.ID() == 25) lhe_prefix = "bb_";
      else if(lhe_p0.ID() * lhe_p1.ID() == 25) lhe_prefix = "bB_";
      else if(lhe_p0.ID() * lhe_p1.ID() == 441) lhe_prefix = "gg_";
      else lhe_prefix = "qQ_";

      double costhetaRecoil_lhe_HSb0 = GetCosThetaRecoil((Particle*)&lhe_l0, (Particle*)&lhe_l1, (Particle*)&b0, (b0.ID() < 0? 1: -1));
      double costhetaRecoil_lhe_HSb1 = GetCosThetaRecoil((Particle*)&lhe_l0, (Particle*)&lhe_l1, (Particle*)&b1, (b1.ID() < 0? 1: -1));
      if(IsNominalRun){
        FillHist("lhe/"+lhe_prefix+"HSb0_idx", (b0 == HSb.at(0)? 0: 1), genweight, 2,0,2);
        FillHist("lhe/"+lhe_prefix+"HSb0_pt", b0.Pt(), genweight, 200,0,200);
        FillHist("lhe/"+lhe_prefix+"HSb0_eta", b0.Eta(), genweight, 200,-10,10);
        FillHist("lhe/"+lhe_prefix+"HSb0Z_pt", ((Particle)b0 + (Particle)lhe_l0 + (Particle)lhe_l1).Pt(), genweight, 200,0,200);
        FillHist("lhe/"+lhe_prefix+"HSb0Z_dR", b0.DeltaR(((Particle)lhe_l0 + (Particle)lhe_l1)), genweight, 200,0,10);
        FillHist("lhe/"+lhe_prefix+"HSb0Z_dPhi", fabs(b0.DeltaPhi(((Particle)lhe_l0 + (Particle)lhe_l1))), genweight, 200,0,10);
        FillHist("lhe/"+lhe_prefix+"HSb0Z_dEta", fabs(b0.Eta() - ((Particle)lhe_l0 + (Particle)lhe_l1).Eta()), genweight, 200,0,10);
        FillHist("lhe/"+lhe_prefix+"HSb1_idx", (b1 == HSb.at(0)? 0: 1), genweight, 2,0,2);
        FillHist("lhe/"+lhe_prefix+"HSb1_pt", b1.Pt(), genweight, 200,0,200);
        FillHist("lhe/"+lhe_prefix+"HSb1_eta", b1.Eta(), genweight, 200,-10,10);
        FillHist("lhe/"+lhe_prefix+"HSb1Z_pt", ((Particle)b1 + (Particle)lhe_l0 + (Particle)lhe_l1).Pt(), genweight, 200,0,200);
        FillHist("lhe/"+lhe_prefix+"HSb1Z_dR", b1.DeltaR(((Particle)lhe_l0 + (Particle)lhe_l1)), genweight, 200,0,10);
        FillHist("lhe/"+lhe_prefix+"HSb1Z_dPhi", fabs(b1.DeltaPhi(((Particle)lhe_l0 + (Particle)lhe_l1))), genweight, 200,0,10);
        FillHist("lhe/"+lhe_prefix+"HSb1Z_dEta", fabs(b1.Eta() - ((Particle)lhe_l0 + (Particle)lhe_l1).Eta()), genweight, 200,0,10);
        FillHist("lhe/"+lhe_prefix+"HSbZ_dR", b0.DeltaR(((Particle)lhe_l0 + (Particle)lhe_l1)), genweight, 200,0,10);
        FillHist("lhe/"+lhe_prefix+"HSbZ_dR", b1.DeltaR(((Particle)lhe_l0 + (Particle)lhe_l1)), genweight, 200,0,10);
        FillHist("lhe/"+lhe_prefix+"HSbZ_dPhi", fabs(b0.DeltaPhi(((Particle)lhe_l0 + (Particle)lhe_l1))), genweight, 200,0,10);
        FillHist("lhe/"+lhe_prefix+"HSbZ_dPhi", fabs(b1.DeltaPhi(((Particle)lhe_l0 + (Particle)lhe_l1))), genweight, 200,0,10);
        FillHist("lhe/"+lhe_prefix+"HSb01_dpt", b0.Pt() - b1.Pt(), genweight, 400,-200,200);
        FillHist("lhe/"+lhe_prefix+"HSb01_dR", b0.DeltaR(b1), genweight, 200,0,10);
        FillHist("lhe/"+lhe_prefix+"HSb01_dPhi", fabs(b0.DeltaPhi(b1)), genweight, 200,0,10);
        FillHist("lhe/"+lhe_prefix+"HSb01_dEta", fabs(b0.Eta() - b1.Eta()), genweight, 200,0,10);
        FillHist("lhe/"+lhe_prefix+"HSb01_dAbsEta", fabs(b0.Eta()) - fabs(b1.Eta()), genweight, 200,-10,10);
        FillHist("lhe/"+lhe_prefix+"costhetaRecoil_HSb0", dimass_lhe, (b0.ID() < 0? 1: -1), costhetaRecoil_lhe_HSb0, genweight, afb_mbinnum,(double*)afb_mbin, afb_chbinnum,(double*)afb_chbin, 20,-1,1);
        FillHist("lhe/"+lhe_prefix+"costhetaRecoil_HSb1", dimass_lhe, (b1.ID() < 0? 1: -1), costhetaRecoil_lhe_HSb1, genweight, afb_mbinnum,(double*)afb_mbin, afb_chbinnum,(double*)afb_chbin, 20,-1,1);
      }
    }

    // qg collisions
    if(lhe_j0_qg.ID() != 0){
      lhe_prefix = "";
      double costhetaRecoil_lhe = GetCosThetaRecoil((Particle*)&lhe_l0, (Particle*)&lhe_l1, (Particle*)&lhe_j0_qg, (lhe_j0_qg.ID() < 0? 1: -1));
      if(IsNominalRun) FillHist("lhe/"+lhe_prefix+"costhetaRecoil_qg", dimass_lhe, (lhe_j0_qg.ID() < 0? 1: -1), costhetaRecoil_lhe, genweight, afb_mbinnum,(double*)afb_mbin, afb_chbinnum,(double*)afb_chbin, 20,-1,1);
      if(fabs(lhe_j0_qg.ID()) == 5) lhe_prefix = "bg_";
      else if(fabs(lhe_j0_qg.ID()) == 4) lhe_prefix = "cg_";
      else if(fabs(lhe_j0_qg.ID()) == 3) lhe_prefix = "sg_";
      else if(fabs(lhe_j0_qg.ID()) == 2) lhe_prefix = "ug_";
      else lhe_prefix = "dg_";

      if(IsNominalRun){
        FillHist("lhe/"+lhe_prefix+"costhetaRecoil_qg", dimass_lhe, (lhe_j0_qg.ID() < 0? 1: -1), costhetaRecoil_lhe, genweight, afb_mbinnum,(double*)afb_mbin, afb_chbinnum,(double*)afb_chbin, 20,-1,1);
	FillHist("lhe/"+lhe_prefix+"mll_qg", dimass_lhe, genweight, 200,40,140);
        FillHist("lhe/"+lhe_prefix+"yll_qg", dirap_lhe, genweight, 100,-5,5);
        FillHist("lhe/"+lhe_prefix+"ptll_qg", dipt_lhe, genweight, 200,0,200);
        FillHist("lhe/"+lhe_prefix+"lpt_qg", lhe_l0.Pt(), genweight, 200,0,200);
        FillHist("lhe/"+lhe_prefix+"lpt_qg", lhe_l1.Pt(), genweight, 200,0,200);
        FillHist("lhe/"+lhe_prefix+"leta_qg", lhe_l0.Eta(), genweight, 100,-5,5);
        FillHist("lhe/"+lhe_prefix+"leta_qg", lhe_l1.Eta(), genweight, 100,-5,5);
        FillHist("lhe/"+lhe_prefix+"jpt_qg", lhe_j0_qg.Pt(), genweight, 200,0,200);
        FillHist("lhe/"+lhe_prefix+"jeta_qg", lhe_j0_qg.Eta(), genweight, 100,-5,5);
      }
    }

    // Signal vs Background defintion
    lhe_prefix = "";
    LHE b = LHE();
    if(HSb.size() == 1){
      lhe_prefix = "sig_";
      b = HSb.at(0);
    }else if(HSb.size() == 2){
      if(HSb.at(0).Pt() < HSb.at(1).Pt()) b = HSb.at(1);
      else b = HSb.at(0);
      if(lhe_p0.ID() * lhe_p1.ID() == 441) lhe_prefix = "sig_";
      else lhe_prefix = "bkg_";
    }else lhe_prefix = "bkg_";

    if(b.ID() != 0){
      double costhetaRecoil_lhe = GetCosThetaRecoil((Particle*)&lhe_l0, (Particle*)&lhe_l1, (Particle*)&b, (b.ID() < 0? 1: -1));
      if(IsNominalRun){
        FillHist("lhe/"+lhe_prefix+"costhetaRecoil", dimass_lhe, (b.ID() < 0? 1: -1), costhetaRecoil_lhe, genweight, afb_mbinnum,(double*)afb_mbin, afb_chbinnum,(double*)afb_chbin, 200,-1,1);
        FillHist("lhe/"+lhe_prefix+"ZbdR", (lhe_l0 + lhe_l1).DeltaR(b), genweight, 100,0,10);
        LHE lm, lp;
        if(lhe_l0.ID() > 0){
          lm = lhe_l0;
          lp = lhe_l1;
        }else{
          lp = lhe_l0;
          lm = lhe_l1;
        }
        FillHist("lhe/"+lhe_prefix+"l0ID", lhe_l0.ID(), genweight, 50,-25,25);
        FillHist("lhe/"+lhe_prefix+"l1ID", lhe_l1.ID(), genweight, 50,-25,25);
        FillHist("lhe/"+lhe_prefix+"lmID", lm.ID(), genweight, 50,-25,25);
        FillHist("lhe/"+lhe_prefix+"lpID", lp.ID(), genweight, 50,-25,25);
        FillHist("lhe/"+lhe_prefix+"lmbdR", lm.DeltaR(b), genweight, 100,0,10);
        FillHist("lhe/"+lhe_prefix+"lpbdR", lp.DeltaR(b), genweight, 100,0,10);
        double mindR = min(lm.DeltaR(b), lp.DeltaR(b));
        FillHist("lhe/"+lhe_prefix+"lbdR", mindR, genweight, 100,0,10);
        FillHist("lhe/"+lhe_prefix+"lbdR_costheta", mindR, costhetaRecoil_lhe, genweight, 100,0,10, 200,-1,1);
      }
    }
    if(IsNominalRun){
      FillHist("lhe/"+lhe_prefix+"mll", dimass_lhe, genweight, 200,40,140);
      FillHist("lhe/"+lhe_prefix+"yll", dirap_lhe, genweight, 100,-5,5);
      FillHist("lhe/"+lhe_prefix+"ptll", dipt_lhe, genweight, 200,0,200);
      FillHist("lhe/"+lhe_prefix+"lpt", lhe_l0.Pt(), genweight, 200,0,200);
      FillHist("lhe/"+lhe_prefix+"lpt", lhe_l1.Pt(), genweight, 200,0,200);
      FillHist("lhe/"+lhe_prefix+"leta", lhe_l0.Eta(), genweight, 100,-5,5);
      FillHist("lhe/"+lhe_prefix+"leta", lhe_l1.Eta(), genweight, 100,-5,5);
      FillHist("lhe/"+lhe_prefix+"jpt", lhe_j0_qg.Pt(), genweight, 200,0,200);
      FillHist("lhe/"+lhe_prefix+"jeta", lhe_j0_qg.Eta(), genweight, 100,-5,5);
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
    }//else gprefix += "tau_";
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

double dybAnalyzer::GetCosThetaRecoil(const Particle *p0, const Particle *p1, Particle *b, const double bcharge, int mode){
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
  }else if(mode == -1){// not looking at bcharge
    return ((*lm - *lp) * (*b)) / ((*lm + *lp) * (*b));
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

void dybAnalyzer::FillHist(TString histname,
                            Double_t value_x, Double_t value_y, Double_t value_u,
                            Double_t weight,
                            Int_t n_binx, const Double_t *xbins,
                            Int_t n_biny, const Double_t *ybins,
                            Int_t n_binu, Double_t u_min, Double_t u_max){

  TH3D *this_hist = GetHist3D(histname);
  if( !this_hist ){
    TAxis uaxis(n_binu,u_min,u_max);
    vector<double> ubins={};
    for(int i=1;i<n_binu+2;i++) ubins.push_back(uaxis.GetBinLowEdge(i));
    this_hist = new TH3D(histname, "", n_binx, xbins, n_biny, ybins, n_binu, &ubins[0]);
    this_hist->SetDirectory(NULL);
    maphist_TH3D[histname] = this_hist;
  }

  this_hist->Fill(value_x, value_y, value_u, weight);

}
void dybAnalyzer::FillHist(TString histname, double value_x, double value_y, double value_u, map<TString,double> weights, int n_binx, const double *xbins, int n_biny, const double *ybins, int n_binu, double u_min, double u_max){
  for(const auto& [suffix,weight]:weights) FillHist(histname+suffix,value_x,value_y,value_u,weight,n_binx,xbins,n_biny,ybins,n_binu,u_min,u_max);
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
  cout<<"[dybAnalyzer::SetupPUJetWeight] setting PUJet Weight"<<endl;
  TString datapath = getenv("DATA_DIR");
  TString pujetpath = datapath+"/"+GetEra()+"/ID/PUJet/PUID_106XTraining_ULRun2_EffSFandUncties_v1.root";
  if(IsExists(pujetpath)){
    cout<<"[dybAnalyzer::SetupPUJetWeight] using file "+pujetpath<<endl;
  }else{
    cout<<"[dybAnalyzer::SetupPUJetWeight] no "+pujetpath<<endl;
    return;
  }
  TFile fPUID(pujetpath);

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

bool dybAnalyzer::PUJetIDPass(const Jet jet, TString ID){
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

void dybAnalyzer::SetupJetVetoMap(){
  cout<<"[dybAnalyzer::SetupJetVetoMap] setting jetveto map"<<endl;
  TString filename = "hotjets-UL16.root";
  TString histname = "h2hot_ul16_plus_hbm2_hbp12_qie11";
  if(DataYear == 2017){
    filename = "hotjets-UL17_v2.root";
    histname = "h2hot_ul17_plus_hep17_plus_hbpw89";
  }else if(DataYear == 2018){
    filename = "hotjets-UL18.root";
    histname = "h2hot_ul18_plus_hem1516_and_hbp2m1";
  }

  TString datapath = getenv("DATA_DIR");
  TString jetvetomappath = datapath+"/"+GetEra()+"/JME/"+filename;
  if(IsExists(jetvetomappath)){
    cout<<"[dybAnalyzer::SetupJetVetoMap] using file "+jetvetomappath<<endl;
  }else{
    cout<<"[dybAnalyzer::SetupJetVetoMap] no "+jetvetomappath<<endl;
    return;
  }
  TFile fvetomap(jetvetomappath);

  hvetomap = (TH2D*)fvetomap.Get(histname);
  if(DataYear == 2016){
    TH2D* additionalmap = (TH2D*)fvetomap.Get("h2hot_mc");
    hvetomap->Add(hvetomap, additionalmap);
  }

  hvetomap->SetDirectory(0);
  fvetomap.Close();
}

double dybAnalyzer::GetbChargeSFWeight(const vector<Jet>& jets, unsigned int mode, int sys, TString bChargeBins){
  double weight = 1.;
  if(IsDATA) return weight;

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

    double Charge = jetCharge(jet); // bjpt > 25->30 GeV // More systs // No Normalization & TTJJ
    double alpha_plus_DATA_eff = 0.63137949; //0.63179105; //0.63181282; //0.63079885;
    double alpha_minus_DATA_eff = 0.61534016; //0.61389421; //0.61395198; //0.61295115;
    double alpha_plus_MC_eff = 0.65458415; //0.65428669; //0.65415066;
    double alpha_minus_MC_eff = 0.63863934; //0.63782927; //0.63769426;
    if(sys > 0){ // asym gets smaller
      alpha_plus_DATA_eff += -0.00135034; //-0.00106574; //-0.00088663; //-0.00086244;
      alpha_minus_DATA_eff += 0.00096904; //0.00085326; //0.00074632; //0.00072373;
    }else if(sys < 0){ // SF gets smaller
      alpha_plus_DATA_eff += 0.00223155; //0.00321587; //0.00250343; //0.00237556;
      alpha_minus_DATA_eff += 0.00310961; //0.0040167; //0.00297406; //0.00283084;
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
        alpha_plus_DATA_eff = 0.52391383; //0.52440445; //0.52466892; //0.52455235;
        alpha_minus_DATA_eff = 0.51844103; //0.51737104; //0.51791587; //0.51779908;
        alpha_plus_MC_eff = 0.52958572; //0.52957027; //0.52968978;
        alpha_minus_MC_eff = 0.52553913; //0.52548965; //0.52540184;
        if(sys > 0 && bChargeBins.Contains("0")){
          alpha_plus_DATA_eff += -0.00264795; //-0.00225534; //-0.00159007; //-0.00157464;
          alpha_minus_DATA_eff += 0.00172429; //0.00212763; //0.0021877; //0.00209931;
        }else if(sys < 0 && bChargeBins.Contains("0")){
          alpha_plus_DATA_eff += 0.00236306; //0.00324108; //0.0029865; //0.00287364;
          alpha_minus_DATA_eff += 0.0036289; //0.00343561; //0.00217065; //0.00215546;
        }
      }else if(fabs(Charge) < afb_chbin[2]){// 0.2
        alpha_plus_DATA_eff = 0.57340855; //0.57499135; //0.57549775; //0.57503119;
        alpha_minus_DATA_eff = 0.56040137; //0.56023694; //0.56120889; //0.56072413;
        alpha_plus_MC_eff = 0.58976601; //0.58955363; //0.58939608;
        alpha_minus_MC_eff = 0.5768458; //0.5764066; //0.57641132;
        if(sys > 0 && bChargeBins.Contains("1")){
          alpha_plus_DATA_eff += -0.00322615; //-0.0028866; //-0.00163046; //-0.00176428;
          alpha_minus_DATA_eff += 0.0011049; //0.00109804; //0.00239976; //0.00222184;
        }else if(sys < 0 && bChargeBins.Contains("1")){
          alpha_plus_DATA_eff += 0.00167391; //0.00204427; //0.00365694; //0.00332861;
          alpha_minus_DATA_eff += 0.00488757; //0.00537413; //0.00248462; //0.00264313;
        }
      }else if(fabs(Charge) < afb_chbin[3]){// 0.6
        alpha_plus_DATA_eff = 0.6685798; //0.668161; //0.66795153; //0.66685555;
        alpha_minus_DATA_eff = 0.64465731; //0.64258252; //0.64183908; //0.6406924;
        alpha_plus_MC_eff = 0.70239881; //0.70205254; //0.7016705;
        alpha_minus_MC_eff = 0.67749213; //0.67617536; //0.67580431;
        if(sys > 0 && bChargeBins.Contains("2")){
          alpha_plus_DATA_eff += -0.00195613; //-0.00257117; //-0.00198435; //-0.00192819;
          alpha_minus_DATA_eff += 0.00206392; //0.00138245; //0.00098174; //0.00094884;
        }else if(sys < 0 && bChargeBins.Contains("2")){
          alpha_plus_DATA_eff += 0.00350697; //0.00373721; //0.00278728; //0.00262054;
          alpha_minus_DATA_eff += 0.00332381; //0.00695071; //0.00563382; //0.00532535;
        }
      }else if(fabs(Charge) < afb_chbin[4]){// 1.0
        alpha_plus_DATA_eff = 0.78123981; //0.78853946; //0.78691713; //0.78271638;
        alpha_minus_DATA_eff = 0.74492833; //0.74319937; //0.74175893; //0.73778292;
        alpha_plus_MC_eff = 0.82792016; //0.8273833; //0.82685472;
        alpha_minus_MC_eff = 0.79454412; //0.79299394; //0.79243187;
        if(sys > 0 && bChargeBins.Contains("3")){
          alpha_plus_DATA_eff += -0.00623389; //-0.00347714; //-0.00277093; //-0.00284074;
          alpha_minus_DATA_eff += 0.00333936; //0.0064814; //0.00567753; //0.00513936;
        }else if(sys < 0 && bChargeBins.Contains("3")){
          alpha_plus_DATA_eff += 0.0066016; //0.01612633; //0.01097872; //0.01040963;
          alpha_minus_DATA_eff += 0.01232379; //0.00865146; //0.0053582; //0.00575384;
        }
      }else if(fabs(Charge) < afb_chbin[5]){// 3.0, soft muons
        alpha_plus_DATA_eff = 0.77488593; //0.77495511; //0.77490976; //0.77182254;
        alpha_minus_DATA_eff = 0.77556002; //0.77649167; //0.77647107; //0.77357429;
        alpha_plus_MC_eff = 0.78471623; //0.78536886; //0.78518778;
        alpha_minus_MC_eff = 0.78474967; //0.78548813; //0.78552888;
        if(sys > 0 && bChargeBins.Contains("4")){
          alpha_plus_DATA_eff += -0.00016469; //-0.00350888; //-0.00930939; //-0.00300315;
          alpha_minus_DATA_eff += 0.01251886; //0.0057526; //0.00629106; //0.01019522;
        }else if(sys < 0 && bChargeBins.Contains("4")){
          alpha_plus_DATA_eff += 0.01046404; //0.01193073; //0.00671166; //0.01005651;
          alpha_minus_DATA_eff += 0.00013766; //0.00727731; //0.00993179; //0.00296229;
        }
      }else{// 5.0, soft electrons
        alpha_plus_DATA_eff = 0.75563201; //0.75874849; //0.75926086; //0.75433632;
        alpha_minus_DATA_eff = 0.75563675; //0.75794022; //0.75920095; //0.75499201;
        alpha_plus_MC_eff = 0.75890386; //0.76086731; //0.76063023;
        alpha_minus_MC_eff = 0.76161921; //0.76353419; //0.76337732;
        if(sys > 0 && bChargeBins.Contains("5")){
          alpha_plus_DATA_eff += -0.00718028; //-0.00587855; //-0.00459789; //-0.00443982;
          alpha_minus_DATA_eff += 0.00669405; //0.00639942; //0.00433558; //0.00398256;
        }else if(sys < 0 && bChargeBins.Contains("5")){
          alpha_plus_DATA_eff += 0.01028731; //0.01097021; //0.00822754; //0.00811232;
          alpha_minus_DATA_eff += 0.01103453; //0.01007731; //0.00872532; //0.00904374;
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

double dybAnalyzer::GetAdhocbChargeSFWeight(float bCharge){
  if(IsDATA) return 1.;

  double sf = 1.;
  if(-5.0 < bCharge && bCharge < -3.0) sf = 18300. / 21068.265;
  else if(bCharge < -1.0) sf = 31046. / 33046.638;
  else if(bCharge < -0.6) sf = 23110. / 26347.320;
  else if(bCharge < -0.2) sf = 229552. / 239533.35;
  else if(bCharge < -0.1) sf = 116808. / 117608.03;
  else if(bCharge < -0.0) sf = 131471. / 132406.60;
  else if(bCharge < 0.1) sf = 132580. / 133977.44;
  else if(bCharge < 0.2) sf = 120521. / 121700.46;
  else if(bCharge < 0.6) sf = 246385. / 253806.95;
  else if(bCharge < 1.0) sf = 25875. / 29350.478;
  else if(bCharge < 3.0) sf = 31988. / 33319.860;
  else if(bCharge < 5.0) sf = 18197. / 20949.792;
  else{
    cout<<"[dybAnalyzer::GetAdhocbChargeSFWeight] |bCharge| >= 5.0?, actual value = "<<bCharge<<endl;
    return sf;
  }

  return sf * 1163115.2 / 1125833.;
}

double dybAnalyzer::GetAdhocZptWeight(float Zpt){
  if(IsDATA) return 1.;

  return 0.94307 + Zpt * 0.000976008;
}

double dybAnalyzer::GetAdhocNjetsWeight(unsigned int njets){
  if(IsDATA) return 1.;

  double sf = 1.;
  if(njets == 1) sf = 0.970012;
  else if(njets == 2) sf = 1.05315;
  else if(njets == 3) sf = 1.08546;
  else if(njets == 4) sf = 1.13618;
  else if(njets == 5) sf = 1.21461;
  else if(njets == 6) sf = 1.33214;
  else if(njets == 7) sf = 1.31229;
  else if(njets == 8) sf = 0.770082;
  else if(njets > 8) sf = 0.754005;
  else{
    cout<<"[dybAnalyzer::GetAdhocNjetsWeight] njets == 0?, actual value = "<<njets<<endl;
    return sf;
  }

  return sf;
}

double dybAnalyzer::GetCFSF(int sys){
  if(IsDATA) return 1.;
  if(!hcfrate_data || !hcfrate_mc) return 1.;

  double sf = 1.;
  double cf_data_l0 = GetBinContentUser(hcfrate_data, lepton0->Eta(), lepton0->Pt(), sys);
  double cf_data_l1 = GetBinContentUser(hcfrate_data, lepton1->Eta(), lepton1->Pt(), sys);
  double cf_mc_l0 = GetBinContentUser(hcfrate_mc, lepton0->Eta(), lepton0->Pt(), -sys);
  double cf_mc_l1 = GetBinContentUser(hcfrate_mc, lepton1->Eta(), lepton1->Pt(), -sys);

  if(lepton0 && !truth_lepton0.IsEmpty()){
    if(lepton0->LeptonFlavour() == Lepton::ELECTRON){
      if(lepton0->Charge() * truth_lepton0.Charge() < 0) sf *= cf_data_l0 / cf_mc_l0;
      else{
        double this_sf = (1 - cf_data_l0) / (1 - cf_mc_l0);
        if(isnormal(this_sf)) sf *= this_sf;
      }
    }
  }

  if(lepton1 && !truth_lepton1.IsEmpty()){
    if(lepton1->LeptonFlavour() == Lepton::ELECTRON){
      if(lepton1->Charge() * truth_lepton1.Charge() < 0) sf *= cf_data_l1 / cf_mc_l1;
      else{
        double this_sf = (1 - cf_data_l1) / (1 - cf_mc_l1);
        if(isnormal(this_sf)) sf *= this_sf;
      }
    }
  }

  return sf;
}

// Hyonsan's NLO Weak Corrections
// From weakWeight_CS.cc, weakWeight_RecoilUDG.cc
double dybAnalyzer::GetDYWeakWeight(double lhe_mass, double lhe_costheta, unsigned int set, unsigned int mem, int lead_pid){
  if(IsDATA || !IsDYSample){
    cout<<"[dybAnalyzer::GetDYWeakWeight] ERROR: IsDATA || !IsDYSample"<<endl;
    exit(1);
  }

  // set = 2 (Old Weak NLO - only mass dependent k-factor)
  if(set == 2) return SMPAnalyzerCore::GetDYWeakWeight(lhe_mass);

  if(lhe_mass == -2. || lhe_costheta == -2.) return 1.;

  int ibin = -1;
  if(lhe_mass < mass_edges[0]) lhe_mass = mass_edges[0] + 1e-4;
  if(lhe_mass >= mass_edges[NBINS]) lhe_mass = mass_edges[NBINS] - 1e-4;
  for(int i=0; i<NBINS; ++i){
    if(lhe_mass >= mass_edges[i] && lhe_mass < mass_edges[i+1]){
      ibin = i;
      break;
    }
  }
  if(ibin < 0) return 1.;

  double x = lhe_costheta;
  if(x > 1.0) x = 1.0;
  if(x < -1.0) x = -1.0;

  int ibin_ang = ibin;
  if(ibin_ang > 24) ibin_ang = 24;

  double A0m = 0.;
  double A4m = 0.;
  double A4v = 0.;
  double w_ang = 1.;
  double kf = 1.;

  if(set == 0){ // Weak NLO based on CS frame
    A0m = A0_minnlo[ibin_ang];
    A4m = A4_minnlo[ibin_ang];
    A4v = A4m + deltaA4[mem][ibin_ang];
    kf = kfactor[mem][ibin];
    if(!(kf > 0.0)) kf = 1.0;
  }else if(set == 1){ // Weak NLO based on Recoil frame
    const int apid = fabs(lead_pid);
    char cat = 'X';
    if(lead_pid == 21) cat = 'G';
    else if(apid >= 1 && apid <= 6){
      if(apid == 2 || apid == 4 || apid == 6) cat = 'U';
      else cat = 'D';
    }else{
      cout<<"[dybAnalyzer::GetDYWeakWeight] ERROR: lead_pid is not a quark/gluon: "<<lead_pid<<endl;
      exit(1);
    }

    const double *A0p = nullptr;
    const double *A4p = nullptr;
    const double (*Kp)[NBINS] = nullptr;
    const double (*Dp)[NBINS] = nullptr;
    if(cat == 'U'){
      A0p = A0_minnlo_U;
      A4p = A4_minnlo_U;
      Kp = kfactor_U;
      Dp = deltaA4_U;
    }else if(cat == 'D'){
      A0p = A0_minnlo_D;
      A4p = A4_minnlo_D;
      Kp = kfactor_D;
      Dp = deltaA4_D;
    }else{
      A0p = A0_minnlo_G;
      A4p = A4_minnlo_G;
      Kp = kfactor_G;
      Dp = deltaA4_G;
    }

    A0m = A0p[ibin_ang];
    A4m = A4p[ibin_ang];
    A4v = A4m + Dp[mem][ibin_ang];

    kf = Kp[mem][ibin];
    if(!(kf > 0.0)) kf = 1.0;
  }else cout<<"[dybAnalyzer::GetDYWeakWeight] set is incorrect"<<endl;

  double f_nom = (1.0 + x * x) + 0.5 * A0m * (1.0 - 3.0 * x * x) + A4m * x;
  double f_var = (1.0 + x * x) + 0.5 * A0m * (1.0 - 3.0 * x * x) + A4v * x;
  if(f_nom != 0.0) w_ang = f_var / f_nom;

  return kf * w_ang;
}
