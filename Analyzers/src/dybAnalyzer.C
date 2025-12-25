#include "dybAnalyzer.h"

void dybAnalyzer::initializeAnalyzer(){
  SMPAnalyzerCore::initializeAnalyzer(); //setup eff zpt roc z0 cf
  IsNominalRun = !HasFlag("SYS") && !HasFlag("PDFSYS") && !HasFlag("LEPSYS");

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

  vector<Jet> lepvetojets = {}, realjets_before_vetomap = {}, realjets = {}, bjets = {}, ajets = {};
  for(const auto& jet:alljets){
    if(lepton0 && jet.DeltaR(*lepton0) < 0.4) continue;
    if(lepton1 && jet.DeltaR(*lepton1) < 0.4) continue;
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
  //map_weight[""] *= zptweight;
  //if(IsNominalLike){
  //  FillHist(prefix+hprefix+"weight_Zpt"+suffix, zptweight, map_weight[""], 200,-5,5);
  //  FillCutflow(prefix+hprefix+"cutflow"+suffix, "Zpt", map_weight[""]);
  //}
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

    if(HasFlag("LEPSYS") && option == ""){
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
    if(channel.Contains("mm"+GetEraShort())){
      muonTriggerSF *= GetDileptonTriggerSF(muonTriggerLeg1SF_key, muonTriggerLeg2SF_key, muonTriggerDZSF_key, leptons, 0,0);
      if(HasFlag("LEPSYS") && option == ""){
        for(unsigned int s=0; s<muonTriggerSF_sys.size(); s++){
          for(unsigned int m=0; m<muonTriggerSF_sys[s].size(); m++){
            muonTriggerSF_sys[s][m] *= GetDileptonTriggerSF(muonTriggerLeg1SF_key, muonTriggerLeg2SF_key, muonTriggerDZSF_key, leptons, s,m);
          }
        }
      }
    }
    if(channel.Contains("ee"+GetEraShort())){
      electronTriggerSF *= GetDileptonTriggerSF(electronTriggerLeg1SF_key, electronTriggerLeg2SF_key, electronTriggerDZSF_key, leptons, 0,0);
      if(HasFlag("LEPSYS") && option == ""){
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

  double weight_noeffSF = map_weight[""];
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

  if(IsNominalRun){
    FillHist(prefix+hprefix+"mll_IncDY"+suffix, dimass, map_weight[""], 80,70,110);
    FillHist(prefix+hprefix+"yll_IncDY"+suffix, dirap, map_weight[""], 50,-2.5,2.5);
    FillHist(prefix+hprefix+"ptll_IncDY"+suffix, dipt, map_weight[""], AFBAnalyzer::unfold_0bjet_ptbinnum_reco,AFBAnalyzer::unfold_0bjet_ptbin_reco);
    FillHist(prefix+hprefix+"lpt_IncDY"+suffix, lepton0->Pt(), map_weight[""], AFBAnalyzer::lptbinnum,AFBAnalyzer::lptbin);
    FillHist(prefix+hprefix+"leta_IncDY"+suffix, lepton0->Eta(), map_weight[""], 50,-2.5,2.5);
    FillHist(prefix+hprefix+"l0pt_IncDY"+suffix, lepton0->Pt(), map_weight[""], AFBAnalyzer::lptbinnum,AFBAnalyzer::lptbin);
    FillHist(prefix+hprefix+"l0eta_IncDY"+suffix, lepton0->Eta(), map_weight[""], 50,-2.5,2.5);
    FillHist(prefix+hprefix+"lpt_IncDY"+suffix, lepton1->Pt(), map_weight[""], AFBAnalyzer::lptbinnum,AFBAnalyzer::lptbin);
    FillHist(prefix+hprefix+"leta_IncDY"+suffix, lepton1->Eta(), map_weight[""], 50,-2.5,2.5);
    FillHist(prefix+hprefix+"l1pt_IncDY"+suffix, lepton1->Pt(), map_weight[""], AFBAnalyzer::lptbinnum,AFBAnalyzer::lptbin);
    FillHist(prefix+hprefix+"l1eta_IncDY"+suffix, lepton1->Eta(), map_weight[""], 50,-2.5,2.5);
    FillHist(prefix+hprefix+"met_IncDY"+suffix, PuppiMET_Type1_pt, map_weight[""], 200,0,200);
    FillHist(prefix+hprefix+"metphi_IncDY"+suffix, PuppiMET_Type1_phi, map_weight[""], 64,-3.2,3.2);
    FillHist(prefix+hprefix+"met_metphi_IncDY"+suffix, PuppiMET_Type1_pt, PuppiMET_Type1_phi, map_weight[""], 40,0,200, 32,-3.2,3.2);
    FillHist(prefix+hprefix+"costhetaCS_IncDY"+suffix, dimass, dirap, dipt, costhetaCS, map_weight[""], afb_mbinnum,(double*)afb_mbin, afb_ybinnum,(double*)afb_ybin, afb_ptbinnum,(double*)afb_ptbin, 20,-1,1);
    FillHist(prefix+hprefix+"nPV_IncDY"+suffix, nPV, map_weight[""], 100,0,100);

    FillHist(prefix+hprefix+"nalljets_IncDY"+suffix, alljets.size(), map_weight[""], 15,0,15);
    FillHist(prefix+hprefix+"nlepvetojets_IncDY"+suffix, lepvetojets.size(), map_weight[""], 15,0,15);
    FillHist(prefix+hprefix+"nrealjets_before_vetomap_IncDY"+suffix, realjets_before_vetomap.size(), map_weight[""], 15,0,15);
    FillHist(prefix+hprefix+"nrealjets_IncDY"+suffix, realjets.size(), map_weight[""], 15,0,15);
    FillHist(prefix+hprefix+"nbjets_IncDY"+suffix, bjets.size(), map_weight[""], 10,0,10);

    // For control plots before/after corrections
    // Without efficiency SF
    FillHist(prefix+hprefix+"mll_IncDY_noeffSF"+suffix, dimass, weight_noeffSF, 80,70,110);
    FillHist(prefix+hprefix+"yll_IncDY_noeffSF"+suffix, dirap, weight_noeffSF, 50,-2.5,2.5);
    FillHist(prefix+hprefix+"ptll_IncDY_noeffSF"+suffix, dipt, weight_noeffSF, AFBAnalyzer::unfold_0bjet_ptbinnum_reco,AFBAnalyzer::unfold_0bjet_ptbin_reco);
    FillHist(prefix+hprefix+"lpt_IncDY_noeffSF"+suffix, lepton0->Pt(), weight_noeffSF, AFBAnalyzer::lptbinnum,AFBAnalyzer::lptbin);
    FillHist(prefix+hprefix+"leta_IncDY_noeffSF"+suffix, lepton0->Eta(), weight_noeffSF, 50,-2.5,2.5);
    FillHist(prefix+hprefix+"lpt_IncDY_noeffSF"+suffix, lepton1->Pt(), weight_noeffSF, AFBAnalyzer::lptbinnum,AFBAnalyzer::lptbin);
    FillHist(prefix+hprefix+"leta_IncDY_noeffSF"+suffix, lepton1->Eta(), weight_noeffSF, 50,-2.5,2.5);

    weight_noeffSF *= muonTrackingSF;
    FillHist(prefix+hprefix+"mll_IncDY_noeffSF_1_Tracking"+suffix, dimass, weight_noeffSF, 80,70,110);
    FillHist(prefix+hprefix+"yll_IncDY_noeffSF_1_Tracking"+suffix, dirap, weight_noeffSF, 50,-2.5,2.5);
    FillHist(prefix+hprefix+"ptll_IncDY_noeffSF_1_Tracking"+suffix, dipt, weight_noeffSF, AFBAnalyzer::unfold_0bjet_ptbinnum_reco,AFBAnalyzer::unfold_0bjet_ptbin_reco);
    FillHist(prefix+hprefix+"lpt_IncDY_noeffSF_1_Tracking"+suffix, lepton0->Pt(), weight_noeffSF, AFBAnalyzer::lptbinnum,AFBAnalyzer::lptbin);
    FillHist(prefix+hprefix+"leta_IncDY_noeffSF_1_Tracking"+suffix, lepton0->Eta(), weight_noeffSF, 50,-2.5,2.5);
    FillHist(prefix+hprefix+"lpt_IncDY_noeffSF_1_Tracking"+suffix, lepton1->Pt(), weight_noeffSF, AFBAnalyzer::lptbinnum,AFBAnalyzer::lptbin);
    FillHist(prefix+hprefix+"leta_IncDY_noeffSF_1_Tracking"+suffix, lepton1->Eta(), weight_noeffSF, 50,-2.5,2.5);

    weight_noeffSF *= muonRECOSF * electronRECOSF;
    FillHist(prefix+hprefix+"mll_IncDY_noeffSF_2_RECO"+suffix, dimass, weight_noeffSF, 80,70,110);
    FillHist(prefix+hprefix+"yll_IncDY_noeffSF_2_RECO"+suffix, dirap, weight_noeffSF, 50,-2.5,2.5);
    FillHist(prefix+hprefix+"ptll_IncDY_noeffSF_2_RECO"+suffix, dipt, weight_noeffSF, AFBAnalyzer::unfold_0bjet_ptbinnum_reco,AFBAnalyzer::unfold_0bjet_ptbin_reco);
    FillHist(prefix+hprefix+"lpt_IncDY_noeffSF_2_RECO"+suffix, lepton0->Pt(), weight_noeffSF, AFBAnalyzer::lptbinnum,AFBAnalyzer::lptbin);
    FillHist(prefix+hprefix+"leta_IncDY_noeffSF_2_RECO"+suffix, lepton0->Eta(), weight_noeffSF, 50,-2.5,2.5);
    FillHist(prefix+hprefix+"lpt_IncDY_noeffSF_2_RECO"+suffix, lepton1->Pt(), weight_noeffSF, AFBAnalyzer::lptbinnum,AFBAnalyzer::lptbin);
    FillHist(prefix+hprefix+"leta_IncDY_noeffSF_2_RECO"+suffix, lepton1->Eta(), weight_noeffSF, 50,-2.5,2.5);

    weight_noeffSF *= muonIDSF * electronIDSF;
    FillHist(prefix+hprefix+"mll_IncDY_noeffSF_3_ID"+suffix, dimass, weight_noeffSF, 80,70,110);
    FillHist(prefix+hprefix+"yll_IncDY_noeffSF_3_ID"+suffix, dirap, weight_noeffSF, 50,-2.5,2.5);
    FillHist(prefix+hprefix+"ptll_IncDY_noeffSF_3_ID"+suffix, dipt, weight_noeffSF, AFBAnalyzer::unfold_0bjet_ptbinnum_reco,AFBAnalyzer::unfold_0bjet_ptbin_reco);
    FillHist(prefix+hprefix+"lpt_IncDY_noeffSF_3_ID"+suffix, lepton0->Pt(), weight_noeffSF, AFBAnalyzer::lptbinnum,AFBAnalyzer::lptbin);
    FillHist(prefix+hprefix+"leta_IncDY_noeffSF_3_ID"+suffix, lepton0->Eta(), weight_noeffSF, 50,-2.5,2.5);
    FillHist(prefix+hprefix+"lpt_IncDY_noeffSF_3_ID"+suffix, lepton1->Pt(), weight_noeffSF, AFBAnalyzer::lptbinnum,AFBAnalyzer::lptbin);
    FillHist(prefix+hprefix+"leta_IncDY_noeffSF_3_ID"+suffix, lepton1->Eta(), weight_noeffSF, 50,-2.5,2.5);

    // Without lepton momentum / energy correction
    Lepton* lepton0_raw = NULL;
    Lepton* lepton1_raw = NULL;
    if(channel.Contains("mm"+GetEraShort())){
      vector<Muon> muons_raw = MuonMomentumCorrection(muons, -1,0, true);
      lepton0_raw = &muons_raw.at(0);
      lepton1_raw = &muons_raw.at(1);
    }else{
      vector<Electron> electrons_raw = ElectronEnergyCorrection(electrons, -1,0, true);
      lepton0_raw = &electrons_raw.at(0);
      lepton1_raw = &electrons_raw.at(1);
    }

    FillHist(prefix+hprefix+"mll_IncDY_noRocAep"+suffix, (*lepton0_raw + *lepton1_raw).M(), map_weight[""], 80,70,110);
    FillHist(prefix+hprefix+"yll_IncDY_noRocAep"+suffix, (*lepton0_raw + *lepton1_raw).Rapidity(), map_weight[""], 50,-2.5,2.5);
    FillHist(prefix+hprefix+"ptll_IncDY_noRocAep"+suffix, (*lepton0_raw + *lepton1_raw).Pt(), map_weight[""], AFBAnalyzer::unfold_0bjet_ptbinnum_reco,AFBAnalyzer::unfold_0bjet_ptbin_reco);
    FillHist(prefix+hprefix+"lpt_IncDY_noRocAep"+suffix, lepton0_raw->Pt(), map_weight[""], AFBAnalyzer::lptbinnum,AFBAnalyzer::lptbin);
    FillHist(prefix+hprefix+"leta_IncDY_noRocAep"+suffix, lepton0_raw->Eta(), map_weight[""], 50,-2.5,2.5);
    FillHist(prefix+hprefix+"lpt_IncDY_noRocAep"+suffix, lepton1_raw->Pt(), map_weight[""], AFBAnalyzer::lptbinnum,AFBAnalyzer::lptbin);
    FillHist(prefix+hprefix+"leta_IncDY_noRocAep"+suffix, lepton1_raw->Eta(), map_weight[""], 50,-2.5,2.5);

    // Without Zpt Correction
    FillHist(prefix+hprefix+"mll_IncDY_noZpt"+suffix, dimass, map_weight[""] / zptweight, 80,70,110);
    FillHist(prefix+hprefix+"yll_IncDY_noZpt"+suffix, dirap, map_weight[""] / zptweight, 50,-2.5,2.5);
    FillHist(prefix+hprefix+"ptll_IncDY_noZpt"+suffix, dipt, map_weight[""] / zptweight, AFBAnalyzer::unfold_0bjet_ptbinnum_reco,AFBAnalyzer::unfold_0bjet_ptbin_reco);
    FillHist(prefix+hprefix+"lpt_IncDY_noZpt"+suffix, lepton0->Pt(), map_weight[""] / zptweight, AFBAnalyzer::lptbinnum,AFBAnalyzer::lptbin);
    FillHist(prefix+hprefix+"leta_IncDY_noZpt"+suffix, lepton0->Eta(), map_weight[""] / zptweight, 50,-2.5,2.5);
    FillHist(prefix+hprefix+"lpt_IncDY_noZpt"+suffix, lepton1->Pt(), map_weight[""] / zptweight, AFBAnalyzer::lptbinnum,AFBAnalyzer::lptbin);
    FillHist(prefix+hprefix+"leta_IncDY_noZpt"+suffix, lepton1->Eta(), map_weight[""] / zptweight, 50,-2.5,2.5);

    // Without Weak Correction
    FillHist(prefix+hprefix+"mll_IncDY_noWeak"+suffix, dimass, map_weight[""] / weakweight, 80,70,110);
    FillHist(prefix+hprefix+"yll_IncDY_noWeak"+suffix, dirap, map_weight[""] / weakweight, 50,-2.5,2.5);
    FillHist(prefix+hprefix+"ptll_IncDY_noWeak"+suffix, dipt, map_weight[""] / weakweight, AFBAnalyzer::unfold_0bjet_ptbinnum_reco,AFBAnalyzer::unfold_0bjet_ptbin_reco);
    FillHist(prefix+hprefix+"lpt_IncDY_noWeak"+suffix, lepton0->Pt(), map_weight[""] / weakweight, AFBAnalyzer::lptbinnum,AFBAnalyzer::lptbin);
    FillHist(prefix+hprefix+"leta_IncDY_noWeak"+suffix, lepton0->Eta(), map_weight[""] / weakweight, 50,-2.5,2.5);
    FillHist(prefix+hprefix+"lpt_IncDY_noWeak"+suffix, lepton1->Pt(), map_weight[""] / weakweight, AFBAnalyzer::lptbinnum,AFBAnalyzer::lptbin);
    FillHist(prefix+hprefix+"leta_IncDY_noWeak"+suffix, lepton1->Eta(), map_weight[""] / weakweight, 50,-2.5,2.5);

    // Without PUweight Correction
    FillHist(prefix+hprefix+"mll_IncDY_noPU"+suffix, dimass, map_weight[""] / PUweight, 80,70,110);
    FillHist(prefix+hprefix+"yll_IncDY_noPU"+suffix, dirap, map_weight[""] / PUweight, 50,-2.5,2.5);
    FillHist(prefix+hprefix+"ptll_IncDY_noPU"+suffix, dipt, map_weight[""] / PUweight, AFBAnalyzer::unfold_0bjet_ptbinnum_reco,AFBAnalyzer::unfold_0bjet_ptbin_reco);
    FillHist(prefix+hprefix+"lpt_IncDY_noPU"+suffix, lepton0->Pt(), map_weight[""] / PUweight, AFBAnalyzer::lptbinnum,AFBAnalyzer::lptbin);
    FillHist(prefix+hprefix+"leta_IncDY_noPU"+suffix, lepton0->Eta(), map_weight[""] / PUweight, 50,-2.5,2.5);
    FillHist(prefix+hprefix+"lpt_IncDY_noPU"+suffix, lepton1->Pt(), map_weight[""] / PUweight, AFBAnalyzer::lptbinnum,AFBAnalyzer::lptbin);
    FillHist(prefix+hprefix+"leta_IncDY_noPU"+suffix, lepton1->Eta(), map_weight[""] / PUweight, 50,-2.5,2.5);
    FillHist(prefix+hprefix+"nPV_IncDY_noPU"+suffix, nPV, map_weight[""] / PUweight, 100,0,100);

    // DY+njets (for Ankita)
    if(realjets.size() == 0) return;
    TString signal = prefix+hprefix;
    if(IsDYSample){
      Gen* gen_b0 = NULL;
      Jet* jet_b0 = NULL;
      jet_b0 = &realjets.at(0);
      //int jet_idx = -1;
      double mindR = 0.4;
      //double maxpT = 0.;
      //for(unsigned int j=0; j<realjets.size(); j++){
      for(unsigned int i=4; i<gens.size(); i++){ // 0,1 - proton, 2,3 - partons
        if(!gens.at(i).isPrompt()) continue;
        if(!gens.at(i).isHardProcess()) continue;
        if(fabs(gens.at(i).PID()) != 5) continue;
        //if(gens.at(i).DeltaR(realjets.at(j)) > 0.4) continue;
        //if(realjets.at(j).Pt() > maxpT) continue;
        if(gens.at(i).DeltaR(realjets.at(0)) > mindR) continue;
        //jet_b0 = &realjets.at(j);
        //jet_idx = j;
        gen_b0 = &gens.at(i);
        mindR = gen_b0->DeltaR(*jet_b0);
      }

      bool isbG = false;
      if(fabs(gens.at(2).PID()) * gens.at(3).PID() == 5 * 21) isbG = true;
      if(gen_b0){
        if(gen_b0->PID() > 0) signal += "dyb_";
        else signal += "dybbar_";
        FillHist(signal+"gen_b0_pt_njets"+suffix, gen_b0->Pt(), map_weight[""], 200,0,200);
        FillHist(signal+"gen_b0_eta_njets"+suffix, gen_b0->Eta(), map_weight[""], 50,-2.5,2.5);
        FillHist(signal+"gen_b0_ptratio_njets"+suffix, gen_b0->Pt() / jet_b0->Pt(), map_weight[""], 200,0,10);
        FillHist(signal+"gen_b0_dR_njets"+suffix, gen_b0->DeltaR(*jet_b0), map_weight[""], 500,0,10);
        //FillHist(signal+"gen_b0_jet_idx_njets"+suffix, jet_idx, map_weight[""], 16,-1,15);

        FillHist(signal+"gen_b0_PID_njets"+suffix, gen_b0->PID(), map_weight[""], 200,-100,100);
        FillHist(signal+"gen_b0_PID_partonflavor_njets"+suffix, gen_b0->PID(), jet_b0->partonFlavour(), map_weight[""], 40,-10,30 , 40,-10,30);
        FillHist(signal+"gen_isbG_njets"+suffix, isbG? 1: 0, map_weight[""], 2,0,2);
      }else{
        FillHist(signal+"gen_b0_PID_njets"+suffix, 0, map_weight[""], 200,-100,100);
        FillHist(signal+"gen_b0_PID_partonflavor_njets"+suffix, 0, jet_b0->partonFlavour(), map_weight[""], 40,-10,30 , 40,-10,30);
        FillHist(signal+"gen_isbG_njets"+suffix, isbG? 1: 0, map_weight[""], 2,0,2);
      }
    }

    FillHist(signal+"mll_njets"+suffix, dimass, map_weight[""], 80,70,110);
    FillHist(signal+"yll_njets"+suffix, dirap, map_weight[""], 50,-2.5,2.5);
    FillHist(signal+"ptll_njets"+suffix, dipt, map_weight[""], 200,0,200);
    FillHist(signal+"lpt_njets"+suffix, lepton0->Pt(), map_weight[""], 200,0,200);
    FillHist(signal+"leta_njets"+suffix, lepton0->Eta(), map_weight[""], 50,-2.5,2.5);
    FillHist(signal+"l0pt_njets"+suffix, lepton0->Pt(), map_weight[""], 200,0,200);
    FillHist(signal+"l0eta_njets"+suffix, lepton0->Eta(), map_weight[""], 50,-2.5,2.5);
    FillHist(signal+"lpt_njets"+suffix, lepton1->Pt(), map_weight[""], 200,0,200);
    FillHist(signal+"leta_njets"+suffix, lepton1->Eta(), map_weight[""], 50,-2.5,2.5);
    FillHist(signal+"l1pt_njets"+suffix, lepton1->Pt(), map_weight[""], 200,0,200);
    FillHist(signal+"l1eta_njets"+suffix, lepton1->Eta(), map_weight[""], 50,-2.5,2.5);
    FillHist(signal+"j0pt_njets"+suffix, realjets.at(0).Pt(), map_weight[""], 200,0,200);
    FillHist(signal+"j0eta_njets"+suffix, realjets.at(0).Eta(), map_weight[""], 50,-2.5,2.5);
    FillHist(signal+"j0phi_njets"+suffix, realjets.at(0).Phi(), map_weight[""], 64,-3.2,3.2);
    FillHist(signal+"met_njets"+suffix, PuppiMET_Type1_pt, map_weight[""], 200,0,200);
    FillHist(signal+"metphi_njets"+suffix, PuppiMET_Type1_phi, map_weight[""], 64,-3.2,3.2);
    FillHist(signal+"met_metphi_njets"+suffix, PuppiMET_Type1_pt, PuppiMET_Type1_phi, map_weight[""], 40,0,200, 32,-3.2,3.2);
  }

  if(bjets.size() == 0) return;
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "Tightnb", map_weight[""]);
  //bjets.at(0) *= bjets.at(0).BJetNNCorrection(); // bJetEnergyCorrection (BjetRegression)?
  jet0 = &bjets.at(0);
  bcharge = jetCharge(*jet0);
  bchargeSF = GetbChargeSFWeight(bjets, 1, 0);
  double costhetaRecoil = GetCosThetaRecoil(lepton0, lepton1, jet0);
  double costhetaRecoil_nobch = GetCosThetaRecoil(lepton0, lepton1, jet0, -1);

  if(IsDYSample){
    vector<Gen> hs_gens = {};
    Gen* gen_b0 = NULL;
    double mindR = 0.4;
    for(unsigned int i=4; i<gens.size(); i++){ // 0,1 - proton, 2,3 - partons
      if(!gens.at(i).isPrompt()) continue;
      if(!gens.at(i).isHardProcess()) continue;
      hs_gens.push_back(gens.at(i));
      if(gens.at(i).DeltaR(*jet0) > mindR) continue;
      gen_b0 = &gens.at(i);
      mindR = gen_b0->DeltaR(*jet0);
    }

    bool isbG = false;
    if(fabs(gens.at(2).PID()) * gens.at(3).PID() == 5 * 21) isbG = true;
    unsigned int nhsb = 0;
    int PID_sum = 0;
    for(unsigned int i=0; i<hs_gens.size(); i++){
      if(fabs(hs_gens.at(i).PID()) == 5){
        if(IsNominalRun) FillHist(prefix+hprefix+"gen_hs_i"+suffix, i, map_weight[""], 10,0,10);
        if(IsNominalRun) FillHist(prefix+hprefix+"gen_hs_index"+suffix, hs_gens.at(i).Index(), map_weight[""], 10,0,10);
        nhsb++;
        PID_sum += hs_gens.at(i).PID();
      }
    }
    if(IsNominalRun){
      FillHist(prefix+hprefix+"gen_nhsb"+suffix, nhsb, map_weight[""], 10,0,10);
      FillHist(prefix+hprefix+"gen_hsb_PIDsum"+suffix, PID_sum, map_weight[""], 40,-20,20);
      if(isbG) FillHist(prefix+hprefix+"gen_nhsb_isbG"+suffix, nhsb, map_weight[""], 10,0,10);
      if(isbG) FillHist(prefix+hprefix+"gen_hsb_PIDsum_isbG"+suffix, PID_sum, map_weight[""], 40,-20,20);
    }

    if(gen_b0){
      if(gen_b0->PID() == 5) prefix += "dyb_";
      else if(gen_b0->PID() == -5) prefix += "dybbar_";
      else if(gen_b0->PID() == 4) prefix += "dyc_";
      else if(gen_b0->PID() == -4) prefix += "dycbar_";

      if(IsNominalRun){
        FillHist(prefix+hprefix+"gen_b0_pt"+suffix, gen_b0->Pt(), map_weight[""], 200,0,200);
        FillHist(prefix+hprefix+"gen_b0_eta"+suffix, gen_b0->Eta(), map_weight[""], 50,-2.5,2.5);
        FillHist(prefix+hprefix+"gen_b0_ptratio"+suffix, gen_b0->Pt() / jet0->Pt(), map_weight[""], 200,0,10);
        FillHist(prefix+hprefix+"gen_b0_dR"+suffix, gen_b0->DeltaR(*jet0), map_weight[""], 500,0,10);

        FillHist(prefix+hprefix+"gen_b0_PID"+suffix, gen_b0->PID(), map_weight[""], 200,-100,100);
        FillHist(prefix+hprefix+"gen_b0_PID_nhsb"+suffix, gen_b0->PID(), nhsb, map_weight[""], 40,-10,30 , 10,0,10);
        FillHist(prefix+hprefix+"gen_b0_PID_partonflavor"+suffix, gen_b0->PID(), jet0->partonFlavour(), map_weight[""], 40,-10,30 , 40,-10,30);
        FillHist(prefix+hprefix+"gen_b0_PID_hadronflavor"+suffix, gen_b0->PID(), jet0->hadronFlavour(), map_weight[""], 40,-10,30 , 40,-10,30);
        FillHist(prefix+hprefix+"gen_b0_PID_genHFflavor"+suffix, gen_b0->PID(), jet0->GenHFHadronMatcherFlavour(), map_weight[""], 40,-10,30 , 40,-10,30);
        FillHist(prefix+hprefix+"gen_b0_PID_genHForigin"+suffix, gen_b0->PID(), jet0->GenHFHadronMatcherOrigin(), map_weight[""], 40,-10,30 , 40,-10,30);
        FillHist(prefix+hprefix+"gen_nhsb_Matching"+suffix, nhsb, map_weight[""], 10,0,10);
        FillHist(prefix+hprefix+"gen_hsb_PIDsum_Matching"+suffix, PID_sum, map_weight[""], 40,-20,20);
        FillHist(prefix+hprefix+"gen_isbG_Matching"+suffix, isbG? 1: 0, map_weight[""], 2,0,2);
      }else{
        FillHist(prefix+hprefix+"gen_b0_PID"+suffix, 0, map_weight[""], 200,-100,100);
        FillHist(prefix+hprefix+"gen_b0_PID_nhsb"+suffix, 0, nhsb, map_weight[""], 40,-10,30 , 10,0,10);
        FillHist(prefix+hprefix+"gen_b0_PID_partonflavor"+suffix, 0, jet0->partonFlavour(), map_weight[""], 40,-10,30 , 40,-10,30);
        FillHist(prefix+hprefix+"gen_b0_PID_hadronflavor"+suffix, 0, jet0->hadronFlavour(), map_weight[""], 40,-10,30 , 40,-10,30);
        FillHist(prefix+hprefix+"gen_b0_PID_genHFflavor"+suffix, 0, jet0->GenHFHadronMatcherFlavour(), map_weight[""], 40,-10,30 , 40,-10,30);
        FillHist(prefix+hprefix+"gen_b0_PID_genHForigin"+suffix, 0, jet0->GenHFHadronMatcherOrigin(), map_weight[""], 40,-10,30 , 40,-10,30);
        FillHist(prefix+hprefix+"gen_nhsb_noMatching"+suffix, nhsb, map_weight[""], 10,0,10);
        FillHist(prefix+hprefix+"gen_hsb_PIDsum_noMatching"+suffix, PID_sum, map_weight[""], 40,-20,20);
        FillHist(prefix+hprefix+"gen_isbG_noMatching"+suffix, isbG? 1: 0, map_weight[""], 2,0,2);
      }
    }
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

  map_weight[""] *= bchargeSF;
  if(IsNominalLike){
    FillHist(prefix+hprefix+"weight_bChargeSF1"+suffix, bchargeSF, map_weight[""], 200,-5,5);
    FillCutflow(prefix+hprefix+"cutflow"+suffix, "bChargeSF1", map_weight[""]);
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
  }else if(!IsDATA && HasFlag("PDFSYS") && option == ""){
    // Zpt Reweight
    //map_weight["_noZpt"] =  map_weight[""] / zptweight;
    //map_weight["_Zpt_gym"] =  map_weight[""] / zptweight * zptweight_gym;
    map_weight["_Zpt"] =  map_weight[""] * zptweight;
    map_weight["_Zpt_gym"] =  map_weight[""] * zptweight_gym;

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

  // Electron charge flip SF
  if(!IsDATA && IsNominalRun){
    map_weight["_noCFSF"] = map_weight[""] / chargeflipSF;
    map_weight["_CFSF_up"] = map_weight[""] / chargeflipSF * GetCFSF(1);
    map_weight["_CFSF_down"] = map_weight[""] / chargeflipSF * GetCFSF(-1);
  }
  if((HasFlag("SYS") || HasFlag("PDFSYS") || HasFlag("LEPSYS")) && option == "") map_weight.erase("");

  // For comparison with Hyonsan's results (8 bins)
  if(jet0->Pt() > 40){
    FillHist(prefix+hprefix+"mll_nbjets"+suffix, dimass, map_weight, 80,70,110);
    FillHist(prefix+hprefix+"yll_nbjets"+suffix, dirap, map_weight, 50,-2.5,2.5);
    FillHist(prefix+hprefix+"ptll_nbjets"+suffix, dipt, map_weight, 200,0,200);
    FillHist(prefix+hprefix+"lpt_nbjets"+suffix, lepton0->Pt(), map_weight, 200,0,200);
    FillHist(prefix+hprefix+"leta_nbjets"+suffix, lepton0->Eta(), map_weight, 50,-2.5,2.5);
    FillHist(prefix+hprefix+"l0pt_nbjets"+suffix, lepton0->Pt(), map_weight, 200,0,200);
    FillHist(prefix+hprefix+"l0eta_nbjets"+suffix, lepton0->Eta(), map_weight, 50,-2.5,2.5);
    FillHist(prefix+hprefix+"lpt_nbjets"+suffix, lepton1->Pt(), map_weight, 200,0,200);
    FillHist(prefix+hprefix+"leta_nbjets"+suffix, lepton1->Eta(), map_weight, 50,-2.5,2.5);
    FillHist(prefix+hprefix+"l1pt_nbjets"+suffix, lepton1->Pt(), map_weight, 200,0,200);
    FillHist(prefix+hprefix+"l1eta_nbjets"+suffix, lepton1->Eta(), map_weight, 50,-2.5,2.5);
    FillHist(prefix+hprefix+"bpt_nbjets"+suffix, jet0->Pt(), map_weight, 200,0,200);
    FillHist(prefix+hprefix+"beta_nbjets"+suffix, jet0->Eta(), map_weight, 50,-2.5,2.5);
    FillHist(prefix+hprefix+"bphi_nbjets"+suffix, jet0->Phi(), map_weight, 64,-3.2,3.2);
    FillHist(prefix+hprefix+"bChargeRaw_nbjets"+suffix, bcharge, map_weight, 200,-5,5);
    FillHist(prefix+hprefix+"bChargeRaw2_nbjets"+suffix, bcharge, map_weight, 100,-5,5);
    FillHist(prefix+hprefix+"bCharge_nbjets"+suffix, (bcharge < 0? -0.5: 0.5), map_weight, 2,-1,1);
    FillHist(prefix+hprefix+"met_nbjets"+suffix, PuppiMET_Type1_pt, map_weight, 200,0,200);
    FillHist(prefix+hprefix+"metphi_nbjets"+suffix, PuppiMET_Type1_phi, map_weight, 64,-3.2,3.2);
    FillHist(prefix+hprefix+"met_metphi_nbjets"+suffix, PuppiMET_Type1_pt, PuppiMET_Type1_phi, map_weight[""], 40,0,200, 32,-3.2,3.2);
    FillHist(prefix+hprefix+"ZbdPhi_nbjets"+suffix, abs((*lepton0 + *lepton1).DeltaPhi(*jet0)), map_weight, 200,0,10);
    FillHist(prefix+hprefix+"Zbpt_nbjets"+suffix, (*lepton0 + *lepton1 + *jet0).Pt(), map_weight, 200,0,200);
    FillHist(prefix+hprefix+"costhetaRecoil_nbjets"+suffix, dimass, fabs(bcharge), jet0->Pt(), costhetaRecoil, map_weight, afb_mbinnum,(double*)afb_mbin, afb_chbinnum,(double*)afb_chbin, afb_ptbinnum,(double*)afb_ptbin, 20,-1,1);
    FillHist(prefix+hprefix+"costhetaRecoil2_nbjets"+suffix, dimass, fabs(bcharge), jet0->Pt(), costhetaRecoil, map_weight, afb_mbinnum_HS,(double*)afb_mbin_HS, afb_chbinnum,(double*)afb_chbin, afb_ptbinnum,(double*)afb_ptbin, 20,-1,1);
  }

  if(bjets.size() != 1) return;
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "Tight1b", map_weight[""]);

  FillHist(prefix+hprefix+"mll_Tight1b"+suffix, dimass, map_weight, 80,70,110);
  FillHist(prefix+hprefix+"yll_Tight1b"+suffix, dirap, map_weight, 50,-2.5,2.5);
  FillHist(prefix+hprefix+"ptll_Tight1b"+suffix, dipt, map_weight, 200,0,200);
  FillHist(prefix+hprefix+"lpt_Tight1b"+suffix, lepton0->Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"leta_Tight1b"+suffix, lepton0->Eta(), map_weight, 50,-2.5,2.5);
  FillHist(prefix+hprefix+"l0pt_Tight1b"+suffix, lepton0->Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"l0eta_Tight1b"+suffix, lepton0->Eta(), map_weight, 50,-2.5,2.5);
  FillHist(prefix+hprefix+"lpt_Tight1b"+suffix, lepton1->Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"leta_Tight1b"+suffix, lepton1->Eta(), map_weight, 50,-2.5,2.5);
  FillHist(prefix+hprefix+"l1pt_Tight1b"+suffix, lepton1->Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"l1eta_Tight1b"+suffix, lepton1->Eta(), map_weight, 50,-2.5,2.5);
  FillHist(prefix+hprefix+"bpt_Tight1b"+suffix, jet0->Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"beta_Tight1b"+suffix, jet0->Eta(), map_weight, 50,-2.5,2.5);
  FillHist(prefix+hprefix+"bphi_Tight1b"+suffix, jet0->Phi(), map_weight, 64,-3.2,3.2);
  FillHist(prefix+hprefix+"bChargeRaw_Tight1b"+suffix, bcharge, map_weight, 200,-5,5);
  FillHist(prefix+hprefix+"bChargeRaw2_Tight1b"+suffix, bcharge, map_weight, 100,-5,5);
  FillHist(prefix+hprefix+"bCharge_Tight1b"+suffix, (bcharge < 0? -0.5: 0.5), map_weight, 2,-1,1);
  FillHist(prefix+hprefix+"met_Tight1b"+suffix, PuppiMET_Type1_pt, map_weight, 200,0,200);
  FillHist(prefix+hprefix+"metphi_Tight1b"+suffix, PuppiMET_Type1_phi, map_weight, 64,-3.2,3.2);
  FillHist(prefix+hprefix+"met_metphi_Tight1b"+suffix, PuppiMET_Type1_pt, PuppiMET_Type1_phi, map_weight[""], 40,0,200, 32,-3.2,3.2);
  FillHist(prefix+hprefix+"ZbdPhi_Tight1b"+suffix, abs((*lepton0 + *lepton1).DeltaPhi(*jet0)), map_weight, 200,0,10);
  FillHist(prefix+hprefix+"Zbpt_Tight1b"+suffix, (*lepton0 + *lepton1 + *jet0).Pt(), map_weight, 200,0,200);
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
  FillHist(prefix+hprefix+"yll_Veto2b"+suffix, dirap, map_weight, 50,-2.5,2.5);
  FillHist(prefix+hprefix+"ptll_Veto2b"+suffix, dipt, map_weight, 200,0,200);
  FillHist(prefix+hprefix+"lpt_Veto2b"+suffix, lepton0->Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"leta_Veto2b"+suffix, lepton0->Eta(), map_weight, 50,-2.5,2.5);
  FillHist(prefix+hprefix+"l0pt_Veto2b"+suffix, lepton0->Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"l0eta_Veto2b"+suffix, lepton0->Eta(), map_weight, 50,-2.5,2.5);
  FillHist(prefix+hprefix+"lpt_Veto2b"+suffix, lepton1->Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"leta_Veto2b"+suffix, lepton1->Eta(), map_weight, 50,-2.5,2.5);
  FillHist(prefix+hprefix+"l1pt_Veto2b"+suffix, lepton1->Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"l1eta_Veto2b"+suffix, lepton1->Eta(), map_weight, 50,-2.5,2.5);
  FillHist(prefix+hprefix+"bpt_Veto2b"+suffix, jet0->Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"beta_Veto2b"+suffix, jet0->Eta(), map_weight, 50,-2.5,2.5);
  FillHist(prefix+hprefix+"bphi_Veto2b"+suffix, jet0->Phi(), map_weight, 64,-3.2,3.2);
  FillHist(prefix+hprefix+"bChargeRaw_Veto2b"+suffix, bcharge, map_weight, 200,-5,5);
  FillHist(prefix+hprefix+"bChargeRaw2_Veto2b"+suffix, bcharge, map_weight, 100,-5,5);
  FillHist(prefix+hprefix+"bCharge_Veto2b"+suffix, (bcharge < 0? -0.5: 0.5), map_weight, 2,-1,1);
  FillHist(prefix+hprefix+"met_Veto2b"+suffix, PuppiMET_Type1_pt, map_weight, 200,0,200);
  FillHist(prefix+hprefix+"metphi_Veto2b"+suffix, PuppiMET_Type1_phi, map_weight, 64,-3.2,3.2);
  FillHist(prefix+hprefix+"met_metphi_Veto2b"+suffix, PuppiMET_Type1_pt, PuppiMET_Type1_phi, map_weight[""], 40,0,200, 32,-3.2,3.2);
  FillHist(prefix+hprefix+"ZbdPhi_Veto2b"+suffix, abs((*lepton0 + *lepton1).DeltaPhi(*jet0)), map_weight, 200,0,10);
  FillHist(prefix+hprefix+"Zbpt_Veto2b"+suffix, (*lepton0 + *lepton1 + *jet0).Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"costhetaRecoil_Veto2b"+suffix, dimass, fabs(bcharge), jet0->Pt(), costhetaRecoil, map_weight, afb_mbinnum,(double*)afb_mbin, afb_chbinnum,(double*)afb_chbin, afb_ptbinnum,(double*)afb_ptbin, 20,-1,1);

  if(PuppiMET_Type1_pt > 60) return;
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "MET60", map_weight[""]);

  FillHist(prefix+hprefix+"mll_met60"+suffix, dimass, map_weight, 80,70,110);
  FillHist(prefix+hprefix+"yll_met60"+suffix, dirap, map_weight, 50,-2.5,2.5);
  FillHist(prefix+hprefix+"ptll_met60"+suffix, dipt, map_weight, 200,0,200);
  FillHist(prefix+hprefix+"lpt_met60"+suffix, lepton0->Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"leta_met60"+suffix, lepton0->Eta(), map_weight, 50,-2.5,2.5);
  FillHist(prefix+hprefix+"l0pt_met60"+suffix, lepton0->Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"l0eta_met60"+suffix, lepton0->Eta(), map_weight, 50,-2.5,2.5);
  FillHist(prefix+hprefix+"lpt_met60"+suffix, lepton1->Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"leta_met60"+suffix, lepton1->Eta(), map_weight, 50,-2.5,2.5);
  FillHist(prefix+hprefix+"l1pt_met60"+suffix, lepton1->Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"l1eta_met60"+suffix, lepton1->Eta(), map_weight, 50,-2.5,2.5);
  FillHist(prefix+hprefix+"bpt_met60"+suffix, jet0->Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"beta_met60"+suffix, jet0->Eta(), map_weight, 50,-2.5,2.5);
  FillHist(prefix+hprefix+"bphi_met60"+suffix, jet0->Phi(), map_weight, 64,-3.2,3.2);
  FillHist(prefix+hprefix+"bChargeRaw_met60"+suffix, bcharge, map_weight, 200,-5,5);
  FillHist(prefix+hprefix+"bChargeRaw2_met60"+suffix, bcharge, map_weight, 100,-5,5);
  FillHist(prefix+hprefix+"bCharge_met60"+suffix, (bcharge < 0? -0.5: 0.5), map_weight, 2,-1,1);
  FillHist(prefix+hprefix+"met_met60"+suffix, PuppiMET_Type1_pt, map_weight, 200,0,200);
  FillHist(prefix+hprefix+"metphi_met60"+suffix, PuppiMET_Type1_phi, map_weight, 64,-3.2,3.2);
  FillHist(prefix+hprefix+"met_metphi_met60"+suffix, PuppiMET_Type1_pt, PuppiMET_Type1_phi, map_weight[""], 40,0,200, 32,-3.2,3.2);
  FillHist(prefix+hprefix+"ZbdPhi_met60"+suffix, abs((*lepton0 + *lepton1).DeltaPhi(*jet0)), map_weight, 200,0,10);
  FillHist(prefix+hprefix+"Zbpt_met60"+suffix, (*lepton0 + *lepton1 + *jet0).Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"costhetaRecoil_met60"+suffix, dimass, fabs(bcharge), jet0->Pt(), costhetaRecoil, map_weight, afb_mbinnum,(double*)afb_mbin, afb_chbinnum,(double*)afb_chbin, afb_ptbinnum,(double*)afb_ptbin, 20,-1,1);

  if((*lepton0 + *lepton1).Pt() < 15) return;
  if(IsNominalLike) FillCutflow(prefix+hprefix+"cutflow"+suffix, "ZpT15", map_weight[""]);

  //if(option.Contains("jet_scale_up")) bjets = ScaleJets(bjets, 1);
  //else if(option.Contains("jet_scale_down")) bjets = ScaleJets(bjets, -1);
  //else if(option.Contains("jet_smear_up")) bjets = SmearJets(bjets, 1);
  //else if(option.Contains("jet_smear_down")) bjets = SmearJets(bjets, -1);
  //jet0 = &bjets.at(0);
  //bcharge = jetCharge(*jet0);

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
  FillHist(prefix+hprefix+"met_metphi"+suffix, PuppiMET_Type1_pt, PuppiMET_Type1_phi, map_weight[""], 40,0,200, 32,-3.2,3.2);
  FillHist(prefix+hprefix+"met2"+suffix, PuppiMET_Type1_PhiCor_pt, map_weight, 200,0,200);
  FillHist(prefix+hprefix+"met2phi"+suffix, PuppiMET_Type1_PhiCor_phi, map_weight, 64,-3.2,3.2);
  FillHist(prefix+hprefix+"met_met2"+suffix, PuppiMET_Type1_pt - PuppiMET_Type1_PhiCor_pt, map_weight, 200,-100,100);
  FillHist(prefix+hprefix+"met2_metphi2"+suffix, PuppiMET_Type1_PhiCor_pt, PuppiMET_Type1_PhiCor_phi, map_weight[""], 40,0,200, 32,-3.2,3.2);
  FillHist(prefix+hprefix+"ZbdPhi"+suffix, abs((*lepton0 + *lepton1).DeltaPhi(*jet0)), map_weight, 200,0,10);
  FillHist(prefix+hprefix+"Zbpt"+suffix, (*lepton0 + *lepton1 + *jet0).Pt(), map_weight, 200,0,200);
  FillHist(prefix+hprefix+"costhetaRecoil"+suffix, dimass, fabs(bcharge), jet0->Pt(), costhetaRecoil, map_weight, afb_mbinnum,(double*)afb_mbin, afb_chbinnum,(double*)afb_chbin, afb_ptbinnum,(double*)afb_ptbin, 20,-1,1);
  FillHist(prefix+hprefix+"costhetaRecoil2"+suffix, dimass, fabs(bcharge), jet0->Pt(), costhetaRecoil, map_weight, afb_mbinnum_HS,(double*)afb_mbin_HS, afb_chbinnum,(double*)afb_chbin, afb_ptbinnum,(double*)afb_ptbin, 20,-1,1);
  FillHist(prefix+hprefix+"costhetaRecoil3"+suffix, dimass, fabs(bcharge), jet0->Pt(), costhetaRecoil, map_weight, afb_mbinnum_original,(double*)afb_mbin_original, afb_chbinnum,(double*)afb_chbin, afb_ptbinnum,(double*)afb_ptbin, 20,-1,1);
  FillHist(prefix+hprefix+"AbscosthetaRecoil"+suffix, dimass, fabs(bcharge), jet0->Pt(), fabs(costhetaRecoil), map_weight, afb_mbinnum,(double*)afb_mbin, afb_chbinnum,(double*)afb_chbin, afb_ptbinnum,(double*)afb_ptbin, 10,0,1);
  FillHist(prefix+hprefix+"costhetaRecoil_nobch"+suffix, dimass, fabs(bcharge), jet0->Pt(), costhetaRecoil_nobch, map_weight, afb_mbinnum,(double*)afb_mbin, afb_chbinnum,(double*)afb_chbin, afb_ptbinnum,(double*)afb_ptbin, 10,0,1);
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
  if(bmuon.size() > 0){
    jetCharge += 2 * bmuon.at(0).Charge();
    if(IsNominalRun){
      FillHist("softmuon/nmu", bmuon.size(), 1, 10,0,10);
      FillHist("softmuon/nel", belectron.size(), 1, 10,0,10);
      FillHist("softmuon/pt", bmuon.at(0).Pt(), 1, 200,0,200);
      FillHist("softmuon/ptRel", bmuon.at(0).P() * sin(bmuon.at(0).Angle(jet.Vect())), 1, 400,0,10);
      FillHist("softmuon/ptRatio", bmuon.at(0).Pt() / jet.Pt(), 1, 200,0,5);
      FillHist("softmuon/jetpt_pt", jet.Pt(), bmuon.at(0).Pt(), 1, 200,0,1000, 80,0,400);
      if(bmuon.at(0).Charge() * jet.Charge() > 0){
        FillHist("softmuon/nmu_SS", bmuon.size(), 1, 10,0,10);
        FillHist("softmuon/nel_SS", belectron.size(), 1, 10,0,10);
        FillHist("softmuon/pt_SS", bmuon.at(0).Pt(), 1, 200,0,200);
        FillHist("softmuon/ptRel_SS", bmuon.at(0).P() * sin(bmuon.at(0).Angle(jet.Vect())), 1, 400,0,10);
        FillHist("softmuon/ptRatio_SS", bmuon.at(0).Pt() / jet.Pt(), 1, 200,0,5);
        FillHist("softmuon/jetpt_pt_SS", jet.Pt(), bmuon.at(0).Pt(), 1, 200,0,1000, 80,0,400);
      }else{
        FillHist("softmuon/nmu_OS", bmuon.size(), 1, 10,0,10);
        FillHist("softmuon/nel_OS", belectron.size(), 1, 10,0,10);
        FillHist("softmuon/pt_OS", bmuon.at(0).Pt(), 1, 200,0,200);
        FillHist("softmuon/ptRel_OS", bmuon.at(0).P() * sin(bmuon.at(0).Angle(jet.Vect())), 1, 400,0,10);
        FillHist("softmuon/ptRatio_OS", bmuon.at(0).Pt() / jet.Pt(), 1, 200,0,5);
        FillHist("softmuon/jetpt_pt_OS", jet.Pt(), bmuon.at(0).Pt(), 1, 200,0,1000, 80,0,400);
      }
    }
  }else if(belectron.size() > 0){
    jetCharge += 4 * belectron.at(0).Charge();
    if(IsNominalRun){
      FillHist("softelectron/nmu", bmuon.size(), 1, 10,0,10);
      FillHist("softelectron/nel", belectron.size(), 1, 10,0,10);
      FillHist("softelectron/pt", belectron.at(0).Pt(), 1, 200,0,200);
      FillHist("softelectron/ptRel", belectron.at(0).P() * sin(belectron.at(0).Angle(jet.Vect())), 1, 400,0,10);
      FillHist("softelectron/ptRatio", belectron.at(0).Pt() / jet.Pt(), 1, 200,0,5);
      FillHist("softelectron/jetpt_pt", jet.Pt(), belectron.at(0).Pt(), 1, 200,0,1000, 80,0,400);
      if(belectron.at(0).Charge() * jet.Charge() > 0){
        FillHist("softelectron/nmu_SS", belectron.size(), 1, 10,0,10);
        FillHist("softelectron/nel_SS", belectron.size(), 1, 10,0,10);
        FillHist("softelectron/pt_SS", belectron.at(0).Pt(), 1, 200,0,200);
        FillHist("softelectron/ptRel_SS", belectron.at(0).P() * sin(belectron.at(0).Angle(jet.Vect())), 1, 400,0,10);
        FillHist("softelectron/ptRatio_SS", belectron.at(0).Pt() / jet.Pt(), 1, 200,0,5);
        FillHist("softelectron/jetpt_pt_SS", jet.Pt(), belectron.at(0).Pt(), 1, 200,0,1000, 80,0,400);
      }else{
        FillHist("softelectron/nmu_OS", belectron.size(), 1, 10,0,10);
        FillHist("softelectron/nel_OS", belectron.size(), 1, 10,0,10);
        FillHist("softelectron/pt_OS", belectron.at(0).Pt(), 1, 200,0,200);
        FillHist("softelectron/ptRel_OS", belectron.at(0).P() * sin(belectron.at(0).Angle(jet.Vect())), 1, 400,0,10);
        FillHist("softelectron/ptRatio_OS", belectron.at(0).Pt() / jet.Pt(), 1, 200,0,5);
        FillHist("softelectron/jetpt_pt_OS", jet.Pt(), belectron.at(0).Pt(), 1, 200,0,1000, 80,0,400);
      }
    }
  }
  return jetCharge;
}

dybAnalyzer::dybAnalyzer(){}
dybAnalyzer::~dybAnalyzer(){}

// From Hyonsan's functions in SMPAnalyzerCore
void dybAnalyzer::executeEventGen(){
  gprefix = "";
  if(IsData) return;

  gens = GetGens();
  if(IsDYSample){
    // LHE Setting
    lhes = GetLHEs();
    LHE lhe_l0, lhe_l1, lhe_p0, lhe_p1, lhe_j0;
    double lhe_j0_pt = 0.1;
    for(int i=0; i<(int)lhes.size(); i++){
      if(lhe_l0.ID() == 0 && (abs(lhes[i].ID()) == 11 || abs(lhes[i].ID()) == 13 || abs(lhes[i].ID()) == 15)) lhe_l0 = lhes[i];
      if(lhe_l0.ID()      && (abs(lhes[i].ID()) == 11 || abs(lhes[i].ID()) == 13 || abs(lhes[i].ID()) == 15)) lhe_l1 = lhes[i];

      if(lhe_p0.ID() == 0 && ((0 < abs(lhes[i].ID()) && abs(lhes[i].ID()) < 7) || lhes[i].ID() == 21)) lhe_p0 = lhes[i];
      if(lhe_p0.ID()      && ((0 < abs(lhes[i].ID()) && abs(lhes[i].ID()) < 7) || lhes[i].ID() == 21)) lhe_p1 = lhes[i];
      if(lhe_p0.ID() && lhe_p1.ID()){
        if(lhe_p0.ID() > lhe_p1.ID()){
          LHE temp = lhe_p0;
          lhe_p0 = lhe_p1;
          lhe_p1 = temp;
        }
        if(lhe_p1.ID() == 21 && (0 < abs(lhe_p0.ID()) && abs(lhe_p0.ID()) < 7) && lhes[i].ID() == lhe_p0.ID() && lhes[i].Pt() > lhe_j0_pt){
          lhe_j0 = lhes[i];
          lhe_j0_pt = lhes[i].Pt();
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

    if(IsNominalRun){
      if(lhe_j0.ID() != 0){ // qG collisions
        bcharge = lhe_j0.ID() < 0? 1: -1;
        double dimass_lhe = ((Particle)lhe_l0 + (Particle)lhe_l1).M();
        double dirap_lhe = ((Particle)lhe_l0 + (Particle)lhe_l1).Rapidity();
        double dipt_lhe = ((Particle)lhe_l0 + (Particle)lhe_l1).Pt();
        double costhetaRecoil_lhe = GetCosThetaRecoil((Particle*)&lhe_l0, (Particle*)&lhe_l1, (Particle*)&lhe_j0);
	double genweight = reductionweight * MCweight() * _event.GetTriggerLumi("Full");
        FillHist(Form("lhe/costhetaRecoil_%dG", abs(lhe_j0.ID())), dimass_lhe, bcharge, dipt_lhe, costhetaRecoil_lhe, genweight, afb_mbinnum,(double*)afb_mbin, afb_chbinnum,(double*)afb_chbin, afb_ptbinnum,(double*)afb_ptbin, 20,-1,1);
        TString lhe_prefix = "";
        if(lhe_j0.ID() == 5) lhe_prefix = "bG_";
        else if(lhe_j0.ID() == -5) lhe_prefix = "BG_";
        else if(lhe_j0.ID() == 4) lhe_prefix = "cG_";
        else if(lhe_j0.ID() == -4) lhe_prefix = "CG_";
        else if(lhe_j0.ID() == 3) lhe_prefix = "sG_";
        else if(lhe_j0.ID() == -3) lhe_prefix = "SG_";
        else if(lhe_j0.ID() == 2) lhe_prefix = "uG_";
        else if(lhe_j0.ID() == -2) lhe_prefix = "UG_";
        else if(lhe_j0.ID() == 1) lhe_prefix = "dG_";
        else if(lhe_j0.ID() == -1) lhe_prefix = "DG_";
        FillHist("lhe/"+lhe_prefix+"costhetaRecoil_qG", dimass_lhe, bcharge, dipt_lhe, costhetaRecoil_lhe, genweight, afb_mbinnum,(double*)afb_mbin, afb_chbinnum,(double*)afb_chbin, afb_ptbinnum,(double*)afb_ptbin, 20,-1,1);
	FillHist("lhe/"+lhe_prefix+"mll_qG", dimass_lhe, genweight, 200,40,140);
        FillHist("lhe/"+lhe_prefix+"yll_qG", dirap_lhe, genweight, 100,-5,5);
        FillHist("lhe/"+lhe_prefix+"ptll_qG", dipt_lhe, genweight, 200,0,200);
        FillHist("lhe/"+lhe_prefix+"lpt_qG", lhe_l0.Pt(), genweight, 200,0,200);
        FillHist("lhe/"+lhe_prefix+"lpt_qG", lhe_l1.Pt(), genweight, 200,0,200);
        FillHist("lhe/"+lhe_prefix+"leta_qG", lhe_l0.Eta(), genweight, 100,-5,5);
        FillHist("lhe/"+lhe_prefix+"leta_qG", lhe_l1.Eta(), genweight, 100,-5,5);
        FillHist("lhe/"+lhe_prefix+"jpt_qG", lhe_j0.Pt(), genweight, 200,0,200);
        FillHist("lhe/"+lhe_prefix+"jeta_qG", lhe_j0.Eta(), genweight, 100,-5,5);
      }
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
