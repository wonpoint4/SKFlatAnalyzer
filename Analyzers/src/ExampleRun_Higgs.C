#include "ExampleRun_Higgs.h"

void ExampleRun_Higgs::initializeAnalyzer(){

}

void ExampleRun_Higgs::executeEvent(){

  AnalyzerParameter param;

  muons = GetMuons("POGTightWithTightIso", 5, 2.4);
  electrons = GetElectrons("passTightID", 5, 2.5);
  photons = GetPhotons("passTightID", 5, 2.5);

  executeEventFromParameter(param);
}

void ExampleRun_Higgs::executeEventFromParameter(AnalyzerParameter param){

  if(!PassMETFilter()) return;

  Event ev = GetEvent();
  bool trig_DoubleMuon = ev.PassTrigger({
    "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v",
    "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",
    "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
    "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
    "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",
    "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"});
  bool trig_EGamma = ev.PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
  bool trig_MuonEG = ev.PassTrigger({
    "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",
    "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
    "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"
    "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"});
  if(DataYear==2017){
    trig_DoubleMuon = ev.PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v");
    trig_EGamma = ev.PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v");
    trig_MuonEG = ev.PassTrigger({
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"});
  }else if(DataYear==2018){
    trig_DoubleMuon = ev.PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v");
    trig_EGamma = ev.PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v");
    trig_MuonEG = ev.PassTrigger({
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"});
  }
  //bool trig_DiPhoton   = ev.PassTrigger({"HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55_v", "HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95_v"});
  bool trig_DoublePhoton = ev.PassTrigger({"HLT_DoublePhoton60_v", "HLT_DoublePhoton70_v", "HLT_DoublePhoton85_v"});

  double weight = 1.;
  if(!IsDATA){
    weight *= MCweight();
    weight *= ev.GetTriggerLumi("Full");
    weight *= GetPrefireWeight(0);
  }

  if(photons.size() >= 2 ){
    if(photons.at(0).Pt() > 30 && photons.at(1).Pt() > 20) FillHist(GetEraShort()+"/2r_HCand_Mass", (photons.at(0)+photons.at(1)).M(), weight, 240, 50., 290.);
    if(photons.at(0).Pt() > 50 && photons.at(1).Pt() > 50) FillHist(GetEraShort()+"/2rHard_HCand_Mass", (photons.at(0)+photons.at(1)).M(), weight, 240, 50., 290.);
    if(photons.at(0).Pt() > 50 && photons.at(1).Pt() > 50 && trig_DoublePhoton) FillHist(GetEraShort()+"/2rHardTrig_HCand_Mass", (photons.at(0)+photons.at(1)).M(), weight, 240, 50., 290.);
  }

  if(muons.size()+electrons.size() < 4) return;

  if((!IsDATA||DataStream.Contains("DoubleMuon")) && (trig_DoubleMuon && !trig_MuonEG && !trig_EGamma)){
    FillHists_4leptons(muons, electrons, weight, "1_RAW");
    if(!Charge_4leptons(muons, electrons)) return;
    FillHists_4leptons(muons, electrons, weight, "2_Charge");
    if(!Mass_4leptons(muons, electrons, 5)) return;
    FillHists_4leptons(muons, electrons, weight, "3_Mass5");
    if(!Mass_4leptons(muons, electrons, 10)) return;
    FillHists_4leptons(muons, electrons, weight, "4_Mass10");
    if(!Mass_4leptons(muons, electrons, 15)) return;
    FillHists_4leptons(muons, electrons, weight, "5_Mass15");
    if(!Mass_4leptons(muons, electrons, 20)) return;
    FillHists_4leptons(muons, electrons, weight, "6_Mass20");

  }else if((!IsDATA||DataStream.Contains("EGamma")||DataStream.Contains("DoubleEG")) && (trig_EGamma && !trig_MuonEG)){
    FillHists_4leptons(muons, electrons, weight, "1_RAW");
    if(!Charge_4leptons(muons, electrons)) return;
    FillHists_4leptons(muons, electrons, weight, "2_Charge");
    if(!Mass_4leptons(muons, electrons, 5)) return;
    FillHists_4leptons(muons, electrons, weight, "3_Mass5");
    if(!Mass_4leptons(muons, electrons, 10)) return;
    FillHists_4leptons(muons, electrons, weight, "4_Mass10");
    if(!Mass_4leptons(muons, electrons, 15)) return;
    FillHists_4leptons(muons, electrons, weight, "5_Mass15");
    if(!Mass_4leptons(muons, electrons, 20)) return;
    FillHists_4leptons(muons, electrons, weight, "6_Mass20");

  }else if((!IsDATA||DataStream.Contains("MuonEG")) && (trig_MuonEG)){
    FillHists_4leptons(muons, electrons, weight, "1_RAW");
    if(!Charge_4leptons(muons, electrons)) return;
    FillHists_4leptons(muons, electrons, weight, "2_Charge");
    if(!Mass_4leptons(muons, electrons, 5)) return;
    FillHists_4leptons(muons, electrons, weight, "3_Mass5");
    if(!Mass_4leptons(muons, electrons, 10)) return;
    FillHists_4leptons(muons, electrons, weight, "4_Mass10");
    if(!Mass_4leptons(muons, electrons, 15)) return;
    FillHists_4leptons(muons, electrons, weight, "5_Mass15");
    if(!Mass_4leptons(muons, electrons, 20)) return;
    FillHists_4leptons(muons, electrons, weight, "6_Mass20");

  }
}

ExampleRun_Higgs::ExampleRun_Higgs(){

}

ExampleRun_Higgs::~ExampleRun_Higgs(){

}

void ExampleRun_Higgs::FillHists_4leptons(const vector<Muon>& muons, const vector<Electron>& electrons, double weight, TString prefix){
  if(prefix != "") prefix = prefix+"_";

  if(muons.size() >= 4){
    FillHist(GetEraShort()+"/"+prefix+"4m_HCand_Mass", (muons.at(0)+muons.at(1)+muons.at(2)+muons.at(3)).M(), weight, 240, 50., 290.);
    FillHist(GetEraShort()+"/"+prefix+"4l_HCand_Mass", (muons.at(0)+muons.at(1)+muons.at(2)+muons.at(3)).M(), weight, 240, 50., 290.);
  }else if(electrons.size() >= 4){
    FillHist(GetEraShort()+"/"+prefix+"4e_HCand_Mass", (electrons.at(0)+electrons.at(1)+electrons.at(2)+electrons.at(3)).M(), weight, 240, 50., 290.);
    FillHist(GetEraShort()+"/"+prefix+"4l_HCand_Mass", (electrons.at(0)+electrons.at(1)+electrons.at(2)+electrons.at(3)).M(), weight, 240, 50., 290.);
  }else if(muons.size() >= 2 && electrons.size() >= 2){
    FillHist(GetEraShort()+"/"+prefix+"2m2e_HCand_Mass", (muons.at(0)+muons.at(1)+electrons.at(0)+electrons.at(1)).M(), weight, 240, 50., 290.);
    FillHist(GetEraShort()+"/"+prefix+"4l_HCand_Mass",   (muons.at(0)+muons.at(1)+electrons.at(0)+electrons.at(1)).M(), weight, 240, 50., 290.);
  }
}

bool ExampleRun_Higgs::Charge_4leptons(const vector<Muon>& muons, const vector<Electron>& electrons){
  if(muons.size() >= 4){
    if(muons.at(0).Charge() + muons.at(1).Charge() + muons.at(2).Charge() + muons.at(3).Charge() == 0) return true;
  }else if(electrons.size() >= 4){
    if(electrons.at(0).Charge() + electrons.at(1).Charge() + electrons.at(2).Charge() + electrons.at(3).Charge() == 0) return true;
  }else if(muons.size() >= 2 && electrons.size() >= 2){
    if(muons.at(0).Charge() + muons.at(1).Charge() == 0 && electrons.at(0).Charge() + electrons.at(1).Charge() == 0) return true;
  }

  return false;
}

bool ExampleRun_Higgs::Mass_4leptons(const vector<Muon>& muons, const vector<Electron>& electrons, double mass_cut){
  if(muons.size() >= 4){
    for(unsigned int i=0; i<4; i++){
      for(unsigned int j=i+1; j<4; j++){
        if(muons.at(i).Charge() + muons.at(j).Charge() == 0){
          if((muons.at(i)+muons.at(j)).M() < mass_cut) return false;
        }
      }
    }
    return true;
  }else if(electrons.size() >= 4){
    for(unsigned int i=0; i<4; i++){
      for(unsigned int j=i+1; j<4; j++){
        if(electrons.at(i).Charge() + electrons.at(j).Charge() == 0){
          if((electrons.at(i)+electrons.at(j)).M() < mass_cut) return false;
        }
      }
    }
    return true;
  }else if(muons.size() >= 2 && electrons.size() >= 2){
    if((muons.at(0)+muons.at(1)).M() > mass_cut && (electrons.at(0)+electrons.at(1)).M() > mass_cut) return true;
  }

  return false;
}
