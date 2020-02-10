#include "ControlplotISR.h"

void ControlplotISR::initializeAnalyzer(){

  if(DataYear==2017) IsoMuTriggerName = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v";
  else if(DataYear==2018) IsoMuTriggerName = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v";
}

void ControlplotISR::executeEvent(){


  AnalyzerParameter param;
  AllMuons = GetAllMuons();
  weight_Prefire = GetPrefireWeight(0);

  param.Muon_Tight_ID = "POGTightWithTightIso";
  param.Muon_ID_SF_Key = "NUM_TightID_DEN_genTracks";
  param.Muon_ISO_SF_Key = "NUM_TightRelIso_DEN_TightIDandIPCut";
  param.Muon_Trigger_SF_Key = "TightID_RelPFIso015";

  executeEventFromParameter(param);

}

void ControlplotISR::executeEventFromParameter(AnalyzerParameter param){

  if(!PassMETFilter()) return;

  Event ev = GetEvent();
  if(! (ev.PassTrigger(IsoMuTriggerName) )) return;

  vector<Muon> muons = SelectMuons(AllMuons, param.Muon_Tight_ID, 10., 2.4);
  std::sort(muons.begin(), muons.end(), PtComparing);


  if(muons.size() != 2) return;
  if(muons.at(0).Charge() * muons.at(1).Charge() > 0.) return;
  if(muons.at(0).Pt() < 20 || muons.at(1).Pt() <10) return;
  Particle ZCand = muons.at(0) + muons.at(1);
  if(ZCand.M() < 60. || ZCand.M() > 180.) return;


  double weight = 1.;

  if(!IsDATA){
    weight *= weight_norm_1invpb*ev.GetTriggerLumi("Full");
    weight *= ev.MCweight();
    weight *= weight_Prefire;
  }

  FillHist("mass", ZCand.M(), weight, 120, 60., 180.);
  FillHist("ZpT", ZCand.Pt(), weight, 120, 0., 120.);
  FillHist("leading_pT", muons.at(0).Pt(), weight, 120, 0., 120.);
  FillHist("leading_eta", muons.at(0).Eta(), weight, 48, -2.4, 2.4);
  FillHist("subleading_pT", muons.at(1).Pt(), weight, 120, 0., 120.);
  FillHist("subleading_eta", muons.at(1).Eta(), weight, 48, -2.4, 2.4);

  if(!IsDATA){
    for(unsigned int i=0; i<muons.size(); i++){
      double this_trigsf = 1.;
      if(i == 0) this_trigsf = mcCorr->MuonTrigger_SF("Mu17_"+param.Muon_Trigger_SF_Key, muons.at(i).Eta(), muons.at(i).MiniAODPt());
      else if(i == 1) this_trigsf = mcCorr->MuonTrigger_SF("Mu8_"+param.Muon_Trigger_SF_Key, muons.at(i).Eta(), muons.at(i).MiniAODPt());
      weight *= this_trigsf;
    }

    FillHist("TrigSF_mass", ZCand.M(), weight, 120, 60., 180.);
    FillHist("TrigSF_ZpT", ZCand.Pt(), weight, 120, 0., 120.);
    FillHist("TrigSF_leading_pT", muons.at(0).Pt(), weight, 120, 0., 120.);
    FillHist("TrigSF_leading_eta", muons.at(0).Eta(), weight, 48, -2.4, 2.4);
    FillHist("TrigSF_subleading_pT", muons.at(1).Pt(), weight, 120, 0., 120.);
    FillHist("TrigSF_subleading_eta", muons.at(1).Eta(), weight, 48, -2.4, 2.4);

    for(unsigned int i=0; i<muons.size(); i++){
      double this_idsf  = mcCorr->MuonID_SF(param.Muon_ID_SF_Key,  muons.at(i).Eta(), muons.at(i).MiniAODPt());
      double this_isosf = mcCorr->MuonISO_SF(param.Muon_ISO_SF_Key,  muons.at(i).Eta(), muons.at(i).MiniAODPt());
      weight *= this_idsf * this_isosf;
    }

    FillHist("AllSF_mass", ZCand.M(), weight, 120, 60., 180.);
    FillHist("AllSF_ZpT", ZCand.Pt(), weight, 120, 0., 120.);
    FillHist("AllSF_leading_pT", muons.at(0).Pt(), weight, 120, 0., 120.);
    FillHist("AllSF_leading_eta", muons.at(0).Eta(), weight, 48, -2.4, 2.4);
    FillHist("AllSF_subleading_pT", muons.at(1).Pt(), weight, 120, 0., 120.);
    FillHist("AllSF_subleading_eta", muons.at(1).Eta(), weight, 48, -2.4, 2.4);

  }

}

ControlplotISR::ControlplotISR(){

}

ControlplotISR::~ControlplotISR(){

}


