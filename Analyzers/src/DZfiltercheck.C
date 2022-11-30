#include "DZfiltercheck.h"

void DZfiltercheck::initializeAnalyzer(){

}

void DZfiltercheck::executeEvent(){


  AnalyzerParameter param;

  muons = SMPGetMuons("POGMediumWithLooseTrkIso", 5, 2.4);
  if(muons.size() < 2) return;
  if(muons.at(0).Charge() * muons.at(1).Charge() > 0) return;
  if(muons.at(0).Pt() < 20 || muons.at(1).Pt() < 10) return;
  //if((muons.at(0)+muons.at(1)).M() < 81 || (muons.at(0)+muons.at(1)).M() > 101) return;

  executeEventFromParameter(param);

}

void DZfiltercheck::executeEventFromParameter(AnalyzerParameter param){

  if(!PassMETFilter()) return;

  Event ev = GetEvent();

  double weight = 1.;
  if(!IsDATA){
    weight *= MCweight();
    weight *= ev.GetTriggerLumi("Full");
    weight *= GetPrefireWeight(0);
  }


  if(DataYear==2016){
    double acceptance_EMTF = 0.;
    if(muons.at(0).Eta() * muons.at(1).Eta() > 0 && fabs(muons.at(0).Eta()) > 1.2 && fabs(muons.at(1).Eta()) > 1.2 && fabs(muons.at(0).DeltaPhi(muons.at(1))) < 3.14159/3) acceptance_EMTF = 1.;

    FillHist("noTrig_acceptance_EMTF", acceptance_EMTF, weight, 2, 0., 2.);
    if(acceptance_EMTF == 1.){
      FillHist("noTrig_acceptance_EMTF_pt0", muons.at(0).Pt(), weight, 200, 0., 200.);
      FillHist("noTrig_acceptance_EMTF_pt1", muons.at(1).Pt(), weight, 200, 0., 200.);
      FillHist("noTrig_acceptance_EMTF_eta0", muons.at(0).Eta(), weight, 60, -3, 3.);
      FillHist("noTrig_acceptance_EMTF_eta1", muons.at(1).Eta(), weight, 60, -3, 3.);
      FillHist("noTrig_acceptance_EMTF_Zmass", (muons.at(0)+muons.at(1)).M(), weight, 200, 0., 200.);
      FillHist("noTrig_acceptance_EMTF_Zpt", (muons.at(0)+muons.at(1)).Pt(), weight, 200, 0., 200.);
      FillHist("noTrig_acceptance_EMTF_Zrap", (muons.at(0)+muons.at(1)).Rapidity(), weight, 100, -5., 5.);
    }

    vector<TString> double_trig = {
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v",
      "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
      "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
      "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",
      "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
    };
    vector<TString> single_trig = {
      "HLT_IsoMu24_v",
      "HLT_IsoTkMu24_v",
    };

    if(ev.PassTrigger(single_trig)) FillHist("singleTrig_acceptance_EMTF", acceptance_EMTF, weight, 2, 0., 2.);
    if(ev.PassTrigger(double_trig)) FillHist("doubleTrig_acceptance_EMTF", acceptance_EMTF, weight, 2, 0., 2.);
    if(ev.PassTrigger(single_trig) || ev.PassTrigger(double_trig)) FillHist("allTrig_acceptance_EMTF", acceptance_EMTF, weight, 2, 0., 2.);

  }else{

    /*
    bool trig_NoDZ = ev.PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");
    bool trig_FullDZ = ev.PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v");
    if(DataYear==2018) trig_FullDZ = ev.PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v");

    if(!trig_NoDZ) return;

    double etabins[NEtaBin+1];
    for(int i=0; i<NEtaBin+1; i++) etabins[i] = vec_etabins.at(i);
    double ptbins[NPtBin+1];
    for(int i=0; i<NPtBin+1; i++) ptbins[i] = vec_ptbins.at(i);

    Muon mu_minus = muons.at(0).Charge() < 0? muons.at(0) : muons.at(1);
    Muon mu_plus  = muons.at(0).Charge() < 0? muons.at(1) : muons.at(0);

    FillHist("eta_denom", mu_minus.Eta(), mu_plus.Eta(), weight, NEtaBin, etabins, NEtaBin, etabins);
    FillHist("pt_denom", mu_minus.Pt(), mu_plus.Pt(), weight, NPtBin, ptbins, NPtBin, ptbins);
    if(mu_minus.Pt() > 20 && mu_plus.Pt() > 20){
      FillHist("eta20_denom", mu_minus.Eta(), mu_plus.Eta(), weight, NEtaBin, etabins, NEtaBin, etabins);
      FillHist("pt20_denom", mu_minus.Pt(), mu_plus.Pt(), weight, NPtBin, ptbins, NPtBin, ptbins);
    }

    if(trig_FullDZ){
      FillHist("eta_num", mu_minus.Eta(), mu_plus.Eta(), weight, NEtaBin, etabins, NEtaBin, etabins);
      FillHist("pt_num", mu_minus.Pt(), mu_plus.Pt(), weight, NPtBin, ptbins, NPtBin, ptbins);
      if(mu_minus.Pt() > 20 && mu_plus.Pt() > 20){
        FillHist("eta20_num", mu_minus.Eta(), mu_plus.Eta(), weight, NEtaBin, etabins, NEtaBin, etabins);
        FillHist("pt20_num", mu_minus.Pt(), mu_plus.Pt(), weight, NPtBin, ptbins, NPtBin, ptbins);
      }
    }
    */
  }
}

DZfiltercheck::DZfiltercheck(){

}

DZfiltercheck::~DZfiltercheck(){

}


