#include "DZfiltercheck.h"

void DZfiltercheck::initializeAnalyzer(){

}

void DZfiltercheck::executeEvent(){


  AnalyzerParameter param;

  muons = SMPGetMuons("POGMediumWithLooseTrkIso", 5, 2.4);
  if(muons.size() < 2) return;
  if(muons.at(0).Charge() * muons.at(1).Charge() > 0) return;
  if(muons.at(0).Pt() < 20 || muons.at(1).Pt() < 10) return;

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
    bool acceptance_EMTF = false;
    if(muons.at(0).Eta() * muons.at(1).Eta() > 0 && fabs(muons.at(0).Eta()) > 1.2 && fabs(muons.at(1).Eta()) > 1.2 && fabs(muons.at(0).DeltaPhi(muons.at(1))) < 3.14159/3) acceptance_EMTF = true;
    bool weak_acceptance_EMTF = false;
    if(muons.at(0).Eta() * muons.at(1).Eta() > 0 && fabs(muons.at(0).Eta()) > 1.2 && fabs(muons.at(1).Eta()) > 1.2) weak_acceptance_EMTF = true;

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

    FillHists(muons, weight, "noT_", acceptance_EMTF);
    if(ev.PassTrigger(single_trig)) FillHists(muons, weight, "singleT_", acceptance_EMTF);
    if(ev.PassTrigger(double_trig)) FillHists(muons, weight, "doubleT_", acceptance_EMTF);
    if(ev.PassTrigger(single_trig) || ev.PassTrigger(double_trig)) FillHists(muons, weight, "allT_", acceptance_EMTF);
    if(ev.PassTrigger(single_trig) && !ev.PassTrigger(double_trig)) FillHists(muons, weight, "onlyST_", acceptance_EMTF);
    if(!ev.PassTrigger(single_trig) && ev.PassTrigger(double_trig)) FillHists(muons, weight, "onlyDT_", acceptance_EMTF);

    if((muons.at(0)+muons.at(1)).M() < 52) return;

    FillHists(muons, weight, "Mass52/noT_", acceptance_EMTF);
    if(ev.PassTrigger(single_trig)) FillHists(muons, weight, "Mass52/singleT_", acceptance_EMTF);
    if(ev.PassTrigger(double_trig)) FillHists(muons, weight, "Mass52/doubleT_", acceptance_EMTF);
    if(ev.PassTrigger(single_trig) || ev.PassTrigger(double_trig)) FillHists(muons, weight, "Mass52/allT_", acceptance_EMTF);
    if(ev.PassTrigger(single_trig) && !ev.PassTrigger(double_trig)) FillHists(muons, weight, "Mass52/onlyST_", acceptance_EMTF);
    if(!ev.PassTrigger(single_trig) && ev.PassTrigger(double_trig)) FillHists(muons, weight, "Mass52/onlyDT_", acceptance_EMTF);

    if((muons.at(0)+muons.at(1)).M() < 81 || (muons.at(0)+muons.at(1)).M() > 101) return;

    FillHists(muons, weight, "Zpeak/noT_", acceptance_EMTF);
    if(ev.PassTrigger(single_trig)) FillHists(muons, weight, "Zpeak/singleT_", acceptance_EMTF);
    if(ev.PassTrigger(double_trig)) FillHists(muons, weight, "Zpeak/doubleT_", acceptance_EMTF);
    if(ev.PassTrigger(single_trig) || ev.PassTrigger(double_trig)) FillHists(muons, weight, "Zpeak/allT_", acceptance_EMTF);
    if(ev.PassTrigger(single_trig) && !ev.PassTrigger(double_trig)) FillHists(muons, weight, "Zpeak/onlyST_", acceptance_EMTF);
    if(!ev.PassTrigger(single_trig) && ev.PassTrigger(double_trig)) FillHists(muons, weight, "Zpeak/onlyDT_", acceptance_EMTF);

    FillHists(muons, weight, "Weak_Zpeak/noT_", weak_acceptance_EMTF);
    if(ev.PassTrigger(single_trig)) FillHists(muons, weight, "Weak_Zpeak/singleT_", weak_acceptance_EMTF);
    if(ev.PassTrigger(double_trig)) FillHists(muons, weight, "Weak_Zpeak/doubleT_", weak_acceptance_EMTF);
    if(ev.PassTrigger(single_trig) || ev.PassTrigger(double_trig)) FillHists(muons, weight, "Weak_Zpeak/allT_", weak_acceptance_EMTF);
    if(ev.PassTrigger(single_trig) && !ev.PassTrigger(double_trig)) FillHists(muons, weight, "Weak_Zpeak/onlyST_", weak_acceptance_EMTF);
    if(!ev.PassTrigger(single_trig) && ev.PassTrigger(double_trig)) FillHists(muons, weight, "Weak_Zpeak/onlyDT_", weak_acceptance_EMTF);

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


void DZfiltercheck::FillHists(vector<Muon> muons, double weight, TString prefix, bool acceptance_EMTF){
  FillHist((prefix+"acceptance_EMTF").ReplaceAll("/","_"), acceptance_EMTF, weight, 2, 0., 2.);
    if(acceptance_EMTF){
      FillHist("1D/"+prefix+"pt0", muons.at(0).Pt(), weight, 200, 0., 200.);
      FillHist("1D/"+prefix+"pt1", muons.at(1).Pt(), weight, 200, 0., 200.);
      FillHist("1D/"+prefix+"eta0", muons.at(0).Eta(), weight, 50, -2.5, 2.5);
      FillHist("1D/"+prefix+"eta1", muons.at(1).Eta(), weight, 50, -2.5, 2.5);
      FillHist("1D/"+prefix+"phi0", muons.at(0).Phi(), weight, 70, -3.5, 3.5);
      FillHist("1D/"+prefix+"phi1", muons.at(1).Phi(), weight, 70, -3.5, 3.5);
      FillHist("1D/"+prefix+"dimass", (muons.at(0)+muons.at(1)).M(), weight, 200, 0., 200.);
      FillHist("1D/"+prefix+"dipt", (muons.at(0)+muons.at(1)).Pt(), weight, 200, 0., 200.);
      FillHist("1D/"+prefix+"dirap", (muons.at(0)+muons.at(1)).Rapidity(), weight, 100, -5., 5.);
      FillHist("2D/"+prefix+"pt", muons.at(0).Pt(), muons.at(1).Pt(), weight, 40, 0., 200., 40, 0., 200.);
      FillHist("2D/"+prefix+"eta", muons.at(0).Eta(), muons.at(1).Eta(), weight, 50, -2.5, 2.5, 50, -2.5, 2.5);
      FillHist("2D/"+prefix+"phi", muons.at(0).Phi(), muons.at(1).Phi(), weight, 70, -3.5, 3.5, 70, -3.5, 3.5);
      FillHist("2D/"+prefix+"phifine", muons.at(0).Phi(), muons.at(1).Phi(), weight, 280, -3.5, 3.5, 280, -3.5, 3.5);
      FillHist("2D/"+prefix+"phifine2", muons.at(0).Phi(), muons.at(1).Phi(), weight, 180, -3.14159, 3.14159, 180, -3.14159, 3.14159);
      FillHist("2D/"+prefix+"phifine3", muons.at(0).Phi(), muons.at(1).Phi(), weight, 360, -3.14159, 3.14159, 360, -3.14159, 3.14159);
      FillHist("2D/"+prefix+"phicoarse", muons.at(0).Phi(), muons.at(1).Phi(), weight, 35, -3.5, 3.5, 35, -3.5, 3.5);
    }
}
