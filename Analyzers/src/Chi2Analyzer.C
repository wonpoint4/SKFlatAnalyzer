#include "Chi2Analyzer.h"

Chi2Analyzer::Chi2Analyzer(){
}
Chi2Analyzer::~Chi2Analyzer(){
}
void Chi2Analyzer::initializeAnalyzer(){
  SMPAnalyzerCore::initializeAnalyzer(); //setup zpt roc z0
  fChain->SetBranchStatus("pfMET_*",false);
  fChain->SetBranchStatus("pfMET_Type1_pt",true);
  //fChain->SetBranchStatus("jet_*",false);
  fChain->SetBranchStatus("fatjet_*",false);
  //fChain->SetBranchStatus("photon_*",false);
}
void Chi2Analyzer::executeEvent(){
  if(!IsDATA||DataStream.Contains("SingleMuon")){
    executeEventWithParameter(MakeParameter("mu")); 
  }
  if(!IsDATA||DataStream.Contains("DoubleMuon")){
    executeEventWithParameter(MakeParameter("mm")); 
  }
  if(!IsDATA||DataStream.Contains("SingleElectron")||DataStream.Contains("EGamma")){
    executeEventWithParameter(MakeParameter("el"));
  }
  if(!IsDATA||DataStream.Contains("DoubleEG")||DataStream.Contains("EGamma")){
    executeEventWithParameter(MakeParameter("ee"));
  }
}
SMPAnalyzerCore::Parameter Chi2Analyzer::MakeParameter(TString key,TString option){
  Parameter p=SMPAnalyzerCore::MakeParameter(key,option);
  p.variationbits|=EfficiencyWeight;
  p.variationbits|=SystematicWeight;
  if(IsDYSample || IsDATA){
    p.variationbits|=PDFWeight;
    p.variationbits|=LeptonCorrection;
  }
  return p;
}
void Chi2Analyzer::EvalDefaultWeight(Parameter& p){  
  SMPAnalyzerCore::EvalDefaultWeight(p);
}
SMPAnalyzerCore::Variations Chi2Analyzer::MakeVariations(const Parameter& p){
  Variations v=SMPAnalyzerCore::MakeVariations(p);
  return v;
}
void Chi2Analyzer::FillHists(Parameter& p){
  TLorentzVector dilepton=(*p.lepton0)+(*p.lepton1);
  double dimass=dilepton.M();
  TString pre=p.prefix+p.hprefix;
  TString suf=p.suffix+p.vsuffix;
  double w=p.weight;
  if(dimass>=54&&dimass<150){
    FillHist(pre+"ym"+suf,fabs(dilepton.Rapidity()),dimass,w,rochester_nybin,rochester_ybins,rochester_nmbin,rochester_mbins);
    if(dimass>=80&&dimass<100){
      for(int i=0;i<(int)p.leptons.size();i++){
	double pt=p.leptons.at(i)->Pt();
	double eta=p.leptons.at(i)->Eta();
	if(p.leptons.at(i)->LeptonFlavour()==Lepton::Flavour::MUON){
	  FillHist(Form("%sl%detapt%s",pre.Data(),i,suf.Data()),eta,pt,w,netabin_muonID,etabins_muonID,nptbin_muonID,ptbins_muonID);
	}else if(p.leptons.at(i)->LeptonFlavour()==Lepton::Flavour::ELECTRON){
	  FillHist(Form("%sl%detapt%s",pre.Data(),i,suf.Data()),eta,pt,w,netabin_electronID,etabins_electronID,nptbin_electronID,ptbins_electronID);
	}
      }
    }
  }
}
