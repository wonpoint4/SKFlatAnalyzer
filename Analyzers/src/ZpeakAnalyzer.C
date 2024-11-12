#include "ZpeakAnalyzer.h"

ZpeakAnalyzer::ZpeakAnalyzer(){
}
ZpeakAnalyzer::~ZpeakAnalyzer(){
}
void ZpeakAnalyzer::initializeAnalyzer(){
  SMPAnalyzerCore::initializeAnalyzer(); //setup zpt roc z0
  fChain->SetBranchStatus("pfMET_*",false);
  fChain->SetBranchStatus("pfMET_Type1_pt",true);
  fChain->SetBranchStatus("jet_*",false);
  fChain->SetBranchStatus("fatjet_*",false);
  fChain->SetBranchStatus("photon_*",false);
}
void ZpeakAnalyzer::executeEvent(){
  if(!IsDATA||DataStream.Contains("DoubleMuon")){
    executeEventWithParameter(MakeParameter("mm"));
    // if(GetEra()=="2016preVFP"||GetEra()=="2016postVFP"){
    //   executeEventWithParameter(MakeParameter("mm","medium_nohip"));
    //   executeEventWithParameter(MakeParameter("mm","tight"));
    // }
    executeEventWithParameter(MakeParameter("MM","fake"));
  }
  if(!IsDATA||DataStream.Contains("SingleMuon")){
    executeEventWithParameter(MakeParameter("mu"));
  }
  if(!IsDATA||DataStream.Contains("DoubleEG")||DataStream.Contains("EGamma")){
    executeEventWithParameter(MakeParameter("ee"));
    executeEventWithParameter(MakeParameter("EE","fake"));
  }
  if(!IsDATA||DataStream.Contains("SingleElectron")||DataStream.Contains("EGamma")){
    executeEventWithParameter(MakeParameter("el"));
  }
}
SMPAnalyzerCore::Parameter ZpeakAnalyzer::MakeParameter(TString key,TString option){
  Parameter p=SMPAnalyzerCore::MakeParameter(key,option);
  p.variationbits|=EfficiencyWeight;
  if(option.Contains("medium_nohip")){
    p.prefix="medium_nohip/"+p.prefix;
    p.SetMuons(MuonMomentumCorrection(SMPGetMuons("POGMediumNoHipWithLooseTrkIso",8.0,2.4),0,0));
    p.option.ReplaceAll("medium_nohip","");
  }
  if(option.Contains("tight")){
    p.prefix="tight/"+p.prefix;
    p.SetMuons(MuonMomentumCorrection(SMPGetMuons("POGTightWithLooseTrkIso",8.0,2.4),0,0));
    p.option.ReplaceAll("tight","");
  }
  return p;
}
void ZpeakAnalyzer::EvalDefaultWeight(Parameter& p){
  p.default_weight=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.CFSF*p.w.btagSF*p.w.bchargeSF*p.w.zptweight*p.w.weakweight*p.w.topptweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.fakeTF;
  // if(!IsDATA){
  //   if(p.channel[0]=='e'){
  //     int ie=p.w.electronIDSF_sys.size()-1;
  //     p.default_weight*=p.w.electronIDSF_sys[ie][0]/p.w.electronIDSF;
  //   }else if(p.channel[0]=='m'){
  //     int im=p.w.muonIDSF_sys.size()-1;
  //     p.default_weight*=p.w.muonIDSF_sys[im][0]/p.w.muonIDSF;
  //   }
  // }
  p.weight=p.default_weight;
}

SMPAnalyzerCore::Variations ZpeakAnalyzer::MakeVariations(const Parameter& p){
  Variations v;
  AddVariationWeight(v,"",p.default_weight);
  if(!IsDATA){
    AddVariationWeight(v,"_z0weight",p.default_weight*p.w.z0weight);
    EvalVariationsCF(p,v);
    if(p.channel[0]=='e'){
      int ie=p.w.electronIDSF_sys.size()-1;
      AddVariationWeight(v,"_efficiency_residual",p.default_weight/p.w.electronIDSF*p.w.electronIDSF_sys[ie][0]);
    }else if(p.channel[0]=='m'){
      int im=p.w.muonIDSF_sys.size()-1;
      AddVariationWeight(v,"_efficiency_residual",p.default_weight/p.w.muonIDSF*p.w.muonIDSF_sys[im][0]);
    }
  }
  if(p.channel[0]=='e'){
    //int ie=p.w.electronIDSF_sys.size()-1;
    //AddVariationWeight(v,"_efficiency_residual",p.default_weight/p.w.electronIDSF*p.w.electronIDSF_sys[ie][0]);      
    AddVariationElectronEnergy(v,"_noroccor",-1,0);
    AddVariationElectronEnergy(v,"_roccor_residual",-2,0);
    AddVariationElectronEnergy(v,"_pogcor",-3,0);
    //AddVariationElectronEnergy(v,"_residual",-2,0,p.default_weight/p.w.electronIDSF*p.w.electronIDSF_sys[ie][0]);
  }else if(p.channel[0]=='m'){
    //int im=p.w.muonIDSF_sys.size()-1;
    //AddVariationWeight(v,"_efficiency_residual",p.default_weight/p.w.muonIDSF*p.w.muonIDSF_sys[im][0]);
    AddVariationMuonMomentum(v,"_noroccor",-1,0);
    AddVariationMuonMomentum(v,"_roccor_residual",-2,0);
    //AddVariationMuonMomentum(v,"_residual",-2,0,p.default_weight/p.w.muonIDSF*p.w.muonIDSF_sys[im][0]);
  }
  return v;
}
void ZpeakAnalyzer::FillHists(Parameter& p){
  TString pre=p.prefix+p.hprefix;
  TString suf=p.suffix+p.vsuffix;
  double weight=p.weight;

  TLorentzVector dilepton=(*p.lepton0)+(*p.lepton1);
  double dimass=dilepton.M();
  if(dimass<52) return;
  FillHist(pre+"dimass"+suf,p.lepton0->Eta(),p.lepton0->Pt(),dimass,weight,netabin,etabins,nptbin,ptbins,nmassbin,massbins);
  FillHist(pre+"dimass"+suf,p.lepton1->Eta(),p.lepton1->Pt(),dimass,weight,netabin,etabins,nptbin,ptbins,nmassbin,massbins);
  Lepton *lm=NULL,*lp=NULL;
  if(p.lepton0->Charge()<0){
    lm=p.lepton0; lp=p.lepton1;
  }else{
    lm=p.lepton0; lp=p.lepton1;
  }
  FillHist(pre+"lmetalpetam2"+suf,fabs(lm->Eta()),fabs(lp->Eta()),dimass*dimass,weight,25,0,2.5,25,0,2.5,170,4500,13000);
  if(!IsDATA){
    if(p.lepton0->Charge()*p.truth_lepton0.Charge()<0){
      FillHist(pre+"dimass_cf"+suf,p.lepton0->Eta(),p.lepton0->Pt(),dimass,weight,netabin,etabins,nptbin,ptbins,nmassbin,massbins);
    }
    if(p.lepton1->Charge()*p.truth_lepton1.Charge()<0){
      FillHist(pre+"dimass_cf"+suf,p.lepton1->Eta(),p.lepton1->Pt(),dimass,weight,netabin,etabins,nptbin,ptbins,nmassbin,massbins);
    }
  }
  FillHist(pre+"ym"+suf,fabs(dilepton.Rapidity()),dimass,weight,rochester_nybin,rochester_ybins,rochester_nmbin,rochester_mbins);
  // if(p.vsuffix=="_roccor_residual"){
  //   double weight_efficiency_residual=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.CFSF*p.w.btagSF*p.w.bchargeSF*p.w.zptweight*p.w.weakweight*p.w.topptweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF;
  //   FillHist(pre+"dimass"+suf+"_efficiency_residual",p.lepton0->Eta(),p.lepton0->Pt(),dimass,weight_efficiency_noresidual,netabin,etabins,nptbin,ptbins,nmassbin,massbins);
  //   FillHist(pre+"dimass"+suf+"_efficiency_residual",p.lepton1->Eta(),p.lepton1->Pt(),dimass,weight_efficiency_noresidual,netabin,etabins,nptbin,ptbins,nmassbin,massbins);
  //   FillHist(pre+"ym"+suf+"_efficiency_residual",fabs(dilepton.Rapidity()),dimass,weight_efficiency_noresidual,rochester_nybin,rochester_ybins,rochester_nmbin,rochester_mbins);
  // }
   
}
