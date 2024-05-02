#include "BBAnalyzer.h"

BBAnalyzer::BBAnalyzer(){
}
BBAnalyzer::~BBAnalyzer(){
}
void BBAnalyzer::initializeAnalyzer(){
  SMPAnalyzerCore::initializeAnalyzer(); //setup zpt roc z0
  fChain->SetBranchStatus("pfMET_*",false);
  fChain->SetBranchStatus("pfMET_Type1_pt",true);
  fChain->SetBranchStatus("fatjet_*",false);
  fChain->SetBranchStatus("photon_*",false);
  IsDileptonSkim=GetSkimName()=="Dilepton";
}
void BBAnalyzer::executeEvent(){
  if(IsDileptonSkim){
    if(!IsDATA||DataStream.Contains("DoubleMuon")){
      executeEventWithParameter(MakeParameter("mm"));
    }
    if(!IsDATA||DataStream.Contains("SingleMuon")){
      executeEventWithParameter(MakeParameter("me"));
    }
    if(!IsDATA||DataStream.Contains("DoubleEG")||DataStream.Contains("EGamma")){
      executeEventWithParameter(MakeParameter("ee"));
    }
  }else{
    if(!IsDATA||DataStream.Contains("SingleMuon")){
      executeEventWithParameter(MakeParameter("mn"));
    }
    if(!IsDATA||DataStream.Contains("SingleElectron")||DataStream.Contains("EGamma")){
      executeEventWithParameter(MakeParameter("en"));
    }
  }
}

void BBAnalyzer::EvalDefaultWeight(Parameter& p){
  p.default_weight=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.CFSF*p.w.btagSF*p.w.zptweight*p.w.weakweight*p.w.topptweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF;
  p.weight=p.default_weight;
}

SMPAnalyzerCore::Variations BBAnalyzer::MakeVariations(const Parameter& p){
  Variations v;
  AddVariationWeight(v,"",p.default_weight);
  if(!IsDATA&&!p.hprefix.Contains("ss_")){
    EvalVariationsPUweight(p,v);
    EvalVariationsPrefireweight(p,v);
    EvalVariationsBtag(p,v);
    EvalVariationsBcharge(p,v);
    EvalVariationsEtc(p,v);
  }
  if(!IsDATA&&p.hprefix.Contains("ss_")){
    EvalVariationsCF(p,v);
  }
  if(!p.hprefix.Contains("ss_")){
    EvalVariationsJetCorrection(p,v);
  }
  if(MCSample.Contains("TTLL")||MCSample.Contains("TTLJ")){
    if(!IsDATA&&p.hprefix==""){
      EvalVariationsPDF(p,v);
    }
  }

  return v;
}

void BBAnalyzer::EvalVariationsBcharge(const Parameter& p,Variations& v){
  AddVariationWeight(v,"_bchargeSF",p.default_weight*p.w.bchargeSF);
  AddVariationWeight(v,"_bchargeSF_s0m0",p.default_weight*p.w.bchargeSF_s0m0);
  AddVariationWeight(v,"_bchargeSF_s0m1",p.default_weight*p.w.bchargeSF_s0m1);
}

void BBAnalyzer::FillHists(Parameter& p){
  TString pre=p.prefix+p.hprefix;
  TString suf=p.suffix+p.vsuffix;
  double weight=p.weight;

  int nbjet=count_if(p.bjets.begin(),p.bjets.end(),[&p](auto& jet){return jet.Pt()>p.c.jetpt;});
  if(nbjet<2) return;

  int bits=0;
  if(p.bjets.at(0).userFloat["AFBCharge"]>0) bits+=1<<0;
  if(p.bjets.at(1).userFloat["AFBCharge"]>0) bits+=1<<1;
  TString csuf=Form("_c%d",bits);
  
  if(IsDileptonSkim){
    TLorentzVector dilepton=(*p.lepton0)+(*p.lepton1);
    double dimass=dilepton.M();
    double dirap=dilepton.Rapidity();
    double dipt=dilepton.Pt();

    if(dimass<52) return;

    FillHist(pre+"dimass"+suf,dimass,weight,nmassbin,massbins);
    if( (p.channel=="ee" || p.channel=="mm") && dimass>76 && dimass<106) return;

    FillHist(pre+"dimass"+csuf+suf,dimass,weight,nmassbin,massbins);
    FillHist(pre+"dirap"+csuf+suf,dirap,weight,netabin,etabins);
    FillHist(pre+"dipt"+csuf+suf,dipt,weight,nptbin,ptbins);
  }

  if(MCSample.Contains("TTLL")||MCSample.Contains("TTLJ")){
    if( 
       (p.bjets.at(0).GenHFHadronMatcherFlavour()==5 && p.bjets.at(1).GenHFHadronMatcherFlavour()==5 )
       && (abs(p.bjets.at(0).GenHFHadronMatcherOrigin())==6 && abs(p.bjets.at(1).GenHFHadronMatcherOrigin())==6)
       && (p.bjets.at(0).GenHFHadronMatcherOrigin()*p.bjets.at(1).GenHFHadronMatcherOrigin()<0) 
	){
      for(int i=0;i<2;i++){
	int correct=p.bjets.at(i).userFloat["AFBCharge"]*p.bjets.at(i).GenHFHadronMatcherOrigin()<0;
	TString scharge=p.bjets.at(i).GenHFHadronMatcherOrigin()<0 ? "p" : "m";
	FillHist(pre+Form("b%d%scorrect",i,scharge.Data())+suf,correct,weight,2,0,2);
	FillHist(pre+Form("b%d%scharge",i,scharge.Data())+suf,p.bjets.at(i).userFloat["AFBCharge"],weight,100,-5,5);
	if(p.vsuffix=""){
	  if(fabs(p.bjets.at(i).userFloat["AFBCharge"])<1){
	    FillHist(pre+Form("b%d%scorrect_type0",i,scharge.Data())+suf,correct,weight,2,0,2);
	  }else if(fabs(p.bjets.at(i).userFloat["AFBCharge"])<3){
	    FillHist(pre+Form("b%d%scorrect_type1",i,scharge.Data())+suf,correct,weight,2,0,2);
	  }else if(fabs(p.bjets.at(i).userFloat["AFBCharge"])<5){
	    FillHist(pre+Form("b%d%scorrect_type2",i,scharge.Data())+suf,correct,weight,2,0,2);	  
	  }
	}
      }
    }else{
      pre+="unmatched_";    
    }
  }

  FillHist(pre+"charge"+suf,bits,weight,4,0,4);
  for(int i=0;i<2;i++){
    const Jet& bjet=p.bjets.at(i);
    double bjet_eta=bjet.Eta();
    double bjet_pt=bjet.Pt();
    FillHist(pre+"bjetapt"+csuf+suf,bjet_eta,bjet_pt,weight,netabin,etabins,nptbin,ptbins);
  }
}
