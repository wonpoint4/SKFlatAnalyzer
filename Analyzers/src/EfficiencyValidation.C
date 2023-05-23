#include "EfficiencyValidation.h"

EfficiencyValidation::EfficiencyValidation(){
}
EfficiencyValidation::~EfficiencyValidation(){
}
void EfficiencyValidation::initializeAnalyzer(){
  SMPAnalyzerCore::initializeAnalyzer(); //setup zpt roc z0
  fChain->SetBranchStatus("pfMET_*",false);
  fChain->SetBranchStatus("pfMET_Type1_pt",true);
  fChain->SetBranchStatus("jet_*",false);
  fChain->SetBranchStatus("fatjet_*",false);
  fChain->SetBranchStatus("photon_*",false);
}
void EfficiencyValidation::executeEvent(){
  //////// nominal channels //////////
  if(!IsDATA||DataStream.Contains("DoubleMuon")) 
    executeEventWithParameter(MakeParameter("mm"));
  if(!IsDATA||DataStream.Contains("SingleMuon"))
    executeEventWithParameter(MakeParameter("mu"));
  if(!IsDATA||DataStream.Contains("DoubleEG")||DataStream.Contains("EGamma")) 
    executeEventWithParameter(MakeParameter("ee"));
  if(!IsDATA||DataStream.Contains("SingleElectron")||DataStream.Contains("EGamma")){
    executeEventWithParameter(MakeParameter("el"));
    executeEventWithParameter(MakeParameter("el","SelQ"));
  }

  //////// testing channels //////////

  if(GetEra()=="2016preVFP"){
    if(!IsDATA||DataStream.Contains("SingleMuon")){
      /*
      {Parameter p=MakeParameter("mu");
      p.prefix="old0/mu2016a/";
      p.SetMuonKeys("Muon_MediumID_trkIsoLoose_old","",{"IsoMu24_MediumID_trkIsoLoose_old"});
      executeEventWithParameter(p);}
      {Parameter p=MakeParameter("mu");
      p.prefix="old1/mu2016a/";
      p.SetMuonKeys("Muon_MediumID_trkIsoLoose","",{"IsoMu24_MediumID_trkIsoLoose_old"});
      executeEventWithParameter(p);}
      {Parameter p=MakeParameter("mu");
      p.prefix="old2/mu2016a/";
      p.SetMuonKeys("Muon_MediumID_trkIsoLoose_old","",{"IsoMu24_MediumID_trkIsoLoose"});
      executeEventWithParameter(p);}
      */
    }
  }else if(GetEra()=="2017"){
    if(!IsDATA||DataStream.Contains("DoubleEG")){
      Parameter p=MakeParameter("ee");
      p.suffix="_noL1";
      p.k.triggerSF={"Ele23Leg1_MediumID_v3_2","Ele12Leg2_MediumID"};
      //executeEventWithParameter(p);
    }
    if(!IsDATA||DataStream.Contains("SingleElectron")){
      Parameter p=MakeParameter("el");
      p.prefix="el201727/";
      p.triggers={"HLT_Ele27_WPTight_Gsf_v"};
      p.k.triggerSF={"Ele27_MediumID"};
      if(!IsDATA) p.w.lumiweight*=31.72/41.54;
      executeEventWithParameter(p);

      p=MakeParameter("el");
      p.prefix="el201732/";
      p.triggers={"HLT_Ele32_WPTight_Gsf_v"};
      p.k.triggerSF={"Ele32_MediumID"};
      p.c.lepton0pt=35;
      if(!IsDATA) p.w.lumiweight*=27.13/41.54;
      executeEventWithParameter(p);
    }
  }else if(GetEra()=="2018"){
    if(!IsDATA||DataStream.Contains("DoubleMuon")){
      Parameter p=MakeParameter("mm");
      p.suffix="_new";
      p.SetMuonKeys("IDISO_SF_MediumID_trkIsoLoose_Q_new","",{"Mu17Leg1_MediumID_trkIsoLoose_Q_new","Mu8Leg2_MediumID_trkIsoLoose_Q_new"});
      //executeEventWithParameter(p);
    }
    if(!IsDATA||DataStream.Contains("SingleMuon")){ 
      Parameter p=MakeParameter("mu");
      p.suffix="_new";
      p.SetMuonKeys("IDISO_SF_MediumID_trkIsoLoose_Q_new","",{"IsoMu24_MediumID_trkIsoLoose_Q_new"});
      //executeEventWithParameter(p);
    }
    if(!IsDATA||DataStream.Contains("EGamma")){
      Parameter p=MakeParameter("ee");
      p.suffix="_noL1";
      p.k.triggerSF={"Ele23Leg1_MediumID_Q_v3","Ele12Leg2_MediumID_Q"};
      //executeEventWithParameter(p);

      p=MakeParameter("el");
      p.prefix="el201828/";
      p.triggers={"HLT_Ele28_WPTight_Gsf_v"};
      p.k.triggerSF={"Ele28_MediumID"};
      if(!IsDATA) p.w.lumiweight*=23687.253/59827.879;
      executeEventWithParameter(p);

      p=MakeParameter("el");
      p.prefix="el201832/";
      p.triggers={"HLT_Ele32_WPTight_Gsf_v"};
      p.k.triggerSF={"Ele32_MediumID"};
      p.c.lepton0pt=35;
      executeEventWithParameter(p);
    }
  }
}
SMPAnalyzerCore::Parameter EfficiencyValidation::MakeParameter(TString key,TString option){
  Parameter p=SMPAnalyzerCore::MakeParameter(key,option);
  p.weightbit|=EfficiencyWeight;
  return p;
}
void EfficiencyValidation::EvalWeights(Parameter& p){
  p.weightmap[""]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
  if(!IsDATA&&p.suffix==""&&!p.hprefix.Contains("ss_")){
    //p.weightmap["_muontrackingSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["muontrackingSF"];
    p.weightmap["_noweight"]=p.w.lumiweight;
    p.weightmap["_noPUweight"]=p.w.lumiweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
    p.weightmap["_noprefireweight"]=p.w.lumiweight*p.w.PUweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
    p.weightmap["_nozptweight"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
    p.weightmap["_noz0weight"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
    p.weightmap["_noweakweight"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
    
    for(int j=0,nj=fEff->nreplica;j<nj;j++){
      double electronRECOSF=p.w.electronRECOSF_sys.size() ? p.w.electronRECOSF_sys[0][j] : 1.;
      double electronIDSF=p.w.electronIDSF_sys.size() ? p.w.electronIDSF_sys[0][j] : 1.;
      double muonIDSF=p.w.muonIDSF_sys.size() ? p.w.muonIDSF_sys[0][j] : 1.;
      double triggerSF=p.w.triggerSF_sys.size() ? p.w.triggerSF_sys[0][j] : 1.;
      p.weightmap[Form("_efficiencySF_stat%d",j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*electronRECOSF*electronIDSF*muonIDSF*p.w.muonISOSF*triggerSF*p.w.CFSF;
    }

    p.weightmap["_noelectronRECOSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
    for(int i=1,ni=p.w.electronRECOSF_sys.size();i<ni;i++){
      for(int j=0,nj=p.w.electronRECOSF_sys[i].size();j<nj;j++){
	p.weightmap[Form("_electronRECOSF_s%d_m%d",i,j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF_sys[i][j]*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
      }
    }	
    
    p.weightmap["_noelectronIDSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
    for(int i=1,ni=p.w.electronIDSF_sys.size();i<ni;i++){
      for(int j=0,nj=p.w.electronIDSF_sys[i].size();j<nj;j++){
	p.weightmap[Form("_electronIDSF_s%d_m%d",i,j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF_sys[i][j]*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
      }
    }	

    p.weightmap["_nomuonIDSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
    for(int i=1,ni=p.w.muonIDSF_sys.size();i<ni;i++){
      for(int j=0,nj=p.w.muonIDSF_sys[i].size();j<nj;j++){
	p.weightmap[Form("_muonIDSF_s%d_m%d",i,j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF_sys[i][j]*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
      }
    }	
    
    p.weightmap["_notriggerSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.CFSF;
    for(int i=1,ni=p.w.triggerSF_sys.size();i<ni;i++){
      for(int j=0,nj=p.w.triggerSF_sys[i].size();j<nj;j++){
	p.weightmap[Form("_triggerSF_s%d_m%d",i,j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF_sys[i][j]*p.w.CFSF;
      }
    }	
    
    p.weightmap["_noefficiencySF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*p.w.CFSF;
    
    p.weightmap["_noCFSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF;
  }
}
void EfficiencyValidation::FillHists(Parameter& p){
  TLorentzVector dilepton=(*p.lepton0)+(*p.lepton1);
  double dimass=dilepton.M();
  if(dimass>=52) FillHist(p.prefix+"m52to3000/"+p.hprefix+"dimass"+p.suffix,dimass,p.weightmap,mbinnum,mbin);
  if(dimass>=52&&dimass<150){
    FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"m52to150",p.weightmap[""]);
    FillHistsEfficiency(p,"m52to150/");
  }
  if(dimass>=80&&dimass<100){
    FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"m80to100",p.weightmap[""]);
    FillHistsEfficiency(p,"m80to100/");
  }
}

void EfficiencyValidation::FillHistsEfficiency(Parameter& p,TString region){
  TString pre=p.prefix+region+p.hprefix;
  for(const auto& [wname,w]:p.weightmap){
    TString suf=p.suffix+wname;
    
    //for leptons
    for(int i=0;i<(int)p.leptons.size();i++){
      double pt=p.leptons.at(i)->Pt();
      double eta=p.leptons.at(i)->Eta();
      TString charge=p.leptons.at(i)->Charge()>0?"p":"m";
      FillHist(Form("%sl%dpt%s",pre.Data(),i,suf.Data()),pt,w,500,0,500);
      FillHist(Form("%sl%deta%s",pre.Data(),i,suf.Data()),eta,w,120,-3,3);
      
      FillHist(Form("%slpt%s",pre.Data(),suf.Data()),pt,w,500,0,500);
      FillHist(Form("%sleta%s",pre.Data(),suf.Data()),eta,w,120,-3,3);

      if(p.leptons.at(i)->LeptonFlavour()==Lepton::Flavour::MUON){
      }else if(p.leptons.at(i)->LeptonFlavour()==Lepton::Flavour::ELECTRON){
	Electron* el=(Electron*)p.leptons.at(i);
	FillHist(Form("%slsceta%s",pre.Data(),suf.Data()),el->scEta(),w,120,-3,3);
	if(el->Pt()<20) FillHist(Form("%slsceta20%s",pre.Data(),suf.Data()),el->scEta(),w,120,-3,3);
      }
    }

    TLorentzVector dilepton=(*p.lepton0)+(*p.lepton1);
    double dimass=dilepton.M();
    double dipt=dilepton.Pt();
    double dirap=dilepton.Rapidity();
    FillHist(pre+"dimass"+suf,dimass,w,196,52,150);
    FillHist(pre+"dipt"+suf,dipt,w,400,0,400);
    FillHist(pre+"dirap"+suf,dirap,w,120,-3,3);
    FillHist(pre+"z0"+suf,vertex_Z,w,100,-20,20);

    if(wname!="") continue;

    //extra; without systematic
    for(int i=0;i<(int)p.leptons.size();i++){
      if(i>1) break;
      double pt=p.leptons.at(i)->Pt();
      double eta=p.leptons.at(i)->Eta();
      TString charge=p.leptons.at(i)->Charge()>0?"p":"m";

      FillHist(Form("%sl%d%spt%s",pre.Data(),i,charge.Data(),suf.Data()),pt,w,500,0,500);
      FillHist(Form("%sl%d%seta%s",pre.Data(),i,charge.Data(),suf.Data()),eta,w,120,-3,3);
      
      FillHist(Form("%sl%spt%s",pre.Data(),charge.Data(),suf.Data()),pt,w,500,0,500);
      FillHist(Form("%sl%seta%s",pre.Data(),charge.Data(),suf.Data()),eta,w,120,-3,3);
      
      FillHist(Form("%slriso%s",pre.Data(),suf.Data()),p.leptons.at(i)->RelIso(),w,30,0,0.3);
      if(p.leptons.at(i)->LeptonFlavour()==Lepton::Flavour::MUON){
	double rtrkiso=((Muon*)p.leptons.at(i))->TrkIso()/pt;
	FillHist(Form("%slrtrkiso%s",pre.Data(),suf.Data()),rtrkiso,w,40,0,0.2);
      }else if(p.leptons.at(i)->LeptonFlavour()==Lepton::Flavour::ELECTRON){
	Electron* el=(Electron*)p.leptons.at(i);
	FillHist(Form("%slrawpt%s",pre.Data(),suf.Data()),el->UncorrPt(),w,500,0,500);
      }
    }
      
    FillHist(pre+"nlepton"+suf,p.muons.size()+p.electrons.size(),w,10,0,10);
    FillHist(pre+"met"+suf,pfMET_Type1_pt,w,100,0,200);
    FillHist(Form("%snPV%s",pre.Data(),suf.Data()),nPV,w,100,0,100);
    if(!IsDATA) FillHist(Form("%snPileUp%s",pre.Data(),suf.Data()),nPileUp,w,100,0,100);

  }
}

bool EfficiencyValidation::PassSelection(Parameter& p){
  if(p.triggers.size()==2){
    if(GetEra()=="2017"&&p.triggers[0]=="HLT_Ele27_WPTight_Gsf_v"&&p.triggers[1]=="HLT_Ele32_WPTight_Gsf_v"){
      //if(!_event.PassTrigger("HLT_Ele27_WPTight_Gsf_v")) p.c.lepton0pt=35;
      if(p.lepton0&&p.lepton0->Pt()<35){
	p.triggers={"HLT_Ele27_WPTight_Gsf_v"};
	if(!_event.PassTrigger(p.triggers)) return false;
	p.k.triggerSF={"Ele27_MediumID"};
	if(!IsDATA) p.w.lumiweight*=31.72/41.54;
      }
    }
    if(GetEra()=="2018"&&p.triggers[0]=="HLT_Ele28_WPTight_Gsf_v"&&p.triggers[1]=="HLT_Ele32_WPTight_Gsf_v"){
      //if(!_event.PassTrigger("HLT_Ele28_WPTight_Gsf_v")) p.c.lepton0pt=35;
      if(p.lepton0&&p.lepton0->Pt()<35){
	p.triggers={"HLT_Ele28_WPTight_Gsf_v"};
	if(!_event.PassTrigger(p.triggers)) return false;
	p.k.triggerSF={"Ele28_MediumID"};
	if(!IsDATA) p.w.lumiweight*=23687.253/59827.879;
      }
    }
    if(GetEra()=="2017"&&p.triggers[0]=="HLT_IsoMu24_v"&&p.triggers[1]=="HLT_IsoMu27_v"){
      //if(!_event.PassTrigger("HLT_IsoMu24_v")) p.c.lepton0pt=30;
      if(p.lepton0&&p.lepton0->Pt()<30){
	p.triggers={"HLT_IsoMu24_v"};
	if(!_event.PassTrigger(p.triggers)) return false;
	p.k.triggerSF={"IsoMu24_MediumID_trkIsoLoose"};
	if(!IsDATA) p.w.lumiweight*=37997.005/41477.878;
      }
    }
  }
  return SMPAnalyzerCore::PassSelection(p);
}
