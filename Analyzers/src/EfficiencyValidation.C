#include "EfficiencyValidation.h"

EfficiencyValidation::EfficiencyValidation(){
}
EfficiencyValidation::~EfficiencyValidation(){
}
void EfficiencyValidation::initializeAnalyzer(){
  SMPAnalyzerCore::initializeAnalyzer(); //setup zpt roc z0
  fChain->SetBranchStatus("pfMET_*",false);
  fChain->SetBranchStatus("pfMET_Type1_pt",true);
  //fChain->SetBranchStatus("jet_*",false);
  fChain->SetBranchStatus("fatjet_*",false);
  //fChain->SetBranchStatus("photon_*",false);
}
void EfficiencyValidation::executeEvent(){
  //////// nominal channels //////////
  if(!IsDATA||DataStream.Contains("DoubleMuon")){
    if(HasFlag("ROCCOR")){
      for(int iset=0;iset<6;iset++){
	int nmem= iset==1 ? 20 : 1;
	for(int imem=0;imem<nmem;imem++)
	  executeEventWithParameter(MakeParameter("mm",Form("roccor_s%dm%d",iset,imem)));	
      } 
    }else{
      executeEventWithParameter(MakeParameter("mm")); 
      //executeEventWithParameter(MakeParameter("mm","mv17"));
    }
  }
  if(!IsDATA||DataStream.Contains("SingleMuon")){
    if(HasFlag("ROCCOR")){
      for(int iset=0;iset<6;iset++){
	int nmem= iset==1 ? 20 : 1;
	for(int imem=0;imem<nmem;imem++)
	  executeEventWithParameter(MakeParameter("mu",Form("roccor_s%dm%d",iset,imem)));	
      } 
    }else{
      executeEventWithParameter(MakeParameter("mu")); 
      //executeEventWithParameter(MakeParameter("mu","mv17"));
    }
  }
  if(!IsDATA||DataStream.Contains("DoubleEG")||DataStream.Contains("EGamma")){
    if(HasFlag("ROCCOR")){
      for(int iset=0;iset<6;iset++){
	int nmem= iset==1 ? 20 : 1;
	for(int imem=0;imem<nmem;imem++)
	  executeEventWithParameter(MakeParameter("ee",Form("roccor_s%dm%d",iset,imem)));	
      } 
    }else{
      executeEventWithParameter(MakeParameter("ee"));
      //executeEventWithParameter(MakeParameter("ee","ev17"));
    }
  }    
  if(!IsDATA||DataStream.Contains("SingleElectron")||DataStream.Contains("EGamma")){
    if(HasFlag("ROCCOR")){
      for(int iset=0;iset<6;iset++){
	int nmem= iset==1 ? 20 : 1;
	for(int imem=0;imem<nmem;imem++)
	  executeEventWithParameter(MakeParameter("el",Form("roccor_s%dm%d",iset,imem)));	
      } 
    }else{
      executeEventWithParameter(MakeParameter("el"));
      executeEventWithParameter(MakeParameter("el","SelQ"));
      //executeEventWithParameter(MakeParameter("el","ev17"));
      //executeEventWithParameter(MakeParameter("el","SelQ ev14"));
    }
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
      if(!IsDATA) p.w.lumiweight*=_event.GetTriggerLumi(p.triggers.at(0))/_event.GetTriggerLumi("Full");
      //executeEventWithParameter(p);

      p=MakeParameter("el");
      p.prefix="v12/el201727/";
      p.triggers={"HLT_Ele27_WPTight_Gsf_v"};
      p.k.triggerSF={"Ele27_MediumID_v12"};
      if(!IsDATA) p.w.lumiweight*=_event.GetTriggerLumi(p.triggers.at(0))/_event.GetTriggerLumi("Full");
      //executeEventWithParameter(p);

      p=MakeParameter("el");
      p.prefix="el201732/";
      p.triggers={"HLT_Ele32_WPTight_Gsf_v"};
      p.k.triggerSF={"Ele32_MediumID"};
      if(!IsDATA) p.w.lumiweight*=_event.GetTriggerLumi(p.triggers.at(0))/_event.GetTriggerLumi("Full");
      //executeEventWithParameter(p);

      p=MakeParameter("el");
      p.prefix="v12/el201732/";
      p.triggers={"HLT_Ele32_WPTight_Gsf_v"};
      p.k.triggerSF={"Ele32_MediumID_v12"};
      if(!IsDATA) p.w.lumiweight*=_event.GetTriggerLumi(p.triggers.at(0))/_event.GetTriggerLumi("Full");
      //executeEventWithParameter(p);
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
      if(!IsDATA) p.w.lumiweight*=_event.GetTriggerLumi(p.triggers.at(0))/_event.GetTriggerLumi("Full");
      //executeEventWithParameter(p);

      p=MakeParameter("el");
      p.prefix="v12/el201828/";
      p.triggers={"HLT_Ele28_WPTight_Gsf_v"};
      p.k.triggerSF={"Ele28_MediumID_v12"};
      p.k.electronIDSF+="_v12";
      if(!IsDATA) p.w.lumiweight*=_event.GetTriggerLumi(p.triggers.at(0))/_event.GetTriggerLumi("Full");
      //executeEventWithParameter(p);

      p=MakeParameter("el");
      p.prefix="el201832/";
      p.triggers={"HLT_Ele32_WPTight_Gsf_v"};
      p.k.triggerSF={"Ele32_MediumID"};
      //executeEventWithParameter(p);

      p=MakeParameter("el");
      p.prefix="v12/el201832/";
      p.triggers={"HLT_Ele32_WPTight_Gsf_v"};
      p.k.triggerSF={"Ele32_MediumID_v12"};
      p.k.electronIDSF+="_v12";
      //executeEventWithParameter(p);

    }
  }
}
SMPAnalyzerCore::Parameter EfficiencyValidation::MakeParameter(TString key,TString option){
  Parameter p=SMPAnalyzerCore::MakeParameter(key,option);
  p.weightbit|=EfficiencyWeight;
  if(option.Contains("ev17")){
    p.prefix="v17/"+p.prefix;
    p.k.electronIDSF+="_v17"; 
    p.k.electronRECOSF+="_v17";
    for(auto& key:p.k.triggerSF)
      key=key+"_v17";
    p.option.ReplaceAll("ev17","");
  }
  if(option.Contains("mv17")){
    p.prefix="v17/"+p.prefix;
    p.k.muonRECOSF+="_v17";
    p.k.muonTrackingSF+="_v17";
    p.k.muonIDSF+="_v17";
    for(auto& key:p.k.triggerSF) key=key+"_v17";
    p.option.ReplaceAll("mv17","");
  }
  if(option.Contains("noroccor")){
    p.suffix="_noroccor";
    if(p.channel[0]=='e'){
      p.SetElectrons(ElectronEnergyCorrection(SMPGetElectrons("passMediumID_SelQ",0.0,2.5),-1,0));
    }
    else if(p.channel[0]=='m'){
      p.SetMuons(MuonMomentumCorrection(SMPGetMuons("POGMediumWithLooseTrkIso",0.0,2.4),0,-1,0));
    }
  }
  return p;
}
void EfficiencyValidation::EvalWeights(Parameter& p){
  if(Userflags.size()==0||HasFlag("ROCCOR"))
    p.weightmap[""]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
  if(Userflags.size()==0){
    if(!IsDATA&&p.suffix==""&&!p.hprefix.Contains("ss_")){
      p.weightmap["_noweight"]=p.w.lumiweight;
      p.weightmap["_noPUweight"]=p.w.lumiweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
      p.weightmap["_PUweight_up"]=p.w.lumiweight*p.w.PUweight_up*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
      p.weightmap["_PUweight_down"]=p.w.lumiweight*p.w.PUweight_down*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
      p.weightmap["_noprefireweight"]=p.w.lumiweight*p.w.PUweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
      p.weightmap["_prefireweight_up"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight_up*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
      p.weightmap["_prefireweight_down"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight_down*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
      
      //// For prefireweight study
      //p.weightmap["_newprefireweight1"]=p.w.lumiweight*p.w.PUweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*GetL1PrefiringWeight(1);
      //p.weightmap["_newprefireweight2"]=p.w.lumiweight*p.w.PUweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*GetL1PrefiringWeight(2);
      p.weightmap["_newprefireweight3"]=p.w.lumiweight*p.w.PUweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*GetL1PrefiringWeight(3);
      p.weightmap["_newprefireweight4"]=p.w.lumiweight*p.w.PUweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*GetL1PrefiringWeight(4);
      p.weightmap["_newprefireweight5"]=p.w.lumiweight*p.w.PUweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*GetL1PrefiringWeight(5);
      p.weightmap["_newprefireweight6"]=p.w.lumiweight*p.w.PUweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*GetL1PrefiringWeight(6);
      
      p.weightmap["_nozptweight"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
      p.weightmap["_z0weight"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.w.z0weight;
      p.weightmap["_noweakweight"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
      
      for(int j=0,nj=fEff->nreplica;j<nj;j++){
	double electronRECOSF=p.w.electronRECOSF_sys.size() ? p.w.electronRECOSF_sys[0][j] : 1.;
	double electronIDSF=p.w.electronIDSF_sys.size() ? p.w.electronIDSF_sys[0][j] : 1.;
	double muonTrackingSF=p.w.muonTrackingSF_sys.size() ? p.w.muonTrackingSF_sys[0][j] : 1.;
	double muonRECOSF=p.w.muonRECOSF_sys.size() ? p.w.muonRECOSF_sys[0][j] : 1.;
	double muonIDSF=p.w.muonIDSF_sys.size() ? p.w.muonIDSF_sys[0][j] : 1.;
	double triggerSF=p.w.triggerSF_sys.size() ? p.w.triggerSF_sys[0][j] : 1.;
	p.weightmap[Form("_efficiencySF_stat%d",j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*electronRECOSF*electronIDSF*muonTrackingSF*muonRECOSF*muonIDSF*p.w.muonISOSF*triggerSF*p.w.CFSF;
      }
      
      p.weightmap["_noelectronRECOSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
      for(int i=1,ni=p.w.electronRECOSF_sys.size();i<ni;i++){
	for(int j=0,nj=p.w.electronRECOSF_sys[i].size();j<nj;j++){
	  p.weightmap[Form("_electronRECOSF_s%d_m%d",i,j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF_sys[i][j]*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
	}
      }	
      
      p.weightmap["_noelectronIDSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
      for(int i=1,ni=p.w.electronIDSF_sys.size();i<ni;i++){
	for(int j=0,nj=p.w.electronIDSF_sys[i].size();j<nj;j++){
	  p.weightmap[Form("_electronIDSF_s%d_m%d",i,j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF_sys[i][j]*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
	}
      }	
      
      p.weightmap["_nomuonTrackingSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
      for(int i=1,ni=p.w.muonTrackingSF_sys.size();i<ni;i++){
	for(int j=0,nj=p.w.muonTrackingSF_sys[i].size();j<nj;j++){
	  p.weightmap[Form("_muonTrackingSF_s%d_m%d",i,j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF_sys[i][j]*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
	}
      }	
      
      p.weightmap["_nomuonRECOSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
      for(int i=1,ni=p.w.muonRECOSF_sys.size();i<ni;i++){
	for(int j=0,nj=p.w.muonRECOSF_sys[i].size();j<nj;j++){
	  p.weightmap[Form("_muonRECOSF_s%d_m%d",i,j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF_sys[i][j]*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
	}
      }	
      
      p.weightmap["_nomuonIDSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
      for(int i=1,ni=p.w.muonIDSF_sys.size();i<ni;i++){
	for(int j=0,nj=p.w.muonIDSF_sys[i].size();j<nj;j++){
	  p.weightmap[Form("_muonIDSF_s%d_m%d",i,j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF_sys[i][j]*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
	}
      }	
      
      p.weightmap["_notriggerSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.CFSF;
      p.weightmap["_triggerSF_mode1"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF_mode1*p.w.CFSF;
      p.weightmap["_triggerSF_interpolation"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF_interpolation*p.w.CFSF;
      for(int i=1,ni=p.w.triggerSF_sys.size();i<ni;i++){
	for(int j=0,nj=p.w.triggerSF_sys[i].size();j<nj;j++){
	  p.weightmap[Form("_triggerSF_s%d_m%d",i,j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF_sys[i][j]*p.w.CFSF;
	}
      }
      
      if(p.channel=="ee"||p.channel=="mm"){
	if(p.k.triggerSF.size()==2&&GetPtThreshold(p.k.triggerSF[0])>GetPtThreshold(p.k.triggerSF[1])){
	  vector<Lepton*> leps;
	  if(p.channel=="ee") leps=MakeLeptonPointerVector(p.electrons);
	  else if(p.channel=="mm") leps=MakeLeptonPointerVector(p.muons);
	  p.weightmap["_noDZSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.CFSF*GetDileptonTriggerSF(p.k.triggerSF[0],p.k.triggerSF[1],"",leps,0,0);
	}
      }
      
      p.weightmap["_noefficiencySF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.CFSF;
      
      p.weightmap["_noCFSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF;
    }
  }
  if(HasFlag("PDFSYS")){
    if(!IsDATA&&p.suffix==""&&!p.hprefix.Contains("ss_")){
      for(unsigned int i=0;i<weight_Scale->size();i++){
	p.weightmap[Form("_scalevariation%d",i)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*weight_Scale->at(i);
      }
      for(unsigned int i=0;i<weight_PDF->size();i++){
	p.weightmap[Form("_pdf%d",i)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*weight_PDF->at(i);
      }
      if(weight_AlphaS->size()==2){
	p.weightmap["_alphaS_down"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*weight_AlphaS->at(0);
	p.weightmap["_alphaS_up"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*weight_AlphaS->at(1);
      }
    }
  }
}
void EfficiencyValidation::FillHists(Parameter& p){
  TLorentzVector dilepton=(*p.lepton0)+(*p.lepton1);
  double dimass=dilepton.M();
  if(dimass>=52) FillHist(p.prefix+"m52to3000/"+p.hprefix+"dimass"+p.suffix,dimass,p.weightmap,mbinnum,mbin);
  if(dimass>=52&&dimass<150){
    if(!Userflags.size()) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"m52to150",p.weightmap[""]);
    FillHistsEfficiency(p,"m52to150/");
    FillHist(p.prefix+p.hprefix+"costym"+p.suffix,fabs(GetCosThetaCS(p.lepton0,p.lepton1)),fabs(dilepton.Rapidity()),dimass,p.weightmap,rochester_ncostbin,rochester_costbins,rochester_nybin,rochester_ybins,rochester_nmbin,rochester_mbins);
  }
  if(dimass>=80&&dimass<100){
    if(!Userflags.size()) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"m80to100",p.weightmap[""]);
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
      //// For prefirieweight check
      if(eta<-2.0) FillHist(Form("%slpt_eta0%s",pre.Data(),suf.Data()),pt,w,100,0,100);
      else if(eta<2.0) FillHist(Form("%slpt_eta1%s",pre.Data(),suf.Data()),pt,w,100,0,100);
      else FillHist(Form("%slpt_eta2%s",pre.Data(),suf.Data()),pt,w,100,0,100);
	
      FillHist(Form("%sleta%s",pre.Data(),suf.Data()),eta,w,120,-3,3);

      if(p.leptons.at(i)->LeptonFlavour()==Lepton::Flavour::MUON){
	FillHist(Form("%sletapt%s",pre.Data(),suf.Data()),eta,pt,w,netabin_muonID,etabins_muonID,nptbin_muonID,ptbins_muonID);
      }else if(p.leptons.at(i)->LeptonFlavour()==Lepton::Flavour::ELECTRON){
	Electron* el=(Electron*)p.leptons.at(i);
	FillHist(Form("%slsceta%s",pre.Data(),suf.Data()),el->scEta(),w,120,-3,3);
	FillHist(Form("%sletapt%s",pre.Data(),suf.Data()),eta,pt,w,netabin_electronID,etabins_electronID,nptbin_electronID,ptbins_electronID);
      }
    }

    TLorentzVector dilepton=(*p.lepton0)+(*p.lepton1);
    double dimass=dilepton.M();
    double dipt=dilepton.Pt();
    double dirap=dilepton.Rapidity();
    FillHist(pre+"dimass"+suf,dimass,w,196,52,150);
    FillHist(pre+"dipt"+suf,dipt,w,400,0,400);
    FillHist(pre+"dirap"+suf,dirap,w,120,-3,3);
    FillHist(pre+"cost"+suf,GetCosThetaCS(p.lepton0,p.lepton1),w,100,-1,1);
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
	//Electron* el=(Electron*)p.leptons.at(i);
	//FillHist(Form("%slrawpt%s",pre.Data(),suf.Data()),el->UncorrPt(),w,500,0,500);
	//FillProfile(pre+"leta_energyscale"+suf,el->scEta(),el->Pt()/el->UncorrPt(),w,60,-3,3);
	//FillProfile(pre+"leta_energyscale2"+suf,el->scEta(),el->Energy()/el->scE(),w,60,-3,3);
      }
      //int iphi=(p.leptons.at(i)->Phi()/TMath::Pi()+1)*4;      
      //FillHist(Form("%sleta%s_phi%d",pre.Data(),suf.Data(),iphi),eta,w,120,-3,3);
      if(vertex_Z<-4){
	//FillHist(Form("%sleta%s_z0",pre.Data(),suf.Data()),eta,w,120,-3,3);
      }else if(vertex_Z<4){
	//FillHist(Form("%sleta%s_z1",pre.Data(),suf.Data()),eta,w,120,-3,3);
      }else{
	//FillHist(Form("%sleta%s_z2",pre.Data(),suf.Data()),eta,w,120,-3,3);
      }
    }
    //if(p.truth_lepton0.Pt()) FillProfile(pre+"leta_deta"+suf,p.lepton0->Eta(),p.lepton0->Eta()-p.truth_lepton0.Eta(),w,60,-3,3);
    //if(p.truth_lepton1.Pt()) FillProfile(pre+"leta_deta"+suf,p.lepton1->Eta(),p.lepton1->Eta()-p.truth_lepton1.Eta(),w,60,-3,3);    

    Lepton *lp=NULL,*lm=NULL;
    if(p.lepton0->Charge()>0){ 
      lp=p.lepton0; 
      lm=p.lepton1; 
    }else{
      lp=p.lepton1; 
      lm=p.lepton0;
    }

    if(lp->LeptonFlavour()==Lepton::Flavour::MUON){
      FillHist(Form("%slpetaptlmetapt%s",pre.Data(),suf.Data()),lp->Eta(),lp->Pt(),lm->Eta(),lm->Pt(),w,netabin_muonID,etabins_muonID,nptbin_muonID,ptbins_muonID,netabin_muonID,etabins_muonID,nptbin_muonID,ptbins_muonID);
    }else if(lp->LeptonFlavour()==Lepton::Flavour::ELECTRON){
      FillHist(Form("%slpetaptlmetapt%s",pre.Data(),suf.Data()),lp->Eta(),lp->Pt(),lm->Eta(),lm->Pt(),w,netabin_electronID,etabins_electronID,nptbin_electronID,ptbins_electronID,netabin_electronID,etabins_electronID,nptbin_electronID,ptbins_electronID);
    }
      
    FillHist(pre+"x0"+suf,vertex_X,w,100,-0.2,0.2);
    FillHist(pre+"y0"+suf,vertex_Y,w,100,-0.2,0.2);
    FillHist(pre+"r0"+suf,sqrt(vertex_X*vertex_X+vertex_Y*vertex_Y),w,100,-0.2,0.2);

    FillHist(pre+"nlepton"+suf,p.muons.size()+p.electrons.size(),w,10,0,10);
    FillHist(pre+"met"+suf,pfMET_Type1_pt,w,100,0,200);
    FillHist(Form("%snPV%s",pre.Data(),suf.Data()),nPV,w,100,0,100);
    if(!IsDATA) FillHist(Form("%snPileUp%s",pre.Data(),suf.Data()),nPileUp,w,100,0,100);

  }
}

bool EfficiencyValidation::PassSelection(Parameter& p){
  /*
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
      //if(p.lepton0&&p.lepton0->Pt()<35){
      //p.triggers={"HLT_Ele28_WPTight_Gsf_v"};
      //if(!_event.PassTrigger(p.triggers)) return false;
      //p.k.triggerSF={"Ele28_MediumID"};
      //if(!IsDATA) p.w.lumiweight*=23687.253/59827.879;
      //}
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
  */
  return SMPAnalyzerCore::PassSelection(p);
}

double EfficiencyValidation::GetLeptonTriggerORSF_old(const Parameter& p,int set,int mem){
  if(IsDATA) return 1;

  vector<Lepton*> leps;
  if(p.channel=="el") leps=MakeLeptonPointerVector(p.electrons);
  else if(p.channel=="mu") leps=MakeLeptonPointerVector(p.muons);
  TString triggerSF_key0=p.k.triggerSF[0];
  TString triggerSF_key1=p.k.triggerSF[1];

  double lumi0=1.; //trigger0 on
  double lumi1=0.; //only trigger1 on
  double lumi2=0.; //both off
  if(DataYear==2017&&triggerSF_key0.Contains("IsoMu24")&&triggerSF_key1.Contains("IsoMu27")){
    lumi0=37997.005;    lumi1=3480.873;    lumi2=0.;
  }else if(DataYear==2017&&triggerSF_key0.Contains("Ele27")&&triggerSF_key1.Contains("Ele32")){
    lumi0=31661.026;    lumi1=9522.208;    lumi2=295.;
  }else if(DataYear==2018&&triggerSF_key0.Contains("Ele28")&&triggerSF_key1.Contains("Ele32")){
    lumi0=23687.253;    lumi1=36140.626;   lumi2=0.;
  }else{
    cout<<"[EfficiencyValidation::GetLeptonTriggerORSF_old] not available combination '"<<triggerSF_key0<<"'||'"<<triggerSF_key1<<"' for "<<DataEra<<endl;
    exit(EXIT_FAILURE);
  }
  
  double data_eff_key0=1.,data_eff_key1=1.,sim_eff=1.;
  for(const auto& lep:leps){
    data_eff_key0*=1-fEff->GetDataEfficiency(triggerSF_key0,lep,set,mem);
    data_eff_key1*=1-fEff->GetDataEfficiency(triggerSF_key1,lep,set,mem);
    sim_eff*=1-fEff->GetSimEfficiency(triggerSF_key0,lep,set,mem);
  }
  data_eff_key0=1-data_eff_key0;
  data_eff_key1=1-data_eff_key1;
  sim_eff=1-sim_eff;
  if(sim_eff==0) return 1.;
  else return (lumi0*data_eff_key0+lumi1*data_eff_key1)/((lumi0+lumi1+lumi2)*sim_eff);
}

double EfficiencyValidation::GetCosThetaCS(const Particle *p0,const Particle *p1,int direction) const {
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
