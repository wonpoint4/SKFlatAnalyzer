#include "EfficiencyValidation.h"

EfficiencyValidation::EfficiencyValidation(){
}
EfficiencyValidation::~EfficiencyValidation(){
}
void EfficiencyValidation::initializeAnalyzer(){
  SMPAnalyzerCore::initializeAnalyzer(); //setup zpt roc z0
  fChain->SetBranchStatus("pfMET_*",false);
  fChain->SetBranchStatus("jet_*",false);
  fChain->SetBranchStatus("fatjet_*",false);
  fChain->SetBranchStatus("photon_*",false);
}
void EfficiencyValidation::executeEvent(){
  event=GetEvent();
  GetEventWeights();

  if(!PassMETFilter()) return;

  map<TString,vector<Muon>> map_muons;
  map_muons["MediumID_trkIsoLoose"]=SMPGetMuons("POGMediumWithLooseTrkIso",0.0,2.4);

  map<TString,vector<Electron>> map_electrons;
  map_electrons["MediumID"]=SMPGetElectrons("passMediumID",0.0,2.4);
  map_electrons["TightID_Selective"]=SMPGetElectrons("passTightID_Selective",0.0,2.4);
  //map_electrons["MediumID_noroccor"]=SMPGetElectrons("passMediumID",0.0,2.4);
  //map_electrons["MediumID"]=ElectronEnergyCorrection(map_electrons["MediumID_noroccor"],0,0);
  //map_electrons["TightID_Selective_noroccor"]=SMPGetElectrons("passTightID_Selective",0.0,2.4);
  //map_electrons["TightID_Selective"]=ElectronEnergyCorrection(map_electrons["TightID_Selective_noroccor"],0,0);

  if(GetEra()=="2016preVFP"){
    vector<TString> doublemuontrigger={
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v",
      "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
      "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
      "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",
      "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
    };
    if(event.PassTrigger(doublemuontrigger))
      if(!IsDATA||DataStream.Contains("DoubleMuon"))
	executeEventWithParameter(Parameter("mm2016a/","",
					    "IDISO_SF_MediumID_trkIsoLoose_Q","",{"Mu17Leg1_MediumID_trkIsoLoose_Q","Mu8Leg2_MediumID_trkIsoLoose_Q"},
					    20,10,MakeLeptonPointerVector(map_muons["MediumID_trkIsoLoose"])));
    if(event.PassTrigger("HLT_IsoMu24_v")||event.PassTrigger("HLT_IsoTkMu24_v"))
      if(!IsDATA||DataStream.Contains("SingleMuon")) 
	executeEventWithParameter(Parameter("mu2016a/","",
					    "IDISO_SF_MediumID_trkIsoLoose_Q","",{"IsoMu24_MediumID_trkIsoLoose_Q"},
					    25,10,MakeLeptonPointerVector(map_muons["MediumID_trkIsoLoose"])));
    if(event.PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"))
      if(!IsDATA||DataStream.Contains("DoubleEG"))
	executeEventWithParameter(Parameter("ee2016a/","",
					    "ID_SF_MediumID_Q",{"Ele23Leg1_MediumID_Q","Ele12Leg2_MediumID_Q"},
					    25,15,MakeLeptonPointerVector(map_electrons["MediumID"])));
    if(event.PassTrigger("HLT_Ele27_WPTight_Gsf_v"))
      if(!IsDATA||DataStream.Contains("SingleElectron")){
	executeEventWithParameter(Parameter("el2016a/","",
					    "ID_SF_MediumID_Q",{"Ele27_MediumID_Q"},
					    30,10,MakeLeptonPointerVector(map_electrons["MediumID"])));
	executeEventWithParameter(Parameter("el2016a/","_ptlt20",
					    "ID_SF_MediumID_Q",{"Ele27_MediumID_Q"},
					    30,10,MakeLeptonPointerVector(map_electrons["MediumID"])));
	executeEventWithParameter(Parameter("el2016a/","_TightID_Selective",
					    "ID_SF_TightID_Selective_Q",{"Ele27_TightID_Selective_Q"},
					    30,10,MakeLeptonPointerVector(map_electrons["TightID_Selective"])));
      }
  }else if(GetEra()=="2016postVFP"){
    vector<TString> doublemuontrigger={
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v",
      "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
      "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
      "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",
      "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
    };
    if(event.PassTrigger(doublemuontrigger))
      if(!IsDATA||DataStream.Contains("DoubleMuon"))
	executeEventWithParameter(Parameter("mm2016b/","",
					    "IDISO_SF_MediumID_trkIsoLoose_Q","",{"Mu17Leg1_MediumID_trkIsoLoose_Q","Mu8Leg2_MediumID_trkIsoLoose_Q"},
					    20,10,MakeLeptonPointerVector(map_muons["MediumID_trkIsoLoose"])));
    if(event.PassTrigger("HLT_IsoMu24_v")||event.PassTrigger("HLT_IsoTkMu24_v"))
      if(!IsDATA||DataStream.Contains("SingleMuon")) 
	executeEventWithParameter(Parameter("mu2016b/","",
					    "IDISO_SF_MediumID_trkIsoLoose_Q","",{"IsoMu24_MediumID_trkIsoLoose_Q"},
					    25,10,MakeLeptonPointerVector(map_muons["MediumID_trkIsoLoose"])));
    if(event.PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"))
      if(!IsDATA||DataStream.Contains("DoubleEG"))
	executeEventWithParameter(Parameter("ee2016b/","",
					    "ID_SF_MediumID_Q",{"Ele23Leg1_MediumID_Q","Ele12Leg2_MediumID_Q"},
					    25,15,MakeLeptonPointerVector(map_electrons["MediumID"])));
    if(event.PassTrigger("HLT_Ele27_WPTight_Gsf_v"))
      if(!IsDATA||DataStream.Contains("SingleElectron")){
	executeEventWithParameter(Parameter("el2016b/","",
					    "ID_SF_MediumID_Q",{"Ele27_MediumID_Q"},
					    30,10,MakeLeptonPointerVector(map_electrons["MediumID"])));
	executeEventWithParameter(Parameter("el2016b/","_ptlt20",
					    "ID_SF_MediumID_Q",{"Ele27_MediumID_Q"},
					    30,10,MakeLeptonPointerVector(map_electrons["MediumID"])));
	executeEventWithParameter(Parameter("el2016b/","_TightID_Selective",
					    "ID_SF_TightID_Selective_Q",{"Ele27_TightID_Selective_Q"},
					    30,10,MakeLeptonPointerVector(map_electrons["TightID_Selective"])));
      }
  }else if(GetEra()=="2017"){
    if(event.PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v"))
      if(!IsDATA||DataStream.Contains("DoubleMuon"))
	executeEventWithParameter(Parameter("mm2017/","",
					    "IDISO_SF_MediumID_trkIsoLoose_Q","",{"Mu17Leg1_MediumID_trkIsoLoose_Q","Mu8Leg2_MediumID_trkIsoLoose_Q"},
					    20,10,MakeLeptonPointerVector(map_muons["MediumID_trkIsoLoose"])));	
    if(event.PassTrigger("HLT_IsoMu24_v")||event.PassTrigger("HLT_IsoMu27_v"))
      if(!IsDATA||DataStream.Contains("SingleMuon"))
	executeEventWithParameter(Parameter("mu2017/","",
					    "IDISO_SF_MediumID_trkIsoLoose_Q","",{"IsoMu24_MediumID_trkIsoLoose_Q","IsoMu27_MediumID_trkIsoLoose_Q"},
					    25,10,MakeLeptonPointerVector(map_muons["MediumID_trkIsoLoose"])));
    if(event.PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"))
      if(!IsDATA||DataStream.Contains("DoubleEG")){
	executeEventWithParameter(Parameter("ee2017/","",
					    "ID_SF_MediumID_Q",{"Ele23Leg1_MediumID_Q","Ele12Leg2_MediumID_Q"},
					    25,15,MakeLeptonPointerVector(map_electrons["MediumID"])));
	executeEventWithParameter(Parameter("ee2017/","_noL1",
					    "ID_SF_MediumID_Q",{"Ele23Leg1_MediumID_Q_v3_2","Ele12Leg2_MediumID_Q"},
					    25,15,MakeLeptonPointerVector(map_electrons["MediumID"])));	
      }	
    if(event.PassTrigger("HLT_Ele27_WPTight_Gsf_v")||event.PassTrigger("HLT_Ele32_WPTight_Gsf_v"))
      if(!IsDATA||DataStream.Contains("SingleElectron")){
	executeEventWithParameter(Parameter("el2017/","",
					    "ID_SF_MediumID_Q",{"Ele27_MediumID_Q","Ele32_MediumID_Q"},
					    30,10,MakeLeptonPointerVector(map_electrons["MediumID"])));
	executeEventWithParameter(Parameter("el2017/","_ptlt20",
					    "ID_SF_MediumID_Q",{"Ele27_MediumID_Q","Ele32_MediumID_Q"},
					    30,10,MakeLeptonPointerVector(map_electrons["MediumID"])));
	executeEventWithParameter(Parameter("el2017/","_TightID_Selective",
					    "ID_SF_TightID_Selective_Q",{"Ele27_TightID_Selective_Q","Ele32_TightID_Selective_Q"},
					    30,10,MakeLeptonPointerVector(map_electrons["TightID_Selective"])));
      }
  }else if(GetEra()=="2018"){
    if(event.PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"))
      if(!IsDATA||DataStream.Contains("DoubleMuon"))
	executeEventWithParameter(Parameter("mm2018/","",
					    "IDISO_SF_MediumID_trkIsoLoose_Q","",{"Mu17Leg1_MediumID_trkIsoLoose_Q","Mu8Leg2_MediumID_trkIsoLoose_Q"},
					    20,10,MakeLeptonPointerVector(map_muons["MediumID_trkIsoLoose"])));
    if(event.PassTrigger("HLT_IsoMu24_v"))
      if(!IsDATA||DataStream.Contains("SingleMuon"))
	executeEventWithParameter(Parameter("mu2018/","",
					    "IDISO_SF_MediumID_trkIsoLoose_Q","",{"IsoMu24_MediumID_trkIsoLoose_Q"},
					    25,10,MakeLeptonPointerVector(map_muons["MediumID_trkIsoLoose"])));
    if(event.PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"))
      if(!IsDATA||DataStream.Contains("EGamma")){
	executeEventWithParameter(Parameter("ee2018/","",
					    "ID_SF_MediumID_Q",{"Ele23Leg1_MediumID_Q","Ele12Leg2_MediumID_Q"},
					    25,15,MakeLeptonPointerVector(map_electrons["MediumID"])));
	executeEventWithParameter(Parameter("ee2018/","_noL1",
					    "ID_SF_MediumID_Q",{"Ele23Leg1_MediumID_Q_v3","Ele12Leg2_MediumID_Q"},
					    25,15,MakeLeptonPointerVector(map_electrons["MediumID"])));
      }
    if(event.PassTrigger("HLT_Ele28_WPTight_Gsf_v")||event.PassTrigger("HLT_Ele32_WPTight_Gsf_v"))
      if(!IsDATA||DataStream.Contains("EGamma")){
	executeEventWithParameter(Parameter("el2018/","",
					    "ID_SF_MediumID_Q",{"Ele28_MediumID_Q","Ele32_MediumID_Q"},
					    30,10,MakeLeptonPointerVector(map_electrons["MediumID"])));
	executeEventWithParameter(Parameter("el2018/","_ptlt20",
					    "ID_SF_MediumID_Q",{"Ele28_MediumID_Q","Ele32_MediumID_Q"},
					    30,10,MakeLeptonPointerVector(map_electrons["MediumID"])));
	executeEventWithParameter(Parameter("el2018/","_TightID_Selective",
					    "ID_SF_TightID_Selective_Q",{"Ele28_TightID_Selective_Q","Ele32_TightID_Selective_Q"},
					    30,10,MakeLeptonPointerVector(map_electrons["TightID_Selective"])));
      }
  }    
}

void EfficiencyValidation::executeEventWithParameter(Parameter p){
  double prescaleweight=1.;

  TString hprefix=tauprefix;
  FillCutflow(p.prefix+hprefix+"cutflow"+p.suffix,"lumi",lumiweight);
  FillCutflow(p.prefix+hprefix+"cutflow"+p.suffix,"PU",lumiweight*PUweight);
  FillCutflow(p.prefix+hprefix+"cutflow"+p.suffix,"prefire",lumiweight*PUweight*prefireweight);
  FillCutflow(p.prefix+hprefix+"cutflow"+p.suffix,"zpt",lumiweight*PUweight*prefireweight*zptweight);
  FillCutflow(p.prefix+hprefix+"cutflow"+p.suffix,"z0",lumiweight*PUweight*prefireweight*zptweight*z0weight);

  double eventweight=lumiweight*PUweight*prefireweight*z0weight*zptweight;
  FillHist(p.prefix+hprefix+"nlepton"+p.suffix,p.leps.size(),eventweight,10,0,10);

  ///////////////////////lepton selection///////////////////////
  if(p.leps.size()>=2){
    if(p.leps.at(0)->Charge()*p.leps.at(1)->Charge()>0)
      hprefix="ss_"+hprefix; //continue;
    FillCutflow(p.prefix+hprefix+"cutflow"+p.suffix,"OS dilepton",eventweight);
    if(p.leps.at(0)->Pt()>p.lep0ptcut&&p.leps.at(1)->Pt()>p.lep1ptcut){
      FillCutflow(p.prefix+hprefix+"cutflow"+p.suffix,"LepPtCut",eventweight);
      /////////////////efficiency scale factors///////////////////
      double IDSF=1.,IDSF_up=1.,IDSF_down=1.;
      double ISOSF=1.,ISOSF_up=1.,ISOSF_down=1.;
      double RECOSF=1.,RECOSF_up=1.,RECOSF_down=1.;
      if(!IsDATA){
	for(const auto& lep:p.leps){	  
	  TString LeptonIDSF_key="";
	  if(lep->LeptonFlavour()==Lepton::ELECTRON){
	    LeptonIDSF_key=p.electronIDSF;

	    double this_pt,this_eta;
	    this_pt=((Electron*)lep)->UncorrPt();
	    this_eta=((Electron*)lep)->scEta();
	    
	    double this_pt_recosf=(!p.suffix.Contains("_ptlt20")&&this_pt<20)?20.1:this_pt;
	    double this_RECOSF=mcCorr->ElectronReco_SF(this_eta,this_pt_recosf,0);
	    double this_RECOSF_up=mcCorr->ElectronReco_SF(this_eta,this_pt_recosf,1);
	    double this_RECOSF_down=mcCorr->ElectronReco_SF(this_eta,this_pt_recosf,-1);
	    RECOSF*=this_RECOSF; RECOSF_up*=this_RECOSF_up; RECOSF_down*=this_RECOSF_down;
	  }else if(lep->LeptonFlavour()==Lepton::MUON){
	    LeptonIDSF_key=p.muonIDSF;

	    double this_ISOSF=Lepton_SF(p.muonISOSF,lep,0);
	    double this_ISOSF_up=Lepton_SF(p.muonISOSF,lep,1);
	    double this_ISOSF_down=Lepton_SF(p.muonISOSF,lep,-1);
	    ISOSF*=this_ISOSF; ISOSF_up*=this_ISOSF_up; ISOSF_down*=this_ISOSF_down;
	  }
	  
	  double this_IDSF=Lepton_SF(LeptonIDSF_key,lep,0);
	  double this_IDSF_up=Lepton_SF(LeptonIDSF_key,lep,1);
	  double this_IDSF_down=Lepton_SF(LeptonIDSF_key,lep,-1);
	  IDSF*=this_IDSF; IDSF_up*=this_IDSF_up; IDSF_down*=this_IDSF_down;
	  
	}
      }
      
      double triggerSF=1.,triggerSF_up=1.,triggerSF_down=1.;
      if(!IsDATA){
	if(p.triggerSF.size()==1){
	  triggerSF*=LeptonTrigger_SF(p.triggerSF[0],p.leps,0);
	  triggerSF_up*=LeptonTrigger_SF(p.triggerSF[0],p.leps,1);
	  triggerSF_down*=LeptonTrigger_SF(p.triggerSF[0],p.leps,-1);
	}else if(p.triggerSF.size()==2){
	  if(p.prefix.Contains("mu")||p.prefix.Contains("el")){
	    triggerSF*=LeptonTrigger_SF_OR(p.triggerSF[0],p.triggerSF[1],p.leps,0);
	    triggerSF_up*=LeptonTrigger_SF_OR(p.triggerSF[0],p.triggerSF[1],p.leps,1);
	    triggerSF_down*=LeptonTrigger_SF_OR(p.triggerSF[0],p.triggerSF[1],p.leps,-1);
	  }else{
	    triggerSF*=DileptonTrigger_SF(p.triggerSF[0],p.triggerSF[1],p.leps,0);
	    triggerSF_up*=DileptonTrigger_SF(p.triggerSF[0],p.triggerSF[1],p.leps,1);
	    triggerSF_down*=DileptonTrigger_SF(p.triggerSF[0],p.triggerSF[1],p.leps,-1);
	  }
	}
      }
      
      FillCutflow(p.prefix+hprefix+"cutflow"+p.suffix,"RECOSF",eventweight*RECOSF);
      FillCutflow(p.prefix+hprefix+"cutflow"+p.suffix,"IDSF",eventweight*RECOSF*IDSF);
      FillCutflow(p.prefix+hprefix+"cutflow"+p.suffix,"ISOSF",eventweight*RECOSF*IDSF*ISOSF);
      FillCutflow(p.prefix+hprefix+"cutflow"+p.suffix,"triggerSF",eventweight*RECOSF*IDSF*ISOSF*triggerSF);
      
      ///////////////////////weight systematics//////////////////
      map<TString,double> map_weight;
      map_weight[""]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF*prescaleweight;
      if(!IsDATA){
	map_weight["_noweight"]=lumiweight;
	map_weight["_noPUweight"]=lumiweight*prefireweight*zptweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF*prescaleweight;
	map_weight["_noprefireweight"]=lumiweight*PUweight*RECOSF*IDSF*ISOSF*triggerSF*zptweight*z0weight*prescaleweight;
	
	map_weight["_noRECOSF"]=lumiweight*PUweight*IDSF*ISOSF*triggerSF*prefireweight*zptweight*z0weight*prescaleweight;
	map_weight["_RECOSF_up"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF_up*IDSF*ISOSF*triggerSF*prescaleweight;
	map_weight["_RECOSF_down"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF_down*IDSF*ISOSF*triggerSF*prescaleweight;
	
	map_weight["_noIDSF"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF*ISOSF*triggerSF*prescaleweight;
	map_weight["_IDSF_up"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF*IDSF_up*ISOSF*triggerSF*prescaleweight;
	map_weight["_IDSF_down"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF*IDSF_down*ISOSF*triggerSF*prescaleweight;
	
	map_weight["_noISOSF"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF*IDSF*triggerSF*prescaleweight;
	map_weight["_ISOSF_up"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF*IDSF*ISOSF_up*triggerSF*prescaleweight;
	map_weight["_ISOSF_down"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF*IDSF*ISOSF_down*triggerSF*prescaleweight;
	
	map_weight["_notriggerSF"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF*IDSF*ISOSF*prescaleweight;
	map_weight["_triggerSF_up"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF_up*prescaleweight;
	map_weight["_triggerSF_down"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF_down*prescaleweight;
	
	map_weight["_noefficiencySF"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*prescaleweight;
	map_weight["_noz0weight"]=lumiweight*PUweight*prefireweight*zptweight*RECOSF*IDSF*ISOSF*triggerSF*prescaleweight;
	map_weight["_nozptweight"]=lumiweight*PUweight*prefireweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF*prescaleweight;
	
	map_weight["_noprescaleweight"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF;
      }
      
      ///////////////////////fill hists///////////////////////
      TLorentzVector dilepton=(*p.leps.at(0))+(*p.leps.at(1));
      double dimass=dilepton.M();
      if(dimass>=60&&dimass<120){
	FillCutflow(p.prefix+hprefix+"cutflow"+p.suffix,"m60to120",map_weight[""]);
	FillHistsEfficiency(p.prefix+"m60to120/"+hprefix,p.suffix,p.leps,map_weight);
	if(dimass>=80&&dimass<100){
	  FillCutflow(p.prefix+hprefix+"cutflow"+p.suffix,"m80to100",map_weight[""]);
	  FillHistsEfficiency(p.prefix+"m80to100/"+hprefix,p.suffix,p.leps,map_weight);
	}
      }
    }
  }
}
void EfficiencyValidation::FillHistsEfficiency(TString pre,TString suffix,const vector<Lepton*>& leps,const map<TString,double>& weights){
  for(const auto& element:weights){
    TString suf=suffix+element.first;
    double w=element.second;
    
    FillHist(Form("%snPV%s",pre.Data(),suf.Data()),nPV,w,100,0,100);
    if(!IsDATA) FillHist(Form("%snPileUp%s",pre.Data(),suf.Data()),nPileUp,w,100,0,100);

    if(leps.at(1)->Pt()<20)
      FillHist(Form("%sl1eta_pt20%s",pre.Data(),suf.Data()),leps.at(0)->Eta(),w,120,-3,3);

    //for leptons
    for(int i=0;i<(int)leps.size();i++){
      if(i>1) break;
      double pt=leps.at(i)->Pt();
      double eta=leps.at(i)->Eta();
      TString charge=leps.at(i)->Charge()>0?"p":"m";
      FillHist(Form("%sl%d%spt%s",pre.Data(),i,charge.Data(),suf.Data()),pt,w,500,0,500);
      FillHist(Form("%sl%d%seta%s",pre.Data(),i,charge.Data(),suf.Data()),eta,w,120,-3,3);
      
      FillHist(Form("%sl%dpt%s",pre.Data(),i,suf.Data()),pt,w,500,0,500);
      FillHist(Form("%sl%deta%s",pre.Data(),i,suf.Data()),eta,w,120,-3,3);
      
      FillHist(Form("%sl%spt%s",pre.Data(),charge.Data(),suf.Data()),pt,w,500,0,500);
      FillHist(Form("%sl%seta%s",pre.Data(),charge.Data(),suf.Data()),eta,w,120,-3,3);
      
      FillHist(Form("%slpt%s",pre.Data(),suf.Data()),pt,w,500,0,500);
      FillHist(Form("%sleta%s",pre.Data(),suf.Data()),eta,w,120,-3,3);
      FillHist(Form("%slriso%s",pre.Data(),suf.Data()),leps.at(i)->RelIso(),w,30,0,0.3);
      if(leps.at(i)->LeptonFlavour()==Lepton::Flavour::MUON){
	double rtrkiso=((Muon*)leps.at(i))->TrkIso()/pt;
	FillHist(Form("%slrtrkiso%s",pre.Data(),suf.Data()),rtrkiso,w,40,0,0.2);
      }else if(leps.at(i)->LeptonFlavour()==Lepton::Flavour::ELECTRON){
	Electron* el=(Electron*)leps.at(i);
	FillHist(Form("%slrawpt%s",pre.Data(),suf.Data()),el->UncorrPt(),w,500,0,500);
	FillHist(Form("%slsceta%s",pre.Data(),suf.Data()),el->scEta(),w,120,-3,3);
      }

      /*
      vector<double> muon_ptbin={10,15,20,25,30,40,50,60,120,500};
      vector<double> muon_etabin={-2.4,-2.1,-1.85,-1.6,-1.4,-1.2,-0.9,-0.6,-0.3,-0.2,0,0.2,0.3,0.6,0.9,1.2,1.4,1.6,1.85,2.1,2.4};
      vector<double> electron_ptbin={10,15,20,25,30,35,40,45,70,100,500};
      vector<double> electron_etabin={-2.4,-2.1,-1.8,-1.57,-1.44,-1,-0.6,-0.3,-0.2,0,0.2,0.3,0.6,1,1.44,1.57,1.8,2.1,2.4};
      vector<double>& ptbin=(leps.at(i)->LeptonFlavour()==Lepton::Flavour::MUON?muon_ptbin:electron_ptbin);
      vector<double>& etabin=(leps.at(i)->LeptonFlavour()==Lepton::Flavour::MUON?muon_etabin:electron_etabin);
      for(unsigned int ib=0;ib<ptbin.size()-1;ib++){
	if(pt>ptbin[ib]&&pt<ptbin[ib+1]){
	  FillHist(Form("%sl%d%seta_pt%.0fto%.0f%s",pre.Data(),i,charge.Data(),ptbin[ib],ptbin[ib+1],suf.Data()),eta,w,120,-3,3);
	  FillHist(Form("%sl%deta_pt%.0fto%.0f%s",pre.Data(),i,ptbin[ib],ptbin[ib+1],suf.Data()),eta,w,120,-3,3);
	  FillHist(Form("%sl%seta_pt%.0fto%.0f%s",pre.Data(),charge.Data(),ptbin[ib],ptbin[ib+1],suf.Data()),eta,w,120,-3,3);
	  FillHist(Form("%sleta_pt%.0fto%.0f%s",pre.Data(),ptbin[ib],ptbin[ib+1],suf.Data()),eta,w,120,-3,3);
	}
      } 
      for(unsigned int ib=0;ib<etabin.size()-1;ib++){
	if(eta>etabin[ib]&&eta<etabin[ib+1]){
	  FillHist(Form("%sl%d%spt_eta%.2fto%.2f%s",pre.Data(),i,charge.Data(),etabin[ib],etabin[ib+1],suf.Data()),pt,w,500,0,500);
	  FillHist(Form("%sl%dpt_eta%.2fto%.2f%s",pre.Data(),i,etabin[ib],etabin[ib+1],suf.Data()),pt,w,500,0,500);
	  FillHist(Form("%sl%spt_eta%.2fto%.2f%s",pre.Data(),charge.Data(),etabin[ib],etabin[ib+1],suf.Data()),pt,w,500,0,500);
	  FillHist(Form("%slpt_eta%.2fto%.2f%s",pre.Data(),etabin[ib],etabin[ib+1],suf.Data()),pt,w,500,0,500);
	}
      } 
      */
    }
      
    TLorentzVector dilepton=(*leps.at(0))+(*leps.at(1));
    double dimass=dilepton.M();
    double dipt=dilepton.Pt();
    double dirap=dilepton.Rapidity();
    FillHist(pre+"dimass"+suf,dimass,w,120,60,120);
    FillHist(pre+"dipt"+suf,dipt,w,400,0,400);
    FillHist(pre+"dirap"+suf,dirap,w,120,-3,3);
    FillHist(pre+"nlepton"+suf,leps.size(),w,10,0,10);
    /*
    double cost=GetCosThetaCS((Particle*)leps.at(0),(Particle*)leps.at(1));
    FillHist(pre+"costhetaCS"+suf,cost,w,40,-1,1);
    FillHist(pre+"abscosthetaCS"+suf,fabs(cost),w,20,0,1);
    double h=0.5*pow(dipt/dimass,2)/(1+pow(dipt/dimass,2))*(1-3*cost*cost);
    double den_weight=0.5*fabs(cost)/pow(1+cost*cost+h,2);
    double num_weight=0.5*cost*cost/pow(1+cost*cost+h,3);
    if(cost>0){
      FillHist(pre+"forward"+suf,dimass,w,mbinnum,(double*)mbin);
      FillHist(pre+"forward_den"+suf,dimass,w*den_weight,mbinnum,(double*)mbin);
      FillHist(pre+"forward_num"+suf,dimass,w*num_weight,mbinnum,(double*)mbin);    
    }else{
      FillHist(pre+"backward"+suf,dimass,w,mbinnum,(double*)mbin);
      FillHist(pre+"backward_den"+suf,dimass,w*den_weight,mbinnum,(double*)mbin);
      FillHist(pre+"backward_num"+suf,dimass,w*num_weight,mbinnum,(double*)mbin);   
    }
    */
  }
}
double EfficiencyValidation::LeptonTrigger_SF_OR(TString triggerSF_key0,TString triggerSF_key1,const vector<Lepton*>& leps,int sys){
  if(IsDATA) return 1;

  double lumi0=1.; //only trigger0 on
  double lumi1=0.; //only trigger1 on
  double lumi2=0.; //both off
  if(DataYear==2017&&triggerSF_key0.Contains("IsoMu24")&&triggerSF_key1.Contains("IsoMu27")){
    lumi0=37997.005;    lumi1=3480.873;    lumi2=0.;
  }else if(DataYear==2017&&triggerSF_key0.Contains("Ele27")&&triggerSF_key1.Contains("Ele32")){
    lumi0=31661.026;    lumi1=9522.208;    lumi2=0.295;
  }else if(DataYear==2018&&triggerSF_key0.Contains("Ele28")&&triggerSF_key1.Contains("Ele32")){
    lumi0=23687.253;    lumi1=36140.626;   lumi2=0.;
  }else{
    cout<<"[EfficiencyValidation::LeptonTrigger_SF_OR] not available combination "<<triggerSF_key0<<"||"<<triggerSF_key1<<" for "<<DataEra<<endl;
    exit(EXIT_FAILURE);
  }

  double data_eff_key0=1.,data_eff_key1=1.,mc_eff=1.;
  for(const auto& lep:leps){
    data_eff_key0*=1-Lepton_SF("Trigger_Eff_DATA_"+triggerSF_key0,lep,sys);
    data_eff_key1*=1-Lepton_SF("Trigger_Eff_DATA_"+triggerSF_key1,lep,sys);
    mc_eff*=1-Lepton_SF("Trigger_Eff_MC_"+triggerSF_key0,lep,-sys);
  }
  data_eff_key0=1-data_eff_key0;
  data_eff_key1=1-data_eff_key1;
  mc_eff=1-mc_eff;
  if(mc_eff==0) return 1.;
  else return (lumi0*data_eff_key0+lumi1*data_eff_key1)/((lumi0+lumi1+lumi2)*mc_eff);
}
