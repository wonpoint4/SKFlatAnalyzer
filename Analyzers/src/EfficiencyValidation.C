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

  FillCutflow(tauprefix+"cutflow","lumi",lumiweight);
  FillCutflow(tauprefix+"cutflow","PU",lumiweight*PUweight);
  FillCutflow(tauprefix+"cutflow","prefire",lumiweight*PUweight*prefireweight);
  FillCutflow(tauprefix+"cutflow","zpt",lumiweight*PUweight*prefireweight*zptweight);
  FillCutflow(tauprefix+"cutflow","z0",lumiweight*PUweight*prefireweight*zptweight*z0weight);

  if(DataYear==2016){
    vector<TString> doublemuontrigger;
    doublemuontrigger.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");
    doublemuontrigger.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v");
    doublemuontrigger.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
    doublemuontrigger.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");
    if(event.PassTrigger(doublemuontrigger))
      if(!IsDATA||DataStream.Contains("DoubleMuon")) executeEventWithChannelName("mm2016");
    if(event.PassTrigger("HLT_IsoMu24_v")||event.PassTrigger("HLT_IsoTkMu24_v"))
      if(!IsDATA||DataStream.Contains("SingleMuon")) executeEventWithChannelName("mu2016");
    if(event.PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"))
      if(!IsDATA||DataStream.Contains("DoubleEG")) executeEventWithChannelName("ee2016");
    if(event.PassTrigger("HLT_Ele27_WPTight_Gsf_v"))
      if(!IsDATA||DataStream.Contains("SingleElectron")) executeEventWithChannelName("el2016");
  }else if(DataYear==2017){
    if(event.PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v"))
      if(!IsDATA||DataStream.Contains("DoubleMuon")) executeEventWithChannelName("mm2017");
    if(event.PassTrigger("HLT_IsoMu27_v"))
      if(!IsDATA||DataStream.Contains("SingleMuon")) executeEventWithChannelName("mu2017");    
    if(event.PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"))
      if(!IsDATA||DataStream.Contains("DoubleEG")) executeEventWithChannelName("ee2017");
    if(event.PassTrigger("HLT_Ele32_WPTight_Gsf_L1DoubleEG_v"))
      if(!IsDATA||DataStream.Contains("SingleElectron")) executeEventWithChannelName("el2017");
  }else if(DataYear==2018){
    if(event.PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"))
      if(!IsDATA||DataStream.Contains("DoubleMuon")) executeEventWithChannelName("mm2018");
    if(event.PassTrigger("HLT_IsoMu24_v"))
      if(!IsDATA||DataStream.Contains("SingleMuon")) executeEventWithChannelName("mu2018");
    if(event.PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"))
      if(!IsDATA||DataStream.Contains("EGamma")) executeEventWithChannelName("ee2018");
    if(event.PassTrigger("HLT_Ele32_WPTight_Gsf_v"))
      if(!IsDATA||DataStream.Contains("EGamma")) executeEventWithChannelName("el2018");
  }    

}

void EfficiencyValidation::executeEventWithChannelName(TString channelname){
  map<TString,vector<Muon>> map_muons;
  map<TString,vector<Electron>> map_electrons;
  map<TString,Parameter> map_parameter;

  if(channelname.Contains(TRegexp("mm20[0-9][0-9]"))){
    map_muons["_MediumID_trkIsoLoose"]=SMPGetMuons("POGMediumWithLooseTrkIso",0.0,2.4);
    map_parameter["_MediumID_trkIsoLoose_Q"]=Parameter("IDISO_SF_MediumID_trkIsoLoose_Q",{"Mu17Leg1_MediumID_trkIsoLoose_Q","Mu8Leg2_MediumID_trkIsoLoose_Q"},20,10,MakeLeptonPointerVector(map_muons["_MediumID_trkIsoLoose"]));
    if(DataYear!=2018){
      map_muons["_TightID_PFIsoTight"]=SMPGetMuons("POGTightWithTightIso",0.0,2.4);
      map_parameter["_TightID_PFIsoTight"]=Parameter("ID_SF_NUM_TightID_DEN_genTracks","ISO_SF_NUM_TightRelIso_DEN_TightIDandIPCut",{"Mu17Leg1_TightID_PFIsoTight","Mu8Leg2_TightID_PFIsoTight"},20,10,MakeLeptonPointerVector(map_muons["_TightID_PFIsoTight"]));
    }
  }else if(channelname.Contains(TRegexp("mu20[0-9][0-9]"))){
    map_muons["_MediumID_trkIsoLoose"]=SMPGetMuons("POGMediumWithLooseTrkIso",0.0,2.4);
    vector<Lepton*> leps=MakeLeptonPointerVector(map_muons["_MediumID_trkIsoLoose"]);
    switch(DataYear){
    case 2016:
      map_parameter["_MediumID_trkIsoLoose_Q"]=Parameter("IDISO_SF_MediumID_trkIsoLoose_Q","",{"IsoMu24_MediumID_trkIsoLoose_Q"},27.,10.,leps);
      break;
    case 2017:
      map_parameter["_MediumID_trkIsoLoose_Q"]=Parameter("IDISO_SF_MediumID_trkIsoLoose_Q","",{"IsoMu27_MediumID_trkIsoLoose_Q"},30.,10.,leps);
      break;
    case 2018:
      map_parameter["_MediumID_trkIsoLoose_Q"]=Parameter("IDISO_SF_MediumID_trkIsoLoose_Q","",{"IsoMu24_MediumID_trkIsoLoose_Q"},27.,10.,leps);
      break;
    default:
      cout<<"[EfficiencyValidation::executeEventFromParameter] wrong year"<<endl;
      exit(EXIT_FAILURE);
    }
  }else if(channelname.Contains(TRegexp("ee20[0-9][0-9]"))){
    Parameter p;
    p.electronIDSF="ID_SF_MediumID_Q";
    p.triggerSF={"Ele23Leg1_MediumID_Q","Ele12Leg2_MediumID_Q"};
    p.lep0ptcut=25.;
    p.lep1ptcut=15.;
    map_electrons["_MediumID_noroccor"]=SMPGetElectrons("passMediumID",0.0,2.4);
    map_electrons["_MediumID"]=ElectronEnergyCorrection(map_electrons["_MediumID_noroccor"],0,0);
    map_parameter["_MediumID_Q"]=p.Clone(MakeLeptonPointerVector(map_electrons["_MediumID"]));
    if(DataYear!=2016){
      p.electronIDSF="ID_SF_MediumID_Q_v1";
      p.triggerSF={"Ele23Leg1_MediumID_Q_v1","Ele12Leg2_MediumID_Q_v1"};
      map_parameter["_MediumID_Q_v1"]=p.Clone(MakeLeptonPointerVector(map_electrons["_MediumID"]));
    }
    p.electronIDSF="ID_SF_MediumID_Q_v1_3";
    p.triggerSF={"Ele23Leg1_MediumID_Q_v1_3","Ele12Leg2_MediumID_Q_v1_3"};
    map_parameter["_MediumID_Q_v1_3"]=p.Clone(MakeLeptonPointerVector(map_electrons["_MediumID"]));
  }else if(channelname.Contains(TRegexp("el20[0-9][0-9]"))){
    map_electrons["_MediumID_noroccor"]=SMPGetElectrons("passMediumID",0.0,2.4);
    map_electrons["_MediumID"]=ElectronEnergyCorrection(map_electrons["_MediumID_noroccor"],0,0);
    map_electrons["_TightID_Selective_noroccor"]=SMPGetElectrons("passMediumID_Selective",0.0,2.4);
    map_electrons["_TightID_Selective"]=ElectronEnergyCorrection(map_electrons["_TightID_Selective_noroccor"],0,0);
    switch(DataYear){
    case 2016:
      map_parameter["_MediumID_Q"]=Parameter("ID_SF_MediumID_Q",{"Ele27_MediumID_Q"},30,10,MakeLeptonPointerVector(map_electrons["_MediumID"])); 
      map_parameter["_MediumID_Q_v1_3"]=Parameter("ID_SF_MediumID_Q_v1_3",{"Ele27_MediumID_Q_v1_3"},30,10,MakeLeptonPointerVector(map_electrons["_MediumID"])); 
      map_parameter["_TightID_Selective_Q"]=Parameter("ID_SF_TightID_Selective_Q",{"Ele27_TightID_Selective_Q"},30,10,MakeLeptonPointerVector(map_electrons["_TightID_Selective"])); 
      map_parameter["_TightID_Selective_Q_v1_3"]=Parameter("ID_SF_TightID_Selective_Q_v1_3",{"Ele27_TightID_Selective_Q_v1_3"},30,10,MakeLeptonPointerVector(map_electrons["_TightID_Selective"]));
      break;
    case 2017:
      map_parameter["_MediumID_Q"]=Parameter("ID_SF_MediumID_Q",{"Ele32_MediumID_Q"},35,10,MakeLeptonPointerVector(map_electrons["_MediumID"])); 
      map_parameter["_MediumID_Q_v1"]=Parameter("ID_SF_MediumID_Q_v1",{"Ele32_MediumID_Q_v1"},35,10,MakeLeptonPointerVector(map_electrons["_MediumID"])); 
      map_parameter["_MediumID_Q_v1_3"]=Parameter("ID_SF_MediumID_Q_v1_3",{"Ele32_MediumID_Q_v1_3"},35,10,MakeLeptonPointerVector(map_electrons["_MediumID"])); 
      map_parameter["_TightID_Selective_Q"]=Parameter("ID_SF_TightID_Selective_Q",{"Ele32_TightID_Selective_Q"},35,10,MakeLeptonPointerVector(map_electrons["_TightID_Selective"])); 
      map_parameter["_TightID_Selective_Q_v1"]=Parameter("ID_SF_TightID_Selective_Q_v1",{"Ele32_TightID_Selective_Q_v1"},35,10,MakeLeptonPointerVector(map_electrons["_TightID_Selective"]));
      map_parameter["_TightID_Selective_Q_v1_3"]=Parameter("ID_SF_TightID_Selective_Q_v1_3",{"Ele32_TightID_Selective_Q_v1_3"},35,10,MakeLeptonPointerVector(map_electrons["_TightID_Selective"]));
      break;
    case 2018: 
      map_parameter["_MediumID_Q"]=Parameter("ID_SF_MediumID_Q",{"Ele32_MediumID_Q"},35,10,MakeLeptonPointerVector(map_electrons["_MediumID"])); 
      map_parameter["_MediumID_Q_v1"]=Parameter("ID_SF_MediumID_Q_v1",{"Ele32_MediumID_Q_v1"},35,10,MakeLeptonPointerVector(map_electrons["_MediumID"])); 
      map_parameter["_MediumID_Q_v1_3"]=Parameter("ID_SF_MediumID_Q_v1_3",{"Ele32_MediumID_Q_v1_3"},35,10,MakeLeptonPointerVector(map_electrons["_MediumID"])); 
      map_parameter["_TightID_Selective_Q"]=Parameter("ID_SF_TightID_Selective_Q",{"Ele32_TightID_Selective_Q"},35,10,MakeLeptonPointerVector(map_electrons["_TightID_Selective"])); 
      map_parameter["_TightID_Selective_Q_v1"]=Parameter("ID_SF_TightID_Selective_Q_v1",{"Ele32_TightID_Selective_Q_v1"},35,10,MakeLeptonPointerVector(map_electrons["_TightID_Selective"]));
      map_parameter["_TightID_Selective_Q_v1_3"]=Parameter("ID_SF_TightID_Selective_Q_v1_3",{"Ele32_TightID_Selective_Q_v1_3"},35,10,MakeLeptonPointerVector(map_electrons["_TightID_Selective"]));
      break;
    default: 
      cout<<"[EfficiencyValidation::executeEventFromParameter] wrong year"<<endl;
      exit(EXIT_FAILURE);
    }
  }else{
    cout<<"[EfficiencyValidation::executeEventFromParameter] wrong channelname"<<endl;
    exit(EXIT_FAILURE);
  }
  
  ///////////////////////lepton selection///////////////////////
  for(const auto& [suffix,p]: map_parameter){
    TString prefix=tauprefix;
    double eventweight=lumiweight*PUweight*prefireweight*z0weight*zptweight;

    FillHist(channelname+"/"+prefix+"nlepton"+suffix,p.leps.size(),eventweight,10,0,10);
    if(p.leps.size()>=2){
      TString prefix=tauprefix;
      if(p.leps.at(0)->Charge()*p.leps.at(1)->Charge()>0)
	//prefix="ss_"+prefix;
	continue;
      FillCutflow(channelname+"/"+prefix+"cutflow"+suffix,"OS dilepton",eventweight);
      if(p.leps.at(0)->Pt()>p.lep0ptcut&&p.leps.at(1)->Pt()>p.lep1ptcut){
	FillCutflow(channelname+"/"+prefix+"cutflow"+suffix,"LepPtCut",eventweight);
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

	      double this_RECOSF=mcCorr->ElectronReco_SF(this_eta,this_pt,0);
	      double this_RECOSF_up=mcCorr->ElectronReco_SF(this_eta,this_pt,1);
	      double this_RECOSF_down=mcCorr->ElectronReco_SF(this_eta,this_pt,-1);
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
            triggerSF*=DileptonTrigger_SF(p.triggerSF[0],p.triggerSF[1],p.leps,0);
            triggerSF_up*=DileptonTrigger_SF(p.triggerSF[0],p.triggerSF[1],p.leps,1);
            triggerSF_down*=DileptonTrigger_SF(p.triggerSF[0],p.triggerSF[1],p.leps,-1);
          }
        }
	
        FillCutflow(channelname+"/"+prefix+"cutflow"+suffix,"RECOSF",eventweight*RECOSF);
        FillCutflow(channelname+"/"+prefix+"cutflow"+suffix,"IDSF",eventweight*RECOSF*IDSF);
        FillCutflow(channelname+"/"+prefix+"cutflow"+suffix,"ISOSF",eventweight*RECOSF*IDSF*ISOSF);
        FillCutflow(channelname+"/"+prefix+"cutflow"+suffix,"triggerSF",eventweight*RECOSF*IDSF*ISOSF*triggerSF);

	///////////////////////weight systematics//////////////////
	map<TString,double> map_weight;
	map_weight[""]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF;
	if(!IsDATA){
	  map_weight["_noweight"]=lumiweight;
	  map_weight["_noPUweight"]=lumiweight*prefireweight*zptweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF;
	  map_weight["_noprefireweight"]=lumiweight*PUweight*RECOSF*IDSF*ISOSF*triggerSF*zptweight*z0weight;
	  
	  map_weight["_noRECOSF"]=lumiweight*PUweight*IDSF*ISOSF*triggerSF*prefireweight*zptweight*z0weight;
	  map_weight["_RECOSF_up"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF_up*IDSF*ISOSF*triggerSF;
	  map_weight["_RECOSF_down"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF_down*IDSF*ISOSF*triggerSF;
	
	  map_weight["_noIDSF"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF*ISOSF*triggerSF;
	  map_weight["_IDSF_up"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF*IDSF_up*ISOSF*triggerSF;
	  map_weight["_IDSF_down"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF*IDSF_down*ISOSF*triggerSF;
	  
	  map_weight["_noISOSF"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF*IDSF*triggerSF;
	  map_weight["_ISOSF_up"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF*IDSF*ISOSF_up*triggerSF;
	  map_weight["_ISOSF_down"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF*IDSF*ISOSF_down*triggerSF;
	
	  map_weight["_notriggerSF"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF*IDSF*ISOSF;
	  map_weight["_triggerSF_up"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF_up;
	  map_weight["_triggerSF_down"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF_down;
	  
	  map_weight["_noefficiencySF"]=lumiweight*PUweight*prefireweight*zptweight*z0weight;
	  map_weight["_noz0weight"]=lumiweight*PUweight*prefireweight*zptweight*RECOSF*IDSF*ISOSF*triggerSF;
	  map_weight["_nozptweight"]=lumiweight*PUweight*prefireweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF;
	}
	
	///////////////////////fill hists///////////////////////
	TLorentzVector dilepton=(*p.leps.at(0))+(*p.leps.at(1));
	double dimass=dilepton.M();
	if(dimass>=60&&dimass<120){
	  FillCutflow(channelname+"/"+prefix+"cutflow"+suffix,"m60to120",map_weight[""]);
	  FillHistsEfficiency(channelname+"/m60to120/"+prefix,suffix,p.leps,map_weight);
	  if(dimass>=80&&dimass<100){
	    FillCutflow(channelname+"/"+prefix+"cutflow"+suffix,"m80to100",map_weight[""]);
	    FillHistsEfficiency(channelname+"/m80to100/"+prefix,suffix,p.leps,map_weight);
	  }
	}
      }
    }
  }
}
void EfficiencyValidation::FillHistsEfficiency(TString pre,TString suffix,const vector<Lepton*>& leps,const map<TString,double>& weights){
  for(const auto& element:weights){
    TString suf=suffix+element.first;
    double w=element.second;
    
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
    FillHist(pre+"dimass"+suf,dimass,w,400,0,400);
    FillHist(pre+"dipt"+suf,dipt,w,400,0,400);
    FillHist(pre+"dirap"+suf,dirap,w,120,-3,3);
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
