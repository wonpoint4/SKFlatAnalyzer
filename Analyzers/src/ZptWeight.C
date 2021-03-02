#include "ZptWeight.h"

ZptWeight::ZptWeight(){
}
ZptWeight::~ZptWeight(){
}
void ZptWeight::executeEvent(){
  event=GetEvent();
  GetEventWeights();

  ////////////////////////check genlevel//////////////////
  if(IsDYSample){
    vector<Gen> gens=GetGens();
    Gen parton0,parton1,l0,l1;
    GetDYGenParticles(gens,parton0,parton1,l0,l1,3);
    if(tauprefix!="tau_"){
      TLorentzVector genZ=(l0+l1);
      FillHist(Form("%s%d/%s",abs(l0.PID())==13?"mm":"ee",DataYear,"/gen_diptdirap"),genZ.Pt(),fabs(genZ.Rapidity()),lumiweight*zptweight,ptbinnum,ptbin,rapbinnum,rapbin);
      FillHist(Form("%s%d/%s",abs(l0.PID())==13?"mm":"ee",DataYear,"/gen_diptdirap_nozptweight"),genZ.Pt(),fabs(genZ.Rapidity()),lumiweight,ptbinnum,ptbin,rapbinnum,rapbin);
    }
  }

  if(!PassMETFilter()) return;

  if(DataYear==2016){
    vector<TString> doublemuontrigger;
    doublemuontrigger.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");
    doublemuontrigger.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v");
    doublemuontrigger.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
    doublemuontrigger.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");
    if(event.PassTrigger(doublemuontrigger))
      if(!IsDATA||DataStream.Contains("DoubleMuon")) executeEventWithChannelName("mm2016");
    
    if(event.PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"))
      if(!IsDATA||DataStream.Contains("DoubleEG")) executeEventWithChannelName("ee2016");
    
  }else if(DataYear==2017){
    if(event.PassTrigger("HLT_IsoMu24_v")||event.PassTrigger("HLT_IsoMu27_v"))
      if(!IsDATA||DataStream.Contains("SingleMuon")) executeEventWithChannelName("mu2017");

    if(event.PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v"))
      if(!IsDATA||DataStream.Contains("DoubleMuon")) executeEventWithChannelName("mm2017");
    
    if(event.PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"))
      if(!IsDATA||DataStream.Contains("DoubleEG")) executeEventWithChannelName("ee2017");
    
  }else if(DataYear==2018){
    if(event.PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"))
      if(!IsDATA||DataStream.Contains("DoubleMuon")) executeEventWithChannelName("mm2018");
    
    if(event.PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"))
      if(!IsDATA||DataStream.Contains("EGamma")) executeEventWithChannelName("ee2018");
    
  }    

}

void ZptWeight::executeEventWithChannelName(TString channelname){
  map<TString,vector<Muon>> map_muons;
  map<TString,vector<Electron>> map_electrons;
  map<TString,Parameter> map_parameter;

  if(channelname.Contains(TRegexp("mm20[0-9][0-9]"))){
    Parameter p;
    p.muonIDSF="IDISO_SF_MediumID_trkIsoLoose_Q";
    p.triggerSF={"Mu17Leg1_MediumID_trkIsoLoose_Q","Mu8Leg2_MediumID_trkIsoLoose_Q"};
    p.lep0ptcut=20.;
    p.lep1ptcut=10.;
    map_muons[""]=MuonMomentumCorrection(SMPGetMuons("POGMediumWithLooseTrkIso",0.0,2.4),0);
    map_parameter[""]=p.Clone(MakeLeptonPointerVector(map_muons[""]));
  }else if(channelname.Contains(TRegexp("mu20[0-9][0-9]"))){
    Parameter p;
    p.muonIDSF="IDISO_SF_MediumID_trkIsoLoose_Q";
    p.triggerSF={"IsoMu2427_MediumID_trkIsoLoose_Q"};
    p.lep0ptcut=25.;
    p.lep1ptcut=10.;
    map_muons[""]=MuonMomentumCorrection(SMPGetMuons("POGMediumWithLooseTrkIso",0.0,2.4),0);
    map_parameter[""]=p.Clone(MakeLeptonPointerVector(map_muons[""]));
  }else if(channelname.Contains(TRegexp("ee20[0-9][0-9]"))){
    Parameter p;
    p.electronIDSF="ID_SF_MediumID_Q";
    p.triggerSF={"Ele23Leg1_MediumID_Q","Ele12Leg2_MediumID_Q"};
    p.lep0ptcut=25.;
    p.lep1ptcut=15.;
    map_electrons[""]=ElectronEnergyCorrection(SMPGetElectrons("passMediumID",0.0,2.4),0,0);
    map_parameter[""]=p.Clone(MakeLeptonPointerVector(map_electrons[""]));
  }else{
    cout<<"[ZptWeight::executeEventFromParameter] wrong channelname"<<endl;
    exit(EXIT_FAILURE);
  }

  FillCutflow(channelname+"/"+tauprefix+"cutflow","lumi",lumiweight);
  FillCutflow(channelname+"/"+tauprefix+"cutflow","PU",lumiweight*PUweight);
  FillCutflow(channelname+"/"+tauprefix+"cutflow","prefire",lumiweight*prefireweight);
  FillCutflow(channelname+"/"+tauprefix+"cutflow","zpt",lumiweight*prefireweight*zptweight);
  FillCutflow(channelname+"/"+tauprefix+"cutflow","z0",lumiweight*prefireweight*zptweight*z0weight);

  ///////////////////////lepton selection///////////////////////
  for(const auto& [suffix,p]:map_parameter){
    TString prefix=tauprefix;
    double eventweight=lumiweight*PUweight*prefireweight*z0weight*zptweight;

    FillHist(channelname+"/"+prefix+"nlepton"+suffix,p.leps.size(),eventweight,10,0,10);
    if(p.leps.size()>=2){
      if(p.leps.at(0)->Charge()*p.leps.at(1)->Charge()>0) return;
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
	  map_weight["_nozptweight"]=lumiweight*PUweight*prefireweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF;
	}
	
	///////////////////////fill hists///////////////////////
	TLorentzVector dilepton=(*p.leps.at(0))+(*p.leps.at(1));
	double dimass=dilepton.M();
	if(dimass>=80&&dimass<100){
	  FillCutflow(channelname+"/"+prefix+"cutflow"+suffix,"m80to100",map_weight[""]);
	  FillHistsZptWeight(channelname+"/m80to100/"+prefix,suffix,p.leps,map_weight);
	}
      }
    }
  }
}
void ZptWeight::FillHistsZptWeight(TString pre,TString suffix,const vector<Lepton*>& leps,const map<TString,double>& weights){
  for(const auto& element:weights){
    TString suf=suffix;
    suf+=element.first;
    double w=element.second;
    
    TLorentzVector dilepton=(*leps.at(0))+(*leps.at(1));
    double dipt=dilepton.Pt();
    double dirap=dilepton.Rapidity();
    FillHist(pre+"diptdirap"+suf,dipt,fabs(dirap),w,ptbinnum,ptbin,rapbinnum,rapbin);
    FillHist(pre+"dipt"+suf,dipt,w,100,0,200);
    for(int i=0;i<rapbinnum;i++){
      if(fabs(dirap)>rapbin[i]&&fabs(dirap)<rapbin[i+1]){
	FillHist(pre+Form("dipt_bin%d",i)+suf,dipt,w,100,0,200);
	break;
      }
    }
  }
}
