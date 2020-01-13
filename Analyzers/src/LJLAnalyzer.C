#include "LJLAnalyzer.h"

LJLAnalyzer::LJLAnalyzer(){
}
LJLAnalyzer::~LJLAnalyzer(){
}
void LJLAnalyzer::executeEvent(){
  tauprefix="";

  if(IsDYSample){
    vector<Gen> gens=GetGens();
    int parton0,parton1,hardl0,hardl1,l0,l1;
    vector<int> photons;
    GetGenIndex(gens,parton0,parton1,hardl0,hardl1,l0,l1,photons);
    if(abs(gens[hardl0].PID())==15) tauprefix="tau_";
  }
  
  if(!PassMETFilter()) return;

  Event* ev=new Event;
  *ev=GetEvent();
  TString channelname,muontrigger,electrontrigger;
  if(DataYear==2016){
    vector<TString> muontriggers;
    muontriggers.push_back("HLT_IsoMu24_v");
    muontriggers.push_back("HLT_IsoTkMu24_v");
    electrontrigger="HLT_Ele27_WPTight_Gsf_v";
    if(ev->PassTrigger(muontriggers)){
      channelname="muon2016";
      if(!IsDATA||DataStream.Contains("SingleMuon")) executeEventFromParameter(channelname,ev);
    }
    if(ev->PassTrigger(electrontrigger)){
      channelname="electron2016";
      if(!IsDATA||DataStream.Contains("SingleEG")) executeEventFromParameter(channelname,ev);
    }
  }else if(DataYear==2017){
    muontrigger="HLT_IsoMu27_v";
    electrontrigger="HLT_Ele32_WPTight_Gsf_L1DoubleEG_v";
    if(ev->PassTrigger(muontrigger)){
      channelname="muon2017";
      if(!IsDATA||DataStream.Contains("SingleMuon")) executeEventFromParameter(channelname,ev);
    }
    if(ev->PassTrigger(electrontrigger)){
      channelname="electron2017";
      if(!IsDATA||DataStream.Contains("SingleEG")) executeEventFromParameter(channelname,ev);
    }
  }else if(DataYear==2018){
    muontrigger="HLT_IsoMu24_v";
    electrontrigger="HLT_Ele32_WPTight_Gsf_v";
    if(ev->PassTrigger(muontrigger)){
      channelname="muon2018";
      if(!IsDATA||DataStream.Contains("SingleMuon")) executeEventFromParameter(channelname,ev);
    }
    if(ev->PassTrigger(electrontrigger)){
      channelname="electron2018";
      if(!IsDATA||DataStream.Contains("EGamma")) executeEventFromParameter(channelname,ev);
    }
  }    

  delete ev;
}

void LJLAnalyzer::executeEventFromParameter(TString channelname,Event* ev){
  map< TString, std::vector<Muon> > map_muons;
  map< TString, std::vector<Electron> > map_electrons;
  //suffix, lepton vector, SF1, SF2, TriggerSF1, TriggerSF2
  map< TString, tuple<vector<Lepton*>,vector<Lepton*>,TString,TString,TString> > map_leps;

  vector<Jet> jets=GetJets("tight",30,2.7);

  double lep0ptcut,lep1ptcut;
  if(channelname.Contains("muon")){
    lep0ptcut=30.;
    lep1ptcut=80.;

    map_muons["iso"]=SMPGetMuons("POGTightWithTightIso",lep0ptcut,2.4);
    map_muons["noniso"]=SMPGetMuons("POGTightWithAntiIso",lep1ptcut,2.4);

    switch(DataYear){
    case 2016: 
      map_leps[""]=make_tuple(MakeLeptonPointerVector(map_muons["iso"]),MakeLeptonPointerVector(map_muons["noniso"]),"ID_SF_NUM_TightID_DEN_genTracks","ISO_SF_NUM_TightRelIso_DEN_TightIDandIPCut","IsoMu24_POGTight");
      break;
    case 2017:  
      map_leps[""]=make_tuple(MakeLeptonPointerVector(map_muons["iso"]),MakeLeptonPointerVector(map_muons["noniso"]),"ID_SF_NUM_TightID_DEN_genTracks","ISO_SF_NUM_TightRelIso_DEN_TightIDandIPCut","IsoMu27_POGTight");
      break;
    case 2018: 
      map_leps[""]=make_tuple(MakeLeptonPointerVector(map_muons["iso"]),MakeLeptonPointerVector(map_muons["noniso"]),"ID_SF_NUM_TightID_DEN_genTracks","ISO_SF_NUM_TightRelIso_DEN_TightIDandIPCut","IsoMu24_POGTight");
      break;
    default: cout<<"[LJLAnalyzer::executeEventFromParameter] wrong year"<<endl;exit(EXIT_FAILURE);break;
    }
 }else if(channelname.Contains("electron")){
    /*
    lep0ptcut=35.;
    lep1ptcut=20.;

    map_electrons[""]=SMPGetElectrons("passMediumID",lep1ptcut,2.5);

    switch(DataYear){
    case 2016: 
      map_leps[""]=make_tuple(MakeLeptonPointerVector(map_electrons[""]),"ID_SF_MediumID_pt10","","TRIGGER"); 
      break;
    case 2017: 
      map_leps[""]=make_tuple(MakeLeptonPointerVector(map_electrons[""]),"ID_SF_MediumID_pt10","","TRIGGER"); 
      break;
    case 2018: 
      map_leps[""]=make_tuple(MakeLeptonPointerVector(map_electrons[""]),"ID_SF_passMediumID","","","TRIGGER"); 
      break;
    default: 
      cout<<"[LJLAnalyzer::executeEventFromParameter] wrong year"<<endl;
      exit(EXIT_FAILURE);
    }
    */
  }else{
    cout<<"[LJLAnalyzer::executeEventFromParameter] wrong channelname"<<endl;
    exit(EXIT_FAILURE);
  }
  
  /////////////////lumi weight///////////////////
  double weight=1.,totalweight=1.;
  if(!IsDATA){
    weight=weight_norm_1invpb*ev->MCweight()*ev->GetTriggerLumi("Full");
  }
  totalweight*=weight;
  FillCutflow(channelname+"/"+tauprefix+"cutflow","lumi",totalweight);

  /////////////////PUreweight///////////////////
  double PUreweight=1.,PUreweight_up=1.,PUreweight_down=1.;
  if(!IsDATA){
    PUreweight=mcCorr->GetPileUpWeight(nPileUp,0);
    PUreweight_up=mcCorr->GetPileUpWeight(nPileUp,1);
    PUreweight_down=mcCorr->GetPileUpWeight(nPileUp,-1);
  }
  totalweight*=PUreweight;
  FillCutflow(channelname+"/"+tauprefix+"cutflow","PU",totalweight);
  
  //////////////////////PrefileWeight////////////////////
  double prefireweight=1.;
  double prefireweight_up=1.;
  double prefireweight_down=1.;
  if(!IsDATA&&DataYear<2018){
    prefireweight=L1PrefireReweight_Central;
    prefireweight_up=L1PrefireReweight_Up;
    prefireweight_down=L1PrefireReweight_Down;
  }
  totalweight*=prefireweight;
  FillCutflow(channelname+"/"+tauprefix+"cutflow","prefire",totalweight);

  ///////////////////////lepton selection///////////////////////
  for(const auto& [suffix,tu]: map_leps){
    auto const& isoleps=get<0>(tu);    
    auto const& nonisoleps=get<1>(tu);    

    if(isoleps.size()&&nonisoleps.size()){
      Lepton* isolep=isoleps[0];
      Lepton* nonisolep=nonisoleps[0];
      TString prefix;
      if(isolep->Charge()>0&&nonisolep->Charge()>0) prefix="pp_"+tauprefix;
      else if(isolep->Charge()<0&&nonisolep->Charge()<0) prefix="mm_"+tauprefix;
      else prefix=tauprefix;      
      FillCutflow(channelname+"/"+tauprefix+"cutflow"+suffix,"LJL",totalweight);
      if((*isolep+*nonisolep).M()<60) return;
      FillCutflow(channelname+"/"+tauprefix+"cutflow"+suffix,"mll>60",totalweight);      
      
      /////////////////efficiency scale factors///////////////////
      double IDSF=1.,IDSF_up=1.,IDSF_down=1.;
      double ISOSF=1.,ISOSF_up=1.,ISOSF_down=1.;
      double RECOSF=1.,RECOSF_up=1.,RECOSF_down=1.;
      if(!IsDATA){
	for(const auto& lep:isoleps){
	  TString LeptonIDSF_key=get<2>(tu);
	  TString LeptonISOSF_key=get<3>(tu);

	  if(isolep->LeptonFlavour()==Lepton::ELECTRON){
	    double this_pt,this_eta;
	    this_pt=lep->Pt();
	    this_eta=((Electron*)lep)->scEta();
	    
	    double this_RECOSF=mcCorr->ElectronReco_SF(this_eta,this_pt,0);
	    double this_RECOSF_up=mcCorr->ElectronReco_SF(this_eta,this_pt,1);
	    double this_RECOSF_down=mcCorr->ElectronReco_SF(this_eta,this_pt,-1);
	    RECOSF*=this_RECOSF; RECOSF_up*=this_RECOSF_up; RECOSF_down*=this_RECOSF_down;
	  }
	  
	  double this_IDSF=Lepton_SF(LeptonIDSF_key,lep,0);
	  double this_IDSF_up=Lepton_SF(LeptonIDSF_key,lep,1);
	  double this_IDSF_down=Lepton_SF(LeptonIDSF_key,lep,-1);
	  IDSF*=this_IDSF; IDSF_up*=this_IDSF_up; IDSF_down*=this_IDSF_down;
	
	  double this_ISOSF=Lepton_SF(LeptonISOSF_key,lep,0);
	  double this_ISOSF_up=Lepton_SF(LeptonISOSF_key,lep,1);
	  double this_ISOSF_down=Lepton_SF(LeptonISOSF_key,lep,-1);
	  ISOSF*=this_ISOSF; ISOSF_up*=this_ISOSF_up; ISOSF_down*=this_ISOSF_down;

	}
	for(const auto& lep:nonisoleps){
	  TString LeptonIDSF_key=get<2>(tu);
	  TString LeptonISOSF_key=get<3>(tu);
	  if(isolep->LeptonFlavour()==Lepton::ELECTRON){
	    double this_pt,this_eta;
	    this_pt=lep->Pt();
	    this_eta=((Electron*)lep)->scEta();
	    
	    double this_RECOSF=mcCorr->ElectronReco_SF(this_eta,this_pt,0);
	    double this_RECOSF_up=mcCorr->ElectronReco_SF(this_eta,this_pt,1);
	    double this_RECOSF_down=mcCorr->ElectronReco_SF(this_eta,this_pt,-1);
	    RECOSF*=this_RECOSF; RECOSF_up*=this_RECOSF_up; RECOSF_down*=this_RECOSF_down;
	  }
	  
	  double this_IDSF=Lepton_SF(LeptonIDSF_key,lep,0);
	  double this_IDSF_up=Lepton_SF(LeptonIDSF_key,lep,1);
	  double this_IDSF_down=Lepton_SF(LeptonIDSF_key,lep,-1);
	  IDSF*=this_IDSF; IDSF_up*=this_IDSF_up; IDSF_down*=this_IDSF_down;
	}
      }
      double triggerSF=1.,triggerSF_up=1.,triggerSF_down=1.;
      TString triggerSF_key=get<4>(tu);
      triggerSF*=LeptonTrigger_SF(triggerSF_key,isoleps,0);
      triggerSF_up*=LeptonTrigger_SF(triggerSF_key,isoleps,1);
      triggerSF_down*=LeptonTrigger_SF(triggerSF_key,isoleps,-1);

      FillCutflow(channelname+"/"+tauprefix+"cutflow"+suffix,"RECO",totalweight*RECOSF);
      FillCutflow(channelname+"/"+tauprefix+"cutflow"+suffix,"ID",totalweight*RECOSF*IDSF);
      FillCutflow(channelname+"/"+tauprefix+"cutflow"+suffix,"ISO",totalweight*RECOSF*IDSF*ISOSF);
      FillCutflow(channelname+"/"+tauprefix+"cutflow"+suffix,"trigger",totalweight*RECOSF*IDSF*ISOSF*triggerSF);

      ///////////////////////weight systematics//////////////////
      map<TString,double> map_weight;
      map_weight[""]=weight*PUreweight*RECOSF*IDSF*ISOSF*triggerSF*prefireweight;
      if(!IsDATA){
	/*
	map_weight["_noefficiencySF"]=weight*PUreweight*prefireweight;
	
	map_weight["_noRECOSF"]=weight*PUreweight*IDSF*ISOSF*triggerSF*prefireweight;
	map_weight["_RECOSF_up"]=weight*PUreweight*RECOSF_up*IDSF*ISOSF*triggerSF*prefireweight;
	map_weight["_RECOSF_down"]=weight*PUreweight*RECOSF_down*IDSF*ISOSF*triggerSF*prefireweight;
	
	map_weight["_noIDSF"]=weight*PUreweight*RECOSF*ISOSF*triggerSF*prefireweight;
	map_weight["_IDSF_up"]=weight*PUreweight*RECOSF*IDSF_up*ISOSF*triggerSF*prefireweight;
	map_weight["_IDSF_down"]=weight*PUreweight*RECOSF*IDSF_down*ISOSF*triggerSF*prefireweight;
	
	map_weight["_noISOSF"]=weight*PUreweight*RECOSF*IDSF*triggerSF*prefireweight;
	map_weight["_ISOSF_up"]=weight*PUreweight*RECOSF*IDSF*ISOSF_up*triggerSF*prefireweight;
	map_weight["_ISOSF_down"]=weight*PUreweight*RECOSF*IDSF*ISOSF_down*triggerSF*prefireweight;
	
	map_weight["_notriggerSF"]=weight*PUreweight*RECOSF*IDSF*ISOSF*prefireweight;
	map_weight["_triggerSF_up"]=weight*PUreweight*RECOSF*IDSF*ISOSF*triggerSF_up*prefireweight;
	map_weight["_triggerSF_down"]=weight*PUreweight*RECOSF*IDSF*ISOSF*triggerSF_down*prefireweight;	
	*/
      }
      
      ///////////////////////fill hists///////////////////////
      FillLJLHists(channelname+"/"+prefix,suffix,isoleps,nonisoleps,jets,map_weight);
      if(jets.size()<5){
	TString region=Form("%djet",int(jets.size()));
	FillLJLHists(channelname+"/"+region+"/"+prefix,suffix,isoleps,nonisoleps,jets,map_weight);
      }
    }
  }
}
void LJLAnalyzer::FillLJLHists(TString pre,TString suffix,const vector<Lepton*>& isoleps,const vector<Lepton*>& nonisoleps,const vector<Jet>& jets,const map<TString,double>& weights){
  if(isoleps.size()&&nonisoleps.size()){
    Lepton* isolep=isoleps[0];
    Lepton* nonisolep=nonisoleps[0];
    TLorentzVector dilepton=(*isolep)+(*nonisolep);
    for(const auto& [weightname,w]:weights){
      TString suf=suffix+weightname;
      FillHist(pre+"jets"+suf,jets.size(),w,10,0,10);
      FillHist(pre+"leps"+suf,isoleps.size()+nonisoleps.size(),w,10,0,10);
      FillHist(pre+"isoleps"+suf,isoleps.size(),w,10,0,10);
      FillHist(pre+"nonisoleps"+suf,nonisoleps.size(),w,10,0,10);
      FillHist(pre+"dimass"+suf,dilepton.M(),w,220,60,500);
      FillHist(pre+"dipt"+suf,dilepton.Pt(),w,150,0,300);
      FillHist(pre+"dirap"+suf,dilepton.Rapidity(),w,100,-5,5);
      FillHist(pre+"l0pt"+suf,isolep->Pt(),w,500,0,1000);
      FillHist(pre+"l0eta"+suf,isolep->Eta(),w,60,-3,3);
      FillHist(pre+"l0riso"+suf,isolep->RelIso(),w,100,0,2);
      FillHist(pre+"l1pt"+suf,nonisolep->Pt(),w,500,0,1000);
      FillHist(pre+"l1eta"+suf,nonisolep->Eta(),w,60,-3,3);
      FillHist(pre+"l1riso"+suf,nonisolep->RelIso(),w,100,0,2);
      FillHist(pre+"lldelr"+suf,nonisolep->DeltaR(*isolep),w,80,0,4);     
      FillHist(pre+"lldelphi"+suf,nonisolep->DeltaPhi(*isolep),w,70,-3.5,3.5);     
      FillHist(pre+"lldeleta"+suf,fabs(nonisolep->Eta()-isolep->Eta()),w,100,-5,5);     	
      if(jets.size()){
	FillHist(pre+"j0pt"+suf,jets[0].Pt(),w,500,0,1000);
	FillHist(pre+"j0eta"+suf,jets[0].Eta(),w,60,-3,3);
	FillHist(pre+"jldelr"+suf,jets[0].DeltaR(*nonisolep),w,80,0,4);
      }
    }
  }    
}
