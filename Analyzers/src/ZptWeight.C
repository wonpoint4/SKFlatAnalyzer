#include "ZptWeight.h"

ZptWeight::ZptWeight(){
}
ZptWeight::~ZptWeight(){
}
void ZptWeight::executeEvent(){
  tauprefix="";
  ////////////////////////check genlevel//////////////////
  zptcor=1.;
  if(IsDYSample){
    vector<Gen> gens=GetGens();
    int parton0,parton1,hardl0,hardl1,l0,l1;
    vector<int> photons;
    GetGenIndex(gens,parton0,parton1,hardl0,hardl1,l0,l1,photons);
    if(abs(gens[hardl0].PID())!=15){
      Gen genhardl0=gens[hardl0],genhardl1=gens[hardl1],genl0=gens[l0],genl1=gens[l1],genphotons;
      for(unsigned int i=0;i<photons.size();i++) genphotons+=gens[photons[i]];
      TLorentzVector genZ=(genl0+genl1+genphotons);
      zptcor*=GetZptWeight(genZ.Pt(),genZ.Rapidity(),abs(genhardl0.PID())==13?Lepton::Flavour::MUON:Lepton::Flavour::ELECTRON);
      FillHist(Form("%s%d/%s",abs(genhardl0.PID())==13?"mm":"ee",DataYear,"/gen_diptdirap"),genZ.Pt(),fabs(genZ.Rapidity()),weight_norm_1invpb*gen_weight*zptcor,ptbinnum,ptbin,rapbinnum,rapbin);
      FillHist(Form("%s%d/%s",abs(genhardl0.PID())==13?"mm":"ee",DataYear,"/gen_diptdirap_nozptcor"),genZ.Pt(),fabs(genZ.Rapidity()),weight_norm_1invpb*gen_weight,ptbinnum,ptbin,rapbinnum,rapbin);
    }else{
      tauprefix="tau_";
    }
  }

  if(!PassMETFilter()) return;

  Event* ev=new Event;
  *ev=GetEvent();
  if(DataYear==2016){
    vector<TString> doublemuontrigger;
    doublemuontrigger.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");
    doublemuontrigger.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v");
    doublemuontrigger.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
    doublemuontrigger.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");
    if(ev->PassTrigger(doublemuontrigger)){
      if(!IsDATA||DataStream.Contains("DoubleMuon")) executeEventFromParameter("mm2016",ev);
    }
    if(ev->PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")){
      if(!IsDATA||DataStream.Contains("DoubleEG")) executeEventFromParameter("ee2016",ev);
    }
  }else if(DataYear==2017){
    if(ev->PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v")){
      if(!IsDATA||DataStream.Contains("DoubleMuon")) executeEventFromParameter("mm2017",ev);
    }
    if(ev->PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v")){
      if(!IsDATA||DataStream.Contains("DoubleEG")) executeEventFromParameter("ee2017",ev);
    }
  }else if(DataYear==2018){
    if(ev->PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"))
      if(!IsDATA||DataStream.Contains("DoubleMuon")) executeEventFromParameter("mm2018",ev);
    
    if(ev->PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"))
      if(!IsDATA||DataStream.Contains("EGamma")) executeEventFromParameter("ee2018",ev);
    
  }    

  delete ev;
}

void ZptWeight::executeEventFromParameter(TString channelname,Event* ev){
  map< TString, std::vector<Muon> > map_muons;
  map< TString, std::vector<Electron> > map_electrons;
  //suffix, lepton vector, SF1, SF2, TriggerSF1, TriggerSF2
  map< TString, tuple<vector<Lepton*>,TString,TString,TString,TString> > map_leps;


  double lep0ptcut,lep1ptcut;
  if(channelname.Contains(TRegexp("^mm"))){
    lep0ptcut=25.;
    lep1ptcut=15.;

    map_muons[""]=SMPGetMuons("POGTightWithTightIso",0.0,2.4);

    switch(DataYear){
    case 2016: 
      map_leps[""]=make_tuple(MakeLeptonPointerVector(map_muons[""]),"ID_SF_NUM_TightID_DEN_genTracks","ISO_SF_NUM_TightRelIso_DEN_TightIDandIPCut","LeadMu17_POGTight","TailMu8_POGTight");
      break;
    case 2017: 
      map_leps[""]=make_tuple(MakeLeptonPointerVector(map_muons[""]),"ID_SF_NUM_TightID_DEN_genTracks","ISO_SF_NUM_TightRelIso_DEN_TightIDandIPCut","LeadMu17_POGTight","TailMu8_POGTight");
      break;
    case 2018: 
      map_leps[""]=make_tuple(MakeLeptonPointerVector(map_muons[""]),"ID_SF_NUM_TightID_DEN_genTracks","ISO_SF_NUM_TightRelIso_DEN_TightIDandIPCut","","");
      break;
    default: cout<<"[ZptWeight::executeEventFromParameter] wrong year"<<endl;exit(EXIT_FAILURE);break;
    }
  }else if(channelname.Contains(TRegexp("^ee"))){
    lep0ptcut=30.;
    lep1ptcut=20.;

    map_electrons[""]=SMPGetElectrons("passMediumID",0.0,2.5);

    switch(DataYear){
    case 2016: 
      map_leps[""]=make_tuple(MakeLeptonPointerVector(map_electrons[""]),"ID_SF_MediumID_pt10","","LeadEle23_MediumID","TailEle12_MediumID"); 
      break;
    case 2017: 
      map_leps[""]=make_tuple(MakeLeptonPointerVector(map_electrons[""]),"ID_SF_MediumID_pt10","","LeadEle23_MediumID","TailEle12_MediumID"); 
      break;
    case 2018: 
      map_leps[""]=make_tuple(MakeLeptonPointerVector(map_electrons[""]),"ID_SF_MediumID_pt10_Q","","Ele23_MediumID_Q","Ele12_MediumID_Q"); 
      break;
    default: 
      cout<<"[ZptWeight::executeEventFromParameter] wrong year"<<endl;
      exit(EXIT_FAILURE);
    }
  }else{
    cout<<"[ZptWeight::executeEventFromParameter] wrong channelname"<<endl;
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

  totalweight*=zptcor;
  FillCutflow(channelname+"/"+tauprefix+"cutflow","zptcor",totalweight);

  ///////////////////////lepton selection///////////////////////
  for(const auto& element_leps: map_leps){
    TString idsuffix=element_leps.first;
    auto const& leps=get<0>(element_leps.second);

    FillHist(channelname+"/"+tauprefix+"nlepton"+idsuffix,leps.size(),totalweight,10,0,10);
    if(leps.size()==2){
      TString prefix=(leps.at(0)->Charge()*leps.at(1)->Charge()>0?"ss_":"")+tauprefix;
      FillCutflow(channelname+"/"+tauprefix+"cutflow"+idsuffix,"dilepton",totalweight);
      if(leps.at(0)->Pt()>lep0ptcut&&leps.at(1)->Pt()>lep1ptcut){
	FillCutflow(channelname+"/"+tauprefix+"cutflow"+idsuffix,"OS",totalweight);
	/////////////////efficiency scale factors///////////////////
	double IDSF=1.,IDSF_up=1.,IDSF_down=1.;
	double ISOSF=1.,ISOSF_up=1.,ISOSF_down=1.;
	double RECOSF=1.,RECOSF_up=1.,RECOSF_down=1.;
	if(!IsDATA){
	  TString LeptonIDSF_key=get<1>(element_leps.second);
	  TString LeptonISOSF_key=get<2>(element_leps.second);
	  for(unsigned int i=0;i<leps.size();i++){	  
            if(leps[i]->LeptonFlavour()==Lepton::ELECTRON){
	      double this_pt,this_eta;
              this_pt=leps.at(i)->Pt();
              this_eta=((Electron*)leps.at(i))->scEta();

	      double this_RECOSF=mcCorr->ElectronReco_SF(this_eta,this_pt,0);
	      double this_RECOSF_up=mcCorr->ElectronReco_SF(this_eta,this_pt,1);
	      double this_RECOSF_down=mcCorr->ElectronReco_SF(this_eta,this_pt,-1);
	      RECOSF*=this_RECOSF; RECOSF_up*=this_RECOSF_up; RECOSF_down*=this_RECOSF_down;
            }
	  
	    double this_IDSF=Lepton_SF(LeptonIDSF_key,leps.at(i),0);
	    double this_IDSF_up=Lepton_SF(LeptonIDSF_key,leps.at(i),1);
	    double this_IDSF_down=Lepton_SF(LeptonIDSF_key,leps.at(i),-1);
	    IDSF*=this_IDSF; IDSF_up*=this_IDSF_up; IDSF_down*=this_IDSF_down;
	  
	    double this_ISOSF=Lepton_SF(LeptonISOSF_key,leps.at(i),0);
	    double this_ISOSF_up=Lepton_SF(LeptonISOSF_key,leps.at(i),1);
	    double this_ISOSF_down=Lepton_SF(LeptonISOSF_key,leps.at(i),-1);
	    ISOSF*=this_ISOSF; ISOSF_up*=this_ISOSF_up; ISOSF_down*=this_ISOSF_down;
	  }
	}
	FillCutflow(channelname+"/"+tauprefix+"cutflow"+idsuffix,"RECOSF",totalweight*RECOSF);
	FillCutflow(channelname+"/"+tauprefix+"cutflow"+idsuffix,"IDSF",totalweight*RECOSF*IDSF);
	FillCutflow(channelname+"/"+tauprefix+"cutflow"+idsuffix,"ISOSF",totalweight*RECOSF*IDSF*ISOSF);
      
	double triggerSF=1.,triggerSF_up=1.,triggerSF_down=1.;
	if(!IsDATA){
	  if(channelname.Contains(TRegexp("^ee"))||channelname.Contains(TRegexp("^mm"))){
	    TString triggerSF_key0=get<3>(element_leps.second);
	    TString triggerSF_key1=get<4>(element_leps.second);
	    triggerSF*=DileptonTrigger_SF(triggerSF_key0,triggerSF_key1,leps,0);
	    triggerSF_up*=DileptonTrigger_SF(triggerSF_key0,triggerSF_key1,leps,1);
	    triggerSF_down*=DileptonTrigger_SF(triggerSF_key0,triggerSF_key1,leps,-1);
	  }else if(channelname.Contains(TRegexp("^el"))||channelname.Contains(TRegexp("^mu"))){
	    TString triggerSF_key=get<3>(element_leps.second);
	    triggerSF*=LeptonTrigger_SF(triggerSF_key,leps,0);
	    triggerSF_up*=LeptonTrigger_SF(triggerSF_key,leps,1);
	    triggerSF_down*=LeptonTrigger_SF(triggerSF_key,leps,-1);
	  }	    
	}
	FillCutflow(channelname+"/"+tauprefix+"cutflow"+idsuffix,"triggerSF",totalweight*RECOSF*IDSF*ISOSF*triggerSF);

	///////////////////////weight systematics//////////////////
	map<TString,double> map_weight;
	map_weight[""]=weight*PUreweight*RECOSF*IDSF*ISOSF*triggerSF*prefireweight*zptcor;
	if(!IsDATA){	  
	  map_weight["_nozptcor"]=weight*PUreweight*RECOSF*IDSF*ISOSF*triggerSF*prefireweight;
	}
	
	///////////////////////fill hists///////////////////////
	TLorentzVector dilepton=(*leps.at(0))+(*leps.at(1));
	double dimass=dilepton.M();
	if(dimass>=80&&dimass<100){
	  FillCutflow(channelname+"/"+tauprefix+"cutflow"+idsuffix,"m80to100",map_weight[""]);
	  FillHistsZptWeight(channelname+"/m80to100/"+prefix,idsuffix,leps,map_weight);
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
