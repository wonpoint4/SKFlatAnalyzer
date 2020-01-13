#include "EfficiencyValidation.h"

EfficiencyValidation::EfficiencyValidation(){
}
EfficiencyValidation::~EfficiencyValidation(){
}
void EfficiencyValidation::executeEvent(){
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

    if(ev->PassTrigger("HLT_Ele32_WPTight_Gsf_v"))
      if(!IsDATA||DataStream.Contains("EGamma")) executeEventFromParameter("el2018",ev);
    
  }    

  delete ev;
}

void EfficiencyValidation::executeEventFromParameter(TString channelname,Event* ev){
  map< TString, std::vector<Muon> > map_muons;
  map< TString, std::vector<Electron> > map_electrons;
  //suffix, lepton vector, SF1, SF2, TriggerSF1, TriggerSF2
  map< TString, tuple<vector<Lepton*>,TString,TString,TString,TString> > map_leps;


  double lep0ptcut,lep1ptcut;
  if(channelname.Contains(TRegexp("^mm"))){
    lep0ptcut=20.;
    lep1ptcut=10.;

    map_muons["POGTight_PFIsoTight"]=SMPGetMuons("POGTightWithTightIso",0.0,2.4);
    map_muons["POGTight_TrkIsoLoose"]=SMPGetMuons("POGTightWithLooseTrkIso",0.0,2.4);

    switch(DataYear){
    case 2016: 
      map_leps["_POGTight_PFIsoTight"]=make_tuple(MakeLeptonPointerVector(map_muons["POGTight_PFIsoTight"]),"ID_SF_NUM_TightID_DEN_genTracks","ISO_SF_NUM_TightRelIso_DEN_TightIDandIPCut","LeadMu17_POGTight","TailMu8_POGTight");
      //map_leps["_POGTight_TrkIsoLoose_Q"]=make_tuple(MakeLeptonPointerVector(map_muons["POGTight_TrkIsoLoose"]),"IDISO_SF_TightID_trkIso_pt20_Q","Default","LeadMu17_TightID_trkIso_Q","TailMu8_TightID_trkIso_Q"); 
      break;
    case 2017: 
      map_leps["_POGTight_PFIsoTight"]=make_tuple(MakeLeptonPointerVector(map_muons["POGTight_PFIsoTight"]),"ID_SF_NUM_TightID_DEN_genTracks","ISO_SF_NUM_TightRelIso_DEN_TightIDandIPCut","LeadMu17_POGTight","TailMu8_POGTight");
      map_leps["_POGTight_TrkIsoLoose_Q"]=make_tuple(MakeLeptonPointerVector(map_muons["POGTight_TrkIsoLoose"]),"IDISO_SF_TightID_trkIso_pt20_Q","Default","LeadMu17_TightID_trkIso_Q","TailMu8_TightID_trkIso_Q");
      break;
    case 2018: 
      map_leps["_POGTight_PFIsoTight"]=make_tuple(MakeLeptonPointerVector(map_muons["POGTight_PFIsoTight"]),"ID_SF_NUM_TightID_DEN_genTracks","ISO_SF_NUM_TightRelIso_DEN_TightIDandIPCut","","");
      map_leps["_POGTight_TrkIsoLoose"]=make_tuple(MakeLeptonPointerVector(map_muons["POGTight_TrkIsoLoose"]),"ID_SF_NUM_TightID_DEN_genTracks","Default","","");
      break;
    default: cout<<"[EfficiencyValidation::executeEventFromParameter] wrong year"<<endl;exit(EXIT_FAILURE);break;
    }
  }else if(channelname.Contains(TRegexp("^ee"))){
    lep0ptcut=25.;
    lep1ptcut=15.;

    map_electrons["MediumID"]=SMPGetElectrons("passMediumID",0.0,2.5);
    map_electrons["MediumID_Selective"]=SMPGetElectrons("passMediumID_Selective",0.0,2.5);

    switch(DataYear){
    case 2016: 
      map_leps["_MediumID"]=make_tuple(MakeLeptonPointerVector(map_electrons["MediumID"]),"ID_SF_MediumID_pt10","","LeadEle23_MediumID","TailEle12_MediumID"); 
      map_leps["_MediumID_selective_Q"]=make_tuple(MakeLeptonPointerVector(map_electrons["MediumID_selective"]),"ID_SF_Selective_MediumID_pt10_Q","","Selective_LeadEle23_MediumID_Q","Selective_TailEle12_MediumID_Q"); 
      break;
    case 2017: 
      map_leps["_MediumID"]=make_tuple(MakeLeptonPointerVector(map_electrons["MediumID"]),"ID_SF_MediumID_pt10","","LeadEle23_MediumID","TailEle12_MediumID"); 
      map_leps["_MediumID_Q"]=make_tuple(MakeLeptonPointerVector(map_electrons["MediumID"]),"ID_SF_MediumID_pt10_Q","","LeadEle23_MediumID_Q","TailEle12_MediumID_Q"); 
      map_leps["_MediumID_selective_Q"]=make_tuple(MakeLeptonPointerVector(map_electrons["MediumID_selective"]),"ID_SF_Selective_MediumID_pt10_Q","","Selective_LeadEle23_MediumID_Q","Selective_TailEle12_MediumID_Q"); 
      break;
    case 2018: 
      map_leps["_MediumID_Q"]=make_tuple(MakeLeptonPointerVector(map_electrons["MediumID"]),"ID_SF_MediumID_pt10_Q","","Ele23_MediumID_Q","Ele12_MediumID_Q"); 
      map_leps["_MediumID_POGSF"]=make_tuple(MakeLeptonPointerVector(map_electrons["MediumID"]),"ID_SF_passMediumID","","Ele23_MediumID_Q","Ele12_MediumID_Q"); 
      break;
    default: 
      cout<<"[EfficiencyValidation::executeEventFromParameter] wrong year"<<endl;
      exit(EXIT_FAILURE);
    }
  }else if(channelname.Contains(TRegexp("^el"))){
    lep0ptcut=35.;
    lep1ptcut=15.;

    map_electrons["MediumID"]=SMPGetElectrons("passMediumID",0.0,2.5);
    map_electrons["TightID_Selective"]=SMPGetElectrons("passTightID_Selective",0.0,2.5);

    switch(DataYear){
    case 2018: 
      map_leps["_TightID_Selective_Q"]=make_tuple(MakeLeptonPointerVector(map_electrons["TightID_Selective"]),"ID_SF_TightID_Selective_pt10_Q","","Ele32_TightID_Selective_Q",""); 
      map_leps["_MediumID_Q"]=make_tuple(MakeLeptonPointerVector(map_electrons["MediumID"]),"ID_SF_MediumID_pt10_Q","","Ele32_MediumID_Q",""); 
      map_leps["_MediumID_POGSF"]=make_tuple(MakeLeptonPointerVector(map_electrons["MediumID"]),"ID_SF_passMediumID","","Ele32_MediumID_Q",""); 
      break;
    default: 
      cout<<"[EfficiencyValidation::executeEventFromParameter] wrong year"<<endl;
      exit(EXIT_FAILURE);
    }
  }else{
    cout<<"[EfficiencyValidation::executeEventFromParameter] wrong channelname"<<endl;
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
	  
	  map_weight["noefficiencySF"]=weight*PUreweight*prefireweight*zptcor;
	
	  map_weight["noRECOSF"]=weight*PUreweight*IDSF*ISOSF*triggerSF*prefireweight*zptcor;
	  map_weight["RECOSF_up"]=weight*PUreweight*RECOSF_up*IDSF*ISOSF*triggerSF*prefireweight*zptcor;
	  map_weight["RECOSF_down"]=weight*PUreweight*RECOSF_down*IDSF*ISOSF*triggerSF*prefireweight*zptcor;
	
	  map_weight["noIDSF"]=weight*PUreweight*RECOSF*ISOSF*triggerSF*prefireweight*zptcor;
	  map_weight["IDSF_up"]=weight*PUreweight*RECOSF*IDSF_up*ISOSF*triggerSF*prefireweight*zptcor;
	  map_weight["IDSF_down"]=weight*PUreweight*RECOSF*IDSF_down*ISOSF*triggerSF*prefireweight*zptcor;
	  
	  map_weight["noISOSF"]=weight*PUreweight*RECOSF*IDSF*triggerSF*prefireweight*zptcor;
	  map_weight["ISOSF_up"]=weight*PUreweight*RECOSF*IDSF*ISOSF_up*triggerSF*prefireweight*zptcor;
	  map_weight["ISOSF_down"]=weight*PUreweight*RECOSF*IDSF*ISOSF_down*triggerSF*prefireweight*zptcor;
	
	  map_weight["notriggerSF"]=weight*PUreweight*RECOSF*IDSF*ISOSF*prefireweight*zptcor;
	  map_weight["triggerSF_up"]=weight*PUreweight*RECOSF*IDSF*ISOSF*triggerSF_up*prefireweight*zptcor;
	  map_weight["triggerSF_down"]=weight*PUreweight*RECOSF*IDSF*ISOSF*triggerSF_down*prefireweight*zptcor;
	  
	  map_weight["nozptcor"]=weight*PUreweight*RECOSF*IDSF*ISOSF*triggerSF*prefireweight;
	}
	
	///////////////////////fill hists///////////////////////
	TLorentzVector dilepton=(*leps.at(0))+(*leps.at(1));
	double dimass=dilepton.M();
	double dipt=dilepton.Pt();
	if(dimass>=80&&dimass<100){
	  FillCutflow(channelname+"/"+tauprefix+"cutflow"+idsuffix,"m80to100",map_weight[""]);
	  FillHistsEfficiency(channelname+"/m80to100/"+prefix,idsuffix,leps,map_weight);
	  /*
	  if(dipt>50){
	    FillHist(channelname+"/"+tauprefix+"cutflow"+idsuffix,11.5,map_weight[""],20,0,20);
	    FillHistsEfficiency(channelname+"/m60to120_pt50/"+prefix,idsuffix,leps,map_weight);
	  }
	  */
	}
	/*
	else if(dimass>=120&&dimass<400){
	  FillHist(channelname+"/"+tauprefix+"cutflow"+idsuffix,12.5,map_weight[""],20,0,20);
	  FillHistsEfficiency(channelname+"/m120to400/"+prefix,idsuffix,leps,map_weight);
	  if(dipt>50){
	    FillHist(channelname+"/"+tauprefix+"cutflow"+idsuffix,13.5,map_weight[""],20,0,20);
	    FillHistsEfficiency(channelname+"/m120to400_pt50/"+prefix,idsuffix,leps,map_weight);
	  }
	}
	*/
      }
    }
  }
}
void EfficiencyValidation::FillHistsEfficiency(TString pre,TString suffix,const vector<Lepton*>& leps,const map<TString,double>& weights){
  for(const auto& element:weights){
    TString suf=suffix;
    if(element.first!="") suf+="_"+element.first;
    double w=element.second;
    
    //for leptons
    for(int i=0;i<(int)leps.size();i++){
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
