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
  tauprefix="";
  ////////////////////////check genlevel//////////////////
  zptcor=1.;
  if(IsDYSample){
    vector<Gen> gens=GetGens();
    Gen parton0,parton1,l0,l1;
    GetDYGenParticles(gens,parton0,parton1,l0,l1,true);
    if(abs(l0.PID())!=15){
      TLorentzVector genZ=(l0+l1);
      zptcor*=GetZptWeight(genZ.Pt(),genZ.Rapidity(),abs(l0.PID())==13?Lepton::Flavour::MUON:Lepton::Flavour::ELECTRON);
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
    if(ev->PassTrigger(doublemuontrigger))
      if(!IsDATA||DataStream.Contains("DoubleMuon")) executeEventFromParameter("mm2016",ev);
    if(ev->PassTrigger("HLT_IsoMu24_v")||ev->PassTrigger("HLT_IsoTkMu24_v"))
      if(!IsDATA||DataStream.Contains("SingleMuon")) executeEventFromParameter("mu2016",ev);
    if(ev->PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"))
      if(!IsDATA||DataStream.Contains("DoubleEG")) executeEventFromParameter("ee2016",ev);
    if(ev->PassTrigger("HLT_Ele27_WPTight_Gsf_v"))
      if(!IsDATA||DataStream.Contains("SingleElectron")) executeEventFromParameter("el2016",ev);
  }else if(DataYear==2017){
    if(ev->PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v"))
      if(!IsDATA||DataStream.Contains("DoubleMuon")) executeEventFromParameter("mm2017",ev);
    if(ev->PassTrigger("HLT_IsoMu27_v"))
      if(!IsDATA||DataStream.Contains("SingleMuon")) executeEventFromParameter("mu2017",ev);    
    if(ev->PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"))
      if(!IsDATA||DataStream.Contains("DoubleEG")) executeEventFromParameter("ee2017",ev);
    if(ev->PassTrigger("HLT_Ele32_WPTight_Gsf_L1DoubleEG_v"))
      if(!IsDATA||DataStream.Contains("SingleElectron")) executeEventFromParameter("el2017",ev);
  }else if(DataYear==2018){
    if(ev->PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"))
      if(!IsDATA||DataStream.Contains("DoubleMuon")) executeEventFromParameter("mm2018",ev);
    if(ev->PassTrigger("HLT_IsoMu24_v"))
      if(!IsDATA||DataStream.Contains("SingleMuon")) executeEventFromParameter("mu2018",ev);    
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
    map_muons["_MediumID_trkIsoLoose"]=SMPGetMuons("POGMediumWithLooseTrkIso",0.0,2.4);
    map_muons["_TightID_PFIsoTight"]=SMPGetMuons("POGTightWithTightIso",0.0,2.4);
    map_leps["_MediumID_trkIsoLoose_Q"]=make_tuple(MakeLeptonPointerVector(map_muons["_MediumID_trkIsoLoose"]),"IDISO_SF_MediumID_trkIsoLoose_Q","","Mu17Leg1_MediumID_trkIsoLoose_Q","Mu8Leg2_MediumID_trkIsoLoose_Q");
    if(DataYear==2016) map_leps["_TightID_PFIsoTight"]=make_tuple(MakeLeptonPointerVector(map_muons["_TightID_PFIsoTight"]),"ID_SF_NUM_TightID_DEN_genTracks","ISO_SF_NUM_TightRelIso_DEN_TightIDandIPCut","Mu17Leg1_POGTight","Mu8Leg2_POGTight");
  }else if(channelname.Contains(TRegexp("^mu"))){
    map_muons["_MediumID_trkIsoLoose"]=SMPGetMuons("POGMediumWithLooseTrkIso",0.0,2.4);
    switch(DataYear){
    case 2016:
      lep0ptcut=27.;
      lep1ptcut=10.;
      map_leps["_MediumID_trkIsoLoose_Q"]=make_tuple(MakeLeptonPointerVector(map_muons["_MediumID_trkIsoLoose"]),"IDISO_SF_MediumID_trkIsoLoose_Q","","IsoMu24_MediumID_trkIsoLoose_Q","");
      break;
    case 2017:
      lep0ptcut=30.;
      lep1ptcut=10.;
      map_leps["_MediumID_trkIsoLoose_Q"]=make_tuple(MakeLeptonPointerVector(map_muons["_MediumID_trkIsoLoose"]),"IDISO_SF_MediumID_trkIsoLoose_Q","","IsoMu27_MediumID_trkIsoLoose_Q","");    
      break;
    case 2018:
      lep0ptcut=27.;
      lep1ptcut=10.;
      map_leps["_MediumID_trkIsoLoose_Q"]=make_tuple(MakeLeptonPointerVector(map_muons["_MediumID_trkIsoLoose"]),"IDISO_SF_MediumID_trkIsoLoose_Q","","IsoMu24_MediumID_trkIsoLoose_Q","");    
      break;
    default: 
      cout<<"[EfficiencyValidation::executeEventFromParameter] wrong year"<<endl;
      exit(EXIT_FAILURE);
    }
  }else if(channelname.Contains(TRegexp("^ee"))){
    lep0ptcut=25.;
    lep1ptcut=15.;
    map_electrons["_MediumID_noroccor"]=SMPGetElectrons("passMediumID",0.0,2.4);
    map_electrons["_MediumID"]=ElectronEnergyCorrection(map_electrons["_MediumID_noroccor"],0,0);
    map_electrons["_MediumID_noEcor"]=ElectronEnergyCorrection(SMPGetElectrons("passMediumID",0.0,2.4),-1,0);
    map_leps["_MediumID_Q"]=make_tuple(MakeLeptonPointerVector(map_electrons["_MediumID"]),"ID_SF_MediumID_Q","","Ele23Leg1_MediumID_Q","Ele12Leg2_MediumID_Q"); 
    map_leps["_MediumID_noEcor_Q"]=make_tuple(MakeLeptonPointerVector(map_electrons["_MediumID_noEcor"]),"ID_SF_MediumID_Q","","Ele23Leg1_MediumID_Q","Ele12Leg2_MediumID_Q"); 

    switch(DataYear){
    case 2016: 
      map_leps["_MediumID_old"]=make_tuple(MakeLeptonPointerVector(map_electrons["_MediumID"]),"ID_SF_MediumID_pt10","","Ele23Leg1_MediumID","Ele12Leg2_MediumID");
      break;
    case 2017: 
      map_leps["_MediumID_old_Q"]=make_tuple(MakeLeptonPointerVector(map_electrons["_MediumID"]),"ID_SF_MediumID_pt10_Q","","Ele23Leg1_MediumID_old_Q","Ele12Leg2_MediumID_old_Q"); 
      break;
    case 2018: 
      break;
    default: 
      cout<<"[EfficiencyValidation::executeEventFromParameter] wrong year"<<endl;
      exit(EXIT_FAILURE);
    }
  }else if(channelname.Contains(TRegexp("^el"))){
    map_electrons["_MediumID_noroccor"]=SMPGetElectrons("passMediumID",0.0,2.4);
    map_electrons["_MediumID"]=ElectronEnergyCorrection(map_electrons["_MediumID_noroccor"],0,0);
    map_electrons["_MediumID_noEcor"]=ElectronEnergyCorrection(SMPGetElectrons("passMediumID",0.0,2.4),-1,0);
    map_electrons["_TightID_Selective_noroccor"]=SMPGetElectrons("passMediumID_Selective",0.0,2.4);
    map_electrons["_TightID_Selective"]=ElectronEnergyCorrection(map_electrons["_TightID_Selective_noroccor"],0,0);

    switch(DataYear){
    case 2016:
      lep0ptcut=30.;
      lep1ptcut=10.;
      map_leps["_MediumID_Q"]=make_tuple(MakeLeptonPointerVector(map_electrons["_MediumID"]),"ID_SF_MediumID_Q","","Ele27_MediumID_Q",""); 
      map_leps["_MediumID_noEcor_Q"]=make_tuple(MakeLeptonPointerVector(map_electrons["_MediumID_noEcor"]),"ID_SF_MediumID_Q","","Ele27_MediumID_Q",""); 
      map_leps["_TightID_Selective_Q"]=make_tuple(MakeLeptonPointerVector(map_electrons["_TightID_Selective"]),"ID_SF_TightID_Selective_Q","","Ele27_TightID_Selective_Q",""); 
      break;
    case 2017:
      lep0ptcut=35.;
      lep1ptcut=10.;
      map_leps["_MediumID_Q"]=make_tuple(MakeLeptonPointerVector(map_electrons["_MediumID"]),"ID_SF_MediumID_Q","","Ele32_MediumID_Q",""); 
      map_leps["_MediumID_noEcor_Q"]=make_tuple(MakeLeptonPointerVector(map_electrons["_MediumID_noEcor"]),"ID_SF_MediumID_Q","","Ele32_MediumID_Q",""); 
      map_leps["_TightID_Selective_Q"]=make_tuple(MakeLeptonPointerVector(map_electrons["_TightID_Selective"]),"ID_SF_TightID_Selective_Q","","Ele32_TightID_Selective_Q",""); 
      break;
    case 2018: 
      lep0ptcut=35.;
      lep1ptcut=10.;
      map_leps["_MediumID_Q"]=make_tuple(MakeLeptonPointerVector(map_electrons["_MediumID"]),"ID_SF_MediumID_Q","","Ele32_MediumID_Q",""); 
      map_leps["_MediumID_noEcor_Q"]=make_tuple(MakeLeptonPointerVector(map_electrons["_MediumID_noEcor"]),"ID_SF_MediumID_Q","","Ele32_MediumID_Q",""); 
      map_leps["_TightID_Selective_Q"]=make_tuple(MakeLeptonPointerVector(map_electrons["_TightID_Selective"]),"ID_SF_TightID_Selective_Q","","Ele32_TightID_Selective_Q",""); 
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

  /////////////////// zpt correction //////////////////////
  totalweight*=zptcor;
  FillCutflow(channelname+"/"+tauprefix+"cutflow","zptcor",totalweight);
  
  //////////////////////Z0 weight///////////////////////
  double z0weight=1.;
  if(!IsDATA){
    z0weight=GetZ0Weight(vertex_Z);
  }
  totalweight*=z0weight;
  FillCutflow(channelname+"/"+tauprefix+"cutflow","z0",totalweight);

  ///////////////////////lepton selection///////////////////////
  for(const auto& element_leps: map_leps){
    TString idsuffix=element_leps.first;
    auto const& leps=get<0>(element_leps.second);

    FillHist(channelname+"/"+tauprefix+"nlepton"+idsuffix,leps.size(),totalweight,10,0,10);
    if(leps.size()==2){
      TString prefix=tauprefix;
      if(leps.at(0)->Charge()*leps.at(1)->Charge()>0)
	//prefix="ss_"+prefix;
	return;
      FillCutflow(channelname+"/"+tauprefix+"cutflow"+idsuffix,"dilepton",totalweight);
      if(leps.at(0)->Pt()>lep0ptcut&&leps.at(1)->Pt()>lep1ptcut){
	FillCutflow(channelname+"/"+tauprefix+"cutflow"+idsuffix,"LepPtCut",totalweight);
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
	map_weight[""]=weight*PUreweight*RECOSF*IDSF*ISOSF*triggerSF*prefireweight*zptcor*z0weight;
	if(!IsDATA){
	  map_weight["_noweight"]=weight;
	  map_weight["_noPUreweight"]=weight*RECOSF*IDSF*ISOSF*triggerSF*prefireweight*zptcor*z0weight;
	  map_weight["_noprefireweight"]=weight*PUreweight*RECOSF*IDSF*ISOSF*triggerSF*zptcor*z0weight;
	  
	  map_weight["_noRECOSF"]=weight*PUreweight*IDSF*ISOSF*triggerSF*prefireweight*zptcor*z0weight;
	  map_weight["_RECOSF_up"]=weight*PUreweight*RECOSF_up*IDSF*ISOSF*triggerSF*prefireweight*zptcor*z0weight;
	  map_weight["_RECOSF_down"]=weight*PUreweight*RECOSF_down*IDSF*ISOSF*triggerSF*prefireweight*zptcor*z0weight;
	
	  map_weight["_noIDSF"]=weight*PUreweight*RECOSF*ISOSF*triggerSF*prefireweight*zptcor*z0weight;
	  map_weight["_IDSF_up"]=weight*PUreweight*RECOSF*IDSF_up*ISOSF*triggerSF*prefireweight*zptcor*z0weight;
	  map_weight["_IDSF_down"]=weight*PUreweight*RECOSF*IDSF_down*ISOSF*triggerSF*prefireweight*zptcor*z0weight;
	  
	  map_weight["_noISOSF"]=weight*PUreweight*RECOSF*IDSF*triggerSF*prefireweight*zptcor*z0weight;
	  map_weight["_ISOSF_up"]=weight*PUreweight*RECOSF*IDSF*ISOSF_up*triggerSF*prefireweight*zptcor*z0weight;
	  map_weight["_ISOSF_down"]=weight*PUreweight*RECOSF*IDSF*ISOSF_down*triggerSF*prefireweight*zptcor*z0weight;
	
	  map_weight["_notriggerSF"]=weight*PUreweight*RECOSF*IDSF*ISOSF*prefireweight*zptcor*z0weight;
	  map_weight["_triggerSF_up"]=weight*PUreweight*RECOSF*IDSF*ISOSF*triggerSF_up*prefireweight*zptcor*z0weight;
	  map_weight["_triggerSF_down"]=weight*PUreweight*RECOSF*IDSF*ISOSF*triggerSF_down*prefireweight*zptcor*z0weight;
	  
	  map_weight["_noefficiencySF"]=weight*PUreweight*prefireweight*zptcor*z0weight;
	  map_weight["_noefficiencySF_nozptcor"]=weight*PUreweight*prefireweight*z0weight;
	  map_weight["_noz0weight"]=weight*PUreweight*RECOSF*IDSF*ISOSF*triggerSF*prefireweight*zptcor;
	  map_weight["_nozptcor"]=weight*PUreweight*RECOSF*IDSF*ISOSF*triggerSF*prefireweight*z0weight;
	}
	
	///////////////////////fill hists///////////////////////
	TLorentzVector dilepton=(*leps.at(0))+(*leps.at(1));
	double dimass=dilepton.M();
	double dipt=dilepton.Pt();
	if(dimass>=60&&dimass<120){
	  FillCutflow(channelname+"/"+tauprefix+"cutflow"+idsuffix,"m60to120",map_weight[""]);
	  FillHistsEfficiency(channelname+"/m60to120/"+prefix,idsuffix,leps,map_weight);
	  if(dimass>=80&&dimass<100){
	    FillCutflow(channelname+"/"+tauprefix+"cutflow"+idsuffix,"m80to100",map_weight[""]);
	    FillHistsEfficiency(channelname+"/m80to100/"+prefix,idsuffix,leps,map_weight);
	    if(leps.at(1)->Pt()<25){
	      FillCutflow(channelname+"/"+tauprefix+"cutflow"+idsuffix,"l1pt<25",map_weight[""]);
	      FillHistsEfficiency(channelname+"/m80to100_l1pt/"+prefix,idsuffix,leps,map_weight);
	    } 
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
