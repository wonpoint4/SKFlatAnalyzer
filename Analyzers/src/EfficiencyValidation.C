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
    if(abs(gens[hardl0].PID())==15){
      tauprefix="tau_";
    }
  }

  if(!PassMETFilter()) return;

  Event* ev=new Event;
  *ev=GetEvent();
  TString channelname,muontrigger,electrontrigger;
  if(DataYear==2016){
    vector<TString> muontrigger;
    muontrigger.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");
    muontrigger.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v");
    muontrigger.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
    muontrigger.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");
    electrontrigger="HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";
    if(ev->PassTrigger(muontrigger)){
      channelname="muon2016";
      if(!IsDATA||DataStream.Contains("DoubleMuon")) executeEventFromParameter(channelname,ev);
    }
    if(ev->PassTrigger(electrontrigger)){
      channelname="electron2016";
      if(!IsDATA||DataStream.Contains("DoubleEG")) executeEventFromParameter(channelname,ev);
    }
  }else if(DataYear==2017){
    muontrigger="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v";
    electrontrigger="HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v";
    if(ev->PassTrigger(muontrigger)){
      channelname="muon2017";
      if(!IsDATA||DataStream.Contains("DoubleMuon")) executeEventFromParameter(channelname,ev);
    }
    if(ev->PassTrigger(electrontrigger)){
      channelname="electron2017";
      if(!IsDATA||DataStream.Contains("DoubleEG")) executeEventFromParameter(channelname,ev);
    }
  }else if(DataYear==2018){
    muontrigger="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v";
    electrontrigger="HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";
    if(ev->PassTrigger(muontrigger)){
      channelname="muon2018";
      if(!IsDATA||DataStream.Contains("DoubleMuon")) executeEventFromParameter(channelname,ev);
    }
    if(ev->PassTrigger(electrontrigger)){
      channelname="electron2018";
      if(!IsDATA||DataStream.Contains("EGamma")) executeEventFromParameter(channelname,ev);
    }
  }    

  delete ev;
}

void EfficiencyValidation::executeEventFromParameter(TString channelname,Event* ev){
  map< TString, std::vector<Muon> > muons;
  map< TString, std::vector<Electron> > electrons;

  double lep0ptcut,lep1ptcut;
  //TString LeptonIDSF_key,LeptonISOSF_key,triggerSF_key0,triggerSF_key1;
  map< TString, tuple<vector<Lepton*>,TString,TString,TString,TString> > map_leps;
  if(channelname.Contains("muon")){
    lep0ptcut=20.;
    lep1ptcut=10.;

    muons["POGTight_PFIsoTight"]=SMPGetMuons("POGTightWithTightIso",0.0,2.4);
    muons["POGTight_TrkIsoLoose"]=SMPGetMuons("POGTightWithLooseTrkIso",0.0,2.4);

    switch(DataYear){
    case 2016: 
      map_leps["POGTight_PFIsoTight"]=make_tuple(MakeLeptonPointerVector(muons["POGTight_PFIsoTight"]),"ID_SF_NUM_TightID_DEN_genTracks","ISO_SF_NUM_TightRelIso_DEN_TightIDandIPCut","LeadMu17_POGTight","TailMu8_POGTight");
      //map_leps["POGTight_TrkIsoLoose_Q"]=make_tuple(MakeLeptonPointerVector(muons["POGTight_TrkIsoLoose"]),"IDISO_SF_TightID_trkIso_pt20_Q","Default","LeadMu17_TightID_trkIso_Q","TailMu8_TightID_trkIso_Q"); 
      break;
    case 2017: 
      map_leps["POGTight_PFIsoTight"]=make_tuple(MakeLeptonPointerVector(muons["POGTight_PFIsoTight"]),"ID_SF_NUM_TightID_DEN_genTracks","ISO_SF_NUM_TightRelIso_DEN_TightIDandIPCut","LeadMu17_POGTight","TailMu8_POGTight");
      map_leps["POGTight_TrkIsoLoose_Q"]=make_tuple(MakeLeptonPointerVector(muons["POGTight_TrkIsoLoose"]),"IDISO_SF_TightID_trkIso_pt20_Q","Default","LeadMu17_TightID_trkIso_Q","TailMu8_TightID_trkIso_Q");
      break;
    case 2018: 
      map_leps["POGTight_PFIsoTight"]=make_tuple(MakeLeptonPointerVector(muons["POGTight_PFIsoTight"]),"ID_SF_NUM_TightID_DEN_genTracks","ISO_SF_NUM_TightRelIso_DEN_TightIDandIPCut","","");
      map_leps["POGTight_TrkIsoLoose"]=make_tuple(MakeLeptonPointerVector(muons["POGTight_TrkIsoLoose"]),"ID_SF_NUM_TightID_DEN_genTracks","Default","","");
      break;
    default: cout<<"[EfficiencyValidation::executeEventFromParameter] wrong year"<<endl;exit(EXIT_FAILURE);break;
    }
 }else if(channelname.Contains("electron")){
    lep0ptcut=25.;
    lep1ptcut=15.;

    electrons["MediumID"]=SMPGetElectrons("passMediumID",0.0,2.5);
    electrons["MediumID_selective"]=SMPGetElectrons("passMediumID_selective",0.0,2.5);

    switch(DataYear){
    case 2016: 
      map_leps["MediumID"]=make_tuple(MakeLeptonPointerVector(electrons["MediumID"]),"ID_SF_MediumID_pt10","","LeadEle23_MediumID","TailEle12_MediumID"); 
      map_leps["MediumID_selective_Q"]=make_tuple(MakeLeptonPointerVector(electrons["MediumID_selective"]),"ID_SF_Selective_MediumID_pt10_Q","","Selective_LeadEle23_MediumID_Q","Selective_TailEle12_MediumID_Q"); 
      break;
    case 2017: 
      map_leps["MediumID"]=make_tuple(MakeLeptonPointerVector(electrons["MediumID"]),"ID_SF_MediumID_pt10","","LeadEle23_MediumID","TailEle12_MediumID"); 
      map_leps["MediumID_Q"]=make_tuple(MakeLeptonPointerVector(electrons["MediumID"]),"ID_SF_MediumID_pt10_Q","","LeadEle23_MediumID_Q","TailEle12_MediumID_Q"); 
      map_leps["MediumID_selective_Q"]=make_tuple(MakeLeptonPointerVector(electrons["MediumID_selective"]),"ID_SF_Selective_MediumID_pt10_Q","","Selective_LeadEle23_MediumID_Q","Selective_TailEle12_MediumID_Q"); 
      break;
    case 2018: 
      map_leps["MediumID"]=make_tuple(MakeLeptonPointerVector(electrons["MediumID"]),"ID_SF_passMediumID","","",""); 
      map_leps["MediumID_selective"]=make_tuple(MakeLeptonPointerVector(electrons["MediumID_selective"]),"","","",""); 
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
  FillHist(channelname+"/"+tauprefix+"cutflow",0.5,totalweight,20,0,20);

  /////////////////PUreweight///////////////////
  double PUreweight=1.,PUreweight_up=1.,PUreweight_down=1.;
  if(!IsDATA&&DataYear==2016){
    PUreweight=mcCorr->GetPileUpWeight(nPileUp,0);
    PUreweight_up=mcCorr->GetPileUpWeight(nPileUp,1);
    PUreweight_down=mcCorr->GetPileUpWeight(nPileUp,-1);
  }
  totalweight*=PUreweight;
  FillHist(channelname+"/"+tauprefix+"cutflow",1.5,totalweight,20,0,20);
  
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
  FillHist(channelname+"/"+tauprefix+"cutflow",2.5,totalweight,20,0,20);

  totalweight*=zptcor;
  FillHist(channelname+"/"+tauprefix+"cutflow",3.5,totalweight,20,0,20);

  ///////////////////////lepton selection///////////////////////
  for(const auto& element_leps: map_leps){
    TString idsuffix="_"+element_leps.first;
    auto const& leps=get<0>(element_leps.second);

    FillHist(channelname+"/"+tauprefix+"nlepton"+idsuffix,leps.size(),totalweight,10,0,10);
    if(leps.size()==2){
      TString prefix=(leps.at(0)->Charge()*leps.at(1)->Charge()>0?"ss_":"")+tauprefix;
      FillHist(channelname+"/"+tauprefix+"cutflow"+idsuffix,4.5,totalweight,20,0,20);
      if(leps.at(0)->Pt()>lep0ptcut&&leps.at(1)->Pt()>lep1ptcut){
	FillHist(channelname+"/"+tauprefix+"cutflow"+idsuffix,5.5,totalweight,20,0,20);
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
	FillHist(channelname+"/"+tauprefix+"cutflow"+idsuffix,6.5,totalweight*RECOSF,20,0,20);
	FillHist(channelname+"/"+tauprefix+"cutflow"+idsuffix,7.5,totalweight*RECOSF*IDSF,20,0,20);
	FillHist(channelname+"/"+tauprefix+"cutflow"+idsuffix,8.5,totalweight*RECOSF*IDSF*ISOSF,20,0,20);
      
	double triggerSF=1.,triggerSF_up=1.,triggerSF_down=1.;
	if(!IsDATA){
	  TString triggerSF_key0=get<3>(element_leps.second);
	  TString triggerSF_key1=get<4>(element_leps.second);
	  triggerSF*=DileptonTrigger_SF(triggerSF_key0,triggerSF_key1,leps,0);
	  triggerSF_up*=DileptonTrigger_SF(triggerSF_key0,triggerSF_key1,leps,1);
	  triggerSF_down*=DileptonTrigger_SF(triggerSF_key0,triggerSF_key1,leps,-1);
	}
	FillHist(channelname+"/"+tauprefix+"cutflow"+idsuffix,9.5,totalweight*RECOSF*IDSF*ISOSF*triggerSF,20,0,20);

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
	if(dimass>=60&&dimass<120){
	  FillHist(channelname+"/"+tauprefix+"cutflow"+idsuffix,10.5,map_weight[""],20,0,20);
	  FillHistsEfficiency(channelname+"/m60to120/"+prefix,idsuffix,leps,map_weight);
	  if(dipt>50){
	    FillHist(channelname+"/"+tauprefix+"cutflow"+idsuffix,11.5,map_weight[""],20,0,20);
	    FillHistsEfficiency(channelname+"/m60to120_pt50/"+prefix,idsuffix,leps,map_weight);
	  }
	}else if(dimass>=120&&dimass<400){
	  FillHist(channelname+"/"+tauprefix+"cutflow"+idsuffix,12.5,map_weight[""],20,0,20);
	  FillHistsEfficiency(channelname+"/m120to400/"+prefix,idsuffix,leps,map_weight);
	  if(dipt>50){
	    FillHist(channelname+"/"+tauprefix+"cutflow"+idsuffix,13.5,map_weight[""],20,0,20);
	    FillHistsEfficiency(channelname+"/m120to400_pt50/"+prefix,idsuffix,leps,map_weight);
	  }
	}	
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
	if(eta>etabin[ib]&&pt<etabin[ib+1]){
	  FillHist(Form("%sl%d%spt_eta%.2fto%.2f%s",pre.Data(),i,charge.Data(),etabin[ib],etabin[ib+1],suf.Data()),pt,w,500,0,500);
	  FillHist(Form("%sl%dpt_eta%.2fto%.2f%s",pre.Data(),i,etabin[ib],etabin[ib+1],suf.Data()),pt,w,500,0,500);
	  FillHist(Form("%sl%spt_eta%.2fto%.2f%s",pre.Data(),charge.Data(),etabin[ib],etabin[ib+1],suf.Data()),pt,w,500,0,500);
	  FillHist(Form("%slpt_eta%.2fto%.2f%s",pre.Data(),etabin[ib],etabin[ib+1],suf.Data()),pt,w,500,0,500);
	}
      } 
    }
    
    TLorentzVector dilepton=(*leps.at(0))+(*leps.at(1));
    double dimass=dilepton.M();
    double dipt=dilepton.Pt();
    double dirap=dilepton.Rapidity();
    FillHist(pre+"dimass"+suf,dimass,w,400,0,400);
    FillHist(pre+"dipt"+suf,dipt,w,400,0,400);
    FillHist(pre+"dirap"+suf,dirap,w,120,-3,3);
    
    double cost=GetCosThetaCS(leps);
    FillHist(pre+"costhetaCS"+suf,cost,w,40,-1,1);
    FillHist(pre+"abscosthetaCS"+suf,fabs(cost),w,20,0,1);
    double h=0.5*pow(dipt/dimass,2)/(1+pow(dipt/dimass,2))*(1-3*cost*cost);
    double den_weight=0.5*fabs(cost)/pow(1+cost*cost+h,2);
    double num_weight=0.5*cost*cost/pow(1+cost*cost+h,3);
    if(cost>0){
      FillHist(pre+"forward"+suf,dimass,w,massbinnum,(double*)massrange);
      FillHist(pre+"forward_den"+suf,dimass,w*den_weight,massbinnum,(double*)massrange);
      FillHist(pre+"forward_num"+suf,dimass,w*num_weight,massbinnum,(double*)massrange);    
    }else{
      FillHist(pre+"backward"+suf,dimass,w,massbinnum,(double*)massrange);
      FillHist(pre+"backward_den"+suf,dimass,w*den_weight,massbinnum,(double*)massrange);
      FillHist(pre+"backward_num"+suf,dimass,w*num_weight,massbinnum,(double*)massrange);   
    }
  }
}
