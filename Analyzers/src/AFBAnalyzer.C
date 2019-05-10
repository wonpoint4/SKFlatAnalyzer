#include "AFBAnalyzer.h"

void AFBAnalyzer::executeEvent(){
  tauprefix="";
  hardprefix="";
  ////////////////////////check genlevel//////////////////
  zptcor=1.;
  if(IsDYSample){
    vector<Gen> gens=GetGens();
    int parton0,parton1,hardl0,hardl1,l0,l1;
    vector<int> photons;
    GetGenIndex(gens,parton0,parton1,hardl0,hardl1,l0,l1,photons);
    int hardj0=0,hardj1=0,hardj2=0;
    Gen genhardj0,genhardj1,genhardj2;
    for(unsigned int i=parton1+1;i<gens.size();i++){
      if(gens[i].isHardProcess()){
	if(gens[i].PID()==21||abs(gens[i].PID())<7){
	  if(!hardj0){
	    hardj0=i;
	    genhardj0=gens[i];
	  }else if(!hardj1){
	    hardj1=i;
	    genhardj1=gens[i];
	  }else if(!hardj2){
	    hardj2=i;
	    genhardj2=gens[i];
	    break;
	  }
	}
      }
    }
    if(abs(gens[hardl0].PID())!=15){
      Gen genparton0=gens[parton0],genparton1=gens[parton1],genhardl0=gens[hardl0],genhardl1=gens[hardl1],genl0=gens[l0],genl1=gens[l1],genphotons;
      genparton0.SetPxPyPzE(0,0,genWeight_X1*6500.,genWeight_X1*6500.);
      genparton1.SetPxPyPzE(0,0,-genWeight_X2*6500.,genWeight_X2*6500.);
      for(unsigned int i=0;i<photons.size();i++) genphotons+=gens[photons[i]];
      TLorentzVector genZ=(genl0+genl1+genphotons);
      zptcor*=GetZPtWeight(genZ.Pt(),genZ.Rapidity(),abs(genhardl0.PID())==13?Lepton::Flavour::MUON:Lepton::Flavour::ELECTRON);
      TString slepton=abs(genhardl0.PID())==13?"muon":"electron";
      if(genparton0.PID()==21){
	if(genparton1.PID()==21) hardprefix="gg_";
	else if(genparton1.PID()>0) hardprefix="gq_";
	else if(genparton1.PID()<0) hardprefix="gqbar_";
      }else if(genparton0.PID()>0){
	if(genparton1.PID()==21) hardprefix="qg_";
	else if(genparton1.PID()>0) hardprefix="qq_";
	else if(genparton1.PID()<0) hardprefix="qqbar_";
      }else if(genparton0.PID()<0){
	if(genparton1.PID()==21) hardprefix="qbarg_";
	else if(genparton1.PID()>0) hardprefix="qbarq_";
	else if(genparton1.PID()<0) hardprefix="qbarqbar_";
      }
      int nhardjet=(hardj0?1:0)+(hardj1?1:0)+(hardj2?1:0);

      ///////////////////////FillGenHist////////////////
      //FillGenHists(Form("%s%dgen/",slepton.Data(),DataYear)+hardprefix,"",genl0,genl1,genphotons,weight_norm_1invpb*gen_weight*zptcor);
      //FillGenHists(Form("%s%dgen/",slepton.Data(),DataYear)+hardprefix,"_nozptcor",genl0,genl1,genphotons,weight_norm_1invpb*gen_weight);
      /*
      FillHist(Form("%s%dgen/",slepton.Data(),DataYear)+hardprefix+"dipty",genZ.Pt(),fabs(genZ.Rapidity()),weight_norm_1invpb*gen_weight*zptcor,zptcor_nptbin,(double*)zptcor_ptbin,zptcor_nybin,(double*)zptcor_ybin);
      FillHist(Form("%s%dgen/",slepton.Data(),DataYear)+hardprefix+"dipty_nozptcor",genZ.Pt(),fabs(genZ.Rapidity()),weight_norm_1invpb*gen_weight,zptcor_nptbin,(double*)zptcor_ptbin,zptcor_nybin,(double*)zptcor_ybin);

      FillHardHists(Form("%s%dgen/",slepton.Data(),DataYear)+hardprefix,"_lab",genparton0,genparton1,genhardl0,genhardl1,genhardj0,weight_norm_1invpb*gen_weight);
      if(nhardjet==0){
	FillHardHists(Form("%s%dgen/",slepton.Data(),DataYear)+hardprefix+"0j_","_lab",genparton0,genparton1,genhardl0,genhardl1,genhardj0,weight_norm_1invpb*gen_weight);
      }else if(nhardjet==1){
	FillHardHists(Form("%s%dgen/",slepton.Data(),DataYear)+hardprefix+"1j_","_lab",genparton0,genparton1,genhardl0,genhardl1,genhardj0,weight_norm_1invpb*gen_weight);
	if(fabs(genhardj0.Eta())<3)
	  FillHardHists(Form("%s%dgen/jeta03/",slepton.Data(),DataYear)+hardprefix+"1j_","_lab",genparton0,genparton1,genhardl0,genhardl1,genhardj0,weight_norm_1invpb*gen_weight);
      }

      TVector3 b=(genhardl0+genhardl1).BoostVector();
      genparton0.Boost(-b);genparton1.Boost(-b);
      genhardl0.Boost(-b);genhardl1.Boost(-b);
      genhardj0.Boost(-b);
      //cout<<"dycm: dy:"<<(genhardl0+genhardl1).P()<<" pp:"<<(genparton0+genparton1).P()<<endl;

      FillHardHists(Form("%s%dgen/",slepton.Data(),DataYear)+hardprefix,"_dycm",genparton0,genparton1,genhardl0,genhardl1,genhardj0,weight_norm_1invpb*gen_weight);
      if(nhardjet==0){
	FillHardHists(Form("%s%dgen/",slepton.Data(),DataYear)+hardprefix+"0j_","_dycm",genparton0,genparton1,genhardl0,genhardl1,genhardj0,weight_norm_1invpb*gen_weight);
      }else if(nhardjet==1){
	FillHardHists(Form("%s%dgen/",slepton.Data(),DataYear)+hardprefix+"1j_","_dycm",genparton0,genparton1,genhardl0,genhardl1,genhardj0,weight_norm_1invpb*gen_weight);
	if(fabs(genhardj0.Eta())<3)
	  FillHardHists(Form("%s%dgen/jeta03/",slepton.Data(),DataYear)+hardprefix+"1j_","_dycm",genparton0,genparton1,genhardl0,genhardl1,genhardj0,weight_norm_1invpb*gen_weight);
      }


      b=(genparton0+genparton1).BoostVector();
      genparton0.Boost(-b);genparton1.Boost(-b);
      genhardl0.Boost(-b);genhardl1.Boost(-b);
      genhardj0.Boost(-b);
      //cout<<"ppcm: dy:"<<(genhardl0+genhardl1).P()<<" pp:"<<(genparton0+genparton1).P()<<endl;

      FillHardHists(Form("%s%dgen/",slepton.Data(),DataYear)+hardprefix,"_ppcm",genparton0,genparton1,genhardl0,genhardl1,genhardj0,weight_norm_1invpb*gen_weight);
      if(nhardjet==0){
	FillHardHists(Form("%s%dgen/",slepton.Data(),DataYear)+hardprefix+"0j_","_ppcm",genparton0,genparton1,genhardl0,genhardl1,genhardj0,weight_norm_1invpb*gen_weight);
      }else if(nhardjet==1){
	FillHardHists(Form("%s%dgen/",slepton.Data(),DataYear)+hardprefix+"1j_","_ppcm",genparton0,genparton1,genhardl0,genhardl1,genhardj0,weight_norm_1invpb*gen_weight);
	if(fabs(genhardj0.Eta())<3)
	  FillHardHists(Form("%s%dgen/jeta03/",slepton.Data(),DataYear)+hardprefix+"1j_","_ppcm",genparton0,genparton1,genhardl0,genhardl1,genhardj0,weight_norm_1invpb*gen_weight);
      }

      if(hardprefix=="qqbar_"||hardprefix=="qbarq_"){
	if(nhardjet==0) hardprefix+="0j_";
	else hardprefix+="nj_";
      }
      */
    }else{
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
  }
  delete ev;
}

void AFBAnalyzer::executeEventFromParameter(TString channelname,Event* ev){
  std::vector<Muon> muons,muons_scale_up,muons_scale_down;
  std::vector<Electron> electrons,electrons_scale_up,electrons_scale_down,electrons_smear_up,electrons_smear_down;

  std::vector<Jet> jets=GetJets("tightLepVeto",30,2.7);
  std::sort(jets.begin(),jets.end(),PtComparing);

  std::vector<Lepton*> leps;
  double lep0ptcut,lep1ptcut;
  double (MCCorrection::*LeptonID_SF)(TString,double,double,int)=NULL;
  double (MCCorrection::*LeptonISO_SF)(TString,double,double,int)=NULL;
  double (MCCorrection::*LeptonReco_SF)(double,double,int)=NULL;
  TString LeptonID_key,LeptonISO_key,LeptonReco_key,triggerSF_key0,triggerSF_key1;
  map<TString,vector<Lepton*>> map_lepton_systematic;
  if(channelname.Contains("muon")){
    muons=GetMuons("POGTightWithTightIso",0.0,2.4);
    std::sort(muons.begin(),muons.end(),PtComparing);
    leps=MakeLeptonPointerVector(muons);
    
    muons_scale_up=ScaleMuons(muons,1);
    std::sort(muons_scale_up.begin(),muons_scale_up.end(),PtComparing);
    map_lepton_systematic["scale_up"]=MakeLeptonPointerVector(muons_scale_up);

    muons_scale_down=ScaleMuons(muons,-1);
    std::sort(muons_scale_down.begin(),muons_scale_down.end(),PtComparing);
    map_lepton_systematic["scale_down"]=MakeLeptonPointerVector(muons_scale_down);
    
    lep0ptcut=20.;
    lep1ptcut=10.;
    LeptonID_SF=&MCCorrection::MuonID_SF;
    LeptonISO_SF=&MCCorrection::MuonISO_SF;
    LeptonID_key="NUM_TightID_DEN_genTracks";
    LeptonISO_key="NUM_TightRelIso_DEN_TightIDandIPCut";
    triggerSF_key0="IsoMu17_POGTight";
    triggerSF_key1="Mu8_POGTight";
 }else if(channelname.Contains("electron")){
    electrons=SMPGetElectrons("passMediumID_selective",0.0,2.5);
    std::sort(electrons.begin(),electrons.end(),PtComparing);
    leps=MakeLeptonPointerVector(electrons);

    electrons_scale_up=ScaleElectrons(electrons,1);
    std::sort(electrons_scale_up.begin(),electrons_scale_up.end(),PtComparing);
    map_lepton_systematic["scale_up"]=MakeLeptonPointerVector(electrons_scale_up);

    electrons_scale_down=ScaleElectrons(electrons,-1);
    std::sort(electrons_scale_down.begin(),electrons_scale_down.end(),PtComparing);
    map_lepton_systematic["scale_down"]=MakeLeptonPointerVector(electrons_scale_down);

    electrons_smear_up=SmearElectrons(electrons,1);
    std::sort(electrons_smear_up.begin(),electrons_smear_up.end(),PtComparing);
    map_lepton_systematic["smear_up"]=MakeLeptonPointerVector(electrons_smear_up);

    electrons_smear_down=SmearElectrons(electrons,-1);
    std::sort(electrons_smear_down.begin(),electrons_smear_down.end(),PtComparing);
    map_lepton_systematic["smear_down"]=MakeLeptonPointerVector(electrons_smear_down);

    lep0ptcut=25.;
    lep1ptcut=15.;
    LeptonID_SF=&MCCorrection::ElectronID_SF;
    LeptonReco_SF=&MCCorrection::ElectronReco_SF;
    LeptonID_key="Strict_MediumID_pt10";
    triggerSF_key0="Strict_LeadEle23_MediumID";
    triggerSF_key1="Strict_TailEle12_MediumID";
  }else{
    cout<<"[AFBAnalyzer::executeEventFromParameter] wrong channelname"<<endl;
    return;
  }
  
  /////////////////lumi weight///////////////////
  double weight=1.,totalweight=1.;
  if(!IsDATA){
    weight=weight_norm_1invpb*ev->MCweight()*ev->GetTriggerLumi("Full");
    totalweight*=weight;
  }
  FillHist(channelname+"/"+tauprefix+"cutflow",0.5,totalweight,20,0,20);

  /////////////////PUreweight///////////////////
  double PUreweight=1.,PUreweight_up=1.,PUreweight_down=1.;
  if(!IsDATA&&!MCSample.Contains("ZToEE")&&!MCSample.Contains("ZToMuMu")){
    PUreweight=GetPileUpWeight(nPileUp,0);
    PUreweight_up=GetPileUpWeight(nPileUp,1);
    PUreweight_down=GetPileUpWeight(nPileUp,-1);
    totalweight*=PUreweight;
  }
  FillHist(channelname+"/"+tauprefix+"cutflow",1.5,totalweight,20,0,20);
  
  //////////////////////PrefileWeight////////////////////
  double prefireweight=1.;
  double prefireweight_up=1.;
  double prefireweight_down=1.;
  if(!IsDATA){
    prefireweight=L1PrefireReweight_Central;
    prefireweight_up=L1PrefireReweight_Up;
    prefireweight_down=L1PrefireReweight_Down;
  }
  totalweight*=prefireweight;
  FillHist(channelname+"/"+tauprefix+"cutflow",2.5,totalweight,20,0,20);

  totalweight*=zptcor;
  FillHist(channelname+"/"+tauprefix+"cutflow",3.5,totalweight,20,0,20);


  ///////////////////////lepton selection///////////////////////
  FillHist(channelname+"/"+tauprefix+"nlepton",leps.size(),totalweight,10,0,10);
  if(leps.size()==2){
    TString prefix=(leps.at(0)->Charge()*leps.at(1)->Charge()>0?"ss_":"")+tauprefix;
    FillHist(channelname+"/"+tauprefix+"cutflow",4.5,totalweight,20,0,20);
    if(leps.at(0)->Pt()>lep0ptcut&&leps.at(1)->Pt()>lep1ptcut){
      FillHist(channelname+"/"+tauprefix+"cutflow",5.5,totalweight,20,0,20);

      /////////////////efficiency scale factors///////////////////
      double IDSF=1.,IDSF_up=1.,IDSF_down=1.;
      double ISOSF=1.,ISOSF_up=1.,ISOSF_down=1.;
      double RECOSF=1.,RECOSF_up=1.,RECOSF_down=1.;
      if(!IsDATA){
	for(unsigned int i=0;i<leps.size();i++){
	  double this_pt,this_eta;
	  TString chargesuffix;
	  if(leps[i]->LeptonFlavour()==Lepton::MUON){
	    this_pt=((Muon*)leps.at(i))->MiniAODPt();
	    this_eta=leps.at(i)->Eta();
	  }else if(leps[i]->LeptonFlavour()==Lepton::ELECTRON){
	    this_pt=leps.at(i)->Pt();
	    this_eta=((Electron*)leps.at(i))->scEta();
	    if(leps[i]->Charge()>0) chargesuffix="_QPlus";
	    else chargesuffix="_QMinus";
	  }else{
	    cout <<"[AFBAnalyzer::executeEventFromParameter] It is not lepton"<<endl;
	    exit(EXIT_FAILURE);
	  }
	  
	  double this_RECOSF=LeptonReco_SF?(mcCorr->*LeptonReco_SF)(this_eta,this_pt,0):1.;
	  double this_RECOSF_up=LeptonReco_SF?(mcCorr->*LeptonReco_SF)(this_eta,this_pt,1):1.;
	  double this_RECOSF_down=LeptonReco_SF?(mcCorr->*LeptonReco_SF)(this_eta,this_pt,-1):1.;
	  RECOSF*=this_RECOSF; RECOSF_up*=this_RECOSF_up; RECOSF_down*=this_RECOSF_down;
	  
	  double this_IDSF=LeptonID_SF?(mcCorr->*LeptonID_SF)(LeptonID_key+chargesuffix,this_eta,this_pt,0):1.;
	  double this_IDSF_up=LeptonID_SF?(mcCorr->*LeptonID_SF)(LeptonID_key+chargesuffix,this_eta,this_pt,1):1.;
	  double this_IDSF_down=LeptonID_SF?(mcCorr->*LeptonID_SF)(LeptonID_key+chargesuffix,this_eta,this_pt,-1):1.;
	  IDSF*=this_IDSF; IDSF_up*=this_IDSF_up; IDSF_down*=this_IDSF_down;
	  
	  double this_ISOSF=LeptonISO_SF?(mcCorr->*LeptonISO_SF)(LeptonISO_key,this_eta,this_pt,0):1.;
	  double this_ISOSF_up=LeptonISO_SF?(mcCorr->*LeptonISO_SF)(LeptonISO_key,this_eta,this_pt,1):1.;
	  double this_ISOSF_down=LeptonISO_SF?(mcCorr->*LeptonISO_SF)(LeptonISO_key,this_eta,this_pt,-1):1.;
	  ISOSF*=this_ISOSF; ISOSF_up*=this_ISOSF_up; ISOSF_down*=this_ISOSF_down;
	}
      }
      totalweight*=RECOSF;
      FillHist(channelname+"/"+tauprefix+"cutflow",6.5,totalweight,20,0,20);
      totalweight*=IDSF;
      FillHist(channelname+"/"+tauprefix+"cutflow",7.5,totalweight,20,0,20);
      totalweight*=ISOSF;
      FillHist(channelname+"/"+tauprefix+"cutflow",8.5,totalweight,20,0,20);
      
      double triggerSF=1.,triggerSF_up=1.,triggerSF_down=1.;
      if(!IsDATA){
	triggerSF*=DileptonTrigger_SF(triggerSF_key0,triggerSF_key1,leps,0);
	triggerSF_up*=DileptonTrigger_SF(triggerSF_key0,triggerSF_key1,leps,1);
	triggerSF_down*=DileptonTrigger_SF(triggerSF_key0,triggerSF_key1,leps,-1);
      }
      totalweight*=triggerSF;
      FillHist(channelname+"/"+tauprefix+"cutflow",9.5,totalweight,20,0,20);

      ///////////////////////weight systematics//////////////////
      map<TString,double> map_weight_systematic;
      if(!IsDATA){
	map_weight_systematic["noPUreweight"]=weight*RECOSF*IDSF*ISOSF*triggerSF*prefireweight*zptcor;
	map_weight_systematic["PUreweight_up"]=weight*PUreweight_up*RECOSF*IDSF*ISOSF*triggerSF*prefireweight*zptcor;
	map_weight_systematic["PUreweight_down"]=weight*PUreweight_down*RECOSF*IDSF*ISOSF*triggerSF*prefireweight*zptcor;
	
	map_weight_systematic["noefficiencySF"]=weight*PUreweight*prefireweight*zptcor;
	
	map_weight_systematic["noRECOSF"]=weight*PUreweight*IDSF*ISOSF*triggerSF*prefireweight*zptcor;
	map_weight_systematic["RECOSF_up"]=weight*PUreweight*RECOSF_up*IDSF*ISOSF*triggerSF*prefireweight*zptcor;
	map_weight_systematic["RECOSF_down"]=weight*PUreweight*RECOSF_down*IDSF*ISOSF*triggerSF*prefireweight*zptcor;
	
	map_weight_systematic["noIDSF"]=weight*PUreweight*RECOSF*ISOSF*triggerSF*prefireweight*zptcor;
	map_weight_systematic["IDSF_up"]=weight*PUreweight*RECOSF*IDSF_up*ISOSF*triggerSF*prefireweight*zptcor;
	map_weight_systematic["IDSF_down"]=weight*PUreweight*RECOSF*IDSF_down*ISOSF*triggerSF*prefireweight*zptcor;
	
	map_weight_systematic["noISOSF"]=weight*PUreweight*RECOSF*IDSF*triggerSF*prefireweight*zptcor;
	map_weight_systematic["ISOSF_up"]=weight*PUreweight*RECOSF*IDSF*ISOSF_up*triggerSF*prefireweight*zptcor;
	map_weight_systematic["ISOSF_down"]=weight*PUreweight*RECOSF*IDSF*ISOSF_down*triggerSF*prefireweight*zptcor;
	
	map_weight_systematic["notriggerSF"]=weight*PUreweight*RECOSF*IDSF*ISOSF*prefireweight*zptcor;
	map_weight_systematic["triggerSF_up"]=weight*PUreweight*RECOSF*IDSF*ISOSF*triggerSF_up*prefireweight*zptcor;
	map_weight_systematic["triggerSF_down"]=weight*PUreweight*RECOSF*IDSF*ISOSF*triggerSF_down*prefireweight*zptcor;
	
	map_weight_systematic["noprefireweight"]=weight*PUreweight*RECOSF*IDSF*ISOSF*triggerSF*zptcor;
	map_weight_systematic["prefireweight_up"]=weight*PUreweight*RECOSF*IDSF*ISOSF*triggerSF*prefireweight_up*zptcor;
	map_weight_systematic["prefireweight_down"]=weight*PUreweight*RECOSF*IDSF*ISOSF*triggerSF*prefireweight_down*zptcor;
	
	map_weight_systematic["nozptcor"]=weight*PUreweight*RECOSF*IDSF*ISOSF*triggerSF*prefireweight;
	
	if(PDFWeights_AlphaS->size()==2){
	  map_weight_systematic["alphaS_up"]=weight*PUreweight*RECOSF*IDSF*ISOSF*triggerSF*prefireweight*zptcor*PDFWeights_AlphaS->at(0);
	  map_weight_systematic["alphaS_down"]=weight*PUreweight*RECOSF*IDSF*ISOSF*triggerSF*prefireweight*zptcor*PDFWeights_AlphaS->at(1);
	}
	
	for(unsigned int i=0;i<PDFWeights_Scale->size();i++){
	  map_weight_systematic[Form("scalevariation%d",i)]=weight*PUreweight*RECOSF*IDSF*ISOSF*triggerSF*prefireweight*zptcor*PDFWeights_Scale->at(i);
	}
	
	for(unsigned int i=0;i<PDFWeights_Error->size();i++){
	  //map_weight_systematic[Form("pdf%d",i)]=weight*PUreweight*RECOSF*IDSF*ISOSF*triggerSF*prefireweight*zptcor*PDFWeights_Error->at(i);
	}	
      }

      ///////////////////////fill hists///////////////////////
      TLorentzVector dilepton=(*leps.at(0))+(*leps.at(1));
      double dimass=dilepton.M();
      double dirap=dilepton.Rapidity();
      double dipt=dilepton.Pt();
      if(dimass>=massrange[0]&&dimass<massrange[massbinnum]){
	FillHist(channelname+"/"+tauprefix+"cutflow",10.5,totalweight,20,0,20);
	FillAFBHists(channelname+Form("/m%.0fto%.0f/",massrange[0],massrange[massbinnum])+prefix,"",leps,jets,totalweight);
	if(IsDYSample&&prefix=="") FillAFBHists(channelname+Form("/m%.0fto%.0f/",massrange[0],massrange[massbinnum])+hardprefix,"",leps,jets,totalweight);
	FillAFBSystematicHists(channelname+Form("/m%.0fto%.0f/",massrange[0],massrange[massbinnum])+prefix,"",leps,jets,map_weight_systematic);
	if(dipt>50){
	  FillAFBHists(channelname+Form("/m%.0fto%.0f_pt50/",massrange[0],massrange[massbinnum])+prefix,"",leps,jets,totalweight);
	  if(IsDYSample&&prefix=="") FillAFBHists(channelname+Form("/m%.0fto%.0f_pt50/",massrange[0],massrange[massbinnum])+hardprefix,"",leps,jets,totalweight);
	  FillAFBSystematicHists(channelname+Form("/m%.0fto%.0f_pt50/",massrange[0],massrange[massbinnum])+prefix,"",leps,jets,map_weight_systematic);	  
	}	
	for(int iy=0;iy<zptcor_nybin;iy++){
	  if(fabs(dirap)>=zptcor_ybin[iy]&&fabs(dirap)<zptcor_ybin[iy+1]){
	    FillAFBHists(channelname+Form("/m%.0fto%.0f_y%.1fto%.1f/",massrange[0],massrange[massbinnum],zptcor_ybin[iy],zptcor_ybin[iy+1])+prefix,"",leps,jets,totalweight);
	    if(IsDYSample&&prefix=="") FillAFBHists(channelname+Form("/m%.0fto%.0f_y%.1fto%.1f/",massrange[0],massrange[massbinnum],zptcor_ybin[iy],zptcor_ybin[iy+1])+hardprefix,"",leps,jets,totalweight);
	    FillAFBSystematicHists(channelname+Form("/m%.0fto%.0f_y%.1fto%.1f/",massrange[0],massrange[massbinnum],zptcor_ybin[iy],zptcor_ybin[iy+1])+prefix,"",leps,jets,map_weight_systematic);
	    if(dipt>50){
	      FillAFBHists(channelname+Form("/m%.0fto%.0f_y%.1fto%.1f_pt50/",massrange[0],massrange[massbinnum],zptcor_ybin[iy],zptcor_ybin[iy+1])+prefix,"",leps,jets,totalweight);
	      if(IsDYSample&&prefix=="") FillAFBHists(channelname+Form("/m%.0fto%.0f_y%.1fto%.1f_pt50/",massrange[0],massrange[massbinnum],zptcor_ybin[iy],zptcor_ybin[iy+1])+hardprefix,"",leps,jets,totalweight);
	      FillAFBSystematicHists(channelname+Form("/m%.0fto%.0f_y%.1fto%.1f_pt50/",massrange[0],massrange[massbinnum],zptcor_ybin[iy],zptcor_ybin[iy+1])+prefix,"",leps,jets,map_weight_systematic);
	    }
	    for(int im=0;im<massbinnum;im++){
	      if(dimass>=massrange[im]&&dimass<massrange[im+1]){
		//FillAFBHists(channelname+Form("/m%.0fto%.0f_y%.1fto%.1f/",massrange[im],massrange[im+1],zptcor_ybin[iy],zptcor_ybin[iy+1])+prefix,"",leps,jets,totalweight);
		//FillAFBSystematicHists(channelname+Form("/m%.0fto%.0f_y%.1fto%.1f/",massrange[im],massrange[im+1],zptcor_ybin[iy],zptcor_ybin[iy+1])+prefix,"",leps,jets,map_weight_systematic);
	      }
	    }
	  }
	}
      }
    }
  }


  ///////////////////////leps systematic fill hists///////////////////////
  for(auto iter=map_lepton_systematic.begin();iter!=map_lepton_systematic.end();iter++){
    totalweight=weight*PUreweight*prefireweight*zptcor;
    FillHist(channelname+"/"+tauprefix+"nlepton"+"_"+iter->first,iter->second.size(),totalweight,10,0,10);
    if(iter->second.size()==2){
      TString prefix=(iter->second.at(0)->Charge()*iter->second.at(1)->Charge()>0?"ss_":"")+tauprefix;
      if(iter->second.at(0)->Pt()>lep0ptcut&&iter->second.at(1)->Pt()>lep1ptcut){

	/////////////////efficiency scale factors///////////////////
	double IDSF=1.;double ISOSF=1.;double RECOSF=1.;
	if(!IsDATA){
	  for(unsigned int i=0;i<iter->second.size();i++){
	    double this_pt,this_eta;
	    TString chargesuffix;
	    if(iter->second.at(i)->LeptonFlavour()==Lepton::MUON){
	      this_pt=((Muon*)iter->second.at(i))->MiniAODPt();
	      this_eta=iter->second.at(i)->Eta();
	    }else if(iter->second.at(i)->LeptonFlavour()==Lepton::ELECTRON){
	      this_pt=iter->second.at(i)->Pt();
	      this_eta=((Electron*)iter->second.at(i))->scEta();
	      if(iter->second.at(i)->Charge()>0) chargesuffix="_QPlus";
	      else chargesuffix="_QMinus";
	    }else{
	      cout <<"[AFBAnalyzer::executeEventFromParameter] It is not lepton"<<endl;
	      exit(EXIT_FAILURE);
	    }
	    
	    double this_RECOSF=LeptonReco_SF?(mcCorr->*LeptonReco_SF)(this_eta,this_pt,0):1.;
	    RECOSF*=this_RECOSF;
	    double this_IDSF=LeptonID_SF?(mcCorr->*LeptonID_SF)(LeptonID_key+chargesuffix,this_eta,this_pt,0):1.;
	    IDSF*=this_IDSF;
	    double this_ISOSF=LeptonISO_SF?(mcCorr->*LeptonISO_SF)(LeptonISO_key,this_eta,this_pt,0):1.;
	    ISOSF*=this_ISOSF;
	  }
	}
	totalweight*=RECOSF;
	totalweight*=IDSF;
	totalweight*=ISOSF;
	
	double triggerSF=1.;
	if(!IsDATA){
	  triggerSF*=DileptonTrigger_SF(triggerSF_key0,triggerSF_key1,iter->second,0);
	}
	totalweight*=triggerSF;

	TLorentzVector dilepton=(*iter->second.at(0))+(*iter->second.at(1));
	double dimass=dilepton.M();
	double dirap=dilepton.Rapidity();
	double dipt=dilepton.Pt();
	if(dimass>=massrange[0]&&dimass<massrange[massbinnum]){
	  FillAFBHists(channelname+Form("/m%.0fto%.0f/",massrange[0],massrange[massbinnum])+prefix,"_"+iter->first,iter->second,jets,totalweight);
	  if(dipt>50){
	    FillAFBHists(channelname+Form("/m%.0fto%.0f_pt50/",massrange[0],massrange[massbinnum])+prefix,"_"+iter->first,iter->second,jets,totalweight);
	  }	
	  for(int iy=0;iy<zptcor_nybin;iy++){
	    if(fabs(dirap)>=zptcor_ybin[iy]&&fabs(dirap)<zptcor_ybin[iy+1]){
	      FillAFBHists(channelname+Form("/m%.0fto%.0f_y%.1fto%.1f/",massrange[0],massrange[massbinnum],zptcor_ybin[iy],zptcor_ybin[iy+1])+prefix,"_"+iter->first,iter->second,jets,totalweight);
	      if(dipt>50){
		FillAFBHists(channelname+Form("/m%.0fto%.0f_y%.1fto%.1f_pt50/",massrange[0],massrange[massbinnum],zptcor_ybin[iy],zptcor_ybin[iy+1])+prefix,"_"+iter->first,iter->second,jets,totalweight);
	      }
	      for(int im=0;im<massbinnum;im++){
		if(dimass>=massrange[im]&&dimass<massrange[im+1]){
		  //FillAFBHists(channelname+Form("/m%.0fto%.0f_y%.1fto%.1f/",massrange[im],massrange[im+1],zptcor_ybin[iy],zptcor_ybin[iy+1])+prefix,"_"+iter->first,iter->second,jets,totalweight);
		}
	      }
	    }
	  }
	}
      }
    }
  }
}
AFBAnalyzer::AFBAnalyzer(){
  random=new TRandom3(4497);
}
AFBAnalyzer::~AFBAnalyzer(){
  if(random) delete random;
}
double AFBAnalyzer::GetCosThetaCS(const vector<Lepton*>& leps){
  TLorentzVector *l0,*l1;
  if(leps.at(0)->Charge()<0&&leps.at(1)->Charge()>0){
    l0=leps.at(0);
    l1=leps.at(1);
  }else if(leps.at(0)->Charge()>0&&leps.at(1)->Charge()<0){
    l0=leps.at(1);
    l1=leps.at(0);
  }else{
    if(random->Rndm()<0.5){
      l0=leps.at(0);
      l1=leps.at(1);
    }else{
      l0=leps.at(1);
      l1=leps.at(0);
    }      
  }

  TLorentzVector dilepton=*l0+*l1;
  double l0pp=(l0->E()+l0->Pz())/sqrt(2);
  double l0pm=(l0->E()-l0->Pz())/sqrt(2);
  double l1pp=(l1->E()+l1->Pz())/sqrt(2);
  double l1pm=(l1->E()-l1->Pz())/sqrt(2);
  double dimass=dilepton.M();
  double dipt=dilepton.Pt();
  double direction=dilepton.Pz()/fabs(dilepton.Pz());
  return direction*2*(l0pp*l1pm-l0pm*l1pp)/sqrt(dimass*dimass*(dimass*dimass+dipt*dipt));
} 
double AFBAnalyzer::GetCosTheta(const vector<Lepton*>& leps,const vector<Jet>& jets,TString option,double fcut){
  TLorentzVector *l0,*l1;
  if(leps.at(0)->Charge()<0&&leps.at(1)->Charge()>0){
    l0=leps.at(0);
    l1=leps.at(1);
  }else if(leps.at(0)->Charge()>0&&leps.at(1)->Charge()<0){
    l0=leps.at(1);
    l1=leps.at(0);
  }else{
    if(random->Rndm()<0.5){
      l0=leps.at(0);
      l1=leps.at(1);
    }else{
      l0=leps.at(1);
      l1=leps.at(0);
    }      
  }
  TLorentzVector dilepton=*l0+*l1;
  TLorentzVector j0;
  if(jets.size()>0){
    if((jets[0]+dilepton).Pt()/dilepton.Pt()<fcut)
      j0=jets[0];
  }
  TLorentzVector system=dilepton+j0;
  double direction=dilepton.Pz()/fabs(dilepton.Pz());
  double systemmass=system.M();
  double systempt=system.Pt();
  double systempz=system.Pz();
  TLorentzVector p0(0,0,(sqrt(systemmass*systemmass+systempz*systempz)+systempz)/2.,(sqrt(systemmass*systemmass+systempz*systempz)+systempz)/2.);
  TLorentzVector p1(0,0,-(sqrt(systemmass*systemmass+systempz*systempz)-systempz)/2.,(sqrt(systemmass*systemmass+systempz*systempz)-systempz)/2.);
  TLorentzVector particle=*l0;
  TVector3 b1=system.BoostVector();
  dilepton.Boost(-b1);particle.Boost(-b1);
  p0.Boost(-b1);p1.Boost(-b1);j0.Boost(-b1);
  if(j0.E()>0){
    if(p0.Angle(j0.Vect())<p1.Angle(j0.Vect())){
      if(option.Contains("C")) p0-=j0;
      if(option.Contains("B")) direction=-1;
    }else{
      if(option.Contains("C")) p1-=j0;
      if(option.Contains("B")) direction=1;
    }
  }
  TVector3 b2=dilepton.BoostVector();
  p0.Boost(-b2);p1.Boost(-b2);particle.Boost(-b2);
  return direction*cos(particle.Angle(p0.Vect().Unit()-p1.Vect().Unit()));
} 
void AFBAnalyzer::FillAFBHists(TString pre,TString suf,const vector<Lepton*>& leps,const vector<Jet>& jets,double w){
  if(leps.at(0)->Pt()>60){
    if(leps.at(0)->Charge()>0) FillHist(pre+"l0eta_pt60_plus"+suf,leps.at(0)->Eta(),w,60,-3,3);
    else FillHist(pre+"l0eta_pt60_minus"+suf,leps.at(0)->Eta(),w,60,-3,3);
  }
  //for leptons
  FillBasicHists(pre,suf,leps,w);
  TLorentzVector dilepton=(*leps.at(0))+(*leps.at(1));
  double dimass=dilepton.M();
  double dipt=dilepton.Pt();

  double cost=GetCosThetaCS(leps);
  FillHist(pre+"costhetaCS"+suf,cost,w,40,-1,1);
  FillHist(pre+"abscosthetaCS"+suf,fabs(cost),w,20,0,1);
  if(pre.Contains(TRegexp("m60to..."))){
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
  /*
  TString options[]={"A","B","C","BC"};
  double fcuts[]={0.5,1,100};
  for(int i=0;i<sizeof(options)/sizeof(TString);i++){
    for(int j=0;j<sizeof(fcuts)/sizeof(double);j++){
      double cost_=GetCosTheta(leps,jets,options[i],fcuts[j]);
      FillHist(pre+"costheta"+options[i]+"_"+Form("%d",int(fcuts[j]*10))+suf,cost_,w,40,-1,1);
      FillHist(pre+"abscostheta"+options[i]+"_"+Form("%d",int(fcuts[j]*10))+suf,fabs(cost_),w,20,0,1);
      if(pre.Contains("m60to120")){
	double h=0.5*pow(dipt/dimass,2)/(1+pow(dipt/dimass,2))*(1-3*cost_*cost_);
	if(cost_>0){
	  FillHist(pre+"forward_"+options[i]+"_"+Form("%d",int(fcuts[j]*10))+suf,dimass,w,massbinnum,(double*)massrange);
	  FillHist(pre+"forward_den_"+options[i]+"_"+Form("%d",int(fcuts[j]*10))+suf,dimass,w*0.5*fabs(cost_)/pow(1+cost_*cost_+h,2),massbinnum,(double*)massrange);
	  FillHist(pre+"forward_num_"+options[i]+"_"+Form("%d",int(fcuts[j]*10))+suf,dimass,w*0.5*cost_*cost_/pow(1+cost_*cost_+h,3),massbinnum,(double*)massrange);    
	}else{
	  FillHist(pre+"backward_"+options[i]+"_"+Form("%d",int(fcuts[j]*10))+suf,dimass,w,massbinnum,(double*)massrange);
	  FillHist(pre+"backward_den_"+options[i]+"_"+Form("%d",int(fcuts[j]*10))+suf,dimass,w*0.5*fabs(cost_)/pow(1+cost_*cost_+h,2),massbinnum,(double*)massrange);
	  FillHist(pre+"backward_num_"+options[i]+"_"+Form("%d",int(fcuts[j]*10))+suf,dimass,w*0.5*cost_*cost_/pow(1+cost_*cost_+h,3),massbinnum,(double*)massrange);   
	}
      }   
    }
  }

  //for jets
  FillHist(pre+"njet"+suf,jets.size(),w,10,0,10);
  for(unsigned int i=0;i<jets.size();i++){
    if(i>2) break;
    FillHist(pre+Form("j%dpt",i)+suf,jets[i].Pt(),w,200,0,200);
    FillHist(pre+Form("j%deta",i)+suf,jets[i].Eta(),w,100,-5,5);
    if(i==0){
      FillHist(pre+"Zj0delphi"+suf,jets[0].DeltaPhi(dilepton),w,80,-4,4);
      FillHist(pre+"j0ptoverZpt"+suf,jets[0].Pt()/dipt,w,200,0,4);
      FillHist(pre+"Zj0pt"+suf,(jets[0]+dilepton).Pt(),w,200,0,200);
      FillHist(pre+"Zj0ptoverZpt"+suf,(jets[0]+dilepton).Pt()/dilepton.Pt(),w,100,0,5);
      FillHist(pre+"Zj0deleta"+suf,dilepton.Eta()-jets[0].Eta(),w,200,-5,5);
      FillHist(pre+"j0mass"+suf,jets[0].M(),w,100,0,100);
    }
  }

  if(jets.size()>0){
    if(fabs(dilepton.Rapidity()<2.0)) return;
    if(dilepton.Pt()<50) return;
    cout<<hardprefix<<endl;
    TLorentzVector j0=jets[0];
    TLorentzVector system=dilepton+j0;
    double direction=dilepton.Pz()/fabs(dilepton.Pz());
    double systemmass=system.M();
    double systempt=system.Pt();
    double systempz=system.Pz();
    TLorentzVector p0(0,0,(sqrt(systemmass*systemmass+systempz*systempz)+systempz)/2.,(sqrt(systemmass*systemmass+systempz*systempz)+systempz)/2.);
    TLorentzVector p1(0,0,-(sqrt(systemmass*systemmass+systempz*systempz)-systempz)/2.,(sqrt(systemmass*systemmass+systempz*systempz)-systempz)/2.);
    //cout<<"###############before boost#########"<<endl;
    //cout<<dilepton.Pt()<<" "<<system.Pt()<<endl;
    //p0.Print();p1.Print();system.Print();dilepton.Print();j0.Print();
    TVector3 b1=system.BoostVector();
    //cout<<"boost: ";(-b1).Print();
    p0.Boost(-b1);p1.Boost(-b1);system.Boost(-b1);dilepton.Boost(-b1);j0.Boost(-b1);
    //cout<<"###############after boost#########"<<endl;
    //p0.Print();p1.Print();system.Print();dilepton.Print();j0.Print();
    FillHist(pre+"p0j0angle_ppcm"+suf,p0.Angle(j0.Vect()),w,200,0,4);
    FillHist(pre+"p1j0angle_ppcm"+suf,p1.Angle(j0.Vect()),w,200,0,4);
    FillHist(pre+"p0j0cos_ppcm"+suf,cos(p0.Angle(j0.Vect())),w,100,-1,1);
    FillHist(pre+"p1j0cos_ppcm"+suf,cos(p1.Angle(j0.Vect())),w,100,-1,1);
    FillHist(pre+"j0direction_ppcm"+suf,p0.Angle(j0.Vect())<p1.Angle(j0.Vect()),w,5,0,5);
    if(p0.Angle(j0.Vect())<p1.Angle(j0.Vect())){
      p0-=j0;      
      cout<<"it is from p0: ";p0.Print();
    }else{
      p1-=j0;
      cout<<"it is from p1: ";p1.Print();
    }
    TVector3 b2=dilepton.BoostVector();
    //cout<<"boost2: ";(-b2).Print();
    p0.Boost(-b2);p1.Boost(-b2);system.Boost(-b2);dilepton.Boost(-b2);j0.Boost(-b2);
    //cout<<"###############after boost2#########"<<endl;
    //p0.Print();p1.Print();system.Print();dilepton.Print();j0.Print();
  }
  */
}
    
void AFBAnalyzer::FillAFBSystematicHists(TString pre,TString suf,const vector<Lepton*>& leps,const vector<Jet>& jets,map<TString,double> map_weight_systematic){
  for(auto iter=map_weight_systematic.begin();iter!=map_weight_systematic.end();iter++){
    FillAFBHists(pre,"_"+iter->first+suf,leps,jets,iter->second);
  }
}

void AFBAnalyzer::FillHardHists(TString pre,TString suf,Gen genparton0,Gen genparton1,Gen genhardl0,Gen genhardl1,Gen genhardj0,double w){
  Gen genhardl=genhardl0.PID()>0?genhardl0:genhardl1;
  TLorentzVector genZ=genhardl0+genhardl1;
  TLorentzVector genpp=genparton0+genparton1;

  FillHist(pre+"cos_l_p0"+suf,cos(genhardl.Angle(genparton0.Vect())),w,100,-1,1);
  FillHist(pre+"cos_l_p0_asym"+suf,cos(genhardl.Angle(genparton0.Vect())),w/2,100,-1,1);
  FillHist(pre+"cos_l_p0_asym"+suf,-1.*cos(genhardl.Angle(genparton0.Vect())),-w/2,100,-1,1);
  if(cos(genhardl.Angle(genparton0.Vect()))>0) FillHist(pre+"cos_l_p0_forward"+suf,genZ.M(),w,massbinnum,(double*)massrange);
  else FillHist(pre+"cos_l_p0_backward"+suf,genZ.M(),w,massbinnum,(double*)massrange);

  FillHist(pre+"cos_l_p1"+suf,cos(genhardl.Angle(genparton1.Vect())),w,100,-1,1);
  FillHist(pre+"cos_l_p1_asym"+suf,cos(genhardl.Angle(genparton1.Vect())),w/2,100,-1,1);
  FillHist(pre+"cos_l_p1_asym"+suf,-1.*cos(genhardl.Angle(genparton1.Vect())),-w/2,100,-1,1);
  if(cos(genhardl.Angle(genparton1.Vect()))>0) FillHist(pre+"cos_l_p1_forward"+suf,genZ.M(),w,massbinnum,(double*)massrange);
  else FillHist(pre+"cos_l_p1_backward"+suf,genZ.M(),w,massbinnum,(double*)massrange);
  
  FillHist(pre+"cos_Z_p0"+suf,cos(genZ.Angle(genparton0.Vect())),w,100,-1,1);
  FillHist(pre+"cos_Z_p0_asym"+suf,cos(genZ.Angle(genparton0.Vect())),w/2,100,-1,1);
  FillHist(pre+"cos_Z_p0_asym"+suf,-1.*cos(genZ.Angle(genparton0.Vect())),-w/2,100,-1,1);
  
  FillHist(pre+"cos_Z_p1"+suf,cos(genZ.Angle(genparton1.Vect())),w,100,-1,1);
  FillHist(pre+"cos_Z_p1_asym"+suf,cos(genZ.Angle(genparton1.Vect())),w/2,100,-1,1);
  FillHist(pre+"cos_Z_p1_asym"+suf,-1.*cos(genZ.Angle(genparton1.Vect())),-w/2,100,-1,1);
  
  
  FillHist(pre+"Zrap"+suf,genZ.Rapidity(),w,100,-5,5);
  FillHist(pre+"Zrap_asym"+suf,genZ.Rapidity(),w/2,100,-5,5);
  FillHist(pre+"Zrap_asym"+suf,-1.*genZ.Rapidity(),-w/2,100,-5,5);
  
  FillHist(pre+"pprap"+suf,genpp.Rapidity(),w,100,-5,5);
  FillHist(pre+"pprap_asym"+suf,genpp.Rapidity(),w/2,100,-5,5);
  FillHist(pre+"pprap_asym"+suf,-1.*genpp.Rapidity(),-w/2,100,-5,5);
  
  if(!genhardj0.IsEmpty()){
    FillHist(pre+"cos_l_j0"+suf,cos(genhardl.Angle(genhardj0.Vect())),w,100,-1,1);
    FillHist(pre+"cos_l_j0_asym"+suf,cos(genhardl.Angle(genhardj0.Vect())),w/2,100,-1,1);
    FillHist(pre+"cos_l_j0_asym"+suf,-1.*cos(genhardl.Angle(genhardj0.Vect())),-w/2,100,-1,1);
    if(cos(genhardl.Angle(genhardj0.Vect()))>0) FillHist(pre+"cos_l_j0_forward"+suf,genZ.M(),w,massbinnum,(double*)massrange);
    else FillHist(pre+"cos_l_j0_backward"+suf,genZ.M(),w,massbinnum,(double*)massrange);

    FillHist(pre+"cos_Z_j0"+suf,cos(genZ.Angle(genhardj0.Vect())),w,100,-1,1);
    FillHist(pre+"cos_Z_j0_asym"+suf,cos(genZ.Angle(genhardj0.Vect())),w/2,100,-1,1);
    FillHist(pre+"cos_Z_j0_asym"+suf,-1.*cos(genZ.Angle(genhardj0.Vect())),-w/2,100,-1,1);
    
    FillHist(pre+"cos_j0_p0"+suf,cos(genhardj0.Angle(genparton0.Vect())),w,100,-1,1);
    FillHist(pre+"cos_j0_p0_asym"+suf,cos(genhardj0.Angle(genparton0.Vect())),w/2,100,-1,1);
    FillHist(pre+"cos_j0_p0_asym"+suf,-1.*cos(genhardj0.Angle(genparton0.Vect())),-w/2,100,-1,1);
    
    FillHist(pre+"cos_j0_p1"+suf,cos(genhardj0.Angle(genparton1.Vect())),w,100,-1,1);
    FillHist(pre+"cos_j0_p1_asym"+suf,cos(genhardj0.Angle(genparton1.Vect())),w/2,100,-1,1);
    FillHist(pre+"cos_j0_p1_asym"+suf,-1.*cos(genhardj0.Angle(genparton1.Vect())),-w/2,100,-1,1);
    
    FillHist(pre+"jeta"+suf,genhardj0.Eta(),w,100,-5,5);
    FillHist(pre+"jeta_asym"+suf,genhardj0.Eta(),w/2,100,-5,5);
    FillHist(pre+"jeta_asym"+suf,-1.*genhardj0.Eta(),-w/2,100,-5,5);
  }
}
