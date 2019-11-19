#include "AFBAnalyzer.h"

void AFBAnalyzer::initializeAnalyzer(){
  SMPAnalyzerCore::initializeAnalyzer();
  if(HasFlag("BJET")){
    vector<Jet::Tagger> vtaggers={Jet::DeepCSV};
    vector<Jet::WP> vwps={Jet::Medium};
    SetupBTagger(vtaggers,vwps,true,true);
  }
}
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
      TLorentzVector gendilepton=genl0+genl1;
      TLorentzVector genharddilepton=genhardl0+genhardl1;
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
      map<TString,double> map_weight;
      map_weight[""]=weight_norm_1invpb*gen_weight*zptcor;
      map_weight["_nozptcor"]=weight_norm_1invpb*gen_weight;
      
      TString channelname=Form("%s%d",slepton.Data(),DataYear);
      if(HasFlag("GEN")) FillHistAll(channelname,"gen_","",&genl0,&genl1,map_weight);
      if(HasFlag("GENHARD")) FillHistAll(channelname,"genhard_","",&genhardl0,&genhardl1,map_weight);
      
    }else{
      tauprefix="tau_";
    }
  }

  if(HasFlag("GEN")||HasFlag("GENHARD")) return;

  if(!PassMETFilter()) return;

  Event* ev=new Event;
  *ev=GetEvent();
  if(DataYear==2016){
    vector<TString> muontrigger={"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v","HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"};
    TString electrontrigger="HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";
    if(ev->PassTrigger(muontrigger)||IsFastSim)
      if(!IsDATA||DataStream.Contains("DoubleMuon")) executeEventFromParameter("muon2016",ev);
    if(ev->PassTrigger(electrontrigger)||IsFastSim)
      if(!IsDATA||DataStream.Contains("DoubleEG")) executeEventFromParameter("electron2016",ev);
  }else if(DataYear==2017){
    TString muontrigger="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v";
    TString electrontrigger="HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v";
    if(ev->PassTrigger(muontrigger))
      if(!IsDATA||DataStream.Contains("DoubleMuon")) executeEventFromParameter("muon2017",ev);
    if(ev->PassTrigger(electrontrigger))
      if(!IsDATA||DataStream.Contains("DoubleEG")) executeEventFromParameter("electron2017",ev);
  }else if(DataYear==2018){
    TString muontrigger="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v";
    TString electrontrigger="HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";
    if(ev->PassTrigger(muontrigger))
      if(!IsDATA||DataStream.Contains("DoubleMuon")) executeEventFromParameter("muon2018",ev);
    if(ev->PassTrigger(electrontrigger))
      if(!IsDATA||DataStream.Contains("EGamma")) executeEventFromParameter("electron2018",ev);
  }    
  
  delete ev;
}

void AFBAnalyzer::executeEventFromParameter(TString channelname,Event* ev){
  map<TString,vector<Muon>> map_muons;
  map<TString,vector<Electron>> map_electrons;
  //suffix, lepton vector, SF1, SF2, TriggerSF1, TriggerSF2
  map<TString,tuple<vector<Lepton*>,TString,TString,TString,TString>> map_leps;

  std::vector<Jet> jets=GetJets("tightLepVeto",30,2.7);
  std::sort(jets.begin(),jets.end(),PtComparing);

  bool IsSYS=HasFlag("SYS")||HasFlag("DYSYS");
  
  double lep0ptcut,lep1ptcut;
  if(channelname.Contains("muon")){
    lep0ptcut=20.;
    lep1ptcut=10.;
    
    map_muons[""]=MuonMomentumCorrection(SMPGetMuons("POGTightWithTightIso",0.0,2.4),0);
    std::sort(map_muons[""].begin(),map_muons[""].end(),PtComparing);
    map_leps[""]=make_tuple(MakeLeptonPointerVector(map_muons[""]),"ID_SF_NUM_TightID_DEN_genTracks","ISO_SF_NUM_TightRelIso_DEN_TightIDandIPCut",DataYear==2018?"":"LeadMu17_POGTight",DataYear==2018?"":"TailMu8_POGTight");
    
    if(HasFlag("SYS")){
      map_muons["_scale_up"]=MuonMomentumCorrection(map_muons[""],+1);
      std::sort(map_muons["_scale_up"].begin(),map_muons["_scale_up"].end(),PtComparing);
      map_leps["_scale_up"]=make_tuple(MakeLeptonPointerVector(map_muons["_scale_up"]),"ID_SF_NUM_TightID_DEN_genTracks","ISO_SF_NUM_TightRelIso_DEN_TightIDandIPCut",DataYear==2018?"":"LeadMu17_POGTight",DataYear==2018?"":"TailMu8_POGTight");
      
      map_muons["_scale_down"]=MuonMomentumCorrection(map_muons[""],-1);
      std::sort(map_muons["_scale_down"].begin(),map_muons["_scale_down"].end(),PtComparing);
      map_leps["_scale_down"]=make_tuple(MakeLeptonPointerVector(map_muons["_scale_down"]),"ID_SF_NUM_TightID_DEN_genTracks","ISO_SF_NUM_TightRelIso_DEN_TightIDandIPCut",DataYear==2018?"":"LeadMu17_POGTight",DataYear==2018?"":"TailMu8_POGTight");

      map_muons["_noroccor"]=MuonMomentumCorrection(map_muons[""],0,-1);
      std::sort(map_muons["_noroccor"].begin(),map_muons["_noroccor"].end(),PtComparing);
      map_leps["_noroccor"]=make_tuple(MakeLeptonPointerVector(map_muons["_noroccor"]),"ID_SF_NUM_TightID_DEN_genTracks","ISO_SF_NUM_TightRelIso_DEN_TightIDandIPCut",DataYear==2018?"":"LeadMu17_POGTight",DataYear==2018?"":"TailMu8_POGTight");
    }
  }else if(channelname.Contains("electron")){
    lep0ptcut=25.;
    lep1ptcut=15.;
    TString LeptonIDSF_key,triggerSF_key0,triggerSF_key1;
    switch(DataYear){
    case 2016: LeptonIDSF_key="ID_SF_MediumID_pt10";triggerSF_key0="LeadEle23_MediumID";triggerSF_key1="TailEle12_MediumID"; break;
    case 2017: LeptonIDSF_key="ID_SF_MediumID_pt10_Q";triggerSF_key0="LeadEle23_MediumID_Q";triggerSF_key1="TailEle12_MediumID_Q"; break;
    case 2018: LeptonIDSF_key="ID_SF_passMediumID";triggerSF_key0="";triggerSF_key1=""; break;
    default: cout<<"[AFBAnalyzer::executeEventFromParameter] wrong year"<<endl;exit(EXIT_FAILURE);break;
    }
    
    map_electrons[""]=SMPGetElectrons("passMediumID",0.0,2.5);
    map_leps[""]=make_tuple(MakeLeptonPointerVector(map_electrons[""]),LeptonIDSF_key,"Default",triggerSF_key0,triggerSF_key1);

    if(HasFlag("SYS")){
      //map_electrons["_roccor0"]=ElectronEnergyCorrection(map_electrons[""],0,0);
      //std::sort(map_electrons["_roccor0"].begin(),map_electrons["_roccor0"].end(),PtComparing);
      //map_leps["_roccor0"]=make_tuple(MakeLeptonPointerVector(map_electrons["_roccor0"]),LeptonIDSF_key,"Default",triggerSF_key0,triggerSF_key1);
      
      map_electrons["_scale_up"]=ScaleElectrons(map_electrons[""],1);
      std::sort(map_electrons["_scale_up"].begin(),map_electrons["_scale_up"].end(),PtComparing);
      map_leps["_scale_up"]=make_tuple(MakeLeptonPointerVector(map_electrons["_scale_up"]),LeptonIDSF_key,"Default",triggerSF_key0,triggerSF_key1);
      
      map_electrons["_scale_down"]=ScaleElectrons(map_electrons[""],-1);
      std::sort(map_electrons["_scale_down"].begin(),map_electrons["_scale_down"].end(),PtComparing);
      map_leps["_scale_down"]=make_tuple(MakeLeptonPointerVector(map_electrons["_scale_down"]),LeptonIDSF_key,"Default",triggerSF_key0,triggerSF_key1);
      
      map_electrons["_smear_up"]=SmearElectrons(map_electrons[""],1);
      std::sort(map_electrons["_smear_up"].begin(),map_electrons["_smear_up"].end(),PtComparing);
      map_leps["_smear_up"]=make_tuple(MakeLeptonPointerVector(map_electrons["_smear_up"]),LeptonIDSF_key,"Default",triggerSF_key0,triggerSF_key1);
      
      map_electrons["_smear_down"]=SmearElectrons(map_electrons[""],-1);
      std::sort(map_electrons["_smear_down"].begin(),map_electrons["_smear_down"].end(),PtComparing);
      map_leps["_smear_down"]=make_tuple(MakeLeptonPointerVector(map_electrons["_smear_down"]),LeptonIDSF_key,"Default",triggerSF_key0,triggerSF_key1);
    }
  }else{
    cout<<"[AFBAnalyzer::executeEventFromParameter] wrong channelname"<<endl;
    return;
  }
  
  /////////////////lumi weight///////////////////
  double weight=1.,totalweight=1.;
  if(!IsDATA){
    weight=weight_norm_1invpb*ev->MCweight()*ev->GetTriggerLumi("Full");
  }
  totalweight*=weight;
  if(!IsSYS) FillCutflow(channelname+"/"+tauprefix+"cutflow","lumi",totalweight);
  
  /////////////////PUreweight///////////////////
  double PUreweight=1.,PUreweight_up=1.,PUreweight_down=1.;
  if(!IsDATA){
    PUreweight=mcCorr->GetPileUpWeight(nPileUp,0);
    PUreweight_up=mcCorr->GetPileUpWeight(nPileUp,1);
    PUreweight_down=mcCorr->GetPileUpWeight(nPileUp,-1);
  }
  totalweight*=PUreweight;
  if(!IsSYS) FillCutflow(channelname+"/"+tauprefix+"cutflow","PU",totalweight);
  
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
  if(!IsSYS) FillCutflow(channelname+"/"+tauprefix+"cutflow","prefire",totalweight);

  totalweight*=zptcor;
  if(!IsSYS) FillCutflow(channelname+"/"+tauprefix+"cutflow","zptcor",totalweight);

  //////////////////////Z0 weight///////////////////////
  double z0weight=1.;
  if(!IsDATA){
    z0weight=GetZ0Weight(vertex_Z);
  }
  totalweight*=z0weight;
  if(!IsSYS) FillCutflow(channelname+"/"+tauprefix+"cutflow","z0",totalweight);

  if(HasFlag("HIGHMET")){
    if(pfMET_Type1_pt<30) return;
    if(!IsSYS) FillCutflow(channelname+"/"+tauprefix+"cutflow","METCut",totalweight);
  }

  int n_bjet=0;
  if(HasFlag("BJET")){
    for(const auto& jet:jets){
      if(IsBTagged(jet,Jet::DeepCSV,Jet::Medium,true,0)) n_bjet++;
    }
    if(!n_bjet) return;
    if(!IsSYS) FillCutflow(channelname+"/"+tauprefix+"cutflow","BJetCut",totalweight);    
  }


  ///////////////////////lepton selection///////////////////////
  for(const auto& element_leps: map_leps){
    TString suffix=element_leps.first;
    bool IsNominal=suffix==""?true:false;

    auto const& leps=get<0>(element_leps.second);

    FillHist(channelname+"/"+tauprefix+"nlepton"+suffix,leps.size(),totalweight,10,0,10);
    if(leps.size()==2){
      TString prefix=tauprefix;
      if(HasFlag("CF")){
	if(leps.at(0)->Charge()>0&&leps.at(1)->Charge()>0) prefix="pp_"+prefix;
	else if(leps.at(0)->Charge()<0&&leps.at(1)->Charge()<0) prefix="mm_"+prefix;
	else return;
      }else{
	if(leps.at(0)->Charge()*leps.at(1)->Charge()>0) prefix="ss_"+prefix;
      }
      if(IsNominal) FillCutflow(channelname+"/"+tauprefix+"cutflow","dilepton",totalweight);
      if(leps.at(0)->Pt()>lep0ptcut&&leps.at(1)->Pt()>lep1ptcut){
	if(IsNominal) FillCutflow(channelname+"/"+tauprefix+"cutflow","ptcut",totalweight);
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
      
	double triggerSF=1.,triggerSF_up=1.,triggerSF_down=1.;
	if(!IsDATA){
	  TString triggerSF_key0=get<3>(element_leps.second);
	  TString triggerSF_key1=get<4>(element_leps.second);
	  triggerSF*=DileptonTrigger_SF(triggerSF_key0,triggerSF_key1,leps,0);
	  triggerSF_up*=DileptonTrigger_SF(triggerSF_key0,triggerSF_key1,leps,1);
	  triggerSF_down*=DileptonTrigger_SF(triggerSF_key0,triggerSF_key1,leps,-1);
	}

	if(IsNominal) FillCutflow(channelname+"/"+tauprefix+"cutflow","RECO",totalweight*RECOSF);
	if(IsNominal) FillCutflow(channelname+"/"+tauprefix+"cutflow","ID",totalweight*RECOSF*IDSF);
	if(IsNominal) FillCutflow(channelname+"/"+tauprefix+"cutflow","ISO",totalweight*RECOSF*IDSF*ISOSF);
	if(IsNominal) FillCutflow(channelname+"/"+tauprefix+"cutflow","trigger",totalweight*RECOSF*IDSF*ISOSF);

	///////////////////////map_weight//////////////////
	map<TString,double> map_weight;
	if(!IsSYS||!IsNominal) map_weight[suffix]=weight*PUreweight*RECOSF*IDSF*ISOSF*triggerSF*prefireweight*zptcor*z0weight;
	if(IsNominal&&!IsDATA){
	  if(HasFlag("SYS")){
	    map_weight["_noefficiencySF"]=weight*PUreweight*prefireweight*zptcor*z0weight;

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
	    
	    map_weight["_nozptcor"]=weight*PUreweight*RECOSF*IDSF*ISOSF*triggerSF*prefireweight*z0weight;

	    map_weight["_noz0weight"]=weight*PUreweight*RECOSF*IDSF*ISOSF*triggerSF*prefireweight*zptcor;
	    
	    map_weight["_noPUreweight"]=weight*RECOSF*IDSF*ISOSF*triggerSF*prefireweight*zptcor*z0weight;
	    map_weight["_PUreweight_up"]=weight*PUreweight_up*RECOSF*IDSF*ISOSF*triggerSF*prefireweight*zptcor*z0weight;
	    map_weight["_PUreweight_down"]=weight*PUreweight_down*RECOSF*IDSF*ISOSF*triggerSF*prefireweight*zptcor*z0weight;
	    
	    map_weight["_noprefireweight"]=weight*PUreweight*RECOSF*IDSF*ISOSF*triggerSF*zptcor*z0weight;
	    map_weight["_prefireweight_up"]=weight*PUreweight*RECOSF*IDSF*ISOSF*triggerSF*prefireweight_up*zptcor*z0weight;
	    map_weight["_prefireweight_down"]=weight*PUreweight*RECOSF*IDSF*ISOSF*triggerSF*prefireweight_down*zptcor*z0weight;
	    
	  }
	  if(HasFlag("DYSYS")){
	    if(IsDYSample){	      
	      for(unsigned int i=0;i<PDFWeights_Scale->size();i++){
		map_weight[Form("_scalevariation%d",i)]=weight*PUreweight*RECOSF*IDSF*ISOSF*triggerSF*prefireweight*zptcor*z0weight*PDFWeights_Scale->at(i);
	      }
	    }
	    for(unsigned int i=0;i<PDFWeights_Error->size();i++){
	      map_weight[Form("_pdf%d",i)]=weight*PUreweight*RECOSF*IDSF*ISOSF*triggerSF*prefireweight*zptcor*z0weight*PDFWeights_Error->at(i);
	    }
	    if(PDFWeights_AlphaS->size()==2){
	      map_weight["_alphaS_up"]=weight*PUreweight*RECOSF*IDSF*ISOSF*triggerSF*prefireweight*zptcor*z0weight*PDFWeights_AlphaS->at(0);
	      map_weight["_alphaS_down"]=weight*PUreweight*RECOSF*IDSF*ISOSF*triggerSF*prefireweight*zptcor*z0weight*PDFWeights_AlphaS->at(1);
	    }
	  }
	}

	///////////////////////fill hists///////////////////////
	FillHistAll(channelname,prefix,"",(Particle*)leps[0],(Particle*)leps[1],map_weight);
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
double AFBAnalyzer::GetCosThetaCS(const Particle *p0,const Particle *p1){
  const TLorentzVector *l0,*l1;
  if(p0->Charge()<0&&p1->Charge()>0){
    l0=p0;
    l1=p1;
  }else if(p0->Charge()>0&&p1->Charge()<0){
    l0=p1;
    l1=p0;
  }else{
    if(random->Rndm()<0.5){
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
void AFBAnalyzer::FillHistAll(TString channelname,TString pre,TString suf,Particle* l0,Particle* l1,map<TString,double> map_weight){
  TLorentzVector dilepton=(*l0)+(*l1);
  double dimass=dilepton.M();
  double dirap=dilepton.Rapidity();
  double dipt=dilepton.Pt();
  
  for(int im=0;im<mbinnum;im++){
    if(dimass>=mbin[im]&&dimass<mbin[im+1]){
      for(int ir=0;ir<ybinnum;ir++){
	if(fabs(dirap)>=ybin[ir]&&fabs(dirap)<ybin[ir+1]){
	  for(int ip=0;ip<ptbinnum;ip++){
	    if(dipt>=ptbin[ip]&&dipt<ptbin[ip+1]){
	      FillAFBSystematicHists(channelname+Form("/m[%.0f,%.0f]/y[%.1f,%.1f]/pt[%.0f,%.0f]/",mbin[im],mbin[im+1],ybin[ir],ybin[ir+1],ptbin[ip],ptbin[ip+1])+pre,suf,l0,l1,map_weight);
	      break;
	    }
	  }
	  break;
	}
      }
      break;
    }
  }  
}
void AFBAnalyzer::FillAFBSystematicHists(TString pre,TString suf,Particle* l0,Particle* l1,map<TString,double> map_weight_systematic){
  for(auto iter=map_weight_systematic.begin();iter!=map_weight_systematic.end();iter++){
    FillAFBHists(pre,iter->first+suf,l0,l1,iter->second);
  }
}
void AFBAnalyzer::FillAFBHists(TString pre,TString suf,Particle* l0,Particle* l1,double w){
  //FillDileptonHists(pre,suf,l0,l1,w);
  
  if(HasFlag("LEP")){
    FillHist(pre+"l0pt"+suf,l0->Pt(),w,lptbinnum,(double*)lptbin);
    FillHist(pre+"l1pt"+suf,l1->Pt(),w,lptbinnum,(double*)lptbin);
    FillHist(pre+"lpt"+suf,l0->Pt(),w,lptbinnum,(double*)lptbin);
    FillHist(pre+"lpt"+suf,l1->Pt(),w,lptbinnum,(double*)lptbin);
    FillHist(pre+"l0eta"+suf,l0->Eta(),w,60,-3,3);
    FillHist(pre+"l1eta"+suf,l1->Eta(),w,60,-3,3);
    FillHist(pre+"leta"+suf,l0->Eta(),w,60,-3,3);
    FillHist(pre+"leta"+suf,l1->Eta(),w,60,-3,3);
  }
  if(HasFlag("ETC")){
    FillHist(pre+"z0"+suf,vertex_Z,w,120,-15,15);
  }
  if(!HasFlag("LEP")&&!HasFlag("ETC")){
    TLorentzVector dilepton=*l0+*l1;
    double dimass=dilepton.M();
    double dirap=dilepton.Rapidity();
    double dipt=dilepton.Pt();
    FillHist(pre+"dimass"+suf,dimass,w,mbinnum_fine,(double*)mbin_fine);
    FillHist(pre+"dirap"+suf,dirap,w,60,-3,3);
    FillHist(pre+"dipt"+suf,dipt,w,ptbinnum_fine,(double*)ptbin_fine);
    
    double cost=GetCosThetaCS(l0,l1);
    double h=0.5*pow(dipt/dimass,2)/(1+pow(dipt/dimass,2))*(1-3*cost*cost);
    double den_weight=0.5*fabs(cost)/pow(1+cost*cost+h,2);
    double num_weight=0.5*cost*cost/pow(1+cost*cost+h,3);
    FillHist(pre+"costhetaCS"+suf,cost,w,40,-1,1);
    //FillHist(pre+"abscosthetaCS"+suf,fabs(cost),w,20,0,1);
    if(cost>0){
      FillHist(pre+"forward"+suf,0,w,3,0,3);
      FillHist(pre+"forward"+suf,1,w*den_weight,3,0,3);
      FillHist(pre+"forward"+suf,2,w*num_weight,3,0,3);    
    }else{
      FillHist(pre+"backward"+suf,0,w,3,0,3);
      FillHist(pre+"backward"+suf,1,w*den_weight,3,0,3);
      FillHist(pre+"backward"+suf,2,w*num_weight,3,0,3);
    }
  }
}





/*
void AFBAnalyzer::FillAFBHists(TString pre,TString suf,const vector<Lepton*>& leps,const vector<Jet>& jets,double w){
  
  if(leps.at(0)->Pt()>60){
    if(leps.at(0)->Charge()>0) FillHist(pre+"l0eta_pt60_plus"+suf,leps.at(0)->Eta(),w,60,-3,3);
    else FillHist(pre+"l0eta_pt60_minus"+suf,leps.at(0)->Eta(),w,60,-3,3);
  }
  //for leptons
  FillBasicHists(pre,suf,leps,w);
  for(int i=0;i<(int)leps.size();i++){
    TString charge=leps.at(i)->Charge()>0?"p":"m";
    FillHist(Form("%sl%d%spt%s",pre.Data(),i,charge.Data(),suf.Data()),leps.at(i)->Pt(),w,1000,0,1000);
    FillHist(Form("%sl%d%seta%s",pre.Data(),i,charge.Data(),suf.Data()),leps.at(i)->Eta(),w,200,-5,5);
    FillHist(Form("%sl%d%sriso%s",pre.Data(),i,charge.Data(),suf.Data()),leps.at(i)->RelIso(),w,30,0,0.3);
    FillHist(Form("%sl%spt%s",pre.Data(),charge.Data(),suf.Data()),leps.at(i)->Pt(),w,1000,0,1000);
    FillHist(Form("%sl%seta%s",pre.Data(),charge.Data(),suf.Data()),leps.at(i)->Eta(),w,200,-5,5);
    FillHist(Form("%sl%sriso%s",pre.Data(),charge.Data(),suf.Data()),leps.at(i)->RelIso(),w,30,0,0.3);
    if(leps.at(i)->LeptonFlavour()==Lepton::Flavour::MUON){
      double rtrkiso=((Muon*)leps.at(i))->TrkIso()/leps.at(i)->Pt();
      FillHist(Form("%sl%drtrkiso%s",pre.Data(),i,suf.Data()),rtrkiso,w,40,0,0.2);
      FillHist(Form("%sl%d%srtrkiso%s",pre.Data(),i,charge.Data(),suf.Data()),rtrkiso,w,40,0,0.2);
      FillHist(Form("%sl%srtrkiso%s",pre.Data(),charge.Data(),suf.Data()),rtrkiso,w,40,0,0.2);
    }
  }
  
  TLorentzVector dilepton=(*leps.at(0))+(*leps.at(1));
  double dimass=dilepton.M();
  double dipt=dilepton.Pt();

  double cost=GetCosThetaCS((Particle*)leps.at(0),(Particle*)leps.at(1));
  FillHist(pre+"costhetaCS"+suf,cost,w,40,-1,1);
  FillHist(pre+"abscosthetaCS"+suf,fabs(cost),w,20,0,1);
  if(pre.Contains(TRegexp("m60to..."))){
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
  }
  
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
	  FillHist(pre+"forward_"+options[i]+"_"+Form("%d",int(fcuts[j]*10))+suf,dimass,w,mbinnum,(double*)mbin);
	  FillHist(pre+"forward_den_"+options[i]+"_"+Form("%d",int(fcuts[j]*10))+suf,dimass,w*0.5*fabs(cost_)/pow(1+cost_*cost_+h,2),mbinnum,(double*)mbin);
	  FillHist(pre+"forward_num_"+options[i]+"_"+Form("%d",int(fcuts[j]*10))+suf,dimass,w*0.5*cost_*cost_/pow(1+cost_*cost_+h,3),mbinnum,(double*)mbin);    
	}else{
	  FillHist(pre+"backward_"+options[i]+"_"+Form("%d",int(fcuts[j]*10))+suf,dimass,w,mbinnum,(double*)mbin);
	  FillHist(pre+"backward_den_"+options[i]+"_"+Form("%d",int(fcuts[j]*10))+suf,dimass,w*0.5*fabs(cost_)/pow(1+cost_*cost_+h,2),mbinnum,(double*)mbin);
	  FillHist(pre+"backward_num_"+options[i]+"_"+Form("%d",int(fcuts[j]*10))+suf,dimass,w*0.5*cost_*cost_/pow(1+cost_*cost_+h,3),mbinnum,(double*)mbin);   
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
  
}
*/    
void AFBAnalyzer::FillHardHists(TString pre,TString suf,const Gen& genparton0,const Gen& genparton1,const Gen& genhardl0,const Gen& genhardl1,const Gen& genhardj0,double w){
  Gen genhardl=genhardl0.PID()>0?genhardl0:genhardl1;
  TLorentzVector genZ=genhardl0+genhardl1;
  TLorentzVector genpp=genparton0+genparton1;

  FillHist(pre+"cos_l_p0"+suf,cos(genhardl.Angle(genparton0.Vect())),w,100,-1,1);
  FillHist(pre+"cos_l_p0_asym"+suf,cos(genhardl.Angle(genparton0.Vect())),w/2,100,-1,1);
  FillHist(pre+"cos_l_p0_asym"+suf,-1.*cos(genhardl.Angle(genparton0.Vect())),-w/2,100,-1,1);
  if(cos(genhardl.Angle(genparton0.Vect()))>0) FillHist(pre+"cos_l_p0_forward"+suf,genZ.M(),w,mbinnum,(double*)mbin);
  else FillHist(pre+"cos_l_p0_backward"+suf,genZ.M(),w,mbinnum,(double*)mbin);

  FillHist(pre+"cos_l_p1"+suf,cos(genhardl.Angle(genparton1.Vect())),w,100,-1,1);
  FillHist(pre+"cos_l_p1_asym"+suf,cos(genhardl.Angle(genparton1.Vect())),w/2,100,-1,1);
  FillHist(pre+"cos_l_p1_asym"+suf,-1.*cos(genhardl.Angle(genparton1.Vect())),-w/2,100,-1,1);
  if(cos(genhardl.Angle(genparton1.Vect()))>0) FillHist(pre+"cos_l_p1_forward"+suf,genZ.M(),w,mbinnum,(double*)mbin);
  else FillHist(pre+"cos_l_p1_backward"+suf,genZ.M(),w,mbinnum,(double*)mbin);
  
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
    if(cos(genhardl.Angle(genhardj0.Vect()))>0) FillHist(pre+"cos_l_j0_forward"+suf,genZ.M(),w,mbinnum,(double*)mbin);
    else FillHist(pre+"cos_l_j0_backward"+suf,genZ.M(),w,mbinnum,(double*)mbin);

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
  
