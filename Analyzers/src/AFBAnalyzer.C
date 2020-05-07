#include "AFBAnalyzer.h"

void AFBAnalyzer::initializeAnalyzer(){
  SMPAnalyzerCore::initializeAnalyzer(); //setup zpt roc z0 
  SetupCosThetaWeight();
  
  if(HasFlag("bjet")||HasFlag("nobjet")){
    vector<JetTagging::Parameters> jtps={JetTagging::Parameters(JetTagging::DeepCSV,JetTagging::Medium,JetTagging::mujets,JetTagging::mujets)};
    mcCorr->SetJetTaggingParameters(jtps);
  }
  SetupToy(100);
  IsNominalRun=!HasFlag("SYS")&&!HasFlag("PDFSYS");
  if(fChain->GetListOfFiles()->GetEntries()){
    TString filename=fChain->GetListOfFiles()->At(0)->GetTitle();
    if(filename.Contains("SkimTree_")) IsSkimmed=true;
    else IsSkimmed=false;
  }else{
    cout<<"[AFBAnalyzer::initializeAnalyzer] no input file"<<endl;
    exit(EXIT_FAILURE);
  }
  if(!IsSkimmed&&!HasFlag("ALL")){
    fChain->SetBranchStatus("pfMET_*",false);
    fChain->SetBranchStatus("HLT_TriggerName",false);
    fChain->SetBranchStatus("jet_*",false);
    fChain->SetBranchStatus("fatjet_*",false);
    fChain->SetBranchStatus("electron_*",false);
    fChain->SetBranchStatus("muon_*",false);
    fChain->SetBranchStatus("photon_*",false);
  }
}
void AFBAnalyzer::executeEvent(){
  GetToyWeight();

  if(!HLT_TriggerName) HLT_TriggerName=new vector<string>;
  event=GetEvent();
  GetEventWeights();

  hardprefix="";
  if(IsDYSample){
    //////////////////////// Check LHE /////////////////////////
    vector<LHE> lhes=GetLHEs();
    LHE lhe_l0,lhe_l1;
    GetDYLHEParticles(lhes,lhe_l0,lhe_l1);
    TString channelname="";
    if(abs(lhe_l0.ID())!=15){
      double l0ptcut,l1ptcut,letacut;
      if(abs(lhe_l0.ID())==11){
	channelname=Form("ee%d",DataYear);
	l0ptcut=25.;
	l1ptcut=15.;
	letacut=2.4;
      }else if(abs(lhe_l0.ID())==13){
	channelname=Form("mm%d",DataYear);
	l0ptcut=20.;
	l1ptcut=10.;
	letacut=2.4;
      }else{
	cout<<"[AFBAnalyzer::executeEvent()] something is wrong l0.ID="<<abs(lhe_l0.ID())<<endl;
	for(auto& lhe:lhes) lhe.Print();
	exit(EXIT_FAILURE);
      }
      
      //////////////////////// GEN /////////////////////////
      vector<Gen> gens=GetGens();
      Gen gen_parton0,gen_parton1,gen_l0,gen_l1;
      GetDYGenParticles(gens,gen_parton0,gen_parton1,gen_l0,gen_l1,true);
      /*
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
    */
      //genparton0.SetPxPyPzE(0,0,genWeight_X1*6500.,genWeight_X1*6500.);
      //genparton1.SetPxPyPzE(0,0,-genWeight_X2*6500.,genWeight_X2*6500.);
      //for(unsigned int i=0;i<photons.size();i++) genphotons+=gens[photons[i]];
      /*
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
      */
     
      map<TString,double> map_weight;
      map_weight[""]=lumiweight*zptweight;
      map_weight["_nozptweight"]=lumiweight;

      //////////////// Fill LHE,Gen hists //////////////////////
      if(IsNominalRun&&!IsSkimmed&&!HasFlag("ALL")){
	FillHists(channelname,"lhe_","",(Particle*)&lhe_l0,(Particle*)&lhe_l1,map_weight);
	FillHists(channelname,"gen_","",(Particle*)&gen_l0,(Particle*)&gen_l1,map_weight);
	if(gen_l0.Pt()>l0ptcut&&gen_l1.Pt()>l1ptcut&&fabs(gen_l0.Eta())<letacut&&fabs(gen_l1.Eta())<letacut){
	  FillHists(channelname,"genfid_","",(Particle*)&gen_l0,(Particle*)&gen_l1,map_weight);
	}
      }
    }
  }

  if(!IsSkimmed&&!HasFlag("ALL")) return;

  if(!PassMETFilter()) return;

  TString prefix="";
  if(HasFlag("bjet")) prefix+="bjet/";
  else if(HasFlag("nobjet")) prefix+="nobjet/";
  if(HasFlag("highmet")) prefix+="highmet/";
			
  ///////////////// cutflow ///////////////////
  if(IsNominalRun){
    FillCutflow(prefix+tauprefix+"cutflow","lumi",lumiweight);
    FillCutflow(prefix+tauprefix+"cutflow","PU",lumiweight*PUweight);
    FillCutflow(prefix+tauprefix+"cutflow","prefire",lumiweight*PUweight*prefireweight);
    FillCutflow(prefix+tauprefix+"cutflow","zpt",lumiweight*PUweight*prefireweight*zptweight);
    FillCutflow(prefix+tauprefix+"cutflow","z0",lumiweight*PUweight*prefireweight*zptweight*z0weight);
  }

  if(HasFlag("highmet")){
    if(pfMET_Type1_pt<60) return;
    if(IsNominalRun) FillCutflow(prefix+tauprefix+"cutflow","METCut",lumiweight*PUweight*prefireweight*zptweight*z0weight);
  }

  int n_bjet=0;
  if(HasFlag("bjet")||HasFlag("nobjet")){
    std::vector<Jet> jets=GetJets("tightLepVeto",30,2.7);
    std::sort(jets.begin(),jets.end(),PtComparing);

    JetTagging::Parameters jtp = JetTagging::Parameters(JetTagging::DeepCSV,JetTagging::Medium,JetTagging::mujets,JetTagging::mujets);
    for(const auto& jet:jets)
      if(mcCorr->IsBTagged_2a(jtp,jet))
	n_bjet++;
    
    if(HasFlag("bjet")&&!n_bjet) return;
    if(HasFlag("nobjet")&&n_bjet) return;
    if(IsNominalRun) FillCutflow(prefix+tauprefix+"cutflow","BJetCut",lumiweight*PUweight*prefireweight*zptweight*z0weight);
  }

  if(DataYear==2016){
    vector<TString> muontrigger={
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v",
      "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
      "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
      "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",
      "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
    };
    vector<TString> emutrigger={
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
    };
    if(!HasFlag("emu") && event.PassTrigger(muontrigger))
      if(!IsDATA||DataStream.Contains("DoubleMuon")) executeEventWithChannelName(prefix+"mm2016");
    if(!HasFlag("emu") && event.PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"))
      if(!IsDATA||DataStream.Contains("DoubleEG")) executeEventWithChannelName(prefix+"ee2016");
    if(HasFlag("emu") && event.PassTrigger(emutrigger))
      if(!IsDATA||DataStream.Contains("MuonEG")) executeEventWithChannelName(prefix+"em2016");
  }else if(DataYear==2017){
    vector<TString> emutrigger={
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
    };
    if(!HasFlag("emu") && event.PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v"))
      if(!IsDATA||DataStream.Contains("DoubleMuon")) executeEventWithChannelName(prefix+"mm2017");
    if(!HasFlag("emu") && event.PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"))
      if(!IsDATA||DataStream.Contains("DoubleEG")) executeEventWithChannelName(prefix+"ee2017");
    if(HasFlag("emu") && event.PassTrigger(emutrigger))
      if(!IsDATA||DataStream.Contains("MuonEG")) executeEventWithChannelName(prefix+"em2017");
  }else if(DataYear==2018){
    vector<TString> emutrigger={
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
    };
    if(!HasFlag("emu") && event.PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"))
      if(!IsDATA||DataStream.Contains("DoubleMuon")) executeEventWithChannelName(prefix+"mm2018");
    if(!HasFlag("emu") && event.PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"))
      if(!IsDATA||DataStream.Contains("EGamma")) executeEventWithChannelName(prefix+"ee2018");
    if(HasFlag("emu")  && event.PassTrigger(emutrigger))
      if(!IsDATA||DataStream.Contains("MuonEG")) executeEventWithChannelName(prefix+"em2018");
  }    

}

void AFBAnalyzer::executeEventWithChannelName(TString channelname){
  map<TString,vector<Muon>> map_muons;
  map<TString,vector<Electron>> map_electrons;
  map<TString,Parameter> map_parameter;
  
  if(channelname.Contains(TRegexp("mm20[0-9][0-9]"))){
    Parameter p("IDISO_SF_MediumID_trkIsoLoose_Q","",{"Mu17Leg1_MediumID_trkIsoLoose_Q","Mu8Leg2_MediumID_trkIsoLoose_Q"},20.,10.);
    
    map_muons[""]=MuonMomentumCorrection(SMPGetMuons("POGMediumWithLooseTrkIso",0.0,2.4),0);
    map_parameter[""]=p.Clone(MakeLeptonPointerVector(map_muons[""]),
			      (IsNominalRun?NominalWeight:0)
			      +(HasFlag("SYS")&&!IsDATA?SystematicWeight:0)
			      +(HasFlag("PDFSYS")&&!IsDATA?PDFWeight:0)
			      );
    
    if(HasFlag("SYS")){
      map_muons["_scale_up"]=MuonMomentumCorrection(map_muons[""],+1);
      map_parameter["_scale_up"]=p.Clone(MakeLeptonPointerVector(map_muons["_scale_up"]));
					 
      map_muons["_scale_down"]=MuonMomentumCorrection(map_muons[""],-1);
      map_parameter["_scale_down"]=p.Clone(MakeLeptonPointerVector(map_muons["_scale_down"]));
      
      map_muons["_noroccor"]=MuonMomentumCorrection(map_muons[""],0,-1);
      map_parameter["_noroccor"]=p.Clone(MakeLeptonPointerVector(map_muons["_noroccor"]));
    }
  }else if(channelname.Contains(TRegexp("ee20[0-9][0-9]"))){
    Parameter p("ID_SF_MediumID_Q",{"Ele23Leg1_MediumID_Q","Ele12Leg2_MediumID_Q"},25.,15.);
    
    map_electrons["_noroccor"]=SMPGetElectrons("passMediumID",0.0,2.4);
    map_electrons[""]=ElectronEnergyCorrection(map_electrons["_noroccor"],0,0);
    map_parameter[""]=p.Clone(MakeLeptonPointerVector(map_electrons[""]),
			      (IsNominalRun?NominalWeight:0)
			      +(HasFlag("SYS")&&!IsDATA?SystematicWeight:0)
			      +(HasFlag("PDFSYS")&&!IsDATA?PDFWeight:0)
			      );

    if(HasFlag("SYS")){
      map_electrons["_scale_up"]=ScaleElectrons(map_electrons[""],1);
      std::sort(map_electrons["_scale_up"].begin(),map_electrons["_scale_up"].end(),PtComparing);
      map_parameter["_scale_up"]=p.Clone(MakeLeptonPointerVector(map_electrons["_scale_up"]));
      
      map_electrons["_scale_down"]=ScaleElectrons(map_electrons[""],-1);
      std::sort(map_electrons["_scale_down"].begin(),map_electrons["_scale_down"].end(),PtComparing);
      map_parameter["_scale_down"]=p.Clone(MakeLeptonPointerVector(map_electrons["_scale_down"]));
      
      map_electrons["_smear_up"]=SmearElectrons(map_electrons[""],1);
      std::sort(map_electrons["_smear_up"].begin(),map_electrons["_smear_up"].end(),PtComparing);
      map_parameter["_smear_up"]=p.Clone(MakeLeptonPointerVector(map_electrons["_smear_up"]));
      
      map_electrons["_smear_down"]=SmearElectrons(map_electrons[""],-1);
      std::sort(map_electrons["_smear_down"].begin(),map_electrons["_smear_down"].end(),PtComparing);
      map_parameter["_smear_down"]=p.Clone(MakeLeptonPointerVector(map_electrons["_smear_down"]));

      map_electrons["_eta2p5"]=ElectronEnergyCorrection(SMPGetElectrons("passMediumID",0.0,2.5),0,0);
      map_parameter["_eta2p5"]=p.Clone(MakeLeptonPointerVector(map_electrons["_eta2p5"]));

      map_parameter["_noroccor"]=p.Clone(MakeLeptonPointerVector(map_electrons["_noroccor"]));
      
      map_electrons["_noEcor"]=ElectronEnergyCorrection(map_electrons[""],-1,0);
      map_parameter["_noEcor"]=p.Clone(MakeLeptonPointerVector(map_electrons["_noEcor"]));

    }
  }else if(channelname.Contains(TRegexp("em20[0-9][0-9]"))){
    Parameter p;
    p.electronIDSF="ID_SF_MediumID_Q";
    p.muonIDSF="IDISO_SF_MediumID_trkIsoLoose_Q";
    p.triggerSF={"",""};
    p.lep0ptcut=25.;
    p.lep1ptcut=15.;

    map_electrons["_noroccor"]=SMPGetElectrons("passMediumID",0.0,2.4);
    map_electrons[""]=ElectronEnergyCorrection(map_electrons["_noroccor"],0,0);
    map_muons[""]=MuonMomentumCorrection(SMPGetMuons("POGMediumWithLooseTrkIso",0.0,2.4),0);

    std::vector<Lepton *> emu = MakeLeptonPointerVector(map_electrons[""]);
    std::vector<Lepton *> mu = MakeLeptonPointerVector(map_muons[""]);
    emu.insert(emu.end(),mu.begin(),mu.end());
    std::sort(emu.begin(),emu.end(),PtComparingPtr);

    map_parameter[""]=p.Clone(emu,
			      (IsNominalRun?NominalWeight:0)
			      +(HasFlag("SYS")&&!IsDATA?SystematicWeight:0)
			      +(HasFlag("PDFSYS")&&!IsDATA?PDFWeight:0)
			      );

  }else{
    cout<<"[AFBAnalyzer::executeEventWithPrefix] wrong channelname"<<endl;
    return;
  }
  
  ///////////////////////lepton selection///////////////////////
  for(const auto& [suffix,p]:map_parameter){
    TString prefix=tauprefix;
    double eventweight=lumiweight*PUweight*prefireweight*z0weight*zptweight;

    if(p.weightbit&NominalWeight) FillHist(channelname+"/"+prefix+"nlepton"+suffix,p.leps.size(),eventweight,10,0,10);
    if(p.leps.size()>=2){
      if(HasFlag("REGION_cf")){
	if(p.leps.at(0)->Charge()>0&&p.leps.at(1)->Charge()>0) prefix="pp_"+prefix;
	else if(p.leps.at(0)->Charge()<0&&p.leps.at(1)->Charge()<0) prefix="mm_"+prefix;
	else return;
      }else{
	if(p.leps.at(0)->Charge()*p.leps.at(1)->Charge()>0) prefix="ss_"+prefix;
      }
      if(channelname.Contains(TRegexp("em20[0-9][0-9]"))){
	if(p.leps.at(0)->LeptonFlavour() == p.leps.at(1)->LeptonFlavour()) return;
      }
      if(p.weightbit&NominalWeight) FillCutflow(channelname+"/"+prefix+"cutflow"+suffix,"dilepton",eventweight);
      if(p.leps.at(0)->Pt()>p.lep0ptcut&&p.leps.at(1)->Pt()>p.lep1ptcut){
	if(p.weightbit&NominalWeight) FillCutflow(channelname+"/"+prefix+"cutflow"+suffix,"ptcut",eventweight);
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

	if(p.weightbit&NominalWeight){
	  FillCutflow(channelname+"/"+tauprefix+"cutflow"+suffix,"RECO",eventweight*RECOSF);
	  FillCutflow(channelname+"/"+tauprefix+"cutflow"+suffix,"ID",eventweight*RECOSF*IDSF);
	  FillCutflow(channelname+"/"+tauprefix+"cutflow"+suffix,"ISO",eventweight*RECOSF*IDSF*ISOSF);
	  FillCutflow(channelname+"/"+tauprefix+"cutflow"+suffix,"trigger",eventweight*RECOSF*IDSF*ISOSF*triggerSF);
	}

	///////////////////////map_weight//////////////////
	map<TString,double> map_weight;
	if(p.weightbit&NominalWeight) map_weight[""]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF;
	if(p.weightbit&SystematicWeight){
	  map_weight["_noefficiencySF"]=lumiweight*PUweight*prefireweight*zptweight*z0weight;
	  
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
	  
	  map_weight["_nozptweight"]=lumiweight*PUweight*prefireweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF;
	  
	  map_weight["_noz0weight"]=lumiweight*PUweight*prefireweight*zptweight*RECOSF*IDSF*ISOSF*triggerSF;
	  
	  map_weight["_noPUweight"]=lumiweight*prefireweight*zptweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF;
	  map_weight["_PUweight_up"]=lumiweight*PUweight_up*prefireweight*zptweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF;
	  map_weight["_PUweight_down"]=lumiweight*PUweight_down*prefireweight*zptweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF;
	  
	  map_weight["_noprefireweight"]=lumiweight*PUweight*zptweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF;
	  map_weight["_prefireweight_up"]=lumiweight*PUweight*prefireweight_up*zptweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF;
	  map_weight["_prefireweight_down"]=lumiweight*PUweight*prefireweight_down*zptweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF;

	}
	if(p.weightbit&PDFWeight){
	  for(unsigned int i=0;i<PDFWeights_Scale->size();i++){
	    map_weight[Form("_scalevariation%d",i)]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF*PDFWeights_Scale->at(i);
	  }
	  for(unsigned int i=0;i<PDFWeights_Error->size();i++){
	    map_weight[Form("_pdf%d",i)]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF*PDFWeights_Error->at(i);
	  }
	  if(PDFWeights_AlphaS->size()==2){
	    map_weight["_alphaS_up"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF*PDFWeights_AlphaS->at(0);
	    map_weight["_alphaS_down"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF*PDFWeights_AlphaS->at(1);
	  }
	}

	///////////////////////fill hists///////////////////////
	if(HasFlag("TOY")) FillHistsToy(channelname,prefix,suffix,(Particle*)p.leps[0],(Particle*)p.leps[1],map_weight);
	else{
	  FillHists(channelname,prefix,suffix,(Particle*)p.leps[0],(Particle*)p.leps[1],map_weight);
	  if(IsDYSample&&prefix==""&&IsNominalRun){
	    vector<Gen> gens=GetGens();
	    Gen truth_l0=GetGenMatchedLepton(*p.leps[0],gens);
	    Gen truth_l1=GetGenMatchedLepton(*p.leps[1],gens);
	    if(!truth_l0.IsEmpty()&&!truth_l1.IsEmpty()) 
		FillHists(channelname,"truth_",suffix,(Particle*)&truth_l0,(Particle*)&truth_l1,map_weight);
	    //else cout<<"no matching"<<endl;
	  }
	}
      }
    }
  }
}
AFBAnalyzer::~AFBAnalyzer(){
  DeleteToy();
  DeleteCosThetaWeight();
}
double AFBAnalyzer::GetCosThetaCS(const Particle *p0,const Particle *p1){
  const TLorentzVector *l0,*l1;
  if(p0->Charge()<0&&p1->Charge()>0){
    l0=p0;
    l1=p1;
  }else if(p0->Charge()>0&&p1->Charge()<0){
    l0=p1;
    l1=p0;
  }else if(strcmp(p0->ClassName(),"LHE")==0){ 
    if(((LHE*)p0)->ID()>0&&((LHE*)p1)->ID()<0){
      l0=p0;
      l1=p1;
    }else if(((LHE*)p0)->ID()<0&&((LHE*)p1)->ID()>0){
      l0=p1;
      l1=p0;
    }else{
      if(gRandom->Rndm()<0.5){
	l0=p0;
	l1=p1;
      }else{
	l0=p1;
	l1=p0;
      }      
    } 
  }else{
    if(gRandom->Rndm()<0.5){
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
    if(gRandom->Rndm()<0.5){
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
  //double systempt=system.Pt();
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
void AFBAnalyzer::FillHistsToy(TString channelname,TString pre,TString suf,Particle* l0,Particle* l1,map<TString,double> map_weight){
  int n_toy=toy_random.size();
  for(int i=0;i<n_toy;i++) FillHists(channelname,pre,suf+Form("_toy%d",i),l0,l1,Multiply(map_weight,toy_weight[i]));
}
void AFBAnalyzer::FillHists(TString channelname,TString pre,TString suf,Particle* l0,Particle* l1,map<TString,double> map_weight){
  TLorentzVector dilepton=(*l0)+(*l1);
  double dimass=dilepton.M();
  double dirap=dilepton.Rapidity();
  double dipt=dilepton.Pt();

  double cost=GetCosThetaCS(l0,l1);
  double h=0.5*pow(dipt/dimass,2)/(1+pow(dipt/dimass,2))*(1-3*cost*cost);
  double den_weight=0.5*fabs(cost)/pow(1+cost*cost+h,2);
  double num_weight=0.5*cost*cost/pow(1+cost*cost+h,3);
  FillHist(channelname+"/"+pre+"costhetaCS"+suf,dimass,dirap,dipt,cost,map_weight,afb_mbinnum,(double*)afb_mbin,afb_ybinnum,(double*)afb_ybin,afb_ptbinnum,(double*)afb_ptbin,20,-1,1);
  FillHist(channelname+"/"+pre+"costhetaCS_den"+suf,dimass,dirap,dipt,cost,Multiply(map_weight,den_weight),afb_mbinnum,(double*)afb_mbin,afb_ybinnum,(double*)afb_ybin,afb_ptbinnum,(double*)afb_ptbin,20,-1,1);
  FillHist(channelname+"/"+pre+"costhetaCS_num"+suf,dimass,dirap,dipt,cost,Multiply(map_weight,num_weight),afb_mbinnum,(double*)afb_mbin,afb_ybinnum,(double*)afb_ybin,afb_ptbinnum,(double*)afb_ptbin,20,-1,1);

  if(!HasFlag("PDFSYS")){
    FillHist(channelname+"/"+pre+"l0pt"+suf,dimass,dirap,dipt,l0->Pt(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,lptbinnum,(double*)lptbin);
    FillHist(channelname+"/"+pre+"l1pt"+suf,dimass,dirap,dipt,l1->Pt(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,lptbinnum,(double*)lptbin);
    FillHist(channelname+"/"+pre+"lpt"+suf,dimass,dirap,dipt,l0->Pt(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,lptbinnum,(double*)lptbin);
    FillHist(channelname+"/"+pre+"lpt"+suf,dimass,dirap,dipt,l1->Pt(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,lptbinnum,(double*)lptbin);
    
    FillHist(channelname+"/"+pre+"l0eta"+suf,dimass,dirap,dipt,l0->Eta(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,60,-3,3);
    FillHist(channelname+"/"+pre+"l1eta"+suf,dimass,dirap,dipt,l1->Eta(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,60,-3,3);
    FillHist(channelname+"/"+pre+"leta"+suf,dimass,dirap,dipt,l0->Eta(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,60,-3,3);
    FillHist(channelname+"/"+pre+"leta"+suf,dimass,dirap,dipt,l1->Eta(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,60,-3,3);

    if(!pre.Contains("gen")&&!pre.Contains("lhe")&&!pre.Contains("truth")){
      FillHist(channelname+"/"+pre+"z0"+suf,dimass,dirap,dipt,vertex_Z,map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,120,-15,15);
      FillHist(channelname+"/"+pre+"met"+suf,dimass,dirap,dipt,pfMET_Type1_pt,map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,100,0,200);  
      FillHist(channelname+"/"+pre+"nPV"+suf,dimass,dirap,dipt,nPV,map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,60,0,60);  
      FillHist(channelname+"/"+pre+"nPileUp"+suf,dimass,dirap,dipt,nPileUp,map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,60,0,60);  
    }
  }
}
void AFBAnalyzer::FillHardHists(TString pre,TString suf,const Gen& genparton0,const Gen& genparton1,const Gen& genhardl0,const Gen& genhardl1,const Gen& genhardj0,double w){
  Gen genhardl=genhardl0.PID()>0?genhardl0:genhardl1;
  TLorentzVector genZ=genhardl0+genhardl1;
  TLorentzVector genpp=genparton0+genparton1;

  FillHist(pre+"cos_l_p0"+suf,cos(genhardl.Angle(genparton0.Vect())),w,100,-1,1);
  FillHist(pre+"cos_l_p0_asym"+suf,cos(genhardl.Angle(genparton0.Vect())),w/2,100,-1,1);
  FillHist(pre+"cos_l_p0_asym"+suf,-1.*cos(genhardl.Angle(genparton0.Vect())),-w/2,100,-1,1);
  if(cos(genhardl.Angle(genparton0.Vect()))>0) FillHist(pre+"cos_l_p0_forward"+suf,genZ.M(),w,fine_mbinnum,(double*)fine_mbin);
  else FillHist(pre+"cos_l_p0_backward"+suf,genZ.M(),w,fine_mbinnum,(double*)fine_mbin);

  FillHist(pre+"cos_l_p1"+suf,cos(genhardl.Angle(genparton1.Vect())),w,100,-1,1);
  FillHist(pre+"cos_l_p1_asym"+suf,cos(genhardl.Angle(genparton1.Vect())),w/2,100,-1,1);
  FillHist(pre+"cos_l_p1_asym"+suf,-1.*cos(genhardl.Angle(genparton1.Vect())),-w/2,100,-1,1);
  if(cos(genhardl.Angle(genparton1.Vect()))>0) FillHist(pre+"cos_l_p1_forward"+suf,genZ.M(),w,fine_mbinnum,(double*)fine_mbin);
  else FillHist(pre+"cos_l_p1_backward"+suf,genZ.M(),w,fine_mbinnum,(double*)fine_mbin);
  
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
    if(cos(genhardl.Angle(genhardj0.Vect()))>0) FillHist(pre+"cos_l_j0_forward"+suf,genZ.M(),w,fine_mbinnum,(double*)fine_mbin);
    else FillHist(pre+"cos_l_j0_backward"+suf,genZ.M(),w,fine_mbinnum,(double*)fine_mbin);

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
void AFBAnalyzer::SetupToy(int n_toy){
  DeleteToy();
  for(int i=0;i<n_toy;i++){
    toy_random.push_back(new TRandom3);
    toy_weight.push_back(-9999.);
  } 
}
void AFBAnalyzer::DeleteToy(){
  int n_toy=toy_random.size();
  for(int i=0;i<n_toy;i++) delete toy_random.at(i);
  toy_random.clear();
  toy_weight.clear();
}
void AFBAnalyzer::GetToyWeight(){
  int n_toy=toy_random.size();
  if(fChain->GetTree()->GetReadEntry()==0){
    TString filename=fChain->GetFile()->GetName();
    for(int i=0;i<n_toy;i++) toy_random[i]->SetSeed((filename+Form("%d",i)).MD5().Hash());
  }
  for(int i=0;i<n_toy;i++)
    toy_weight[i]=toy_random[i]->PoissonD(1.);
}

void AFBAnalyzer::FillHistToy(TString histname, double value, double weight, int n_bin, double x_min, double x_max){
  int n_toy=toy_random.size();
  for(int i=0;i<n_toy;i++) FillHist(histname+Form("_toy%d",i),value,weight*toy_weight[i],n_bin,x_min,x_max);
}
void AFBAnalyzer::FillHistToy(TString histname, double value, map<TString,double> weights, int n_bin, double x_min, double x_max){
  for(const auto& [suffix,weight]:weights) FillHistToy(histname+suffix,value,weight,n_bin,x_min,x_max);
}
void AFBAnalyzer::FillHistToy(TString histname, double value, double weight, int n_bin, double *xbins){ 
  int n_toy=toy_random.size();
  for(int i=0;i<n_toy;i++) FillHist(histname+Form("_toy%d",i),value,weight*toy_weight[i],n_bin,xbins);
}
void AFBAnalyzer::FillHistToy(TString histname, double value, map<TString,double> weights, int n_bin, double *xbins){ 
  for(const auto& [suffix,weight]:weights) FillHistToy(histname+suffix,value,weight,n_bin,xbins);
}
void AFBAnalyzer::FillHistToy(TString histname, double value_x, double value_y, double value_z, double weight, int n_binx, double *xbins, int n_biny, double *ybins, int n_binz, double *zbins){
  int n_toy=toy_random.size();
  for(int i=0;i<n_toy;i++) FillHist(histname+Form("_toy%d",i),value_x,value_y,value_z,weight,n_binx,xbins,n_biny,ybins,n_binz,zbins);
}
void AFBAnalyzer::FillHistToy(TString histname, double value_x, double value_y, double value_z, map<TString,double> weights, int n_binx, double *xbins, int n_biny, double *ybins, int n_binz, double *zbins){
  for(const auto& [suffix,weight]:weights) FillHistToy(histname+suffix,value_x,value_y,value_z,weight,n_binx,xbins,n_biny,ybins,n_binz,zbins);
}
void AFBAnalyzer::SetupCosThetaWeight(){
}
void AFBAnalyzer::DeleteCosThetaWeight(){
}
