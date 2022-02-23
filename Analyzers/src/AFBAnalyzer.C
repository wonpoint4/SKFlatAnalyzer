#include "AFBAnalyzer.h"

void AFBAnalyzer::initializeAnalyzer(){
  SMPAnalyzerCore::initializeAnalyzer(); //setup zpt roc z0 PUJet 
  //SetupCosThetaWeight();
  IsNominalRun=!HasFlag("SYS")&&!HasFlag("PDFSYS");

  vector<JetTagging::Parameters> jtps={JetTagging::Parameters(JetTagging::DeepJet,JetTagging::Tight,JetTagging::incl,JetTagging::comb),
                                       JetTagging::Parameters(JetTagging::DeepJet,JetTagging::Loose,JetTagging::incl,JetTagging::comb)};
  mcCorr->SetJetTaggingParameters(jtps);

  if(fChain->GetListOfFiles()->GetEntries()){
    TString filename=fChain->GetListOfFiles()->At(0)->GetTitle();
    if(filename.Contains("SkimTree_")) IsSkimmed=true;
    else IsSkimmed=false;
  }else{
    cout<<"[AFBAnalyzer::initializeAnalyzer] no input file"<<endl;
    exit(EXIT_FAILURE);
  }
}
void AFBAnalyzer::executeEvent(){
  //// FIXME some events of DYJets has nan PDF weights. I don't know why...
  if(MCSample=="DYJets"&&!isnormal(weight_Scale->at(0))) return;

  ///////////////// GEN level /////////////////////
  executeEventGen();

  ///////////////// RECO level /////////////////////
  if(!IsDATA||DataStream.Contains("SingleMuon")){
    //executeEventWithParameter(MakeParameter("me"));
    //executeEventWithParameter(MakeParameter("mM"));
    //executeEventWithParameter(MakeParameter("mu"));
  }
  if(!IsDATA||DataStream.Contains("DoubleMuon")){
    executeEventWithParameter(MakeParameter("mm"));
    //executeEventWithParameter(MakeParameter("MM"));
  }
  if(!IsDATA||DataStream.Contains("SingleElectron")||DataStream.Contains("EGamma")){
    //executeEventWithParameter(MakeParameter("em"));
    //executeEventWithParameter(MakeParameter("eE"));
    //executeEventWithParameter(MakeParameter("el"));
  }
  if(!IsDATA||DataStream.Contains("DoubleEG")||DataStream.Contains("EGamma")){
    executeEventWithParameter(MakeParameter("ee"));
    //executeEventWithParameter(MakeParameter("EE"));
  }
}
SMPAnalyzerCore::Parameter AFBAnalyzer::MakeParameter(TString key){
  Parameter p=SMPAnalyzerCore::MakeParameter(key);

  p.weightbit=0;
  if(IsNominalRun) p.weightbit|=NominalWeight;
  if(HasFlag("SYS")&&!IsDATA&&(p.channel=="ee"||p.channel=="mm")) p.weightbit|=SystematicWeight|EfficiencyWeight;
  if(HasFlag("PDFSYS")&&!IsDATA&&(p.channel=="ee"||p.channel=="mm")) p.weightbit|=PDFWeight;

  if(HasFlag("nbjet")) p.prefix+="nbjet/";
  else if(HasFlag("0bjet")) p.prefix+="0bjet/";
  if(HasFlag("highmet")) p.prefix+="highmet/";

  return p;
}
bool AFBAnalyzer::PassSelection(Parameter& p){
  if(p.prefix.Contains("highmet")){
    if(pfMET_Type1_pt<60) return false;
    if(IsNominalRun) FillCutflow(p.prefix+p.hprefix+"cutflow","METCut",p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight);
  }

  if(!SMPAnalyzerCore::PassSelection(p)) return false;

  jets=GetJets("tightLepVeto",20,2.4);
  std::sort(jets.begin(),jets.end(),PtComparing);
  p.doublemap["PUJetSF"] = GetPUJetWeight(jets, "Loose", 0);

  realjets.clear();
  for(unsigned int l=0; l<jets.size();l++){
    //jets passing PileUp ID Loose
    if(!PUJetIDPass(jets.at(l), "Loose")) continue;

    //jets spliting from (sub)leading electrons, and muons - jet cleaning
    if(p.lepton0->DeltaR(jets.at(l)) <0.4) continue;
    if(p.lepton1->DeltaR(jets.at(l)) <0.4) continue;

    realjets.push_back(jets.at(l));
  }

  //bjet setting
  n_loosebjet=0;
  n_loose30bjet=0;
  bjets.clear();
  JetTagging::Parameters jtpT = JetTagging::Parameters(JetTagging::DeepJet,JetTagging::Tight,JetTagging::incl,JetTagging::comb);
  JetTagging::Parameters jtpL = JetTagging::Parameters(JetTagging::DeepJet,JetTagging::Loose,JetTagging::incl,JetTagging::comb);

  for(const auto& jet:realjets){
    // For leading b, ID:Tight, pT>30
    if(jet.GetTaggerResult(jtpT.j_Tagger) > mcCorr->GetJetTaggingCutValue(jtpT.j_Tagger, jtpT.j_WP) && jet.Pt()>30){
      bjets.push_back(jet);
    }
    // For 2nd b-veto, ID:Loose, pT>20
    else if(jet.GetTaggerResult(jtpL.j_Tagger) > mcCorr->GetJetTaggingCutValue(jtpL.j_Tagger, jtpL.j_WP)){
      n_loosebjet++;
      if(jet.Pt()>30) n_loose30bjet++;
    }
  }
  p.doublemap["btagSF"] = GetBTaggingReweight_1a_2WP(realjets, jtpT, jtpL, "central");

  //if(IsNominalRun) FillCutflow(p.prefix+p.hprefix+"cutflow","BJetCut",p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight);
  //if(IsNominalRun) FillCutflow(p.prefix+p.hprefix+"cutflow","PUJetSF",p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.doublemap["PUJetSF"]);
  //if(IsNominalRun) FillCutflow(p.prefix+p.hprefix+"cutflow","btagSF",p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.doublemap["PUJetSF"]*p.doublemap["btagSF"]);

  return true;
}
void AFBAnalyzer::executeEventGen(){
  costhetaweight=1.;
  costhetaweight_up=1.;
  costhetaweight_down=1.;
  if(IsDYSample && abs(lhe_l0.ID())!=15){
    //////////////////////// Check LHE /////////////////////////
    Parameter p;
    double letacut=2.4;
    if(abs(lhe_l0.ID())==11){
      p=MakeParameter("ee");
      p.c.lepton0pt=25;
      p.c.lepton1pt=15;
    }else if(abs(lhe_l0.ID())==13){
      p=MakeParameter("mm");
      p.c.lepton0pt=20;
      p.c.lepton1pt=10;
    }else{
      cout<<"[AFBAnalyzer::executeEvent()] something is wrong l0.ID="<<abs(lhe_l0.ID())<<endl;
      vector<LHE> lhes=GetLHEs();
      for(auto& lhe:lhes) lhe.Print();
      exit(EXIT_FAILURE);
    }

    //////////////////////// GEN /////////////////////////
    TLorentzVector gen_Z=gen_l0+gen_l1;
    double gen_Zmass=gen_Z.M();
    double gen_Zrap=gen_Z.Rapidity();
    double gen_Zpt=gen_Z.Pt();
    double gen_cost_correct=-999;
    if(gen_p0.PID()==21){
      if(gen_p1.PID()==21) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,0);
      else if(gen_p1.PID()>0) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,-1);
      else if(gen_p1.PID()<0) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,1);
    }else if(gen_p0.PID()>0){
      if(gen_p1.PID()==21) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,1);
      else if(gen_p1.PID()>0) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,0);
      else if(gen_p1.PID()<0) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,1);
    }else if(gen_p0.PID()<0){
      if(gen_p1.PID()==21) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,-1);
      else if(gen_p1.PID()>0) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,-1);
      else if(gen_p1.PID()<0) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,0);
    }
    if(gen_cost_correct==-999){
      cout<<"wrong pid for parton"<<endl;
      exit(EXIT_FAILURE);
    }
    //costhetaweight=GetCosThetaWeight(gen_Zmass,gen_Zpt,gen_cost_correct,"_pdg");
    //costhetaweight_up=GetCosThetaWeight(gen_Zmass,gen_Zpt,gen_cost_correct,"_up");
    //costhetaweight_down=GetCosThetaWeight(gen_Zmass,gen_Zpt,gen_cost_correct,"_down");

    map<TString,double> map_weight;
    map_weight[""]=p.w.lumiweight;
    //map_weight[""]=p.w.lumiweight*p.w.zptweight*costhetaweight;
    //map_weight["_noweight"]=p.w.lumiweight;
    //map_weight["_nozptweight"]=p.w.lumiweight*costhetaweight;
    //map_weight["_nocosthetaweight"]=p.w.lumiweight*p.w.zptweight;

    /// Matching Study
    /*
    if(bjets.size()>0){
      int nparton=0;
      for(unsigned int i=0; i<gens.size(); i++){
        if(!gens.at(i).isHardProcess()) continue;
        if(abs(gens.at(i).PID())>=11 && abs(gens.at(i).PID())<=16) continue; // No Lepton
        if(gens.at(i).PID()==22 || gens.at(i).PID()==23) continue; // No gamma, Z
        if(gens.at(i).DeltaR(bjets.at(0)) > 0.4) continue; //dR 0.4 matching

        if(nparton==0) gen_j0=gens.at(i);
        else{
          if(gen_j0.PID()+gens.at(i).PID()==0) p.prefix += ""; //"Dyg_";
          else if(gen_j0.PID()==gens.at(i).PID()){ nparton--;}
          else if(gen_j0.PID()==21){ gen_j0=gens.at(i); nparton--;}
          else if(gens.at(i).PID()==21){ nparton--;}
          else gen_j0=(gens.at(i).Pt()>gen_j0.Pt()?gens.at(i):gen_j0);
        }
        nparton++;
      }

      if(nparton==1){
        if(gen_j0.PID()==5) p.prefix += "Dyb_";
        else if(gen_j0.PID()==-5) p.prefix += "Dybbar_";
        else if(gen_j0.PID()==4) p.prefix += "Dyc_";
        else if(gen_j0.PID()==-4) p.prefix += "Dycbar_";
      }
    }
    */

    //////////////// Fill LHE,Gen hists //////////////////////
    if(IsNominalRun){
      TLorentzVector lhe_Z=lhe_l0+lhe_l1;
      double lhe_Zmass=lhe_Z.M();
      double lhe_Zrap=lhe_Z.Rapidity();
      double lhe_Zpt=lhe_Z.Pt();
      FillHistsAFB(p.prefix,"lhe_","",(Particle*)&lhe_l0,(Particle*)&lhe_l1,map_weight);
      FillHist(p.prefix+"lhe_costhetaCS",lhe_Zmass,lhe_Zrap,lhe_Zpt,GetCosThetaCS(&lhe_l0,&lhe_l1,0),map_weight,afb_mbinnum,(double*)afb_mbin,afb_ybinnum,(double*)afb_ybin,afb_ptbinnum,(double*)afb_ptbin,20,-1,1);
      if(lhe_j0.Pt()){
	FillHist(p.prefix+"lhe_costhetaR",lhe_Zmass,lhe_Zrap,lhe_Zpt,GetCosThetaR(&lhe_l0,&lhe_l1,&lhe_j0,0),map_weight,afb_mbinnum,(double*)afb_mbin,afb_ybinnum,(double*)afb_ybin,afb_ptbinnum,(double*)afb_ptbin,20,-1,1);
	FillHist(p.prefix+"lhe_costhetaT",lhe_Zmass,lhe_Zrap,lhe_Zpt,GetCosThetaT(&lhe_l0,&lhe_l1,&lhe_j0,0),map_weight,afb_mbinnum,(double*)afb_mbin,afb_ybinnum,(double*)afb_ybin,afb_ptbinnum,(double*)afb_ptbin,20,-1,1);
      }
      if(lhe_l0.Pt()>p.c.lepton0pt&&lhe_l1.Pt()>p.c.lepton1pt&&fabs(lhe_l0.Eta())<letacut&&fabs(lhe_l1.Eta())<letacut){
	FillHistsAFB(p.prefix,"lhefid_","",(Particle*)&lhe_l0,(Particle*)&lhe_l1,map_weight);
      }
      FillHistsAFB(p.prefix,"gen_","",(Particle*)&gen_l0,(Particle*)&gen_l1,map_weight);
      FillHistsAFB(p.prefix,"gen_","_dressed",(Particle*)&gen_l0_dressed,(Particle*)&gen_l1_dressed,map_weight);
      FillHistsAFB(p.prefix,"gen_","_bare",(Particle*)&gen_l0_bare,(Particle*)&gen_l1_bare,map_weight);
      if(gen_l0.Pt()>p.c.lepton0pt&&gen_l1.Pt()>p.c.lepton1pt&&fabs(gen_l0.Eta())<letacut&&fabs(gen_l1.Eta())<letacut){
        FillHistsAFB(p.prefix,"genfid_","",(Particle*)&gen_l0,(Particle*)&gen_l1,map_weight);
      }
      if(gen_l0_dressed.Pt()>p.c.lepton0pt&&gen_l1_dressed.Pt()>p.c.lepton1pt&&fabs(gen_l0_dressed.Eta())<letacut&&fabs(gen_l1_dressed.Eta())<letacut){
        FillHistsAFB(p.prefix,"genfid_","_dressed",(Particle*)&gen_l0_dressed,(Particle*)&gen_l1_dressed,map_weight);
      }
      if(gen_l0_bare.Pt()>p.c.lepton0pt&&gen_l1_bare.Pt()>p.c.lepton1pt&&fabs(gen_l0_bare.Eta())<letacut&&fabs(gen_l1_bare.Eta())<letacut){
        FillHistsAFB(p.prefix,"genfid_","_bare",(Particle*)&gen_l0_bare,(Particle*)&gen_l1_bare,map_weight);
      }
      FillHist(p.prefix+"gen_costhetaCS_correct",gen_Zmass,gen_Zrap,gen_Zpt,gen_cost_correct,map_weight,afb_mbinnum,(double*)afb_mbin,afb_ybinnum,(double*)afb_ybin,afb_ptbinnum,(double*)afb_ptbin,20,-1,1);
      FillHist(p.prefix+"gen_nPU_noPUweight",gen_Zmass,gen_Zrap,gen_Zpt,nPileUp,p.w.lumiweight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,100,0,100);
      FillHist(p.prefix+"gen_nPU",gen_Zmass,gen_Zrap,gen_Zpt,nPileUp,p.w.lumiweight*p.w.PUweight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,100,0,100);
      FillHist(p.prefix+"gen_nPU_PUweight_up",gen_Zmass,gen_Zrap,gen_Zpt,nPileUp,p.w.lumiweight*p.w.PUweight_up,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,100,0,100);
      FillHist(p.prefix+"gen_nPU_PUweight_down",gen_Zmass,gen_Zrap,gen_Zpt,nPileUp,p.w.lumiweight*p.w.PUweight_down,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,100,0,100);
    }
  }
}

void AFBAnalyzer::FillHists(Parameter& p){
  ///////////////////////map_weight//////////////////
  map<TString,double> map_weight;
  if(p.weightbit&NominalWeight){
    map_weight[""]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["PUJetSF"]*p.doublemap["btagSF"];
    map_weight["_nbtag"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["PUJetSF"];
    map_weight["_nopujet"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"];

    if(MCSample.Contains("MiNNLO")){
      map_weight["_sthw3"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["PUJetSF"]*p.doublemap["btagSF"]*weight_sthw2->at(3);
      map_weight["_sthw4"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["PUJetSF"]*p.doublemap["btagSF"]*weight_sthw2->at(4);
      map_weight["_sthw10"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["PUJetSF"]*p.doublemap["btagSF"]*weight_sthw2->at(10);
    }
  }
  if(p.weightbit&SystematicWeight){
    if(!IsDATA){
      //map_weight["_noPUweight"]=p.w.lumiweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"];
      map_weight["_PUweight_up"]=p.w.lumiweight*p.w.PUweight_up*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"];
      map_weight["_PUweight_down"]=p.w.lumiweight*p.w.PUweight_down*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"];
      
      map_weight["_noprefireweight"]=p.w.lumiweight*p.w.PUweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"];
      map_weight["_prefireweight_up"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight_up*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"];
      map_weight["_prefireweight_down"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight_down*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"];
      
      map_weight["_nozptweight"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"];
      
      map_weight["_nobtagSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;

      //map_weight["_noz0weight"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"];
      
      //map_weight["_nocosthetaweight"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"];
      //map_weight["_costhetaweight_up"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight_up*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"];
      //map_weight["_costhetaweight_down"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight_down*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"];
      
      //map_weight["_noefficiencySF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.CFSF*p.doublemap["btagSF"];
      
      //map_weight["_noelectronRECOSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"];
      //map_weight["_noIDSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.electronRECOSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"];
      //map_weight["_nomuonISOSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"];      
      //map_weight["_notriggerSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.CFSF*p.doublemap["btagSF"];
      for(int j=0,nj=fEff->nreplica;j<nj;j++){
	double electronRECOSF=p.w.electronRECOSF_sys.size() ? p.w.electronRECOSF_sys[0][j] : 1.;
	double electronIDSF=p.w.electronIDSF_sys.size() ? p.w.electronIDSF_sys[0][j] : 1.;
	double muonIDSF=p.w.muonIDSF_sys.size() ? p.w.muonIDSF_sys[0][j] : 1.;
	double triggerSF=p.w.triggerSF_sys.size() ? p.w.triggerSF_sys[0][j] : 1.;
	map_weight[Form("_efficiencySF_stat%d",j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*electronRECOSF*electronIDSF*muonIDSF*p.w.muonISOSF*triggerSF*p.w.CFSF*p.doublemap["btagSF"];
      }

      for(int i=1,ni=p.w.electronRECOSF_sys.size();i<ni;i++){
	for(int j=0,nj=p.w.electronRECOSF_sys[i].size();j<nj;j++){
	  map_weight[Form("_electronRECOSF_s%d_m%d",i,j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.electronRECOSF_sys[i][j]*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"];
	}
      }

      for(int i=1,ni=p.w.electronIDSF_sys.size();i<ni;i++){
	for(int j=0,nj=p.w.electronIDSF_sys[i].size();j<nj;j++){
	  map_weight[Form("_electronIDSF_s%d_m%d",i,j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF_sys[i][j]*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"];
	}
      }

      for(int i=1,ni=p.w.muonIDSF_sys.size();i<ni;i++){
	for(int j=0,nj=p.w.muonIDSF_sys[i].size();j<nj;j++){
	  map_weight[Form("_muonIDSF_s%d_m%d",i,j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF_sys[i][j]*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"];
	}
      }

      for(int i=1,ni=p.w.triggerSF_sys.size();i<ni;i++){
	for(int j=0,nj=p.w.triggerSF_sys[i].size();j<nj;j++){
	  map_weight[Form("_triggerSF_s%d_m%d",i,j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF_sys[i][j]*p.w.CFSF*p.doublemap["btagSF"];
	}
      }

      //map_weight["_noCFSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"];
      map_weight["_CFSF_up"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF_up*p.doublemap["btagSF"];
      map_weight["_CFSF_down"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF_down*p.doublemap["btagSF"];
    }
  }
  if(p.weightbit&PDFWeight){
    for(unsigned int i=0;i<weight_Scale->size();i++){
      map_weight[Form("_scalevariation%d",i)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"]*weight_Scale->at(i);
    }
    for(unsigned int i=0;i<weight_PDF->size();i++){
      map_weight[Form("_pdf%d",i)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"]*weight_PDF->at(i);
    }
    if(weight_AlphaS->size()==2){
      map_weight["_alphaS_down"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"]*weight_AlphaS->at(0);
      map_weight["_alphaS_up"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"]*weight_AlphaS->at(1);
    }

    if(MCSample.Contains("MiNNLO")){
      //map_weight["_sthw2_down"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"]*weight_sthw2->at(0);
      //map_weight["_sthw2_up"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"]*weight_sthw2->at(2);
      map_weight["_largeptscales"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"]*weight_largeptscales->at(0);
      map_weight["_q0_up"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"]*weight_q0->at(0);
      map_weight["_q0_down"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"]*weight_q0->at(2);      
    }
  }

  ///////////////////////DY+b Specific Event Selection///////////////////////
  double eventweight = p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight;

  eventweight *= p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
  double nopujetweight = eventweight;
  double nobtagweight = eventweight;

  eventweight *= p.doublemap["PUJetSF"];
  nobtagweight *= p.doublemap["PUJetSF"];
  FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"PUJetSF",eventweight);
  eventweight *= p.doublemap["btagSF"];;
  nopujetweight *= p.doublemap["btagSF"];
  FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"btagSF",eventweight);

  TLorentzVector *l0,*l1;
  l0=p.lepton0;
  l1=p.lepton1;

  if(IsDYSample && abs(lhe_l0.ID())!=15){
    if(bjets.size()>0){
      int nparton=0;
      for(unsigned int i=0; i<gens.size(); i++){
        if(!gens.at(i).isHardProcess()) continue;
        if(abs(gens.at(i).PID())>=11 && abs(gens.at(i).PID())<=16) continue; // No Lepton
        if(gens.at(i).PID()==22 || gens.at(i).PID()==23) continue; // No gamma, Z
        if(gens.at(i).DeltaR(bjets.at(0)) > 0.4) continue; //dR 0.4 matching

        if(nparton==0) gen_j0=gens.at(i);
        else{
          if(gen_j0.PID()+gens.at(i).PID()==0) p.prefix += ""; //"Dyg_";
          else if(gen_j0.PID()==gens.at(i).PID()){ nparton--;}
          else if(gen_j0.PID()==21){ gen_j0=gens.at(i); nparton--;}
          else if(gens.at(i).PID()==21){ nparton--;}
          else gen_j0=(gens.at(i).Pt()>gen_j0.Pt()?gens.at(i):gen_j0);
        }
        nparton++;
      }

      /*
      if(nparton==0) p.prefix += "DyNo_";
      else if(nparton==1){
        if(gen_j0.PID()==5) p.prefix += "Dyb_";
        else if(gen_j0.PID()==-5) p.prefix += "Dybbar_";
        else if(gen_j0.PID()==4) p.prefix += "Dyc_";
        else if(gen_j0.PID()==-4) p.prefix += "Dycbar_";
        else if(gen_j0.PID()<4) p.prefix += "Dyuds_";
        else p.prefix += "Dyg_";
      */
      if(nparton==0) FillHist(p.prefix+p.hprefix+"MatchedParton_PID"+p.suffix,0,eventweight,20,-10,10);
      else if(nparton==1) FillHist(p.prefix+p.hprefix+"MatchedParton_PID"+p.suffix,gen_j0.PID(),eventweight,20,-10,10);
      else if(nparton==2) FillHist(p.prefix+p.hprefix+"MatchedParton_PID"+p.suffix,8,eventweight,20,-10,10);
      else FillHist(p.prefix+p.hprefix+"MatchedParton_PID"+p.suffix,9,eventweight,20,-10,10);
    }
  }

  if((*l0+*l1).M() <52 ) return;
  if(p.weightbit&NominalWeight) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"mass52Cut",eventweight);

  FillHist(p.prefix+p.hprefix+"OneTightb_bef_mll"+p.suffix,               (*l0+*l1).M(),eventweight,100,50,150);
  FillHist(p.prefix+p.hprefix+"OneTightb_bef_yll"+p.suffix,               (*l0+*l1).Rapidity(),eventweight,60,-3,3);
  FillHist(p.prefix+p.hprefix+"OneTightb_bef_pTll"+p.suffix,              (*l0+*l1).Pt(),eventweight,100,0,100);
  FillHist(p.prefix+p.hprefix+"OneTightb_bef_bjet"+p.suffix,              bjets.size()+n_loosebjet,eventweight,5,0,5);
  FillHist(p.prefix+p.hprefix+"OneTightb_bef_tightbjet"+p.suffix,         bjets.size(),eventweight,5,0,5);
  FillHist(p.prefix+p.hprefix+"OneTightb_bef_loosebjet"+p.suffix,         n_loosebjet,eventweight,5,0,5);
  FillHist(p.prefix+p.hprefix+"OneTightb_bef_loose30bjet"+p.suffix,       n_loose30bjet,eventweight,5,0,5);
  FillHist(p.prefix+p.hprefix+"OneTightb_bef_jet"+p.suffix,               realjets.size(),eventweight,10,0,10);
  FillHist(p.prefix+p.hprefix+"OneTightb_bef_normjet"+p.suffix,           jets.size(),eventweight,10,0,10);
  FillHist(p.prefix+p.hprefix+"OneTightb_bef_nobtagSF"+p.suffix,          bjets.size(),nobtagweight,5,0,5);
  FillHist(p.prefix+p.hprefix+"OneTightb_bef_jet_nobtagSF"+p.suffix,      realjets.size(),nobtagweight,10,0,10);
  FillHist(p.prefix+p.hprefix+"OneTightb_bef_normjet_nobtagSF"+p.suffix,  jets.size(),nobtagweight,10,0,10);
  FillHist(p.prefix+p.hprefix+"OneTightb_bef_noPUJetSF"+p.suffix,         bjets.size(),nopujetweight,5,0,5);
  FillHist(p.prefix+p.hprefix+"OneTightb_bef_jet_noPUJetSF"+p.suffix,     realjets.size(),nopujetweight,10,0,10);
  FillHist(p.prefix+p.hprefix+"OneTightb_bef_normjet_noPUJetSF"+p.suffix, jets.size(),nopujetweight,10,0,10);
  if(bjets.size() !=1) return;
  if(p.weightbit&NominalWeight) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"oneTightb",eventweight);
  bjet = bjets.at(0);
  bjet_charge = bjetCharge(bjets.at(0),0,p.prefix+p.hprefix);

  FillHist(p.prefix+p.hprefix+"2bveto_bef_mll"+p.suffix,        (*l0+*l1).M(),eventweight,100,50,150);
  FillHist(p.prefix+p.hprefix+"2bveto_bef_yll"+p.suffix,        (*l0+*l1).Rapidity(),eventweight,60,-3,3);
  FillHist(p.prefix+p.hprefix+"2bveto_bef_pTll"+p.suffix,       (*l0+*l1).Pt(),eventweight,100,0,100);
  FillHist(p.prefix+p.hprefix+"2bveto_bef"+p.suffix,            n_loosebjet,eventweight,5,0,5);
  FillHist(p.prefix+p.hprefix+"2bveto_bef_loose30bjet"+p.suffix,n_loose30bjet,eventweight,5,0,5);
  if(n_loosebjet !=0) return;
  if(p.weightbit&NominalWeight) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"2bveto",eventweight);

  FillHist(p.prefix+p.hprefix+"2jveto_bef_mll"+p.suffix,     (*l0+*l1).M(),eventweight,100,50,150);
  FillHist(p.prefix+p.hprefix+"2jveto_bef_yll"+p.suffix,     (*l0+*l1).Rapidity(),eventweight,60,-3,3);
  FillHist(p.prefix+p.hprefix+"2jveto_bef_pTll"+p.suffix,    (*l0+*l1).Pt(),eventweight,100,0,100);
  FillHist(p.prefix+p.hprefix+"2jveto_bef"+p.suffix,         realjets.size(),eventweight,7,0,7);
  FillHist(p.prefix+p.hprefix+"2jveto_bef_normjet"+p.suffix, jets.size(),eventweight,10,0,10);
  int n_powerjet = 0;
  for(unsigned int k=0; k<realjets.size(); k++){
    if(realjets.at(k).Pt() > 30) n_powerjet += 1;
  }
  FillHist(p.prefix+p.hprefix+"2jveto_bef_powerjet"+p.suffix,n_powerjet,eventweight,7,0,7);
  //if(realjets.size() >1) return;
  if(n_powerjet >1) return;
  if(p.weightbit&NominalWeight) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"2jveto",eventweight);

  FillHist(p.prefix+p.hprefix+"MET90Cut_bef_mll"+p.suffix,  (*l0+*l1).M(),eventweight,100,50,150);
  FillHist(p.prefix+p.hprefix+"MET90Cut_bef_yll"+p.suffix,  (*l0+*l1).Rapidity(),eventweight,60,-3,3);
  FillHist(p.prefix+p.hprefix+"MET90Cut_bef_pTll"+p.suffix, (*l0+*l1).Pt(),eventweight,100,0,100);
  FillHist(p.prefix+p.hprefix+"MET90Cut_bef"+p.suffix,      pfMET_Type1_pt,eventweight,150,0,150);
  if(pfMET_Type1_pt >90) return;
  if(p.weightbit&NominalWeight) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"MET90Cut",eventweight);

  FillHist(p.prefix+p.hprefix+"ZbdPhiCut_bef"+p.suffix,(*l0+*l1).DeltaPhi(bjets.at(0)),eventweight,140,-3.5,3.5);
  if(abs((*l0+*l1).DeltaPhi(bjets.at(0))) <2) return;
  if(p.weightbit&NominalWeight) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"ZbdPhiCut",eventweight);

  FillHist(p.prefix+p.hprefix+"ZbpT50Cut_bef"+p.suffix,(*l0+*l1+bjets.at(0)).Pt(),eventweight,200,0,200);
  if((*l0+*l1+bjets.at(0)).Pt() >50) return;
  if(p.weightbit&NominalWeight) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"ZbpT50Cut",eventweight);

  FillHist(p.prefix+p.hprefix+"ZpT20Cut_bef"+p.suffix,(*l0+*l1).Pt(),eventweight,200,0,200);
  if((*l0+*l1).Pt() <20) return;
  if(p.weightbit&NominalWeight) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"ZpT20Cut",eventweight);

  bjet_charge = bjetCharge(bjets.at(0), 1, p.prefix+p.hprefix, eventweight);

  if(fabs(bjet_charge) > 1.) FillHist(p.prefix+p.hprefix+"bjetCharge_raw_SL",bjet_charge,eventweight,600,-6,6);
  else FillHist(p.prefix+p.hprefix+"bjetCharge_raw_noSL",bjet_charge,eventweight,600,-6,6);
  if(fabs(bjet_charge) < 0.2) return;
  if(p.weightbit&NominalWeight) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"bjetcharge020Cut",eventweight);

  ///////////////////////fill hists///////////////////////
  FillHist(p.prefix+p.hprefix+"bjetpT",bjets.at(0).Pt(),map_weight,200,0,1000);
  FillHist(p.prefix+p.hprefix+"bjeteta",bjets.at(0).Eta(),map_weight,100,-5,5);
  FillHist(p.prefix+p.hprefix+"bjetM",bjets.at(0).M(),map_weight,200,0,20);
  FillHist(p.prefix+p.hprefix+"bjetCharge",bjet_charge,map_weight,600,-6,6);
  FillHist(p.prefix+p.hprefix+"bjetPUID",bjets.at(0).PileupJetId(),map_weight,200,-2,2);
  //FillHist(p.prefix+p.hprefix+"bjet_DeepJet",bjets.at(0).DeepJet(),map_weight,100,0,1);
  //FillHist(p.prefix+p.hprefix+"bjet_DeepCSV",bjets.at(0).DeepCSV(),map_weight,100,0,1);
  //FillHist(p.prefix+p.hprefix+"bjet_DeepJet_b",bjets.at(0).DeepJet_b(),map_weight,100,0,1);
  //FillHist(p.prefix+p.hprefix+"bjet_DeepJet_bb",bjets.at(0).DeepJet_bb(),map_weight,100,0,1);
  //FillHist(p.prefix+p.hprefix+"bjet_DeepJet_lepb",bjets.at(0).DeepJet_lepb(),map_weight,100,0,1);
  //FillHist(p.prefix+p.hprefix+"bjet_DeepJet_bonlyratio",(bjets.at(0).DeepJet_b()+bjets.at(0).DeepJet_lepb())/bjets.at(0).DeepJet(),map_weight,100,0,1);

  FillHist(p.prefix+p.hprefix+"yZ",(*l0+*l1).Rapidity(),map_weight,60,-3,3);
  FillHist(p.prefix+p.hprefix+"yb",bjets.at(0).Rapidity(),map_weight,60,-3,3);
  FillHist(p.prefix+p.hprefix+"Zb_y",(*l0+*l1+bjets.at(0)).Rapidity(),map_weight,100,-5,5);
  FillHist(p.prefix+p.hprefix+"Zb_pT",(*l0+*l1+bjets.at(0)).Pt(),map_weight,100,0,100);
  FillHist(p.prefix+p.hprefix+"Zb_dy",(*l0+*l1).Rapidity()-bjets.at(0).Rapidity(),map_weight,100,-5,5);
  FillHist(p.prefix+p.hprefix+"Zb_dphi",(*l0+*l1).DeltaPhi(bjets.at(0)),map_weight,100,-5,5);
  FillHist(p.prefix+p.hprefix+"Zb_y2D",(*l0+*l1).Rapidity(),bjets.at(0).Rapidity(),map_weight,30,-3,3,30,-3,3);

  FillHist(p.prefix+p.hprefix+"mll",(*l0+*l1).M(),map_weight,100,50,150);
  FillHist(p.prefix+p.hprefix+"yll",(*l0+*l1).Rapidity(),map_weight,60,-3,3);
  FillHist(p.prefix+p.hprefix+"pTll",(*l0+*l1).Pt(),map_weight,100,0,100);

  TLorentzVector dilepton=*p.lepton0+*p.lepton1;
  double dimass=dilepton.M();
  double dirap=dilepton.Rapidity();
  double dipt=dilepton.Pt();
  FillHistsAFB(p.prefix,p.hprefix,p.suffix,(Particle*)p.lepton0,(Particle*)p.lepton1,map_weight);
  FillHist(p.prefix+p.hprefix+"nbjet"+p.suffix,dimass,dirap,dipt,p.intmap["nbjet"],map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,10,0,10);
  FillHist(p.prefix+p.hprefix+"z0"+p.suffix,dimass,dirap,dipt,vertex_Z,SelectWeights(map_weight,{""}),grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,120,-15,15);
  map<TString,double> map_PUweight=SelectWeights(map_weight,{"","_noPUweight","_PUweight_up","_PUweight_down"});
  FillHist(p.prefix+p.hprefix+"nPV"+p.suffix,dimass,dirap,dipt,nPV,map_PUweight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,100,0,100);
  FillHist(p.prefix+p.hprefix+"rho"+p.suffix,dimass,dirap,dipt,Rho,map_PUweight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,50,0,50);
  FillHist(p.prefix+p.hprefix+"met"+p.suffix,dimass,dirap,dipt,pfMET_Type1_pt,map_PUweight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,100,0,200);
  if(IsDYSample&&p.hprefix==""&&IsNominalRun){
    vector<Gen> gens=GetGens();
    Gen truth_l0=GetGenMatchedLepton(*p.lepton0,gens);
    Gen truth_l1=GetGenMatchedLepton(*p.lepton1,gens);
    if(!truth_l0.IsEmpty()&&!truth_l1.IsEmpty())  FillHistsAFB(p.prefix,"truth_",p.suffix,(Particle*)&truth_l0,(Particle*)&truth_l1,map_weight);
    //else cout<<"no matching"<<endl;
  }
}

AFBAnalyzer::AFBAnalyzer(){}
AFBAnalyzer::~AFBAnalyzer(){
  //DeleteCosThetaWeight();
}
double AFBAnalyzer::GetCosThetaCS(const Particle *p0,const Particle *p1,int direction){
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
      gRandom->SetSeed((run<<15)+(lumi<<10)+(event<<5)+p0->Eta()*100);
      if(gRandom->Rndm()<0.5){
        l0=p0;
        l1=p1;
      }else{
        l0=p1;
        l1=p0;
      }      
    } 
  }else{
    gRandom->SetSeed((run<<15)+(lumi<<10)+(event<<5)+p0->Eta()*100);
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
  if(bjet*bjet==0){
    if(direction==0) direction=(dilepton.Pz()>0?1:-1);
  }else{
    if(direction==0) direction=((0.75*dilepton.Rapidity()-bjet.Rapidity())>0?1:-1);
    if(bjet_charge>0) direction *= -1.;
  }
  return direction*2*(l0pp*l1pm-l0pm*l1pp)/sqrt(dimass*dimass*(dimass*dimass+dipt*dipt));
}

double AFBAnalyzer::GetCosThetaR(const Particle *l0,const Particle *l1,const Particle *j0,int direction){
  const Particle *lm=NULL,*lp=NULL;
  if(l0->Charge()<0&&l1->Charge()>0){
    lm=l0; lp=l1;
  }else if(l0->Charge()>0&&l1->Charge()<0){
    lm=l1; lp=l0;
  }else{
    gRandom->SetSeed((run<<15)+(lumi<<10)+(event<<5)+l0->Eta()*100);
    if(gRandom->Rndm()<0.5){
      lm=l0; lp=l1;
    }else{
      lm=l1; lp=l0;
    }
  }
  if(j0->E()){
    int jid=0;
    if(j0->InheritsFrom("LHE")) jid=((LHE*)j0)->ID();
    if(jid==22){
      TLorentzVector dilepton=*lm+*lp;
      TLorentzVector jet=*j0;
      TLorentzVector system=dilepton+jet;
      TLorentzVector p0(0,0,0.5*system.M()*exp(system.Rapidity()),0.5*system.M()*exp(system.Rapidity()));
      TLorentzVector p1(0,0,-0.5*system.M()*exp(-system.Rapidity()),0.5*system.M()*exp(-system.Rapidity()));
      TLorentzVector lepton=*lm;
      TVector3 b1=system.BoostVector();
      dilepton.Boost(-b1);lepton.Boost(-b1);
      p0.Boost(-b1);p1.Boost(-b1);jet.Boost(-b1);
      if(p0.Angle(jet.Vect())<p1.Angle(jet.Vect())){
        p0-=jet;
      }else{
        p1-=jet;
      }
      TVector3 b2=dilepton.BoostVector();
      p0.Boost(-b2);p1.Boost(-b2);lepton.Boost(-b2);
      if(direction==0) direction=system.Pz()/fabs(system.Pz());
      return direction*cos(lepton.Angle(p0.Vect().Unit()-p1.Vect().Unit()));
    }else if(0<jid&&jid<=6){
      return ((*lm-*lp)*(*j0))/((*lm+*lp)*(*j0));
    }else if(-6<=jid&&jid<0){
      return -1*((*lm-*lp)*(*j0))/((*lm+*lp)*(*j0));
    }
  }
  return GetCosThetaCS(l0,l1,direction);
}
double AFBAnalyzer::GetCosThetaT(const Particle *l0,const Particle *l1,const Particle *j0,int direction){
  const Particle *lm=NULL,*lp=NULL;
  if(l0->Charge()<0&&l1->Charge()>0){
    lm=l0; lp=l1;
  }else if(l0->Charge()>0&&l1->Charge()<0){
    lm=l1; lp=l0;
  }else{
    gRandom->SetSeed((run<<15)+(lumi<<10)+(event<<5)+l0->Eta()*100);
    if(gRandom->Rndm()<0.5){
      lm=l0; lp=l1;
    }else{
      lm=l1; lp=l0;
    }
  }
  if(j0->E()){
    int jid=0;
    if(j0->InheritsFrom("LHE")) jid=((LHE*)j0)->ID();
    TLorentzVector dilepton=*lm+*lp;
    TLorentzVector jet=*j0;
    TLorentzVector system=dilepton+jet;
    TLorentzVector p0(0,0,0.5*system.M()*exp(system.Rapidity()),0.5*system.M()*exp(system.Rapidity()));
    TLorentzVector p1(0,0,-0.5*system.M()*exp(-system.Rapidity()),0.5*system.M()*exp(-system.Rapidity()));
    TLorentzVector lepton=*lm;
    TVector3 b1=system.BoostVector();
    dilepton.Boost(-b1);lepton.Boost(-b1);
    p0.Boost(-b1);p1.Boost(-b1);jet.Boost(-b1);
    int temp_direction=direction;
    if(p0.Angle(jet.Vect())<p1.Angle(jet.Vect())){
      p0-=jet;
      if(temp_direction==0) direction=-1;
    }else{
      p1-=jet;
      if(temp_direction==0) direction=+1;
    }
    TVector3 b2=dilepton.BoostVector();
    p0.Boost(-b2);p1.Boost(-b2);lepton.Boost(-b2);
    if(direction==0){
      if(jid==22){
        direction=system.Pz()/fabs(system.Pz());
      }else if(0<jid&&jid<=6){
        direction=temp_direction;
      }else if(-6<=jid&&jid<0){
        direction=-temp_direction;
      }
      return direction*cos(lepton.Angle(p0.Vect().Unit()-p1.Vect().Unit()));
    }
  }
  return GetCosThetaCS(l0,l1,direction);
}

double AFBAnalyzer::GetCosThetaRecoil(const Particle *p0,const Particle *p1,int direction){
  const TLorentzVector *lm,*lp;
  if(p0->Charge()<0&&p1->Charge()>0){
    lm=p0;
    lp=p1;
  }else if(p0->Charge()>0&&p1->Charge()<0){
    lm=p1;
    lp=p0;
  }else if(strcmp(p0->ClassName(),"LHE")==0){ 
    if(((LHE*)p0)->ID()>0&&((LHE*)p1)->ID()<0){
      lm=p0;
      lp=p1;
    }else if(((LHE*)p0)->ID()<0&&((LHE*)p1)->ID()>0){
      lm=p1;
      lp=p0;
    }else{
      if(gRandom->Rndm()<0.5){
        lm=p0;
        lp=p1;
      }else{
      	lm=p1;
      	lp=p0;
      }      
    } 
  }else{
    if(gRandom->Rndm()<0.5){
      lm=p0;
      lp=p1;
    }else{
      lm=p1;
      lp=p0;
    }      
  }
  if(bjet_charge>0) direction *= -1.;
  return direction*((*lm-*lp)*bjet)/((*lm+*lp)*bjet);
}

void AFBAnalyzer::FillHistsAFB(TString pre,TString hpre,TString suf,Particle* l0,Particle* l1,map<TString,double> map_weight){
  TLorentzVector dilepton=(*l0)+(*l1);
  double dimass=dilepton.M();
  double dirap=dilepton.Rapidity();
  double dipt=dilepton.Pt();

  double cost=GetCosThetaCS(l0,l1);
  FillHist(pre+hpre+"costhetaCS"+suf,dimass,dirap,dipt,cost,map_weight,afb_mbinnum,(double*)afb_mbin,afb_ybinnum,(double*)afb_ybin,afb_ptbinnum,(double*)afb_ptbin,20,-1,1);
  //double h=0.5*pow(dipt/dimass,2)/(1+pow(dipt/dimass,2))*(1-3*cost*cost);
  //double den_weight=0.5*fabs(cost)/pow(1+cost*cost+h,2);
  //double num_weight=0.5*cost*cost/pow(1+cost*cost+h,3);
  //FillHist(pre+hpre+"costhetaCS_den"+suf,dimass,dirap,dipt,cost,Multiply(map_weight,den_weight),afb_mbinnum,(double*)afb_mbin,afb_ybinnum,(double*)afb_ybin,afb_ptbinnum,(double*)afb_ptbin,20,-1,1);
  //FillHist(pre+hpre+"costhetaCS_num"+suf,dimass,dirap,dipt,cost,Multiply(map_weight,num_weight),afb_mbinnum,(double*)afb_mbin,afb_ybinnum,(double*)afb_ybin,afb_ptbinnum,(double*)afb_ptbin,20,-1,1);
  if(bjet*bjet!=0){
    double cosr=GetCosThetaRecoil(l0,l1);
    FillHist(pre+hpre+"costhetaRecoil"+suf,dimass,dirap,dipt,cosr,map_weight,afb_mbinnum,(double*)afb_mbin,afb_ybinnum,(double*)afb_ybin,afb_ptbinnum,(double*)afb_ptbin,20,-1,1);
  }

  FillHist(pre+hpre+"l0pt"+suf,dimass,dirap,dipt,l0->Pt(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,lptbinnum,(double*)lptbin);
  FillHist(pre+hpre+"l1pt"+suf,dimass,dirap,dipt,l1->Pt(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,lptbinnum,(double*)lptbin);
  FillHist(pre+hpre+"lpt"+suf,dimass,dirap,dipt,l0->Pt(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,lptbinnum,(double*)lptbin);
  FillHist(pre+hpre+"lpt"+suf,dimass,dirap,dipt,l1->Pt(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,lptbinnum,(double*)lptbin);
  
  FillHist(pre+hpre+"l0eta"+suf,dimass,dirap,dipt,l0->Eta(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,60,-3,3);
  FillHist(pre+hpre+"l1eta"+suf,dimass,dirap,dipt,l1->Eta(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,60,-3,3);
  FillHist(pre+hpre+"leta"+suf,dimass,dirap,dipt,l0->Eta(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,60,-3,3);
  FillHist(pre+hpre+"leta"+suf,dimass,dirap,dipt,l1->Eta(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,60,-3,3);

  if(!hpre.Contains("gen")&&!hpre.Contains("lhe")&&!hpre.Contains("truth")){
    FillHist(pre+hpre+"z0"+suf,dimass,dirap,dipt,vertex_Z,map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,120,-15,15);
    FillHist(pre+hpre+"met"+suf,dimass,dirap,dipt,pfMET_Type1_pt,map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,100,0,200);  
    FillHist(pre+hpre+"bjetCh"+suf,dimass,dirap,dipt,bjet_charge,map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,20,-5,5);
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
void AFBAnalyzer::SetupCosThetaWeight(){
  cout<<"[AFBAnalyzer::SetupCosThetaWeight] Setup"<<endl;
  TString datapath=getenv("DATA_DIR");
  ifstream file_check(datapath+"/"+GetEra()+"/SMP/CosThetaWeight.root");
  bool isexist=file_check.is_open();
  file_check.close();
  if(!isexist){
    cout<<"[AFBAnalyzer::SetupCosThetaWeight] no CosThetaWeight.root"<<endl;
    return;
  }
  TFile fcost(datapath+"/"+GetEra()+"/SMP/CosThetaWeight.root");
  for(const auto&& key:*(fcost.GetListOfKeys())){
    TObject* obj=((TKey*)key)->ReadObj();
    if(!obj->InheritsFrom("TH3D")) continue;
    TH3D* hist=(TH3D*)obj;
    cout<<"[AFBAnalyzer::SetupCosThetaWeight] get "<<hist->GetName()<<endl;
    map_hist_cost[hist->GetName()]=hist;
    hist->SetDirectory(0);
  }
}
void AFBAnalyzer::DeleteCosThetaWeight(){
  for(auto& iter:map_hist_cost)
    if(iter.second) delete iter.second;
}
double AFBAnalyzer::GetCosThetaWeight(double mass,double pt,double cost,TString suffix){
  double val=1.;
  if(!IsDYSample) return val;
  TString MCName=MCSample;
  if(MCName.Contains(TRegexp("^DY[0-9]Jets$"))) MCName="DYJets";
  if(MCName.Contains(TRegexp("^DYJets_Pt-[0-9]*To[0-9Inf]*$"))) MCName="DYJets";
  if(MCName.Contains(TRegexp("^DYJets_M-[0-9]*to[0-9Inf]*$"))) MCName="DYJets";
  TString hname=MCName+suffix;
  auto it=map_hist_cost.find(hname);
  if(it!=map_hist_cost.end())
    val*=GetBinContentUser(it->second,mass,pt,cost,0);
  if(val==0) val=1.;
  return val;
}
void AFBAnalyzer::test(){
  vector<LHE> lhes=GetLHEs();
  for(auto lhe:lhes) lhe.Print();
  if(lhe_j0.E()){
    cout<<"Jet event"<<endl;
    TLorentzVector z,j,system;
    system=lhe_l0+lhe_l1+lhe_j0;
    cout<<"system:";system.Print();
    TLorentzVector p0(0,0,0.5*system.M()*exp(system.Rapidity()),0.5*system.M()*exp(system.Rapidity()));
    TLorentzVector p1(0,0,-0.5*system.M()*exp(-system.Rapidity()),0.5*system.M()*exp(-system.Rapidity()));
    j=lhe_j0;
    TVector3 b=system.BoostVector();
    p0.Boost(-b);
    p1.Boost(-b);
    j.Boost(-b);
    (-b).Print();
    cout<<"After Boost"<<endl;
    p0.Print();
    p1.Print();
    j.Print();

    cout<<"JetID: "<<lhe_j0.ID()<<" Angle0: "<<p0.Angle(j.Vect())<<" Angle1:"<<p1.Angle(j.Vect())<<endl;
  }else{
    cout<<"No Jet"<<endl;
    TLorentzVector z,system,l;
    system=lhe_l0+lhe_l1;
    cout<<"system:";system.Print();
    TLorentzVector p0(0,0,0.5*system.M()*exp(system.Rapidity()),0.5*system.M()*exp(system.Rapidity()));
    TLorentzVector p1(0,0,-0.5*system.M()*exp(-system.Rapidity()),0.5*system.M()*exp(-system.Rapidity()));
    l=lhe_l0;
    TVector3 b=system.BoostVector();
    p0.Boost(-b);
    p1.Boost(-b);
    l.Boost(-b);
    cout<<"After Boost"<<endl;
    p0.Print();
    p1.Print();
    l.Print();
    cout<<"CosTheta: "<<cos(l.Angle(p0.Vect().Unit()-p1.Vect().Unit()))<<" CosThetaCS: "<<GetCosThetaCS(&lhe_l0,&lhe_l1)<<endl;
  }
}
