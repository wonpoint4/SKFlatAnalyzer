#include "AFBAnalyzer.h"

void AFBAnalyzer::initializeAnalyzer(){
  SMPAnalyzerCore::initializeAnalyzer(); //setup zpt roc z0 
  //SetupCosThetaWeight();
  
  vector<JetTagging::Parameters> jtps={JetTagging::Parameters(JetTagging::DeepJet,JetTagging::Tight,JetTagging::incl,JetTagging::comb)};
  mcCorr->SetJetTaggingParameters(jtps);

  if(fChain->GetListOfFiles()->GetEntries()){
    TString filename=fChain->GetListOfFiles()->At(0)->GetTitle();
    if(filename.Contains("SkimTree_")) IsSkimmed=true;
    else IsSkimmed=false;
  }else{
    cout<<"[AFBAnalyzer::initializeAnalyzer] no input file"<<endl;
    exit(EXIT_FAILURE);
  }
  IsNominalRun=!HasFlag("SYS")&&!HasFlag("PDFSYS")&&IsSkimmed;
}
void AFBAnalyzer::executeEvent(){
  //// FIXME some events of DYJets has nan PDF weights. I don't know why...
  if(MCSample=="DYJets"&&!isnormal(weight_Scale->at(0))) return;

  ///////////////// GEN level /////////////////////
  executeEventGen();

  ///////////////// RECO level /////////////////////
  if(!IsDATA||DataStream.Contains("SingleMuon")){
    executeEventWithParameter(MakeParameter("me"));
    executeEventWithParameter(MakeParameter("mu"));
  }
  if(!IsDATA||DataStream.Contains("DoubleMuon")){
    executeEventWithParameter(MakeParameter("mm"));
    executeEventWithParameter(MakeParameter("mM"));
    executeEventWithParameter(MakeParameter("MM"));
  }
  if(!IsDATA||DataStream.Contains("SingleElectron")||DataStream.Contains("EGamma")){
    executeEventWithParameter(MakeParameter("em"));
    executeEventWithParameter(MakeParameter("el"));
  }
  if(!IsDATA||DataStream.Contains("DoubleEG")||DataStream.Contains("EGamma")){
    executeEventWithParameter(MakeParameter("ee"));
    executeEventWithParameter(MakeParameter("eE"));
    executeEventWithParameter(MakeParameter("EE"));
  }
}
SMPAnalyzerCore::Parameter AFBAnalyzer::MakeParameter(TString key){
  Parameter p=SMPAnalyzerCore::MakeParameter(key);

  p.weightbit=0;
  if(IsSkimmed){
    if(IsNominalRun) p.weightbit|=NominalWeight;
    if(HasFlag("SYS")&&!IsDATA&&(p.channel=="ee"||p.channel=="mm")) p.weightbit|=SystematicWeight|EfficiencyWeight;
    if(HasFlag("PDFSYS")&&!IsDATA&&(p.channel=="ee"||p.channel=="mm")) p.weightbit|=PDFWeight;
  }else p.weightbit|=NominalWeight|SystematicWeight|EfficiencyWeight|PDFWeight;

  if(HasFlag("nbjet")) p.prefix+="nbjet/";
  else if(HasFlag("0bjet")) p.prefix+="0bjet/";
  if(HasFlag("highmet")) p.prefix+="highmet/";

  return p;
}
bool AFBAnalyzer::PassSelection(Parameter& p){
  if(p.prefix.Contains("highmet")){
    if(pfMET_Type1_pt<60) return false;
    if(IsNominalRun) FillCutflow(p.prefix+p.hprefix+"cutflow","METCut",p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight);
  }

  int n_bjet=0;
  std::vector<Jet> jets=GetJets("tightLepVeto",40,2.4);
  std::sort(jets.begin(),jets.end(),PtComparing);
  JetTagging::Parameters jtp = JetTagging::Parameters(JetTagging::DeepJet,JetTagging::Tight,JetTagging::incl,JetTagging::comb);
  for(const auto& jet:jets)
    if(jet.GetTaggerResult(jtp.j_Tagger) > mcCorr->GetJetTaggingCutValue(jtp.j_Tagger, jtp.j_WP))
      n_bjet++;
  p.intmap["nbjet"]=n_bjet;
  p.doublemap["btagSF"]=mcCorr->GetBTaggingReweight_1a(jets,jtp);
  p.doublemap["btagSF_hup"]=mcCorr->GetBTaggingReweight_1a(jets,jtp,"SystUpHTag");
  p.doublemap["btagSF_hdown"]=mcCorr->GetBTaggingReweight_1a(jets,jtp,"SystDownHTag");
  p.doublemap["btagSF_lup"]=mcCorr->GetBTaggingReweight_1a(jets,jtp,"SystUpLTag");
  p.doublemap["btagSF_ldown"]=mcCorr->GetBTaggingReweight_1a(jets,jtp,"SystDownLTag");

  if(p.prefix.Contains("nbjet")&&!n_bjet) return false;
  if(p.prefix.Contains("0bjet")&&n_bjet) return false;
  if(IsNominalRun) FillCutflow(p.prefix+p.hprefix+"cutflow","BJetCut",p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight);
  if(IsNominalRun) FillCutflow(p.prefix+p.hprefix+"cutflow","BJetCutSF",p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*p.doublemap["btagSF"]);

  if(!SMPAnalyzerCore::PassSelection(p)) return false;  
  return true;
}
void AFBAnalyzer::executeEventGen(){
  costhetaweight=1.;
  costhetaweight_up=1.;
  costhetaweight_down=1.;
  if(IsDYSample||MCSample.Contains("GamGamToLL")||MCSample.Contains("TTLL")){
    //////////////////////// Check LHE /////////////////////////
    if(abs(lhe_l0.ID())!=15&&abs(lhe_l1.ID())!=15){
      Parameter p;
      double letacut=2.4;
      if( (abs(lhe_l0.ID())==11&&abs(lhe_l1.ID())==11) || (!lhes.size()&&abs(gen_l0.PID())==11&&abs(gen_l1.PID())==11) ){
	p=MakeParameter("ee");
	p.c.lepton0pt=25;
	p.c.lepton1pt=15;
      }else if( (abs(lhe_l0.ID())==13&&abs(lhe_l1.ID())==13) || (!lhes.size()&&abs(gen_l0.PID())==13&&abs(gen_l1.PID())==13) ){
	p=MakeParameter("mm");
	p.c.lepton0pt=25; //sync with electron
	p.c.lepton1pt=15; //sync with electron
      }else{
	if(IsDYSample||MCSample.Contains("GamGamToLL")){
	  cout<<"[AFBAnalyzer::executeEvent()] something is wrong l0.ID="<<abs(lhe_l0.ID())<<endl;
	  vector<LHE> lhes=GetLHEs();
	  for(auto& lhe:lhes) lhe.Print();
	  exit(EXIT_FAILURE);
	}else if(MCSample.Contains("TTLL")){
	  return;
	}
      }
      
      //////////////////////// GEN /////////////////////////
      TLorentzVector gen_Z=gen_l0+gen_l1;
      double gen_Zmass=gen_Z.M();
      double gen_Zrap=gen_Z.Rapidity();
      double gen_Zpt=gen_Z.Pt();
      double gen_cost_correct=-999;
      if(gen_p0.PID()==21||gen_p0.PID()==22){
	if(gen_p1.PID()==21||gen_p1.PID()==22) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,0);
	else if(gen_p1.PID()>0) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,-1);
	else if(gen_p1.PID()<0) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,1);
      }else if(gen_p0.PID()>0){
	if(gen_p1.PID()==21||gen_p1.PID()==22) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,1);
	else if(gen_p1.PID()>0) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,0);
	else if(gen_p1.PID()<0) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,1);
      }else if(gen_p0.PID()<0){
	if(gen_p1.PID()==21||gen_p1.PID()==22) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,-1);
	else if(gen_p1.PID()>0) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,-1);
	else if(gen_p1.PID()<0) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,0);
      }
      if(gen_cost_correct==-999){
	cout<<"wrong pid for parton: "<<gen_p0.PID()<<" "<<gen_p1.PID()<<endl;
	exit(EXIT_FAILURE);
      }
      //costhetaweight=GetCosThetaWeight(gen_Zmass,gen_Zpt,gen_cost_correct,"_pdg");
      //costhetaweight_up=GetCosThetaWeight(gen_Zmass,gen_Zpt,gen_cost_correct,"_up");
      //costhetaweight_down=GetCosThetaWeight(gen_Zmass,gen_Zpt,gen_cost_correct,"_down");
      
      map<TString,double> map_weight;
      map_weight[""]=p.w.lumiweight*p.w.zptweight*costhetaweight;
      map_weight["_noweight"]=p.w.lumiweight;
      map_weight["_nozptweight"]=p.w.lumiweight*costhetaweight;
      map_weight["_nocosthetaweight"]=p.w.lumiweight*p.w.zptweight;

      //////////////// Fill LHE,Gen hists //////////////////////
      if(!IsSkimmed){
	TLorentzVector lhe_Z=lhe_l0+lhe_l1;
	double lhe_Zmass=lhe_Z.M();
	double lhe_Zrap=lhe_Z.Rapidity();
	double lhe_Zpt=lhe_Z.Pt();
	FillHistsAFB(p.prefix,"lhe_","",(Particle*)&lhe_l0,(Particle*)&lhe_l1,map_weight);
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
}
int AFBAnalyzer::GetUnfoldBin(int nbin,const double* bins,double value,double cost){
  int i;
  int forward=cost>0?1:0;
  if(value<bins[0]) return 0;
  value=TMath::Min(value,bins[nbin]-0.1);
  i=TMath::BinarySearch(nbin+1,bins,value);
  return forward*nbin+i+1;
}
void AFBAnalyzer::executeEventWithParameter(Parameter& p){
  SMPAnalyzerCore::executeEventWithParameter(p);
  // response matrix for unfolding
  if(IsSkimmed) return;
  if(p.channel=="ee"){
    if(abs(lhe_l0.ID())!=11||abs(lhe_l1.ID())!=11) return;
  }else if(p.channel=="mm"){
    if(abs(lhe_l0.ID())!=13||abs(lhe_l1.ID())!=13) return;
  }else return;
  TLorentzVector gen_ll_dressed=gen_l0_dressed+gen_l1_dressed;
  double genm=-1;
  double geny=-100;
  double genpt=-1;
  double gencost=0;
  if(gen_l0_dressed.Pt()>25||gen_l1_dressed.Pt()>25){
    if(gen_l0_dressed.Pt()>15&&gen_l1_dressed.Pt()>15){
      if(fabs(gen_l0_dressed.Eta())<2.4&&fabs(gen_l1_dressed.Eta())<2.4){
	if(gen_ll_dressed.M()>52){
	  genm=gen_ll_dressed.M();
	  geny=gen_ll_dressed.Rapidity();
	  genpt=gen_ll_dressed.Pt();
	  gencost=GetCosThetaCS(&gen_l0_dressed,&gen_l1_dressed);
	}
      }
    }
  }
  Parameter pgen=p;
  ResetRecoWeights(pgen);
  EvalWeights(pgen);
  for(auto [wname,genweight]:pgen.weightmap){
    double recoweight=0;
    double recom=-1;
    double recoy=-100;
    double recopt=-1;
    double recocost=0;
    if(p.weightmap.find(wname)!=p.weightmap.end()){
      recoweight=p.weightmap[wname];
      recom=(*p.lepton0+*p.lepton1).M();
      recoy=(*p.lepton0+*p.lepton1).Rapidity();
      recopt=(*p.lepton0+*p.lepton1).Pt();
      recocost=GetCosThetaCS(p.lepton0,p.lepton1);
    }
    //cout<<"wname:"<<wname<<" genweight:"<<genweight<<" recoweight:"<<recoweight<<" recom:"<<recom<<endl;
    int imbin=GetUnfoldBin(afb_mbinnum,afb_mbin,genm,gencost);
    int jmbin=GetUnfoldBin(afb_mbinnum,afb_mbin,recom,recocost);
    int iybin=GetUnfoldBin(afb_ybinnum,afb_ybin,geny,gencost);
    int jybin=GetUnfoldBin(afb_ybinnum,afb_ybin,recoy,recocost);
    int iptbin=GetUnfoldBin(afb_ptbinnum,afb_ptbin,genpt,gencost);
    int jptbin=GetUnfoldBin(afb_ptbinnum,afb_ptbin,recopt,recocost);

    FillHist(p.prefix+p.hprefix+"response_afbm"+p.suffix+wname,imbin,jmbin,recoweight,2*afb_mbinnum,1,2*afb_mbinnum+1,2*afb_mbinnum,1,2*afb_mbinnum+1);
    FillHist(p.prefix+p.hprefix+"response_afbm"+p.suffix+wname,imbin,0,genweight-recoweight,2*afb_mbinnum,1,2*afb_mbinnum+1,2*afb_mbinnum,1,2*afb_mbinnum+1);
    FillHist(p.prefix+p.hprefix+"response_afby"+p.suffix+wname,iybin,jybin,recoweight,2*afb_ybinnum,1,2*afb_ybinnum+1,2*afb_ybinnum,1,2*afb_ybinnum+1);
    FillHist(p.prefix+p.hprefix+"response_afby"+p.suffix+wname,iybin,0,genweight-recoweight,2*afb_ybinnum,1,2*afb_ybinnum+1,2*afb_ybinnum,1,2*afb_ybinnum+1);
    FillHist(p.prefix+p.hprefix+"response_afbpt"+p.suffix+wname,iptbin,jptbin,recoweight,2*afb_ptbinnum,1,2*afb_ptbinnum+1,2*afb_ptbinnum,1,2*afb_ptbinnum+1);
    FillHist(p.prefix+p.hprefix+"response_afbpt"+p.suffix+wname,iptbin,0,genweight-recoweight,2*afb_ptbinnum,1,2*afb_ptbinnum+1,2*afb_ptbinnum,1,2*afb_ptbinnum+1);
    double massregion[]={52,77,106,280,3000};
    for(int k=0;k<4;k++){
      int ibin=(genm>=massregion[k]&&genm<massregion[k+1]) ? iybin : 0;
      int jbin=(recom>=massregion[k]&&recom<massregion[k+1]) ? jybin : 0;
      FillHist(p.prefix+p.hprefix+Form("response_afby_m%d",k)+p.suffix+wname,ibin,jbin,recoweight,2*afb_ybinnum,1,2*afb_ybinnum+1,2*afb_ybinnum,1,2*afb_ybinnum+1);
      FillHist(p.prefix+p.hprefix+Form("response_afby_m%d",k)+p.suffix+wname,ibin,0,genweight-recoweight,2*afb_ybinnum,1,2*afb_ybinnum+1,2*afb_ybinnum,1,2*afb_ybinnum+1);

      ibin=(genm>=massregion[k]&&genm<massregion[k+1]) ? iptbin : 0;
      jbin=(recom>=massregion[k]&&recom<massregion[k+1]) ? jptbin : 0;
      FillHist(p.prefix+p.hprefix+Form("response_afbpt_m%d",k)+p.suffix+wname,ibin,jbin,recoweight,2*afb_ptbinnum,1,2*afb_ptbinnum+1,2*afb_ptbinnum,1,2*afb_ptbinnum+1);
      FillHist(p.prefix+p.hprefix+Form("response_afbpt_m%d",k)+p.suffix+wname,ibin,0,genweight-recoweight,2*afb_ptbinnum,1,2*afb_ptbinnum+1,2*afb_ptbinnum,1,2*afb_ptbinnum+1);
    }
  }
}
void AFBAnalyzer::EvalWeights(Parameter& p){
  if(p.weightbit&NominalWeight){
    p.weightmap[""]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"];
  }
  if(p.weightbit&SystematicWeight){
    if(!IsDATA){
      p.weightmap["_noPUweight"]=p.w.lumiweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"]; //need for AN
      p.weightmap["_PUweight_up"]=p.w.lumiweight*p.w.PUweight_up*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"];
      p.weightmap["_PUweight_down"]=p.w.lumiweight*p.w.PUweight_down*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"];
      
      p.weightmap["_noprefireweight"]=p.w.lumiweight*p.w.PUweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"];
      p.weightmap["_prefireweight_up"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight_up*p.w.zptweight*p.w.z0weight*p.w.weakweight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"];
      p.weightmap["_prefireweight_down"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight_down*p.w.zptweight*p.w.z0weight*p.w.weakweight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"];
      
      p.weightmap["_nozptweight"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"];
      p.weightmap["_noz0weight"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"];
      p.weightmap["_noweakweight"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"];
      
      p.weightmap["_nobtagSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
      p.weightmap["_btagSF_hup"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF_hup"];
      p.weightmap["_btagSF_hdown"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF_hdown"];
      p.weightmap["_btagSF_lup"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF_lup"];
      p.weightmap["_btagSF_ldown"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF_ldown"];

      //p.weightmap["_nocosthetaweight"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"];
      //p.weightmap["_costhetaweight_up"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*costhetaweight_up*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"];
      //p.weightmap["_costhetaweight_down"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*costhetaweight_down*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"];
      
      //p.weightmap["_noefficiencySF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*costhetaweight*p.w.CFSF*p.doublemap["btagSF"];
      
      //p.weightmap["_noelectronRECOSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*costhetaweight*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"];
      //p.weightmap["_noIDSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*costhetaweight*p.w.electronRECOSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"];
      //p.weightmap["_nomuonISOSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"];      
      //p.weightmap["_notriggerSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.CFSF*p.doublemap["btagSF"];
      for(int j=0,nj=fEff->nreplica;j<nj;j++){
	double electronRECOSF=p.w.electronRECOSF_sys.size() ? p.w.electronRECOSF_sys[0][j] : 1.;
	double electronIDSF=p.w.electronIDSF_sys.size() ? p.w.electronIDSF_sys[0][j] : 1.;
	double muonIDSF=p.w.muonIDSF_sys.size() ? p.w.muonIDSF_sys[0][j] : 1.;
	double triggerSF=p.w.triggerSF_sys.size() ? p.w.triggerSF_sys[0][j] : 1.;
	p.weightmap[Form("_efficiencySF_stat%d",j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*costhetaweight*electronRECOSF*electronIDSF*muonIDSF*p.w.muonISOSF*triggerSF*p.w.CFSF*p.doublemap["btagSF"];
      }

      for(int i=1,ni=p.w.electronRECOSF_sys.size();i<ni;i++){
	for(int j=0,nj=p.w.electronRECOSF_sys[i].size();j<nj;j++){
	  p.weightmap[Form("_electronRECOSF_s%d_m%d",i,j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*costhetaweight*p.w.electronRECOSF_sys[i][j]*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"];
	}
      }

      for(int i=1,ni=p.w.electronIDSF_sys.size();i<ni;i++){
	for(int j=0,nj=p.w.electronIDSF_sys[i].size();j<nj;j++){
	  p.weightmap[Form("_electronIDSF_s%d_m%d",i,j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF_sys[i][j]*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"];
	}
      }

      for(int i=1,ni=p.w.muonIDSF_sys.size();i<ni;i++){
	for(int j=0,nj=p.w.muonIDSF_sys[i].size();j<nj;j++){
	  p.weightmap[Form("_muonIDSF_s%d_m%d",i,j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF_sys[i][j]*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"];
	}
      }

      for(int i=1,ni=p.w.triggerSF_sys.size();i<ni;i++){
	for(int j=0,nj=p.w.triggerSF_sys[i].size();j<nj;j++){
	  p.weightmap[Form("_triggerSF_s%d_m%d",i,j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF_sys[i][j]*p.w.CFSF*p.doublemap["btagSF"];
	}
      }

      
      //p.weightmap["_noCFSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"];
      p.weightmap["_CFSF_up"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF_up*p.doublemap["btagSF"];
      p.weightmap["_CFSF_down"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF_down*p.doublemap["btagSF"];

    }
  }
  if(p.weightbit&PDFWeight){
    for(unsigned int i=0;i<weight_Scale->size();i++){
      p.weightmap[Form("_scalevariation%d",i)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"]*weight_Scale->at(i);
    }
    for(unsigned int i=0;i<weight_PDF->size();i++){
      p.weightmap[Form("_pdf%d",i)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"]*weight_PDF->at(i);
    }
    if(weight_AlphaS->size()==2){
      p.weightmap["_alphaS_down"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"]*weight_AlphaS->at(0);
      p.weightmap["_alphaS_up"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"]*weight_AlphaS->at(1);
    }

    if(MCSample.Contains("MiNNLO")){
      p.weightmap["_sthw2_down"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"]*weight_sthw2->at(0);
      p.weightmap["_sthw2_up"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"]*weight_sthw2->at(2);
      p.weightmap["_largeptscales"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"]*weight_largeptscales->at(0);
      p.weightmap["_q0_up"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"]*weight_q0->at(0);
      p.weightmap["_q0_down"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*costhetaweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.doublemap["btagSF"]*weight_q0->at(2);      
    }
  }
  return;
}
void AFBAnalyzer::ResetRecoWeights(Parameter& p){
  p.w.prefireweight=1.; p.w.prefireweight_up=1.; p.w.prefireweight_down=1.;
  p.w.z0weight=1.;
  p.w.electronRECOSF=1.;
  p.w.electronRECOSF_sys=fEff->GetStructure(p.k.electronRECOSF);
  p.w.electronIDSF=1.;
  p.w.electronIDSF_sys=fEff->GetStructure(p.k.electronIDSF);
  p.w.muonIDSF=1.;
  p.w.muonIDSF_sys=fEff->GetStructure(p.k.muonIDSF);
  p.w.muonISOSF=1.;
  p.w.muonISOSF_sys=fEff->GetStructure(p.k.muonISOSF);
  p.w.triggerSF=1.;
  p.w.triggerSF_sys=fEff->GetStructure(p.k.triggerSF[0]);
  p.doublemap["btagSF"]=1.;
}
void AFBAnalyzer::FillHists(Parameter& p){
  if(!IsSkimmed) return;
  TLorentzVector dilepton=*p.lepton0+*p.lepton1;
  double dimass=dilepton.M();
  double dirap=dilepton.Rapidity();
  double dipt=dilepton.Pt();
  FillHistsAFB(p.prefix,p.hprefix,p.suffix,(Particle*)p.lepton0,(Particle*)p.lepton1,p.weightmap);
  FillHist(p.prefix+p.hprefix+"nbjet"+p.suffix,dimass,dirap,dipt,p.intmap["nbjet"],p.weightmap,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,10,0,10);
  FillHist(p.prefix+p.hprefix+"z0"+p.suffix,dimass,dirap,dipt,vertex_Z,SelectWeights(p.weightmap,{"","_noz0weight"}),grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,120,-15,15);
  map<TString,double> map_PUweight=SelectWeights(p.weightmap,{"","_noPUweight","_PUweight_up","_PUweight_down"});
  FillHist(p.prefix+p.hprefix+"nPV"+p.suffix,dimass,dirap,dipt,nPV,map_PUweight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,100,0,100);
  FillHist(p.prefix+p.hprefix+"rho"+p.suffix,dimass,dirap,dipt,Rho,map_PUweight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,50,0,50);
  FillHist(p.prefix+p.hprefix+"met"+p.suffix,dimass,dirap,dipt,pfMET_Type1_pt,map_PUweight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,100,0,200);
  if(IsDYSample&&p.hprefix==""&&IsNominalRun){
    vector<Gen> gens=GetGens();
    Gen truth_l0=GetGenMatchedLepton(*p.lepton0,gens);
    Gen truth_l1=GetGenMatchedLepton(*p.lepton1,gens);
    if(!truth_l0.IsEmpty()&&!truth_l1.IsEmpty()) 
      FillHistsAFB(p.prefix,"truth_",p.suffix,(Particle*)&truth_l0,(Particle*)&truth_l1,p.weightmap);
    //else cout<<"no matching"<<endl;
  }
  // fill fake hists
  /*
  if(p.channel=="EE"&&p.prefix.Contains("EE")){
    for(int i=0,n=p.aelectrons.size();i<n;i++){
      for(int j=i+1,n=p.aelectrons.size();j<n;j++){
	Parameter this_p=p;
	this_p.prefix.ReplaceAll("EE","ee");
	this_p.hprefix="fake_";
	this_p.lepton0=&this_p.aelectrons.at(i);
	this_p.lepton1=&this_p.aelectrons.at(j);
	this_p.w.lumiweight*=GetFakeRate(&this_p.aelectrons.at(i))*GetFakeRate(&this_p.aelectrons.at(j));
	for(int k=j+1,n=p.aelectrons.size();k<n;k++) this_p.w.lumiweight*=1+GetFakeRate(&this_p.aelectrons.at(k));
	{
	  double pt=this_p.aelectrons.at(i).Pt();
	  double riso=this_p.aelectrons.at(i).RelIso();
	  double f=TMath::Max(0.,TMath::Min(1.,(pt-30)/30));
	  //if(fabs(this_p.aelectrons.at(i).Eta())<1.479) this_p.aelectrons.at(i)*=(1+f*riso-f*0.506/pt)/(1+f*0.0478);
	  //else this_p.aelectrons.at(i)*=(1+f*riso-f*0.963/pt)/(1+f*0.0658);
	}
	{
	  double pt=this_p.aelectrons.at(j).Pt();
	  double riso=this_p.aelectrons.at(j).RelIso();
	  double f=TMath::Max(0.,TMath::Min(1.,(pt-30)/30));
	  //if(fabs(this_p.aelectrons.at(j).Eta())<1.479) this_p.aelectrons.at(j)*=(1+f*riso-f*0.506/pt)/(1+f*0.0478);
	  //else this_p.aelectrons.at(j)*=(1+f*riso-f*0.963/pt)/(1+f*0.0658);
	}
	if(PassSelection(this_p)) FillHists(this_p);
      }
    }
  }
  if(p.channel=="MM"&&p.prefix.Contains("MM")){
    for(int i=0,n=p.amuons.size();i<n;i++){
      for(int j=i+1,n=p.amuons.size();j<n;j++){
	Parameter this_p=p;
	this_p.prefix.ReplaceAll("MM","mm");
	this_p.hprefix="fake_";
	this_p.lepton0=&this_p.amuons.at(i);
	this_p.lepton1=&this_p.amuons.at(j);
	this_p.w.lumiweight*=GetFakeRate(&this_p.amuons.at(i))*GetFakeRate(&this_p.amuons.at(j));
	for(int k=j+1,n=p.amuons.size();k<n;k++) this_p.w.lumiweight*=1+GetFakeRate(&this_p.amuons.at(k));
	{
	  double pt=this_p.amuons.at(i).Pt();
	  double riso=this_p.amuons.at(i).RelIso();
	  double f=TMath::Max(0.,TMath::Min(1.,(pt-30)/30));
	  //this_p.amuons.at(i)*=(1+f*riso)/(1+f*0.1);
	}
	{
	  double pt=this_p.amuons.at(j).Pt();
	  double riso=this_p.amuons.at(j).RelIso();
	  double f=TMath::Max(0.,TMath::Min(1.,(pt-30)/30));
	  //this_p.amuons.at(j)*=(1+f*riso)/(1+f*0.1);
	}
	if(PassSelection(this_p)) FillHists(this_p);
      }
    }
  }
  */
  if(p.channel=="EE"&&p.prefix.Contains("EE")){
    Parameter this_p=p;
    this_p.prefix.ReplaceAll("EE","ee");
    this_p.hprefix="fake_"+this_p.hprefix;
    this_p.lepton0=&this_p.aelectrons.at(0);
    this_p.lepton1=&this_p.aelectrons.at(1);
    this_p.w.lumiweight*=GetFakeRate(this_p.lepton0)*GetFakeRate(this_p.lepton1);
    FillHists(this_p);
  }
  if(p.channel=="MM"&&p.prefix.Contains("MM")){
    Parameter this_p=p;
    this_p.prefix.ReplaceAll("MM","mm");
    this_p.hprefix="fake_"+this_p.hprefix;
    this_p.lepton0=&this_p.amuons.at(0);
    this_p.lepton1=&this_p.amuons.at(1);
    this_p.w.lumiweight*=GetFakeRate(this_p.lepton0)*GetFakeRate(this_p.lepton1);
    FillHists(this_p);
  }
}

AFBAnalyzer::AFBAnalyzer(){}
AFBAnalyzer::~AFBAnalyzer(){
  //DeleteCosThetaWeight();
}
double AFBAnalyzer::GetCosThetaCS(const Particle *p0,const Particle *p1,int direction){
  if(!p0||!p1) return 0.;
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
  if(direction==0) direction=dilepton.Pz()>0?1:-1;
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

  FillHist(pre+hpre+"l0pt"+suf,dimass,dirap,dipt,l0->Pt(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,lptbinnum,(double*)lptbin);
  FillHist(pre+hpre+"l1pt"+suf,dimass,dirap,dipt,l1->Pt(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,lptbinnum,(double*)lptbin);
  FillHist(pre+hpre+"lpt"+suf,dimass,dirap,dipt,l0->Pt(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,lptbinnum,(double*)lptbin);
  FillHist(pre+hpre+"lpt"+suf,dimass,dirap,dipt,l1->Pt(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,lptbinnum,(double*)lptbin);
  
  FillHist(pre+hpre+"l0eta"+suf,dimass,dirap,dipt,l0->Eta(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,60,-3,3);
  FillHist(pre+hpre+"l1eta"+suf,dimass,dirap,dipt,l1->Eta(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,60,-3,3);
  FillHist(pre+hpre+"leta"+suf,dimass,dirap,dipt,l0->Eta(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,60,-3,3);
  FillHist(pre+hpre+"leta"+suf,dimass,dirap,dipt,l1->Eta(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,60,-3,3);
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
