#include "AFBAnalyzer.h"

void AFBAnalyzer::initializeAnalyzer(){
  SMPAnalyzerCore::initializeAnalyzer(); //setup zpt roc z0 
  
  IsSkimmed= GetSkimName()!="" ? true : false;
  IsNominalRun=!HasFlag("SYS")&&!HasFlag("PDFSYS")&&IsSkimmed;

  PDFbase=LHAPDF::mkPDF(306000);
  PDFnf4=LHAPDF::mkPDF(325500);

}
void AFBAnalyzer::executeEvent(){
  //cout<<"Event:"<<event<<endl;
  //// FIXME some events of DYJets has nan PDF weights. I don't know why...
  if(MCSample=="DYJets"&&!isnormal(weight_Scale->at(0))) return;

  ///////////////// GEN level /////////////////////
  genfid_b0=NULL;
  for(auto& lhe:lhes){
    if(lhe.Status()!=1) continue;
    if(abs(lhe.ID())!=5) continue;
    if(lhe.Pt()<40) continue;
    if(fabs(lhe.Eta())>2.4) continue;
    if(genfid_b0&&genfid_b0->Pt()>lhe.Pt()) continue;
    genfid_b0=&lhe;
  }
  executeEventGen();

  ///////////////// RECO level /////////////////////
  if(!IsDATA||DataStream.Contains("SingleMuon")){
    if(IsNominalRun) executeEventWithParameter(MakeParameter("mu"));
  }
  if(!IsDATA||DataStream.Contains("DoubleMuon")){
    executeEventWithParameter(MakeParameter("mm"));
    if(HasFlag("SYS")){
      for(TString syst:{"jet_scale_up","jet_scale_down","jet_smear_up","jet_smear_down"}){
	executeEventWithParameter(MakeParameter("mm",syst));
      }
    }
    //For the fake estimation with transfer-factor method
    //if(IsNominalRun) executeEventWithParameter(MakeParameter("mM"));
    //if(IsNominalRun) executeEventWithParameter(MakeParameter("MM"));
  }
  if(!IsDATA||DataStream.Contains("SingleElectron")||DataStream.Contains("EGamma")){
    if(IsNominalRun) executeEventWithParameter(MakeParameter("el"));
  }
  if(!IsDATA||DataStream.Contains("DoubleEG")||DataStream.Contains("EGamma")){
    executeEventWithParameter(MakeParameter("ee"));
    if(HasFlag("SYS")){
      for(TString syst:{"jet_scale_up","jet_scale_down","jet_smear_up","jet_smear_down"}){
	executeEventWithParameter(MakeParameter("ee",syst));
      }
    }
    //For the fake estimation with transfer-factor method
    //if(IsNominalRun) executeEventWithParameter(MakeParameter("eE"));
    //if(IsNominalRun) executeEventWithParameter(MakeParameter("EE"));
  }
}
SMPAnalyzerCore::Parameter AFBAnalyzer::MakeParameter(TString key,TString option){
  Parameter p=SMPAnalyzerCore::MakeParameter(key,option);

  p.variationbits=0;
  if(IsSkimmed){
    if(IsNominalRun) p.variationbits|=NominalWeight;
    if(HasFlag("SYS")&&!IsDATA&&(p.channel=="ee"||p.channel=="mm")){
      if(p.suffix==""){
	p.variationbits|=SystematicWeight|EfficiencyWeight;
      }else{
	p.variationbits|=NominalWeight;
      }
    }
    if(HasFlag("PDFSYS")&&!IsDATA&&(p.channel=="ee"||p.channel=="mm")) p.variationbits|=PDFWeight;
  }else p.variationbits|=NominalWeight|SystematicWeight|EfficiencyWeight|PDFWeight;
  if(p.option.Contains("DeepCSV")) p.prefix="DeepCSV/"+p.prefix;
  else if(p.option.Contains("DeepJet::Tight::mujets")) p.prefix="mujets/"+p.prefix;

  if(HasFlag("nbjet")) p.prefix+="nbjet/";
  else if(HasFlag("0bjet")) p.prefix+="0bjet/";
  if(HasFlag("highmet")) p.prefix+="highmet/";

  if(IsDYSample&&p.hprefix==""){
    if(abs(genWeight_id1)==5||abs(genWeight_id2)==5){
      //p.hprefix="bx_";
    }
  }

  return p;
}
bool AFBAnalyzer::PassSelection(Parameter& p,bool cutflow){
  if(p.prefix.Contains("highmet")){
    if(pfMET_Type1_pt<60) return false;
    if(IsNominalRun) FillCutflow(p.prefix+p.hprefix+"cutflow","METCut",p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight);
  }

  if(p.prefix.Contains("nbjet")){
    if(!p.bjets.size()) return false;
    if(p.c.jetpt>0&&p.bjets[0].Pt()<p.c.jetpt) return false;
  }
  if(p.prefix.Contains("0bjet")&&p.bjets.size()&&p.bjets[0].Pt()>p.c.jetpt) return false;
  if(IsNominalRun) FillCutflow(p.prefix+p.hprefix+"cutflow","BJetCut",p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight);
  if(IsNominalRun) FillCutflow(p.prefix+p.hprefix+"cutflow","BJetCutSF",p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*p.w.btagSF);

  if(!SMPAnalyzerCore::PassSelection(p,cutflow)) return false;  
  return true;
}
void AFBAnalyzer::executeEventGen(){
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
      
      map<TString,double> map_weight;

      map_weight[""]=p.w.lumiweight*p.w.PUweight*p.w.zptweight*p.w.weakweight*p.w.topptweight;
      map_weight["_nozptweight"]=p.w.lumiweight*p.w.PUweight*p.w.weakweight*p.w.topptweight;

      //////////////// Fill LHE,Gen hists //////////////////////
      if(!IsSkimmed){
	for(auto& [vsuf,weight]:map_weight){
	  TLorentzVector lhe_Z=lhe_l0+lhe_l1;
	  double lhe_Zmass=lhe_Z.M();
	  double lhe_Zrap=lhe_Z.Rapidity();
	  double lhe_Zpt=lhe_Z.Pt();
	  FillHistsAFB(p.prefix,"lhe_",p.suffix+vsuf,(Particle*)&lhe_l0,(Particle*)&lhe_l1,weight);
	  if(lhe_j0.Pt()){
	    //FillHist(p.prefix+"lhe_costhetaR"+p.suffix+vsuf,lhe_Zmass,lhe_Zrap,lhe_Zpt,GetCosThetaR(&lhe_l0,&lhe_l1,&lhe_j0,0),weight,afb_mbinnum,(double*)afb_mbin,afb_ybinnum,(double*)afb_ybin,afb_ptbinnum,(double*)afb_ptbin,20,-1,1);
	    //FillHist(p.prefix+"lhe_costhetaT"+p.suffix+vsuf,lhe_Zmass,lhe_Zrap,lhe_Zpt,GetCosThetaT(&lhe_l0,&lhe_l1,&lhe_j0,0),weight,afb_mbinnum,(double*)afb_mbin,afb_ybinnum,(double*)afb_ybin,afb_ptbinnum,(double*)afb_ptbin,20,-1,1);
	  }
	  if(lhe_l0.Pt()>p.c.lepton0pt&&lhe_l1.Pt()>p.c.lepton1pt&&fabs(lhe_l0.Eta())<letacut&&fabs(lhe_l1.Eta())<letacut){
	    FillHistsAFB(p.prefix,"lhefid_",p.suffix+vsuf,(Particle*)&lhe_l0,(Particle*)&lhe_l1,weight);
	  }
	  FillHistsAFB(p.prefix,"gen_",p.suffix+vsuf,(Particle*)&gen_l0,(Particle*)&gen_l1,weight);
	  FillHistsAFB(p.prefix,"gen_","_dressed"+p.suffix+vsuf,(Particle*)&gen_l0_dressed,(Particle*)&gen_l1_dressed,weight);
	  FillHistsAFB(p.prefix,"gen_","_bare"+p.suffix+vsuf,(Particle*)&gen_l0_bare,(Particle*)&gen_l1_bare,weight);
	  if(TMath::Max(gen_l0.Pt(),gen_l1.Pt())>p.c.lepton0pt&&TMath::Min(gen_l0.Pt(),gen_l1.Pt())>p.c.lepton1pt&&fabs(gen_l0.Eta())<letacut&&fabs(gen_l1.Eta())<letacut){
	    FillHistsAFB(p.prefix,"genfid_",p.suffix+vsuf,(Particle*)&gen_l0,(Particle*)&gen_l1,weight);
	  }
	  if(TMath::Max(gen_l0_dressed.Pt(),gen_l1_dressed.Pt())>p.c.lepton0pt&&TMath::Min(gen_l0_dressed.Pt(),gen_l1_dressed.Pt())>p.c.lepton1pt&&fabs(gen_l0_dressed.Eta())<letacut&&fabs(gen_l1_dressed.Eta())<letacut){
	    FillHistsAFB(p.prefix,"genfid_","_dressed"+p.suffix+vsuf,(Particle*)&gen_l0_dressed,(Particle*)&gen_l1_dressed,weight);
	    if(genfid_b0){
	      //cout<<"Fill genfid_dressed "<<p.prefix<<" "<<map_weight[""]<<endl;
	      FillHistsRecoil(p.prefix,"genfid_","_dressed"+p.suffix+vsuf,(Particle*)&gen_l0_dressed,(Particle*)&gen_l1_dressed,genfid_b0,weight);
	    }
	  }
	  if(TMath::Max(gen_l0_bare.Pt(),gen_l1_bare.Pt())>p.c.lepton0pt&&TMath::Min(gen_l0_bare.Pt(),gen_l1_bare.Pt())>p.c.lepton1pt&&fabs(gen_l0_bare.Eta())<letacut&&fabs(gen_l1_bare.Eta())<letacut){
	    FillHistsAFB(p.prefix,"genfid_","_bare"+p.suffix+vsuf,(Particle*)&gen_l0_bare,(Particle*)&gen_l1_bare,weight);
	  }
	  //FillHist(p.prefix+"gen_costhetaCS_correct"+p.suffix+vsuf,gen_Zmass,gen_Zrap,gen_Zpt,gen_cost_correct,weight,afb_mbinnum,(double*)afb_mbin,afb_ybinnum,(double*)afb_ybin,afb_ptbinnum,(double*)afb_ptbin,20,-1,1);
	  if(vsuf==""){
	    FillHist(p.prefix+"gen_nPU"+p.suffix+vsuf,gen_Zmass,gen_Zrap,gen_Zpt,nPileUp,p.w.lumiweight*p.w.PUweight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,100,0,100);
	    FillHist(p.prefix+"gen_nPU_noPUweight"+p.suffix+vsuf,gen_Zmass,gen_Zrap,gen_Zpt,nPileUp,p.w.lumiweight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,100,0,100);
	    FillHist(p.prefix+"gen_nPU_PUweight_down"+p.suffix+vsuf,gen_Zmass,gen_Zrap,gen_Zpt,nPileUp,p.w.lumiweight*p.w.PUweight_down,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,100,0,100);
	    FillHist(p.prefix+"gen_nPU_PUweight_up"+p.suffix+vsuf,gen_Zmass,gen_Zrap,gen_Zpt,nPileUp,p.w.lumiweight*p.w.PUweight_up,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,100,0,100);
	  }
	}
      }
    }
  }
}
int AFBAnalyzer::GetUnfoldBin(int nbin,const double* bins,double value,double cost){
  int i;
  int forward=cost>0?1:0;
  if(value<bins[0]) return 0;
  if(value>=bins[nbin]) return 0;
  //value=TMath::Min(value,bins[nbin]-0.1);
  i=TMath::BinarySearch(nbin+1,bins,value);
  return forward*nbin+i+1;
}
void AFBAnalyzer::executeEventWithParameter(Parameter& p){
  SMPAnalyzerCore::executeEventWithParameter(p);
  // response matrix for unfolding
  if(IsSkimmed) return;
  if(p.channel!="ee"&&p.channel!="mm") return;

  TString gen_channel="";
  if(abs(lhe_l0.ID())==11&&abs(lhe_l1.ID())==11) gen_channel="ee";
  else if(abs(lhe_l0.ID())==13&&abs(lhe_l1.ID())==13) gen_channel="mm";

  TLorentzVector gen_ll_dressed=gen_l0_dressed+gen_l1_dressed;
  double genm=-1;
  double geny=-100;
  double genpt=-1;
  double gencost=0;
  if(gen_channel==p.channel){
    if(TMath::Max(gen_l0_dressed.Pt(),gen_l1_dressed.Pt())>25){
      if(TMath::Min(gen_l0_dressed.Pt(),gen_l1_dressed.Pt())>15){
	if(fabs(gen_l0_dressed.Eta())<2.4&&fabs(gen_l1_dressed.Eta())<2.4){
	  if(gen_ll_dressed.M()>52){
	    if(p.prefix.Contains("nbjet")){
	      if(genfid_b0){
		genm=gen_ll_dressed.M();
		geny=gen_ll_dressed.Rapidity();
		genpt=gen_ll_dressed.Pt();
		gencost=GetCosThetaRecoil(&gen_l0_dressed,&gen_l1_dressed,genfid_b0);
	      }
	    }else{
	      genm=gen_ll_dressed.M();
	      geny=gen_ll_dressed.Rapidity();
	      genpt=gen_ll_dressed.Pt();
	      gencost=GetCosThetaCS(&gen_l0_dressed,&gen_l1_dressed);
	    }
	  }
	}
      }
    }
  }
  Parameter preco=p.Clone();
  Variations recovariations=MakeVariations(preco);

  Parameter pgen=p.Clone();
  ResetRecoWeights(pgen);
  Variations genvariations=MakeVariations(pgen);

  for(auto& [vsuf,genvariation]:genvariations){
    Apply(pgen,vsuf,genvariation);
    Apply(preco,vsuf,recovariations[vsuf]);
    double genweight=pgen.weight;

    double recoweight=0;
    double recom=-1;
    double recoy=-100;
    double recopt=-1;
    double recocost=0;
    if(PassSelection(preco)){
      if(p.prefix.Contains("nbjet")){
	recoweight=preco.weight;
	recom=(*preco.lepton0+*preco.lepton1).M();
	recoy=(*preco.lepton0+*preco.lepton1).Rapidity();
	recopt=(*preco.lepton0+*preco.lepton1).Pt();
	recocost=GetCosThetaRecoil(preco.lepton0,preco.lepton1,&preco.bjets.at(0));
      }else{
	recoweight=preco.weight;
	recom=(*preco.lepton0+*preco.lepton1).M();
	recoy=(*preco.lepton0+*preco.lepton1).Rapidity();
	recopt=(*preco.lepton0+*preco.lepton1).Pt();
	recocost=GetCosThetaCS(preco.lepton0,preco.lepton1);
      }
    }
    //cout<<"vsuf:"<<vsuf<<" genweight:"<<genweight<<" recoweight:"<<recoweight<<" recom:"<<recom<<endl;
    int imbin=GetUnfoldBin(afb_mbinnum,afb_mbin,genm,gencost);
    int jmbin=GetUnfoldBin(afb_mbinnum,afb_mbin,recom,recocost);
    int iybin=GetUnfoldBin(afb_ybinnum,afb_ybin,geny,gencost);
    int jybin=GetUnfoldBin(afb_ybinnum,afb_ybin,recoy,recocost);
    int iptbin=GetUnfoldBin(afb_ptbinnum,afb_ptbin,genpt,gencost);
    int jptbin=GetUnfoldBin(afb_ptbinnum,afb_ptbin,recopt,recocost);

    FillHist(p.prefix+p.hprefix+"response_afbm"+p.suffix+vsuf,imbin,jmbin,recoweight,2*afb_mbinnum,1,2*afb_mbinnum+1,2*afb_mbinnum,1,2*afb_mbinnum+1);
    FillHist(p.prefix+p.hprefix+"response_afbm"+p.suffix+vsuf,imbin,0,genweight-recoweight,2*afb_mbinnum,1,2*afb_mbinnum+1,2*afb_mbinnum,1,2*afb_mbinnum+1);
    FillHist(p.prefix+p.hprefix+"response_afby"+p.suffix+vsuf,iybin,jybin,recoweight,2*afb_ybinnum,1,2*afb_ybinnum+1,2*afb_ybinnum,1,2*afb_ybinnum+1);
    FillHist(p.prefix+p.hprefix+"response_afby"+p.suffix+vsuf,iybin,0,genweight-recoweight,2*afb_ybinnum,1,2*afb_ybinnum+1,2*afb_ybinnum,1,2*afb_ybinnum+1);
    FillHist(p.prefix+p.hprefix+"response_afbpt"+p.suffix+vsuf,iptbin,jptbin,recoweight,2*afb_ptbinnum,1,2*afb_ptbinnum+1,2*afb_ptbinnum,1,2*afb_ptbinnum+1);
    FillHist(p.prefix+p.hprefix+"response_afbpt"+p.suffix+vsuf,iptbin,0,genweight-recoweight,2*afb_ptbinnum,1,2*afb_ptbinnum+1,2*afb_ptbinnum,1,2*afb_ptbinnum+1);
    double massregion[]={52,77,106,280,3000};
    for(int k=0;k<4;k++){
      int ibin=(genm>=massregion[k]&&genm<massregion[k+1]) ? iybin : 0;
      int jbin=(recom>=massregion[k]&&recom<massregion[k+1]) ? jybin : 0;

      //if(vsuf==""&&ibin>0) cout<<"Fill response "<<p.prefix<<" "<<genweight<<endl;
      FillHist(p.prefix+p.hprefix+Form("response_afby_m%d",k)+p.suffix+vsuf,ibin,jbin,recoweight,2*afb_ybinnum,1,2*afb_ybinnum+1,2*afb_ybinnum,1,2*afb_ybinnum+1);
      FillHist(p.prefix+p.hprefix+Form("response_afby_m%d",k)+p.suffix+vsuf,ibin,0,genweight-recoweight,2*afb_ybinnum,1,2*afb_ybinnum+1,2*afb_ybinnum,1,2*afb_ybinnum+1);

      ibin=(genm>=massregion[k]&&genm<massregion[k+1]) ? iptbin : 0;
      jbin=(recom>=massregion[k]&&recom<massregion[k+1]) ? jptbin : 0;
      FillHist(p.prefix+p.hprefix+Form("response_afbpt_m%d",k)+p.suffix+vsuf,ibin,jbin,recoweight,2*afb_ptbinnum,1,2*afb_ptbinnum+1,2*afb_ptbinnum,1,2*afb_ptbinnum+1);
      FillHist(p.prefix+p.hprefix+Form("response_afbpt_m%d",k)+p.suffix+vsuf,ibin,0,genweight-recoweight,2*afb_ptbinnum,1,2*afb_ptbinnum+1,2*afb_ptbinnum,1,2*afb_ptbinnum+1);
    }

    if(vsuf==""){
      FillHist(p.prefix+p.hprefix+"genfid_myptcost_dressed_check"+p.suffix,genm,geny,genpt,gencost,genweight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,grid_costbinnum,(double*)grid_costbin);
      FillHist(p.prefix+p.hprefix+"myptcost_check"+p.suffix,recom,recoy,recopt,recocost,recoweight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,grid_costbinnum,(double*)grid_costbin);
    }
  }
}
AFBAnalyzer::Variations AFBAnalyzer::MakeVariations(const Parameter& p){
  Variations v=SMPAnalyzerCore::MakeVariations(p);
  if(p.variationbits&PDFWeight){
    if(MCSample.Contains("MiNNLO")){
      double pdfreweight=LHAPDF::weightxxQ(genWeight_id1,genWeight_id2,genWeight_X1,genWeight_X2,genWeight_Q,PDFbase,PDFnf4,-1);
      if(!isnormal(pdfreweight)&&pdfreweight!=0) pdfreweight=1.;
      if(pdfreweight>5) pdfreweight=5;
      if(pdfreweight<-5) pdfreweight=-5;
      AddVariationWeight(v,"_nf4",p.default_weight*pdfreweight);
    }
  }
  return v;
}
void AFBAnalyzer::ResetRecoWeights(Parameter& p){
  p.w.prefireweight=1.; p.w.prefireweight_up=1.; p.w.prefireweight_down=1.;
  p.w.z0weight=1.;
  p.w.electronRECOSF=1.;
  p.w.electronRECOSF_sys=Make2DWeights(fEff->GetStructure(p.k.electronRECOSF));
  p.w.electronIDSF=1.;
  p.w.electronIDSF_sys=Make2DWeights(fEff->GetStructure(p.k.electronIDSF));
  p.w.muonTrackingSF=1.;
  p.w.muonTrackingSF_sys=Make2DWeights(fEff->GetStructure(p.k.muonTrackingSF));
  p.w.muonRECOSF=1.;
  p.w.muonRECOSF_sys=Make2DWeights(fEff->GetStructure(p.k.muonRECOSF));
  p.w.muonIDSF=1.;
  p.w.muonIDSF_sys=Make2DWeights(fEff->GetStructure(p.k.muonIDSF));
  p.w.muonISOSF=1.;
  p.w.muonISOSF_sys=Make2DWeights(fEff->GetStructure(p.k.muonISOSF));
  p.w.triggerSF=1.;
  if(p.k.triggerSF.size())
    p.w.triggerSF_sys=Make2DWeights(fEff->GetStructure(p.k.triggerSF[0]));
  p.w.CFSF=1.; p.w.CFSF_up=1.; p.w.CFSF_down=1.;
  p.w.btagSF=1.; p.w.btagSF_hup=1.; p.w.btagSF_hdown=1.; p.w.btagSF_lup=1.; p.w.btagSF_ldown=1.;
  p.w.bchargeSF=1.; p.w.bchargeSF_s0m0=1.; p.w.bchargeSF_s0m1=1.;
}
void AFBAnalyzer::FillHistsSyst(Parameter p,Variations& vs){
  if(!IsSkimmed) return;
  for(auto& [vsuf,variation]:vs){
    Apply(p,vsuf,variation);
    if(!PassSelection(p)) continue;
    TLorentzVector dilepton=*p.lepton0+*p.lepton1;
    double dimass=dilepton.M();
    double dirap=dilepton.Rapidity();
    double dipt=dilepton.Pt();
    FillHistsAFB(p.prefix,p.hprefix,p.suffix+vsuf,(Particle*)p.lepton0,(Particle*)p.lepton1,p.weight);
    //FillHist(p.prefix+p.hprefix+"zpmass"+p.suffix+vsuf,dimass,p.weight,100,600,800);
    FillHist(p.prefix+p.hprefix+"jets"+p.suffix+vsuf,dimass,dirap,dipt,p.jets.size(),p.weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,10,0,10);
    if(p.jets.size()){
      FillHist(p.prefix+p.hprefix+"j0pt"+p.suffix+vsuf,dimass,dirap,dipt,p.jets.at(0).Pt(),p.weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,100,0,200);
    }
    FillHist(p.prefix+p.hprefix+"bjets"+p.suffix+vsuf,dimass,dirap,dipt,p.bjets.size(),p.weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,10,0,10);
    if(p.bjets.size()){
      FillHistsRecoil(p.prefix,p.hprefix,p.suffix+vsuf,(Particle*)p.lepton0,(Particle*)p.lepton1,(Particle*)&p.bjets[0],p.weight);
      //FillHist(p.prefix+p.hprefix+"costhetaRecoil2"+p.suffix+vsuf,dimass,dirap,dipt,GetCosThetaRecoil(p.lepton0,p.lepton1,&p.bjets[0],1),p.weight,afb_mbinnum,(double*)afb_mbin,afb_ybinnum,(double*)afb_ybin,afb_ptbinnum,(double*)afb_ptbin,20,-1,1);
      FillHist(p.prefix+p.hprefix+"b0pt"+p.suffix+vsuf,dimass,dirap,dipt,p.bjets.at(0).Pt(),p.weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,100,0,200);
      FillHist(p.prefix+p.hprefix+"b0charge"+p.suffix+vsuf,dimass,dirap,dipt,p.bjets.at(0).userFloat["AFBCharge"],p.weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,100,-5,5);
      if(vsuf=="")
	FillHist(p.prefix+p.hprefix+"zb0dphi"+p.suffix+vsuf,dilepton.DeltaPhi(p.bjets.at(0)),p.weight,100,-5,5);
    }
    if(vsuf==""||vsuf.Contains("z0weight"))
      FillHist(p.prefix+p.hprefix+"z0"+p.suffix+vsuf,dimass,dirap,dipt,vertex_Z,p.weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,120,-15,15);
    if(vsuf.Contains("PUWeight")){
      FillHist(p.prefix+p.hprefix+"nPV"+p.suffix+vsuf,dimass,dirap,dipt,nPV,p.weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,100,0,100);
      FillHist(p.prefix+p.hprefix+"rho"+p.suffix+vsuf,dimass,dirap,dipt,Rho,p.weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,50,0,50);
      FillHist(p.prefix+p.hprefix+"met"+p.suffix+vsuf,dimass,dirap,dipt,pfMET_Type1_pt,p.weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,100,0,200);
    }
    if(IsDYSample&&p.hprefix==""&&IsNominalRun){
      vector<Gen> gens=GetGens();
      Gen truth_l0=GetGenMatchedLepton(*p.lepton0,gens);
      Gen truth_l1=GetGenMatchedLepton(*p.lepton1,gens);
      if(!truth_l0.IsEmpty()&&!truth_l1.IsEmpty()) 
	FillHistsAFB(p.prefix,"truth_",p.suffix+vsuf,(Particle*)&truth_l0,(Particle*)&truth_l1,p.weight);
      //else cout<<"no matching"<<endl;
    }
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
  /*
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
  */
}

AFBAnalyzer::AFBAnalyzer(){}
AFBAnalyzer::~AFBAnalyzer(){}

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
double AFBAnalyzer::GetCosThetaRecoil(const Particle *p0,const Particle *p1,Particle *b,int mode){
  if(!p0||!p1) return 0.;
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
  int direction=0;
  if(b->InheritsFrom("LHE")||b->InheritsFrom("Gen")){
    if(b->Charge()>0) direction = -1;
    else direction=1;    
  }else if(b->userFloat.find("AFBCharge")!=b->userFloat.end()){
    if(b->userFloat["AFBCharge"]>0) direction = -1;
    else direction=1;
  }
  if(mode==0){
    return direction*((*lm-*lp)*(*b))/((*lm+*lp)*(*b));
  }else{
    TLorentzVector b_m0;
    b_m0.SetPtEtaPhiM(b->Pt(),b->Eta(),b->Phi(),0);
    b_m0*=b->E()/b_m0.E();
    return direction*((*lm-*lp)*(b_m0))/((*lm+*lp)*(b_m0));
  }
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
void AFBAnalyzer::FillHistsAFB(TString pre,TString hpre,TString suf,Particle* l0,Particle* l1,double weight){
  TLorentzVector dilepton=(*l0)+(*l1);
  double dimass=dilepton.M();
  double dirap=dilepton.Rapidity();
  double dipt=dilepton.Pt();

  double cost=GetCosThetaCS(l0,l1);
  FillHist(pre+hpre+"myptcostCS"+suf,dimass,dirap,dipt,cost,weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,grid_costbinnum,(double*)grid_costbin);
  FillHist(pre+hpre+"dimassCS"+suf,dimass,dirap,dipt,cost,weight,afb_mbinnum,(double*)afb_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,grid_costbinnum,(double*)grid_costbin);
  FillHist(pre+hpre+"dirapCS"+suf,dimass,dirap,dipt,cost,weight,grid_mbinnum,(double*)grid_mbin,afb_ybinnum,(double*)afb_ybin,grid_ptbinnum,(double*)grid_ptbin,grid_costbinnum,(double*)grid_costbin);
  FillHist(pre+hpre+"diptCS"+suf,dimass,dirap,dipt,cost,weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,afb_ptbinnum,(double*)afb_ptbin,grid_costbinnum,(double*)grid_costbin);
  FillHist(pre+hpre+"costhetaCS"+suf,dimass,dirap,dipt,cost,weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,afb_costbinnum,(double*)afb_costbin);

  //double h=0.5*pow(dipt/dimass,2)/(1+pow(dipt/dimass,2))*(1-3*cost*cost);
  //double den_weight=0.5*fabs(cost)/pow(1+cost*cost+h,2);
  //double num_weight=0.5*cost*cost/pow(1+cost*cost+h,3);
  //FillHist(pre+hpre+"costhetaCS_den"+suf,dimass,dirap,dipt,cost,weight*den_weight,afb_mbinnum,(double*)afb_mbin,afb_ybinnum,(double*)afb_ybin,afb_ptbinnum,(double*)afb_ptbin,20,-1,1);
  //FillHist(pre+hpre+"costhetaCS_num"+suf,dimass,dirap,dipt,cost,weight*num_weight,afb_mbinnum,(double*)afb_mbin,afb_ybinnum,(double*)afb_ybin,afb_ptbinnum,(double*)afb_ptbin,20,-1,1);

  FillHist(pre+hpre+"l0pt"+suf,dimass,dirap,dipt,l0->Pt(),weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,lptbinnum,(double*)lptbin);
  FillHist(pre+hpre+"l1pt"+suf,dimass,dirap,dipt,l1->Pt(),weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,lptbinnum,(double*)lptbin);
  FillHist(pre+hpre+"lpt"+suf,dimass,dirap,dipt,l0->Pt(),weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,lptbinnum,(double*)lptbin);
  FillHist(pre+hpre+"lpt"+suf,dimass,dirap,dipt,l1->Pt(),weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,lptbinnum,(double*)lptbin);
  
  FillHist(pre+hpre+"l0eta"+suf,dimass,dirap,dipt,l0->Eta(),weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,60,-3,3);
  FillHist(pre+hpre+"l1eta"+suf,dimass,dirap,dipt,l1->Eta(),weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,60,-3,3);
  FillHist(pre+hpre+"leta"+suf,dimass,dirap,dipt,l0->Eta(),weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,60,-3,3);
  FillHist(pre+hpre+"leta"+suf,dimass,dirap,dipt,l1->Eta(),weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,60,-3,3);
}
void AFBAnalyzer::FillHistsRecoil(TString pre,TString hpre,TString suf,Particle* l0,Particle* l1,Particle* b,double weight){
  TLorentzVector dilepton=(*l0)+(*l1);
  double dimass=dilepton.M();
  double dirap=dilepton.Rapidity();
  double dipt=dilepton.Pt();

  double cost=GetCosThetaRecoil(l0,l1,b);
  FillHist(pre+hpre+"myptcostRecoil"+suf,dimass,dirap,dipt,cost,weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,grid_costbinnum,(double*)grid_costbin);
  FillHist(pre+hpre+"dimassRecoil"+suf,dimass,dirap,dipt,cost,weight,afb_mbinnum,(double*)afb_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,grid_costbinnum,(double*)grid_costbin);
  FillHist(pre+hpre+"dirapRecoil"+suf,dimass,dirap,dipt,cost,weight,grid_mbinnum,(double*)grid_mbin,afb_ybinnum,(double*)afb_ybin,grid_ptbinnum,(double*)grid_ptbin,grid_costbinnum,(double*)grid_costbin);
  FillHist(pre+hpre+"diptRecoil"+suf,dimass,dirap,dipt,cost,weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,afb_ptbinnum,(double*)afb_ptbin,grid_costbinnum,(double*)grid_costbin);
  FillHist(pre+hpre+"costhetaRecoil"+suf,dimass,dirap,dipt,cost,weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,afb_costbinnum,(double*)afb_costbin);

}
void AFBAnalyzer::FillHardHists(TString pre,TString suf,const Gen& genparton0,const Gen& genparton1,const Gen& genhardl0,const Gen& genhardl1,const Gen& genhardj0,double w){
  Gen genhardl=genhardl0.PID()>0?genhardl0:genhardl1;
  TLorentzVector genZ=genhardl0+genhardl1;
  TLorentzVector genpp=genparton0+genparton1;

  FillHist(pre+"cos_l_p0"+suf,cos(genhardl.Angle(genparton0.Vect())),w,100,-1,1);
  FillHist(pre+"cos_l_p0_asym"+suf,cos(genhardl.Angle(genparton0.Vect())),w/2,100,-1,1);
  FillHist(pre+"cos_l_p0_asym"+suf,-1.*cos(genhardl.Angle(genparton0.Vect())),-w/2,100,-1,1);
  if(cos(genhardl.Angle(genparton0.Vect()))>0) FillHist(pre+"cos_l_p0_forward"+suf,genZ.M(),w,afb_mbinnum,(double*)afb_mbin);
  else FillHist(pre+"cos_l_p0_backward"+suf,genZ.M(),w,afb_mbinnum,(double*)afb_mbin);

  FillHist(pre+"cos_l_p1"+suf,cos(genhardl.Angle(genparton1.Vect())),w,100,-1,1);
  FillHist(pre+"cos_l_p1_asym"+suf,cos(genhardl.Angle(genparton1.Vect())),w/2,100,-1,1);
  FillHist(pre+"cos_l_p1_asym"+suf,-1.*cos(genhardl.Angle(genparton1.Vect())),-w/2,100,-1,1);
  if(cos(genhardl.Angle(genparton1.Vect()))>0) FillHist(pre+"cos_l_p1_forward"+suf,genZ.M(),w,afb_mbinnum,(double*)afb_mbin);
  else FillHist(pre+"cos_l_p1_backward"+suf,genZ.M(),w,afb_mbinnum,(double*)afb_mbin);
  
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
    if(cos(genhardl.Angle(genhardj0.Vect()))>0) FillHist(pre+"cos_l_j0_forward"+suf,genZ.M(),w,afb_mbinnum,(double*)afb_mbin);
    else FillHist(pre+"cos_l_j0_backward"+suf,genZ.M(),w,afb_mbinnum,(double*)afb_mbin);

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
