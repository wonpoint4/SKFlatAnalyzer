#include "LTAnalyzer.h"

LTAnalyzer::LTAnalyzer(){
}
LTAnalyzer::~LTAnalyzer(){
}
void LTAnalyzer::initializeAnalyzer(){
  SMPAnalyzerCore::initializeAnalyzer(); //setup zpt roc z0
  fChain->SetBranchStatus("pfMET_*",false);
  fChain->SetBranchStatus("pfMET_Type1_pt",true);
  fChain->SetBranchStatus("fatjet_*",false);
}
void LTAnalyzer::executeEvent(){
  Parameter p;
  if(IsDYSample){
    if(abs(lhe_l0.ID())==11&&abs(lhe_l1.ID())==11){
      p=MakeParameter("ee");
    }else if(abs(lhe_l0.ID())==13&&abs(lhe_l1.ID())==13){
      p=MakeParameter("mm");
    }
  }
  if(p.channel!=""){
    int njet=0;
    vector<Jet> jets=GetAllJets();
    for(auto& jet:jets){
      if(!jet.IsGenMatched()) continue;
      if(jet.DeltaR(gen_l0_dressed)<0.4) continue;
      if(jet.DeltaR(gen_l1_dressed)<0.4) continue;
      if(jet.Pt()<20) continue;
      njet++;
    }

    TLorentzVector dilepton=gen_l0+gen_l1;
    double m=dilepton.M();
    double pt=dilepton.Pt();
    pair<double,double> costphi=GetCostAndPhiCS(&gen_l0,&gen_l1);
    double cost=costphi.first;
    cost=fabs(cost);
    double phi=costphi.second;
    phi=TMath::Pi()/2-fabs(TMath::Pi()/2-fabs(phi));
    if(81<m&&m<101){
      //FillHist(p.prefix+p.hprefix+"gen_ptcostphi",pt,cost,phi,p.w.lumiweight,nptbin,ptbins,ncostbin,costbins,nphibin,phibins);
      int ptbin=0;
      for(int i=0;i<nptbin+1;i++)
	if(pt>ptbins[i]) 
	  ptbin=i+1;
      FillHist(p.prefix+p.hprefix+Form("gen_pt%d_costphi",ptbin),cost,phi,p.w.lumiweight,ncostbin,costbins,nphibin,phibins);
      if(njet<2){
	FillHist(p.prefix+p.hprefix+Form("jet0/gen_pt%d_costphi",ptbin),cost,phi,p.w.lumiweight,ncostbin,costbins,nphibin,phibins);
      }else{
	FillHist(p.prefix+p.hprefix+Form("jet1/gen_pt%d_costphi",ptbin),cost,phi,p.w.lumiweight,ncostbin,costbins,nphibin,phibins);
      }
    }
    int igen=GetUnfoldBin(njet,m,pt,cost,phi);
    FillHist(p.prefix+p.hprefix+"gen"+p.suffix,igen,p.w.lumiweight,nresponsebin,0,nresponsebin);
  }
  //////// nominal channels //////////
  if(!IsDATA||DataStream.Contains("DoubleMuon")){
    executeEventWithParameter(MakeParameter("mm")); 
  }
  if(!IsDATA||DataStream.Contains("DoubleEG")||DataStream.Contains("EGamma")){
    executeEventWithParameter(MakeParameter("ee"));
  }    
}
void LTAnalyzer::ResetRecoWeights(Parameter& p){
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
  p.w.btagSF=1.; p.w.btagSF_hup=1.; p.w.btagSF_hdown=1.; p.w.btagSF_lup=1.; p.w.btagSF_ldown=1.;
}
int LTAnalyzer::GetUnfoldBin(int njet,double mass,double pt,double cost,double phi){
  if(mass<81||mass>101) return -1;
  int ijet=njet<2 ? 0 : 1;
  int ipt=TMath::BinarySearch(nptbin,ptbins,pt);
  if(ipt<0) ipt=0;
  cost=fabs(cost);
  int icost=int(cost*ncostbin);
  if(icost<0) icost=0;
  if(icost>=ncostbin) icost=ncostbin-1;
  phi=TMath::Pi()/2-fabs(TMath::Pi()/2-fabs(phi));
  int iphi=int(phi/TMath::Pi()*2*nphibin);
  if(iphi<0) iphi=0;
  if(iphi>=nphibin) iphi=nphibin-1;

  return nphibin*(ncostbin*(nptbin*ijet+ipt)+icost)+iphi;
}
void LTAnalyzer::executeEventWithParameter(Parameter& p){
  SMPAnalyzerCore::executeEventWithParameter(p);
  if(!IsDYSample) return;
  TLorentzVector gen_dilepton=gen_l0+gen_l1;
  double gen_mass=gen_dilepton.M();
  double gen_pt=gen_dilepton.Pt();
  pair<double,double> gen_costphi=GetCostAndPhiCS(&gen_l0,&gen_l1);
  double gen_cost=gen_costphi.first;
  gen_cost=fabs(gen_cost);
  double gen_phi=gen_costphi.second;
  gen_phi=TMath::Pi()/2-fabs(TMath::Pi()/2-fabs(gen_phi));
  int gen_njet=0;
  vector<Jet> alljets=GetAllJets();
  for(auto& jet:alljets){
    if(!jet.IsGenMatched()) continue;
    if(jet.DeltaR(gen_l0_dressed)<0.4) continue;
    if(jet.DeltaR(gen_l1_dressed)<0.4) continue;
    if(jet.Pt()<20) continue;
    gen_njet++;
  }

  int igen=-1;
  if(p.channel=="ee"){
    if(abs(lhe_l0.ID())==11&&abs(lhe_l1.ID())==11)
      igen=GetUnfoldBin(gen_njet,gen_mass,gen_pt,gen_cost,gen_phi);
  }else if(p.channel=="mm"){
    if(abs(lhe_l0.ID())==13&&abs(lhe_l1.ID())==13)
      igen=GetUnfoldBin(gen_njet,gen_mass,gen_pt,gen_cost,gen_phi);
  }

  Parameter pgen=p;
  ResetRecoWeights(pgen);
  EvalWeights(pgen);

  TLorentzVector dilepton;
  pair<double,double> costphi=make_pair(0.,0.);
  if(p.lepton0&&p.lepton1){
    dilepton=*p.lepton0+*p.lepton1;
    costphi=GetCostAndPhiCS(p.lepton0,p.lepton1);
  }
  double mass=dilepton.M();
  double pt=dilepton.Pt();
  double cost=costphi.first;
  cost=fabs(cost);
  double phi=costphi.second;
  phi=TMath::Pi()/2-fabs(TMath::Pi()/2-fabs(phi));
  int njet=0;
  vector<Jet> jets=GetJets("tightLepVeto",20,5.0);
  for(int i=0,n=jets.size();i<n;i++){
    if(p.lepton0&&jets[i].DeltaR(*p.lepton0)<0.4) continue;
    if(p.lepton1&&jets[i].DeltaR(*p.lepton1)<0.4) continue;
    njet++;
  }
  for(auto [wname,genweight]:pgen.weightmap){
    double recoweight=0;
    int ireco=-1;
    if(p.weightmap.find(wname)!=p.weightmap.end()){
      recoweight=p.weightmap[wname];
      ireco=GetUnfoldBin(njet,mass,pt,cost,phi);
    }
    if(igen>=0){
      FillHist(p.prefix+p.hprefix+"response"+p.suffix+wname,igen,ireco,recoweight,nresponsebin,0,nresponsebin,nresponsebin,0,nresponsebin);
      FillHist(p.prefix+p.hprefix+"response"+p.suffix+wname,igen,-1,genweight-recoweight,nresponsebin,0,nresponsebin,nresponsebin,0,nresponsebin);
    }else if(ireco>=0){
      FillHist(p.prefix+p.hprefix+"response"+p.suffix+wname,igen,ireco,recoweight,nresponsebin,0,nresponsebin,nresponsebin,0,nresponsebin);      
    }
  } 
}


pair<double,double> LTAnalyzer::GetCostAndPhiCS(Particle* l0,Particle* l1){
  if(!l0||!l1) return make_pair(0.,0.);
  const TLorentzVector *lm,*lp;
  if(l0->Charge()<0&&l1->Charge()>0){
    lm=l0;
    lp=l1;
  }else if(l0->Charge()>0&&l1->Charge()<0){
    lm=l1;
    lp=l0;
  }else if(strcmp(l0->ClassName(),"LHE")==0){
    if(((LHE*)l0)->ID()>0&&((LHE*)l1)->ID()<0){
      lm=l0;
      lp=l1;
    }else if(((LHE*)l0)->ID()<0&&((LHE*)l1)->ID()>0){
      lm=l1;
      lp=l0;
    }else{
      if(gRandom->Rndm()<0.5){
        lm=l0;
        lp=l1;
      }else{
        lm=l1;
        lp=l0;
      }
    }
  }else{
    if(gRandom->Rndm()<0.5){
      lm=l0;
      lp=l1;
    }else{
      lm=l1;
      lp=l0;
    }
  }
  TLorentzVector dilepton=*lm+*lp;
  TLorentzVector lmcs=*lm,lpcs=*lp,p0(0,0,1,1),p1(0,0,-1,1);
  TVector3 b=dilepton.BoostVector();
  lmcs.Boost(-b);
  p0.Boost(-b);
  p1.Boost(-b);
  TVector3 v0=p0.Vect().Unit(),v1=p1.Vect().Unit();
  TVector3 z=(v0-v1).Unit();
  TVector3 x=(-v0-v1).Unit();
  if(dilepton.Pz()<0){
    z=-z;
    x=-x;
  }
  TVector3 y=z.Cross(x);
  TRotation rot;
  rot.RotateAxes(x,y,z);
  //x.Print();y.Print();z.Print();
  rot.Invert();
  lmcs.Transform(rot);
  
  /*
  p0.Transform(rot);
  p1.Transform(rot);
  double cost=0;
  {
    double lmpp=(lm->E()+lm->Pz())/sqrt(2);
    double lmpm=(lm->E()-lm->Pz())/sqrt(2);
    double lppp=(lp->E()+lp->Pz())/sqrt(2);
    double lppm=(lp->E()-lp->Pz())/sqrt(2);
    double dimass=dilepton.M();
    double dipt=dilepton.Pt();
    int direction=dilepton.Pz()>0?1:-1;
    cost=direction*2*(lmpp*lppm-lmpm*lppp)/sqrt(dimass*dimass*(dimass*dimass+dipt*dipt));
  }
  cout<<"cost: "<<lmcs.CosTheta()<<" "<<cost<<endl;
  cout<<"phi: "<<lmcs.Phi()<<endl;
  cout<<"p0 phi: "<<p0.Phi()<<" p1 phi: "<<p1.Phi()<<endl;
  */
  return make_pair(lmcs.CosTheta(),lmcs.Phi());
}


SMPAnalyzerCore::Parameter LTAnalyzer::MakeParameter(TString key,TString option){
  Parameter p=SMPAnalyzerCore::MakeParameter(key,option);
  p.weightbit|=EfficiencyWeight;
  return p;
}
void LTAnalyzer::EvalWeights(Parameter& p){
  p.weightmap[""]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
  if(!IsDATA&&p.suffix==""&&!p.hprefix.Contains("ss_")){
    p.weightmap["_noweight"]=p.w.lumiweight;
    p.weightmap["_noPUweight"]=p.w.lumiweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
    p.weightmap["_noprefireweight"]=p.w.lumiweight*p.w.PUweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
    p.weightmap["_nozptweight"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
    p.weightmap["_z0weight"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*p.w.z0weight;
    p.weightmap["_noweakweight"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
    /*
    for(int j=0,nj=fEff->nreplica;j<nj;j++){
      double electronRECOSF=p.w.electronRECOSF_sys.size() ? p.w.electronRECOSF_sys[0][j] : 1.;
      double electronIDSF=p.w.electronIDSF_sys.size() ? p.w.electronIDSF_sys[0][j] : 1.;
      double muonTrackingSF=p.w.muonTrackingSF_sys.size() ? p.w.muonTrackingSF_sys[0][j] : 1.;
      double muonRECOSF=p.w.muonRECOSF_sys.size() ? p.w.muonRECOSF_sys[0][j] : 1.;
      double muonIDSF=p.w.muonIDSF_sys.size() ? p.w.muonIDSF_sys[0][j] : 1.;
      double triggerSF=p.w.triggerSF_sys.size() ? p.w.triggerSF_sys[0][j] : 1.;
      p.weightmap[Form("_efficiencySF_stat%d",j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*electronRECOSF*electronIDSF*muonTrackingSF*muonRECOSF*muonIDSF*p.w.muonISOSF*triggerSF*p.w.CFSF;
    }

    p.weightmap["_noelectronRECOSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
    for(int i=1,ni=p.w.electronRECOSF_sys.size();i<ni;i++){
      for(int j=0,nj=p.w.electronRECOSF_sys[i].size();j<nj;j++){
	p.weightmap[Form("_electronRECOSF_s%d_m%d",i,j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF_sys[i][j]*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
      }
    }	
    
    p.weightmap["_noelectronIDSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
    for(int i=1,ni=p.w.electronIDSF_sys.size();i<ni;i++){
      for(int j=0,nj=p.w.electronIDSF_sys[i].size();j<nj;j++){
	p.weightmap[Form("_electronIDSF_s%d_m%d",i,j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF_sys[i][j]*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
      }
    }	

    p.weightmap["_nomuonTrackingSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
    for(int i=1,ni=p.w.muonTrackingSF_sys.size();i<ni;i++){
      for(int j=0,nj=p.w.muonTrackingSF_sys[i].size();j<nj;j++){
	p.weightmap[Form("_muonTrackingSF_s%d_m%d",i,j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF_sys[i][j]*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
      }
    }	

    p.weightmap["_nomuonRECOSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
    for(int i=1,ni=p.w.muonRECOSF_sys.size();i<ni;i++){
      for(int j=0,nj=p.w.muonRECOSF_sys[i].size();j<nj;j++){
	p.weightmap[Form("_muonRECOSF_s%d_m%d",i,j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF_sys[i][j]*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
      }
    }	

    p.weightmap["_nomuonIDSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
    for(int i=1,ni=p.w.muonIDSF_sys.size();i<ni;i++){
      for(int j=0,nj=p.w.muonIDSF_sys[i].size();j<nj;j++){
	p.weightmap[Form("_muonIDSF_s%d_m%d",i,j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF_sys[i][j]*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
      }
    }	
    
    p.weightmap["_notriggerSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.CFSF;
    for(int i=1,ni=p.w.triggerSF_sys.size();i<ni;i++){
      for(int j=0,nj=p.w.triggerSF_sys[i].size();j<nj;j++){
	p.weightmap[Form("_triggerSF_s%d_m%d",i,j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF_sys[i][j]*p.w.CFSF;
      }
    }
        
    p.weightmap["_noefficiencySF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.CFSF;
    */
    p.weightmap["_noCFSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF;
  }
}
void LTAnalyzer::FillHists(Parameter& p){
  TLorentzVector dilepton=(*p.lepton0)+(*p.lepton1);
  double mass=dilepton.M();
  double pt=dilepton.Pt();
  pair<double,double> costphi=GetCostAndPhiCS(p.lepton0,p.lepton1);
  double cost=costphi.first;
  double phi=costphi.second;
  int njet=0;
  vector<Jet> jets=GetJets("tightLepVeto",20,5.0);
  for(int i=0,n=jets.size();i<n;i++){
    if(jets[i].DeltaR(*p.lepton0)<0.4) continue;
    if(jets[i].DeltaR(*p.lepton1)<0.4) continue;
    njet++;
  }
  int ireco=GetUnfoldBin(njet,mass,pt,cost,phi);
  FillHist(p.prefix+p.hprefix+"reco"+p.suffix,ireco,p.weightmap,nresponsebin,0,nresponsebin);
}
