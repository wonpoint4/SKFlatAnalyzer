#include "SMPAnalyzerCore.h"

SMPAnalyzerCore::SMPAnalyzerCore(){}
SMPAnalyzerCore::~SMPAnalyzerCore(){
  if(roc) delete roc;
  if(rocele) delete rocele;
  for(auto& [_,h]:fRoccorResidual){
    if(h) delete h;
  }
  for(std::map< TString, TH4D* >::iterator mapit = maphist_TH4D.begin(); mapit!=maphist_TH4D.end(); mapit++){
    delete mapit->second;
  }
  maphist_TH4D.clear();
  DeleteEfficiency();
  DeleteZptWeight();
  DeleteCFRate();
  DeleteFakeRate();
  DeleteL1PrefiringWeight();
}

void SMPAnalyzerCore::initializeAnalyzer(){
  if(MaxEvent>0) reductionweight=1.*fChain->GetEntries()/MaxEvent;
  else reductionweight=1.;
  SetupEfficiency();
  SetupRoccoR();
  SetupCFRate();
  SetupFakeRate();
  SetupL1PrefiringWeight();
  IsDYSample=false;
  IsTTSample=false;
  IsTTLLSample=false;
  if(MCSample.Contains("DYJets")||MCSample.Contains("ZToEE")||MCSample.Contains("ZToMuMu")||MCSample.Contains(TRegexp("DY[0-9]Jets"))) IsDYSample=true;
  if(IsDYSample) SetupZptWeight();
  if(MCSample.Contains(TRegexp("TT[LJ][LJ]"))) IsTTSample=true;
  if(MCSample.Contains("TTLL")) IsTTLLSample=true;
  mcCorr->SetJetTaggingParameters({JetTagging::Parameters(JetTagging::DeepJet,JetTagging::Tight,JetTagging::incl,JetTagging::comb),
	JetTagging::Parameters(JetTagging::DeepJet,JetTagging::Medium,JetTagging::incl,JetTagging::comb),
	JetTagging::Parameters(JetTagging::DeepJet,JetTagging::Tight,JetTagging::incl,JetTagging::mujets),
	JetTagging::Parameters(JetTagging::DeepCSV,JetTagging::Medium,JetTagging::incl,JetTagging::comb),
	JetTagging::Parameters(JetTagging::DeepCSV,JetTagging::Tight,JetTagging::incl,JetTagging::comb)});
}
void SMPAnalyzerCore::beginEvent(){
  _event=GetEvent();
  if(!IsDATA){
    lhes=GetLHEs();
    gens=GetGens();
    if(IsDYSample||MCSample.Contains("GamGamToLL")||MCSample.Contains("TTLL")){
      GetAFBLHEParticles(lhes,lhe_p0,lhe_p1,lhe_l0,lhe_l1,lhe_j0);
      GetAFBGenParticles(gens,gen_p0,gen_p1,gen_l0,gen_l1,3);
      GetAFBGenParticles(gens,gen_p0,gen_p1,gen_l0_dressed,gen_l1_dressed,1);
      GetAFBGenParticles(gens,gen_p0,gen_p1,gen_l0_bare,gen_l1_bare,0);
    }
  }
}
void SMPAnalyzerCore::executeEventWithParameter(Parameter& p){
  p.SetLeptons();
  
  ///////////////// efficiency scale factors ///////////////////
  EvalIDSF(p);
  EvalTriggerSF(p);

  ///////////// define variations ///////////////
  EvalDefaultWeight(p);
  Variations v=MakeVariations(p);

  ////// Fill histograms //////////
  FillHistsSyst(p,v);
}
void SMPAnalyzerCore::EvalIDSF(Parameter& p){
  //p.doublemap["muontrackingSF"]=1.;
  if(!IsDATA){
    if(p.variationbits&EfficiencyWeight){
      p.w.electronRECOSF_sys=Make2DWeights(fEff->GetStructure(p.k.electronRECOSF));
      p.w.electronIDSF_sys=Make2DWeights(fEff->GetStructure(p.k.electronIDSF));
      p.w.muonTrackingSF_sys=Make2DWeights(fEff->GetStructure(p.k.muonTrackingSF));
      p.w.muonRECOSF_sys=Make2DWeights(fEff->GetStructure(p.k.muonRECOSF));
      p.w.muonIDSF_sys=Make2DWeights(fEff->GetStructure(p.k.muonIDSF));
      p.w.muonISOSF_sys=Make2DWeights(fEff->GetStructure(p.k.muonISOSF));
    }
    for(const Lepton* lepton:p.leptons){
      if(!lepton) continue;
      if(lepton->LeptonFlavour()==Lepton::ELECTRON){
	const Electron* electron=(const Electron*)lepton;
	p.w.electronRECOSF*=fEff->GetEfficiencySF(p.k.electronRECOSF,electron,0,0);
	if(p.variationbits&EfficiencyWeight){
	  int nset=p.w.electronRECOSF_sys.size();
	  for(int s=0;s<nset;s++){
	    int nmem=p.w.electronRECOSF_sys[s].size();
	    for(int m=0;m<nmem;m++){
	      p.w.electronRECOSF_sys[s][m]*=fEff->GetEfficiencySF(p.k.electronRECOSF,electron,s,m);
	    }
	  }
	}
	p.w.electronIDSF*=fEff->GetEfficiencySF(p.k.electronIDSF,electron,0,0);
	p.w.electronIDSF*=fEff->GetEfficiencySF(p.k.electronIDSF2,electron,0,0);
	if(p.variationbits&EfficiencyWeight){
	  int nset=p.w.electronIDSF_sys.size();
	  for(int s=0;s<nset;s++){
	    int nmem=p.w.electronIDSF_sys[s].size();
	    for(int m=0;m<nmem;m++){
	      p.w.electronIDSF_sys[s][m]*=fEff->GetEfficiencySF(p.k.electronIDSF,electron,s,m);
	      // electronIDSF2 is for the selective charge ID
	      if(fEff->Get(p.k.electronIDSF2)){
		if(s<(int)fEff->Get(p.k.electronIDSF2)->fDataPlus.size())
		  p.w.electronIDSF_sys[s][m]*=fEff->GetEfficiencySF(p.k.electronIDSF2,electron,s,m);
		else
		  p.w.electronIDSF_sys[s][m]*=fEff->GetEfficiencySF(p.k.electronIDSF2,electron,0,0);
	      }
	    }
	  }
	}
      }else if(lepton->LeptonFlavour()==Lepton::MUON){
	const Muon* muon=(const Muon*)lepton;
	//p.doublemap["muontrackingSF"]*=GetMuonTrackingSF(muon.Eta());
	p.w.muonTrackingSF*=fEff->GetEfficiencySF(p.k.muonTrackingSF,muon,0,0);
	if(p.variationbits&EfficiencyWeight){
	  int nset=p.w.muonTrackingSF_sys.size();
	  for(int s=0;s<nset;s++){
	    int nmem=p.w.muonTrackingSF_sys[s].size();
	    for(int m=0;m<nmem;m++){
	      p.w.muonTrackingSF_sys[s][m]*=fEff->GetEfficiencySF(p.k.muonTrackingSF,muon,s,m);
	    }
	  }
	}
	p.w.muonRECOSF*=fEff->GetEfficiencySF(p.k.muonRECOSF,muon,0,0);
	if(p.variationbits&EfficiencyWeight){
	  int nset=p.w.muonRECOSF_sys.size();
	  for(int s=0;s<nset;s++){
	    int nmem=p.w.muonRECOSF_sys[s].size();
	    for(int m=0;m<nmem;m++){
	      p.w.muonRECOSF_sys[s][m]*=fEff->GetEfficiencySF(p.k.muonRECOSF,muon,s,m);
	    }
	  }
	}
	p.w.muonIDSF*=fEff->GetEfficiencySF(p.k.muonIDSF,muon,0,0);
	if(p.variationbits&EfficiencyWeight){
	  int nset=p.w.muonIDSF_sys.size();
	  for(int s=0;s<nset;s++){
	    int nmem=p.w.muonIDSF_sys[s].size();
	    for(int m=0;m<nmem;m++){
	      p.w.muonIDSF_sys[s][m]*=fEff->GetEfficiencySF(p.k.muonIDSF,muon,s,m);
	    }
	  }
	}
	p.w.muonISOSF*=fEff->GetEfficiencySF(p.k.muonISOSF,muon,0,0);
	if(p.variationbits&EfficiencyWeight){
	  int nset=p.w.muonISOSF_sys.size();
	  for(int s=0;s<nset;s++){
	    int nmem=p.w.muonISOSF_sys[s].size();
	    for(int m=0;m<nmem;m++){
	      p.w.muonISOSF_sys[s][m]*=fEff->GetEfficiencySF(p.k.muonISOSF,muon,s,m);
	    }
	  }
	}
      }
    }
  }
  p.w.CFSF=GetCFSF(p,0);
  p.w.CFSF_up=GetCFSF(p,1);
  p.w.CFSF_down=GetCFSF(p,-1);
}
void SMPAnalyzerCore::EvalTriggerSF(Parameter& p){
  if(!IsDATA){
    vector<Lepton*> triggerables;
    vector<Lepton*> triggerables_mode1=p.leptons;
    if(p.k.triggerSF.size()){
      if(p.variationbits&EfficiencyWeight){
	p.w.triggerSF_sys=Make2DWeights(fEff->GetStructure(p.k.triggerSF[0]));
      }
      if(p.k.triggerSF.at(0).Contains("Mu")) triggerables=MakeLeptonPointerVector(p.muons);
      else if(p.k.triggerSF.at(0).Contains("Ele")) triggerables=MakeLeptonPointerVector(p.electrons);
    }
    if(p.k.triggerSF.size()==1){
      p.w.triggerSF*=GetLeptonTriggerSF(p.k.triggerSF[0],triggerables,0,0);
      p.w.triggerSF_mode1*=GetLeptonTriggerSF(p.k.triggerSF[0],triggerables_mode1,0,0);
      p.w.triggerSF_interpolation*=GetLeptonTriggerSF(p.k.triggerSF[0],triggerables,0,0,"interpolation");
      if(p.variationbits&EfficiencyWeight){
	int nset=p.w.triggerSF_sys.size();
	for(int s=0;s<nset;s++){
	  int nmem=p.w.triggerSF_sys[s].size();
	  for(int m=0;m<nmem;m++){
	    p.w.triggerSF_sys[s][m]*=GetLeptonTriggerSF(p.k.triggerSF[0],triggerables,s,m);
	  }
	}
      }
    }else if(p.k.triggerSF.size()==2){
      if(GetPtThreshold(p.k.triggerSF[0])<GetPtThreshold(p.k.triggerSF[1])){
	p.w.triggerSF*=GetLeptonTriggerORSF(p,triggerables,0,0);
	p.w.triggerSF_mode1*=GetLeptonTriggerORSF(p,triggerables_mode1,0,0);
	p.w.triggerSF_interpolation*=GetLeptonTriggerORSF(p,triggerables,0,0,"interpolation");
	if(p.variationbits&EfficiencyWeight){
	  int nset=p.w.triggerSF_sys.size();
	  for(int s=0;s<nset;s++){
	    int nmem=p.w.triggerSF_sys[s].size();
	    for(int m=0;m<nmem;m++){
	      p.w.triggerSF_sys[s][m]*=GetLeptonTriggerORSF(p,triggerables,s,m);
	    }
	  }
	}
      }else{
	p.w.triggerSF*=GetDileptonTriggerSF(p.k.triggerSF[0],p.k.triggerSF[1],p.k.DZSF,triggerables,0,0);
	p.w.triggerSF_mode1*=GetDileptonTriggerSF(p.k.triggerSF[0],p.k.triggerSF[1],p.k.DZSF,triggerables_mode1,0,0);
	p.w.triggerSF_interpolation*=GetDileptonTriggerSF(p.k.triggerSF[0],p.k.triggerSF[1],p.k.DZSF,triggerables_mode1,0,0,"interpolation");
	if(p.variationbits&EfficiencyWeight){
	  int nset=p.w.triggerSF_sys.size();
	  for(int s=0;s<nset;s++){
	    int nmem=p.w.triggerSF_sys[s].size();
	    for(int m=0;m<nmem;m++){
	      p.w.triggerSF_sys[s][m]*=GetDileptonTriggerSF(p.k.triggerSF[0],p.k.triggerSF[1],p.k.DZSF,triggerables,s,m);
	    }
	  }
	}
      }
    }
  }
}
void SMPAnalyzerCore::EvalDefaultWeight(Parameter& p){
  p.default_weight=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.CFSF*p.w.btagSF*p.w.bchargeSF*p.w.zptweight*p.w.weakweight*p.w.topptweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonTrackingSF*p.w.muonRECOSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF;
  p.weight=p.default_weight;
}
SMPAnalyzerCore::Variations SMPAnalyzerCore::MakeVariations(const Parameter& p){
  Variations v;
  if(p.variationbits&NominalWeight){
    AddVariationWeight(v,"",p.default_weight);
  }
  if(p.variationbits&SystematicWeight){
    if(!IsDATA&&!p.hprefix.Contains("ss_")){
      EvalVariationsPUweight(p,v);
      EvalVariationsPrefireweight(p,v);
      EvalVariationsBtag(p,v);
      EvalVariationsBcharge(p,v);
      EvalVariationsEtc(p,v);
    }
    if(!IsDATA&&p.hprefix.Contains("ss_")){
      EvalVariationsCF(p,v);
    }
    if(!p.hprefix.Contains("ss_")){
      EvalVariationsJetCorrection(p,v);
    }
  }
  if(p.variationbits&EfficiencyWeight){
    if(!IsDATA&&!p.hprefix.Contains("ss_")){
      EvalVariationsEfficiency(p,v);
    }
  }
  if(p.variationbits&PDFWeight){
    if(!IsDATA&&p.hprefix==""){
      EvalVariationsPDF(p,v);
    }
  }
  if(p.variationbits&LeptonCorrection){
    if(p.hprefix==""){
      if(p.channel=="mm"||p.channel=="mu"){
	EvalVariationsMuonMomentum(p,v);
      }else if(p.channel=="ee"||p.channel=="el"){
	EvalVariationsElectronEnergy(p,v);
      }
    }
  }
  return v;
}

void SMPAnalyzerCore::EvalVariationsPUweight(const Parameter& p,Variations& v){
  AddVariationWeight(v,"_noPUweight",p.default_weight/p.w.PUweight);
  AddVariationWeight(v,"_PUweight_down",p.default_weight/p.w.PUweight*p.w.PUweight_down);
  AddVariationWeight(v,"_PUweight_up",p.default_weight/p.w.PUweight*p.w.PUweight_up);
}
void SMPAnalyzerCore::EvalVariationsPrefireweight(const Parameter& p,Variations& v){
  AddVariationWeight(v,"_noprefireweight",p.default_weight/p.w.prefireweight);
  AddVariationWeight(v,"_prefireweight_down",p.default_weight/p.w.prefireweight*p.w.prefireweight_down);
  AddVariationWeight(v,"_prefireweight_up",p.default_weight/p.w.prefireweight*p.w.prefireweight_up);
}
void SMPAnalyzerCore::EvalVariationsCF(const Parameter& p,Variations& v){
  AddVariationWeight(v,"_noCFSF",p.default_weight/p.w.CFSF);
  AddVariationWeight(v,"_CFSF_down",p.default_weight/p.w.CFSF*p.w.CFSF_down);
  AddVariationWeight(v,"_CFSF_up",p.default_weight/p.w.CFSF*p.w.CFSF_up);
}
void SMPAnalyzerCore::EvalVariationsBtag(const Parameter& p,Variations& v){
  AddVariationWeight(v,"_nobtagSF",p.default_weight/p.w.btagSF);
  AddVariationWeight(v,"_btagSF_hdown",p.default_weight/p.w.btagSF*p.w.btagSF_hdown);
  AddVariationWeight(v,"_btagSF_hup",p.default_weight/p.w.btagSF*p.w.btagSF_hup);
  AddVariationWeight(v,"_btagSF_hcorr",p.default_weight/p.w.btagSF*p.w.btagSF_hcorr);
  AddVariationWeight(v,"_btagSF_huncorr",p.default_weight/p.w.btagSF*p.w.btagSF_huncorr);
  AddVariationWeight(v,"_btagSF_ldown",p.default_weight/p.w.btagSF*p.w.btagSF_ldown);
  AddVariationWeight(v,"_btagSF_lup",p.default_weight/p.w.btagSF*p.w.btagSF_lup);
  AddVariationWeight(v,"_btagSF_lcorr",p.default_weight/p.w.btagSF*p.w.btagSF_lcorr);
  AddVariationWeight(v,"_btagSF_luncorr",p.default_weight/p.w.btagSF*p.w.btagSF_luncorr);
}
void SMPAnalyzerCore::EvalVariationsBcharge(const Parameter& p,Variations& v){
  AddVariationWeight(v,"_nobchargeSF",p.default_weight/p.w.bchargeSF);
  AddVariationWeight(v,"_bchargeSF_s0m0",p.default_weight/p.w.bchargeSF*p.w.bchargeSF_s0m0);
  AddVariationWeight(v,"_bchargeSF_s0m1",p.default_weight/p.w.bchargeSF*p.w.bchargeSF_s0m1);
}
void SMPAnalyzerCore::EvalVariationsEtc(const Parameter& p,Variations& v){
  AddVariationWeight(v,"_z0weight",p.default_weight*p.w.z0weight);
  AddVariationWeight(v,"_zptweight_gym",p.default_weight/p.w.zptweight*p.w.zptweight_gym);
  AddVariationWeight(v,"_zptweight_gy",p.default_weight/p.w.zptweight*p.w.zptweight_gy);
  AddVariationWeight(v,"_zptweight_g",p.default_weight/p.w.zptweight*p.w.zptweight_g);
  AddVariationWeight(v,"_nozptweight",p.default_weight/p.w.zptweight);
  AddVariationWeight(v,"_noweakweight",p.default_weight/p.w.weakweight);
  AddVariationWeight(v,"_notopptweight",p.default_weight/p.w.topptweight);
}
void SMPAnalyzerCore::EvalVariationsEfficiency(const Parameter& p,Variations& v){
  Weight noeffSF=p.default_weight/p.w.electronRECOSF/p.w.electronIDSF/p.w.muonTrackingSF/p.w.muonRECOSF/p.w.muonIDSF/p.w.muonISOSF/p.w.triggerSF;
  AddVariationWeight(v,"_noefficiencySF",noeffSF);
  for(int j=0,nj=fEff->nreplica;j<nj;j++){
    double electronRECOSF=p.w.electronRECOSF_sys.size() ? (double)p.w.electronRECOSF_sys[0][j] : 1.;
    double electronIDSF=p.w.electronIDSF_sys.size() ? (double)p.w.electronIDSF_sys[0][j] : 1.;
    double muonTrackingSF=p.w.muonTrackingSF_sys.size() ? (double)p.w.muonTrackingSF_sys[0][j] : 1.;
    double muonRECOSF=p.w.muonRECOSF_sys.size() ? (double)p.w.muonRECOSF_sys[0][j] : 1.;
    double muonIDSF=p.w.muonIDSF_sys.size() ? (double)p.w.muonIDSF_sys[0][j] : 1.;
    double triggerSF=p.w.triggerSF_sys.size() ? (double)p.w.triggerSF_sys[0][j] : 1.;
    AddVariationWeight(v,Form("_efficiencySF_stat%d",j),noeffSF*electronRECOSF*electronIDSF*muonTrackingSF*muonRECOSF*muonIDSF*p.w.muonISOSF*triggerSF);
  }
  
  AddVariationWeight(v,"_noelectronRECOSF",p.default_weight/p.w.electronRECOSF);
  for(int i=1,ni=p.w.electronRECOSF_sys.size();i<ni;i++){
    for(int j=0,nj=p.w.electronRECOSF_sys[i].size();j<nj;j++){
      AddVariationWeight(v,Form("_electronRECOSF_s%dm%d",i,j),p.default_weight/p.w.electronRECOSF*p.w.electronRECOSF_sys[i][j]);
    }
  }

  AddVariationWeight(v,"_noelectronIDSF",p.default_weight/p.w.electronIDSF);
  for(int i=1,ni=p.w.electronIDSF_sys.size();i<ni;i++){
    for(int j=0,nj=p.w.electronIDSF_sys[i].size();j<nj;j++){
      AddVariationWeight(v,Form("_electronIDSF_s%dm%d",i,j),p.default_weight/p.w.electronIDSF*p.w.electronIDSF_sys[i][j]);
    }
  }

  AddVariationWeight(v,"_nomuonTrackingSF",p.default_weight/p.w.muonTrackingSF);
  for(int i=1,ni=p.w.muonTrackingSF_sys.size();i<ni;i++){
    for(int j=0,nj=p.w.muonTrackingSF_sys[i].size();j<nj;j++){
      AddVariationWeight(v,Form("_muonTrackingSF_s%dm%d",i,j),p.default_weight/p.w.muonTrackingSF*p.w.muonTrackingSF_sys[i][j]);
    }
  }
    
  AddVariationWeight(v,"_nomuonRECOSF",p.default_weight/p.w.muonRECOSF);
  for(int i=1,ni=p.w.muonRECOSF_sys.size();i<ni;i++){
    for(int j=0,nj=p.w.muonRECOSF_sys[i].size();j<nj;j++){
      AddVariationWeight(v,Form("_muonRECOSF_s%dm%d",i,j),p.default_weight/p.w.muonRECOSF*p.w.muonRECOSF_sys[i][j]);
    }
  }

  AddVariationWeight(v,"_nomuonIDSF",p.default_weight/p.w.muonIDSF);
  for(int i=1,ni=p.w.muonIDSF_sys.size();i<ni;i++){
    for(int j=0,nj=p.w.muonIDSF_sys[i].size();j<nj;j++){
      AddVariationWeight(v,Form("_muonIDSF_s%dm%d",i,j),p.default_weight/p.w.muonIDSF*p.w.muonIDSF_sys[i][j]);
    }
  }
  
  AddVariationWeight(v,"_notriggerSF",p.default_weight/p.w.triggerSF);
  AddVariationWeight(v,"_triggerSF_mode1",p.default_weight/p.w.triggerSF*p.w.triggerSF_mode1);
  AddVariationWeight(v,"_triggerSF_interpolation",p.default_weight/p.w.triggerSF*p.w.triggerSF_interpolation);
  for(int i=1,ni=p.w.triggerSF_sys.size();i<ni;i++){
    for(int j=0,nj=p.w.triggerSF_sys[i].size();j<nj;j++){
      AddVariationWeight(v,Form("_triggerSF_s%dm%d",i,j),p.default_weight/p.w.triggerSF*p.w.triggerSF_sys[i][j]);
    }
  }  
}
void SMPAnalyzerCore::EvalVariationsPDF(const Parameter& p,Variations& v){
  for(unsigned int i=0;i<weight_Scale->size();i++){
    AddVariationWeight(v,Form("_scalevariation%d",i),p.default_weight*TMath::Range(-10.,10.,weight_Scale->at(i)));
  }
  for(unsigned int i=0;i<weight_PDF->size();i++){
    AddVariationWeight(v,Form("_pdf%d",i),p.default_weight*TMath::Range(-10.,10.,weight_PDF->at(i)));
  }
  if(weight_AlphaS->size()==2){
    AddVariationWeight(v,"_alphaS_down",p.default_weight*weight_AlphaS->at(0));
    AddVariationWeight(v,"_alphaS_up",p.default_weight*weight_AlphaS->at(1));
  }
  
  if(weight_PSSyst->size()){
    AddVariationWeight(v,"_FSR_down",p.default_weight*TMath::Range(-5.,5.,weight_PSSyst->at(0)));
    AddVariationWeight(v,"_FSR_up",p.default_weight*TMath::Range(-5.,5.,weight_PSSyst->at(1)));
    AddVariationWeight(v,"_ISR_down",p.default_weight*TMath::Range(-5.,5.,weight_PSSyst->at(2)));
    AddVariationWeight(v,"_ISR_up",p.default_weight*TMath::Range(-5.,5.,weight_PSSyst->at(3)));
  }
  
  if(MCSample.Contains("MiNNLO")){
    AddVariationWeight(v,"_sthw2_down",p.default_weight*weight_sthw2->at(0));
    AddVariationWeight(v,"_sthw2_up",p.default_weight*weight_sthw2->at(2));
    AddVariationWeight(v,"_largeptscales",p.default_weight*weight_largeptscales->at(0));
    AddVariationWeight(v,"_q0_up",p.default_weight*weight_q0->at(0));
    AddVariationWeight(v,"_q0_down",p.default_weight*weight_q0->at(2));
  }
}
void SMPAnalyzerCore::EvalVariationsMuonMomentum(const Parameter& p,Variations& v){
  const vector<double> nmem={1,40,1,1,1,1};
  for(int i=0,ni=nmem.size();i<ni;i++){
    for(int j=0,nj=nmem[i];j<nj;j++){
      AddVariationMuonMomentum(v,Form("_muonmomentum_s%dm%d",i,j),i,j);
    }
  }
  AddVariationMuonMomentum(v,"_muonmomentum_residual",-2,0);
}
void SMPAnalyzerCore::EvalVariationsElectronEnergy(const Parameter& p,Variations& v){
  const vector<double> nmem={1,40,1,1,1,1,1,1,1};
  for(int i=0,ni=nmem.size();i<ni;i++){
    for(int j=0,nj=nmem[i];j<nj;j++){
      AddVariationElectronEnergy(v,Form("_electronenergy_s%dm%d",i,j),i,j);
    }
  }
  AddVariationElectronEnergy(v,"_electronenergy_residual",-2,0);
}
void SMPAnalyzerCore::EvalVariationsJetCorrection(const Parameter& p,Variations& v){
  AddVariationJES(v,"_jet_scale_down",-1);
  AddVariationJES(v,"_jet_scale_up",1);
  if(!IsDATA){
    AddVariationJER(v,"_jet_smear_down",-1);
    AddVariationJER(v,"_jet_smear_up",1);
  }
}

bool SMPAnalyzerCore::PassSelection(Parameter& p,bool cutflow){
  double weight=p.w.lumiweight;
  if(cutflow) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"lumi",weight);

  if(!_event.PassTrigger(p.triggers)) return false;
  if(cutflow) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"trigger",weight);

  if(!PassMETFilter()) return false;
  if(cutflow) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"METfilter",weight);  

  if(p.c.nelectronmax>=0&&(int)p.electrons.size()>p.c.nelectronmax) return false;
  if(p.c.nmuonmax>=0&&(int)p.muons.size()>p.c.nmuonmax) return false;
    
  if(p.c.nleptonmin>=2)
    if(!p.lepton0||!p.lepton1) return false;
  if(p.c.nleptonmin==1)
    if(!p.lepton0) return false;
    
  if(cutflow) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"Dilepton",weight);

  if(p.c.lepton0pt>0){
    if(p.lepton0->Pt()<p.c.lepton0pt) return false;
  }
  if(p.c.lepton1pt>0){
    if(p.lepton1->Pt()<p.c.lepton1pt) return false;
  }
  if(p.c.muon0pt>0){
    if(p.muons.size()<1) return false;
    if(p.muons.at(0).Pt()<p.c.muon0pt) return false;
  }
  if(p.c.muon1pt>0){
    if(p.muons.size()<2) return false;
    if(p.muons.at(1).Pt()<p.c.muon1pt) return false;
  }
  if(p.c.amuon0pt>0){
    if(p.amuons.size()<1) return false;
    if(p.amuons.at(0).Pt()<p.c.amuon0pt) return false;
  }
  if(p.c.amuon1pt>0){
    if(p.amuons.size()<2) return false;
    if(p.amuons.at(1).Pt()<p.c.amuon1pt) return false;
  }
  if(p.c.electron0pt>0){
    if(p.electrons.size()<1) return false;
    if(p.electrons.at(0).Pt()<p.c.electron0pt) return false;
  }
  if(p.c.electron1pt>0){
    if(p.electrons.size()<2) return false;
    if(p.electrons.at(1).Pt()<p.c.electron1pt) return false;
  }
  if(p.c.aelectron0pt>0){
    if(p.aelectrons.size()<1) return false;
    if(p.aelectrons.at(0).Pt()<p.c.aelectron0pt) return false;
  }
  if(p.c.aelectron1pt>0){
    if(p.aelectrons.size()<2) return false;
    if(p.aelectrons.at(1).Pt()<p.c.aelectron1pt) return false;
  }
  if(cutflow) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"LepPtCut",weight);

  if(cutflow) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"charge",weight);

  if(p.option.Contains("triggermatching")){
    if(p.triggers.size()){
      if(p.triggers.at(0).Contains(TRegexp("HLT_[Tk]*Mu17_TrkIsoVVL_[Tk]*Mu8_TrkIsoVVL"))){
	if(p.lepton0->LeptonFlavour()!=Lepton::MUON) return false;
	if(p.lepton1->LeptonFlavour()!=Lepton::MUON) return false;
	Muon *muon0=(Muon*)p.lepton0,*muon1=(Muon*)p.lepton1;
	if(!PassDLT1(muon0)&&!PassDLT1(muon1)) return false;
	if(!PassDLT2(muon0)||!PassDLT2(muon1)) return false;
      }else if(p.triggers.at(0).Contains("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL")){
	if(p.lepton0->LeptonFlavour()!=Lepton::ELECTRON) return false;
	if(p.lepton1->LeptonFlavour()!=Lepton::ELECTRON) return false;
	Electron *electron0=(Electron*)p.lepton0,*electron1=(Electron*)p.lepton1;
	if(!PassDLT1(electron0)&&!PassDLT1(electron1)) return false;
	if(!PassDLT2(electron0)||!PassDLT2(electron1)) return false;
      }else if(p.channel=="mj"||p.channel=="Mj"){
	bool passor=false;
	for(auto trigger:p.triggers){
	  if(((Muon*)p.lepton0)->PassPath(trigger)) passor=true;
	}
	if(!passor) return false;
      }else{
	cout<<"[SMPAnalyzerCore::PassSelection] trigger matching for "<<p.triggers.at(0)<<" is not implemented"<<endl;
	exit(EXIT_FAILURE);
      }
    }else{
      cout<<"[SMPAnalyzerCore::PassSelection] no trigger for trigger matching"<<endl;
      exit(EXIT_FAILURE);
    }
    if(cutflow) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"TriggerMatching",weight);
  }
  if(p.option.Contains("strictorder")){
    vector<Lepton*> leptons;
    vector<Lepton*> electrons=MakeLeptonPointerVector(p.electrons);
    vector<Lepton*> aelectrons=MakeLeptonPointerVector(p.aelectrons);
    vector<Lepton*> muons=MakeLeptonPointerVector(p.muons);
    vector<Lepton*> amuons=MakeLeptonPointerVector(p.amuons);
    leptons.insert(leptons.end(),electrons.begin(),electrons.end());
    leptons.insert(leptons.end(),aelectrons.begin(),aelectrons.end());
    leptons.insert(leptons.end(),muons.begin(),muons.end());
    leptons.insert(leptons.end(),amuons.begin(),amuons.end());
    std::sort(leptons.begin(),leptons.end(),PtComparingPtr);
    if(p.c.nleptonmin>=1&&leptons.at(0)!=p.lepton0) return false;
    if(p.c.nleptonmin>=2&&leptons.at(1)!=p.lepton1) return false;
    if(cutflow) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"StrictOrder",weight);
  }

  if(cutflow){
    weight*=p.w.PUweight;
    FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"PU",weight);
    weight*=p.w.prefireweight;
    FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"prefire",weight);
    weight*=p.w.zptweight;
    FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"zpt",weight);
    weight*=p.w.weakweight;
    FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"weak",weight);
    weight*=p.w.electronRECOSF*p.w.muonRECOSF;
    FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"RECOSF",weight);
    weight*=p.w.electronIDSF*p.w.muonIDSF;
    FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"IDSF",weight);
    weight*=p.w.triggerSF;
    FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"triggerSF",weight);
    weight*=p.w.CFSF;
    FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"CFSF",weight);
  }
  return true;
}

TH4D* SMPAnalyzerCore::GetHist4D(TString histname){
  TH4D *h = NULL;
  std::map<TString, TH4D*>::iterator mapit = maphist_TH4D.find(histname);
  if(mapit != maphist_TH4D.end()) return mapit->second;
  return h;
}

void SMPAnalyzerCore::FillProfile(TString histname,
				  Double_t value_x, Double_t value_y, Double_t weight,
				  Int_t n_binx, Double_t x_min, Double_t x_max){
  TProfile *this_hist = NULL;
  auto mapit = maphist_TH1D.find(histname);
  if(mapit == maphist_TH1D.end()){
    this_hist = new TProfile(histname, "", n_binx, x_min, x_max);
    this_hist->SetDirectory(NULL);
    maphist_TH1D[histname] = this_hist;
  }else{
    this_hist = (TProfile*)mapit->second;
  }
  this_hist->Fill(value_x, value_y, weight);

}

void SMPAnalyzerCore::FillHist(TString histname,
                            Double_t value_x, Double_t value_y, Double_t value_z, Double_t value_u,
                            Double_t weight,
                            Int_t n_binx, Double_t x_min, Double_t x_max,
                            Int_t n_biny, Double_t y_min, Double_t y_max,
                            Int_t n_binz, Double_t z_min, Double_t z_max,
                            Int_t n_binu, Double_t u_min, Double_t u_max){

  TH4D *this_hist = GetHist4D(histname);
  if( !this_hist ){
    this_hist = new TH4D(histname, "", n_binx, x_min, x_max, n_biny, y_min, y_max, n_binz, z_min, z_max, n_binu, u_min, u_max);
    this_hist->SetDirectory(NULL);
    maphist_TH4D[histname] = this_hist;
  }

  this_hist->Fill(value_x, value_y, value_z, value_u, weight);

}

void SMPAnalyzerCore::FillHist(TString histname,
                            Double_t value_x, Double_t value_y, Double_t value_z, Double_t value_u,
                            Double_t weight,
                            Int_t n_binx, const Double_t *xbins,
                            Int_t n_biny, const Double_t *ybins,
                            Int_t n_binz, const Double_t *zbins,
                            Int_t n_binu, const Double_t *ubins){

  TH4D *this_hist = GetHist4D(histname);
  if( !this_hist ){
    this_hist = new TH4D(histname, "", n_binx, xbins, n_biny, ybins, n_binz, zbins, n_binu, ubins);
    this_hist->SetDirectory(NULL);
    maphist_TH4D[histname] = this_hist;
  }

  this_hist->Fill(value_x, value_y, value_z, value_u, weight);

}

void SMPAnalyzerCore::FillHist(TString histname,
                            Double_t value_x, Double_t value_y, Double_t value_z, Double_t value_u,
                            Double_t weight,
                            Int_t n_binx, const Double_t *xbins,
                            Int_t n_biny, const Double_t *ybins,
                            Int_t n_binz, const Double_t *zbins,
			    Int_t n_binu, Double_t u_min, Double_t u_max){

  TH4D *this_hist = GetHist4D(histname);
  if( !this_hist ){
    TAxis uaxis(n_binu,u_min,u_max);
    vector<double> ubins={};
    for(int i=1;i<n_binu+2;i++) ubins.push_back(uaxis.GetBinLowEdge(i));
    this_hist = new TH4D(histname, "", n_binx, xbins, n_biny, ybins, n_binz, zbins, n_binu, &ubins[0]);
    this_hist->SetDirectory(NULL);
    maphist_TH4D[histname] = this_hist;
  }

  this_hist->Fill(value_x, value_y, value_z, value_u, weight);

}
void SMPAnalyzerCore::FillHistsSyst(Parameter p,Variations& v){
  for(auto& [vsuf,variation]:v){
    Apply(p,vsuf,variation);
    if(PassSelection(p,vsuf=="")){
      FillHists(p);
    }
  }
}
void SMPAnalyzerCore::FillHists(Parameter& p){
  TLorentzVector dilepton=(*p.lepton0)+(*p.lepton1);
  double dimass=dilepton.M();
  if(dimass>=60&&dimass<120){
    if(p.vsuffix=="") FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix+p.vsuffix,"m60to120",p.weight);
    FillHist(p.prefix+"m60to120/"+p.hprefix+"dimass"+p.suffix+p.vsuffix,dimass,p.weight,60,60,120);
    if(dimass>=80&&dimass<100){
      if(p.vsuffix=="") FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix+p.vsuffix,"m80to100",p.weight);
      FillHist(p.prefix+"m80to100/"+p.hprefix+"dimass"+p.suffix+p.vsuffix,dimass,p.weight,40,80,100);
    }
  }
}

void SMPAnalyzerCore::WriteHist(){
  AnalyzerCore::WriteHist();
  outfile->cd();
  for(std::map< TString, TH4D* >::iterator mapit = maphist_TH4D.begin(); mapit!=maphist_TH4D.end(); mapit++){
    TString this_fullname=mapit->second->GetName();
    TString this_name=this_fullname(this_fullname.Last('/')+1,this_fullname.Length());
    TString this_suffix=this_fullname(0,this_fullname.Last('/'));
    TDirectory *dir = outfile->GetDirectory(this_suffix);
    if(!dir){
      outfile->mkdir(this_suffix);
    }
    outfile->cd(this_suffix);
    mapit->second->Write(this_name);
    outfile->cd();
  }
}

void SMPAnalyzerCore::FillHist(TString histname, double value, map<TString,double> weights, int n_bin, double x_min, double x_max){
  for(const auto& [suffix,weight]:weights) FillHist(histname+suffix,value,weight,n_bin,x_min,x_max);
}
void SMPAnalyzerCore::FillHist(TString histname, double value, map<TString,double> weights, int n_bin, const double *xbins){
  for(const auto& [suffix,weight]:weights) FillHist(histname+suffix,value,weight,n_bin,xbins);
}
void SMPAnalyzerCore::FillHist(TString histname, double value_x, double value_y, map<TString,double> weights, int n_binx, double x_min, double x_max, int n_biny, double y_min, double y_max){
  for(const auto& [suffix,weight]:weights) FillHist(histname+suffix,value_x,value_y,weight,n_binx,x_min,x_max,n_biny,y_min,y_max);
}
void SMPAnalyzerCore::FillHist(TString histname, double value_x, double value_y, map<TString,double> weights, int n_binx, const double *xbins, int n_biny, const double *ybins){
  for(const auto& [suffix,weight]:weights) FillHist(histname+suffix,value_x,value_y,weight,n_binx,xbins,n_biny,ybins);
}
void SMPAnalyzerCore::FillHist(TString histname, double value_x, double value_y, double value_z, map<TString,double> weights, int n_binx, double x_min, double x_max, int n_biny, double y_min, double y_max, int n_binz, double z_min, double z_max){
  for(const auto& [suffix,weight]:weights) FillHist(histname+suffix,value_x,value_y,value_z,weight,n_binx,x_min,x_max,n_biny,y_min,y_max,n_binz,z_min,z_max);
}
void SMPAnalyzerCore::FillHist(TString histname, double value_x, double value_y, double value_z, map<TString,double> weights, int n_binx, const double *xbins, int n_biny, const double *ybins, int n_binz, const double *zbins){
  for(const auto& [suffix,weight]:weights) FillHist(histname+suffix,value_x,value_y,value_z,weight,n_binx,xbins,n_biny,ybins,n_binz,zbins);
}
void SMPAnalyzerCore::FillHist(TString histname, double value_x, double value_y, double value_z, double value_u, map<TString,double> weights, int n_binx, double x_min, double x_max, int n_biny, double y_min, double y_max, int n_binz, double z_min, double z_max, int n_binu, double u_min, double u_max){
  for(const auto& [suffix,weight]:weights) FillHist(histname+suffix,value_x,value_y,value_z,value_u,weight,n_binx,x_min,x_max,n_biny,y_min,y_max,n_binz,z_min,z_max,n_binu,u_min,u_max);
}
void SMPAnalyzerCore::FillHist(TString histname, double value_x, double value_y, double value_z, double value_u, map<TString,double> weights, int n_binx, const double *xbins, int n_biny, const double *ybins, int n_binz, const double *zbins, int n_binu, const double *ubins){
  for(const auto& [suffix,weight]:weights) FillHist(histname+suffix,value_x,value_y,value_z,value_u,weight,n_binx,xbins,n_biny,ybins,n_binz,zbins,n_binu,ubins);
}
void SMPAnalyzerCore::FillHist(TString histname, double value_x, double value_y, double value_z, double value_u, map<TString,double> weights, int n_binx, const double *xbins, int n_biny, const double *ybins, int n_binz, const double *zbins, int n_binu, double u_min, double u_max){
  for(const auto& [suffix,weight]:weights) FillHist(histname+suffix,value_x,value_y,value_z,value_u,weight,n_binx,xbins,n_biny,ybins,n_binz,zbins,n_binu,u_min,u_max);
}

void SMPAnalyzerCore::FillDileptonHists(TString pre,TString suf,Particle *l0,Particle *l1,double w){
  TLorentzVector dilepton=*l0+*l1;
  double dimass=dilepton.M();
  double dipt=dilepton.Pt();
  double dirap=dilepton.Rapidity();
  FillHist(pre+"dimass"+suf,dimass,w,200,0,400);
  FillHist(pre+"dipt"+suf,dipt,w,200,0,400);
  FillHist(pre+"dirap"+suf,dirap,w,120,-6,6);
  vector<Particle*> leps;
  if(l0->Pt()>l1->Pt()) leps={l0,l1};
  else leps={l1,l0};
  for(int i=0;i<(int)leps.size();i++){
    FillHist(Form("%sl%dpt%s",pre.Data(),i,suf.Data()),leps.at(i)->Pt(),w,200,0,400);
    FillHist(Form("%sl%deta%s",pre.Data(),i,suf.Data()),leps.at(i)->Eta(),w,100,-5,5);
    FillHist(Form("%sl%dphi%s",pre.Data(),i,suf.Data()),leps.at(i)->Phi(),w,80,-4,4);
  }
  FillHist(pre+"lldelR"+suf,l0->DeltaR(*l1),w,70,0,7);  
  FillHist(pre+"lldelphi"+suf,l0->DeltaPhi(*l1),w,80,-4,4);
  FillHist(pre+"lldeleta"+suf,fabs(l0->Eta()-l1->Eta()),w,100,-5,5);
}
double SMPAnalyzerCore::GetPtThreshold(TString path){
  TString str=path(TRegexp("[0-9]+"));
  if(str.IsFloat()) return str.Atof();
  else return -1.;
}
bool SMPAnalyzerCore::IsExists(TString filepath){
  ifstream fcheck(filepath);
  return fcheck.good();
}
void SMPAnalyzerCore::SetupEfficiency(){
  TString configpath=getenv("DATA_DIR")+TString("/")+GetEra()+"/ID/eff.conf";
  if(IsExists(configpath)){
    fEff=new EfficiencyTool(configpath);
  }
}
void SMPAnalyzerCore::DeleteEfficiency(){
  if(fEff) delete fEff;
}
double SMPAnalyzerCore::GetLeptonTriggerSF(TString triggerSF_key,const vector<Lepton*>& leps,int set,int mem,TString option){
  if(IsDATA) return 1;
  if(triggerSF_key=="") return 1;
  if(triggerSF_key=="Default") return 1;

  double data_eff=1.,sim_eff=1.;
  for(const auto& lep:leps){
    if(!lep) continue;
    data_eff*=1-fEff->GetDataEfficiency(triggerSF_key,lep,set,mem,option);
    sim_eff*=1-fEff->GetSimEfficiency(triggerSF_key,lep,set,mem,option);
  }
  data_eff=1-data_eff;
  sim_eff=1-sim_eff;
  if(sim_eff==0) return 1.;
  else return data_eff/sim_eff;
}
double SMPAnalyzerCore::GetLeptonTriggerORSF(const Parameter& p,const vector<Lepton*>& leps,int set,int mem,TString option){
  if(IsDATA) return 1;
  if(p.triggers.size()!=2){
    cout<<"[SMPAnalyzerCore::LeptonTriggerOR_SF] p.triggers.size()= "<<p.triggers.size()<<endl;
    exit(EXIT_FAILURE);
  }
  if(p.k.triggerSF.size()!=2){
    cout<<"[SMPAnalyzerCore::LeptonTriggerOR_SF] p.k.triggerSF.size()= "<<p.k.triggerSF.size()<<endl;
    exit(EXIT_FAILURE);
  }
  double lumi=_event.GetTriggerLumi("Full");
  double lumi0,lumi1,lumi01;
  if(DataYear==2017&&p.k.triggerSF[0].Contains("IsoMu24")&&p.k.triggerSF[1].Contains("IsoMu27")){
    lumi0=_event.GetTriggerLumi(p.triggers[0]); lumi1=_event.GetTriggerLumi(p.triggers[1]); lumi01=lumi0;
  }else if(DataYear==2017&&p.k.triggerSF[0].Contains("Ele27")&&p.k.triggerSF[1].Contains("Ele32")){
    lumi0=_event.GetTriggerLumi(p.triggers[0]); lumi1=_event.GetTriggerLumi(p.triggers[1]); lumi01=17599.732185;
  }else if(DataYear==2018&&p.k.triggerSF[0].Contains("Ele28")&&p.k.triggerSF[1].Contains("Ele32")){
    lumi0=_event.GetTriggerLumi(p.triggers[0]); lumi1=_event.GetTriggerLumi(p.triggers[1]); lumi01=lumi0;
  }else{
    cout<<"[SMPAnalyzerCore::GetLeptonTriggerORSF] not available combination '"<<p.k.triggerSF[0]<<"'||'"<<p.k.triggerSF[1]<<"' for "<<DataEra<<endl;
    exit(EXIT_FAILURE);
  }

  bool newflag=set<0; //temp
  if(newflag) set=0; //temp
  double data_eff0=1.,sim_eff0=1.;
  double data_eff1=1.,sim_eff1=1.;
  for(const auto& lep:leps){
    if(!lep) continue;
    data_eff0*=1-fEff->GetDataEfficiency(p.k.triggerSF[0],lep,set,mem,option);
    sim_eff0*=1-fEff->GetSimEfficiency(p.k.triggerSF[0],lep,set,mem,option);
    data_eff1*=1-fEff->GetDataEfficiency(p.k.triggerSF[1],lep,set,mem,option);
    sim_eff1*=1-fEff->GetSimEfficiency(p.k.triggerSF[1],lep,set,mem,option);
  }
  data_eff0=1-data_eff0;
  sim_eff0=1-sim_eff0;
  data_eff1=1-data_eff1;
  sim_eff1=1-sim_eff1;
  double sf=0.;
  if(_event.PassTrigger(p.triggers[1])){
    double this_sf=(lumi1-lumi01)/lumi;
    if(sim_eff1) this_sf*=data_eff1/sim_eff1;
    sf+=this_sf;
  }
  if(_event.PassTrigger(p.triggers[0])){
    double this_sf=(lumi0-lumi01)/lumi;
    if(sim_eff0) this_sf*=data_eff0/sim_eff0;
    sf+=this_sf;
  }
  //overlap region
  if(!newflag){
    if(lumi0>lumi1){
      if(_event.PassTrigger(p.triggers[0])){
	double this_sf=lumi01/lumi;
	if(sim_eff0) this_sf*=data_eff0/sim_eff0;
	sf+=this_sf;
      }
    }else{
      if(_event.PassTrigger(p.triggers[1])){
	double this_sf=lumi01/lumi;
	if(sim_eff1) this_sf*=data_eff1/sim_eff1;
	sf+=this_sf;
      }else if(_event.PassTrigger(p.triggers[0])){
	double this_sf=lumi01/lumi;
	if(sim_eff0) this_sf*=data_eff0/sim_eff0;
	sf+=this_sf;
      }
    }    
  }else{
    if(_event.PassTrigger(p.triggers[0])){
      double this_sf=lumi01/2/lumi;
      if(sim_eff0) this_sf*=data_eff0/sim_eff0;
      sf+=this_sf;
    }
    if(_event.PassTrigger(p.triggers[1])){
      double this_sf=lumi01/2/lumi;
      if(sim_eff1) this_sf*=data_eff1/sim_eff1;
      sf+=this_sf;
    }else if(_event.PassTrigger(p.triggers[0])){
      double this_sf=lumi01/2/lumi;
      if(sim_eff0) this_sf*=data_eff0/sim_eff0;
      sf+=this_sf;
    }    
  }    
  return sf;
}
double SMPAnalyzerCore::GetDileptonTriggerSF(TString triggerSF_key0,TString triggerSF_key1,TString DZSF,const vector<Lepton*>& leps,int set,int mem,TString option){
  if(IsDATA) return 1;
  if((triggerSF_key0==""||triggerSF_key0=="Default")&&(triggerSF_key1==""||triggerSF_key1=="Default")) return 1;
  int nlep=leps.size();
  if(nlep<2){
    //cout<<"[SMPAnalyzerCore::DileptonTrigger_SF] nlep < 2. return 1."<<endl;
    return 1.;
  }
  double data_noleg1=1.,sim_noleg1=1.;
  vector<double> data_oneleg1_noleg2(nlep,1.);
  vector<double> sim_oneleg1_noleg2(nlep,1.);
  for(int i=0;i<nlep;i++){
    if(!leps.at(i)) return 1.;
    double data_eff_leg1=fEff->GetDataEfficiency(triggerSF_key0,leps.at(i),set,mem,option);
    double data_eff_leg2=fEff->GetDataEfficiency(triggerSF_key1,leps.at(i),set,mem,option);
    double sim_eff_leg1=fEff->GetSimEfficiency(triggerSF_key0,leps.at(i),set,mem,option);
    double sim_eff_leg2=fEff->GetSimEfficiency(triggerSF_key1,leps.at(i),set,mem,option);
    if(DZSF!=""){
      double data_eff_dz=fEff->GetDataEfficiency(DZSF,leps.at(i),0,0,option);
      double sim_eff_dz=fEff->GetSimEfficiency(DZSF,leps.at(i),0,0,option);
      data_eff_leg1*=data_eff_dz;
      data_eff_leg2*=data_eff_dz;
      sim_eff_leg1*=sim_eff_dz;
      sim_eff_leg2*=sim_eff_dz;
    } 
    data_noleg1*=(1-data_eff_leg1);
    sim_noleg1*=(1-sim_eff_leg1);
    for(int j=0;j<nlep;j++){
      if(i==j){
	data_oneleg1_noleg2[j]*=data_eff_leg1;
	sim_oneleg1_noleg2[j]*=sim_eff_leg1;
      }else{
	data_oneleg1_noleg2[j]*=(1-data_eff_leg2);
	sim_oneleg1_noleg2[j]*=(1-sim_eff_leg2);
      }
    }
  }
  double data_eff=1.-data_noleg1;
  double sim_eff=1.-sim_noleg1;
  for(int i=0;i<nlep;i++){
    data_eff-=data_oneleg1_noleg2[i];
    sim_eff-=sim_oneleg1_noleg2[i];
  }
  double sf=1.;
  if(sim_eff==0) return sf=1.;
  else sf=data_eff/sim_eff;
  return sf;
}

// ZptWeight
void SMPAnalyzerCore::SetupZptWeight(){
  TString _MCSample=MCSample;
  if(MCSample.Contains("MiNNLO")) _MCSample="MiNNLO";
  TString zptpath=(TString)getenv("SKFlat_WD")+"/external/ZptCorrection/ZptWeight_"+_MCSample+".root";
  if(IsExists(zptpath)){
    cout<<"[SMPAnalyzerCore::SetupZptWeight] using file "+zptpath<<endl;
  }else{
    cout<<"[SMPAnalyzerCore::SetupZptWeight] no "+zptpath<<endl;
    return;
  }
  DeleteZptWeight();
  fZptCorrection=new ZptCorrection(zptpath);
}
void SMPAnalyzerCore::DeleteZptWeight(){
  if(fZptCorrection) delete fZptCorrection;
}

double SMPAnalyzerCore::GetTopPtReweight2(const std::vector<Gen>& gens){
  //==== ref: https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting2017
  //==== Only top quarks in SM ttbar events must be reweighted,
  //==== not single tops or tops from BSM production mechanisms.
  if(!MCSample.Contains("TT") || !MCSample.Contains("powheg")){
    return 1.;
  }
  //==== initialize with large number                                                                                                                                                                                                                                                                                                                           
  double toppt1=10000, toppt2=10000;
  bool found_top = false, found_atop = false;
  
  for(vector<Gen>::const_iterator genit=gens.begin(); genit!=gens.end(); genit++){
    
    if(genit->Status() == 22){
      if(genit->PID() == 6){
        toppt1= genit->Pt();
        found_top = true;
      }
      else if(genit->PID() == -6){
        toppt2= genit->Pt();
        found_atop = true;
      }
    }
    //==== after we found top pair, break the loop
    if(found_top && found_atop) break;
  }
  double pt_reweight = 1.;
  //==== if top pair is not found, return 1.
  pt_reweight*=0.103*exp(-0.0118*toppt1)-0.000134*toppt1+0.973;
  pt_reweight*=0.103*exp(-0.0118*toppt2)-0.000134*toppt2+0.973;
  pt_reweight = sqrt(pt_reweight);
  return pt_reweight;
}


void SMPAnalyzerCore::SetupRoccoR(){
  cout<<"[SMPAnalyzerCore::SetupRoccoR] setting Rocheseter Correction"<<endl;
  TString erashort=GetEraShort();

  //TString rocpath=datapath+"/"+GetEra()+"/RoccoR/RoccoR"+GetEraShort()+"UL.txt"; //central roccor for amc
  TString rocpath=TString(getenv("SKFlat_WD"))+"/external/Aepcor/u_"+erashort(2,3)+"UL_1.txt"; //roccor for minnlo
  if(IsExists(rocpath)) roc=new RoccoR(rocpath.Data());
  else cout<<"[SMPAnalyzerCore::SetupRoccoR] no "+rocpath<<endl;

  TString rocelepath=TString(getenv("SKFlat_WD"))+"/external/Aepcor/e_"+erashort(2,3)+"UL_1.txt";
  if(IsExists(rocelepath)){
    rocele=new Aepcor;
    rocele->init(rocelepath.Data(),Aepres::CB);
  }
  else cout<<"[SMPAnalyzerCore::SetupRoccoR] no "+rocelepath<<endl;

  TString rocresidualpath=TString(getenv("DATA_DIR"))+"/"+GetEra()+"/SMP/RoccorResidual.root";
  if(IsExists(rocresidualpath)){
    TFile f(rocresidualpath);
    for(TString suffix:{"_scale","_centralRes","_leftFrac","_leftRes","_rightFrac","_rightRes"}){
      fRoccorResidual["electron"+suffix]=(TH1*)f.Get("ee"+GetEraShort()+suffix);
      fRoccorResidual["muon"+suffix]=(TH1*)f.Get("mm"+GetEraShort()+suffix);
    }
    for(auto& [_,h]:fRoccorResidual){
      if(h){
	h->SetDirectory(0);
      }
    }
  }
  else cout<<"[SMPAnalyzerCore::SetupRoccoR] no "+rocresidualpath<<endl;
}
double SMPAnalyzerCore::GetZ0Weight(double valx){
  if(IsDATA) return 1.;
  double rt=1.;
  if(GetEra()=="2016preVFP"){
    double data_val=TMath::Gaus(valx,2.46312e-01,3.50458e+00,true);
    double mc_val=TMath::Gaus(valx,9.28612e-01,3.65203e+00,true);
    rt=data_val/mc_val;
  }else if(GetEra()=="2016postVFP"){
    double data_val=TMath::Gaus(valx,2.41640e-01,3.63717e+00,true);
    double mc_val=TMath::Gaus(valx,9.30108e-01,3.65454e+00,true);
    rt=data_val/mc_val;
  }else if(GetEra()=="2017"){
    double data_val=TMath::Gaus(valx,3.81830e-01,3.67614e+00,true);
    double mc_val=TMath::Gaus(valx,8.19642e-01,3.50992e+00,true);
    rt=data_val/mc_val;
  }else if(GetEra()=="2018"){
    double data_val=TMath::Gaus(valx,-1.36030e-01,3.41464e+00,true);
    double mc_val=TMath::Gaus(valx,3.58575e-02,3.50953e+00,true);
    rt=data_val/mc_val;
  } 
  if(rt>2) rt=2;
  return rt;
}
void SMPAnalyzerCore::SetupCFRate(){
  cout<<"[SMPAnalyzerCore::SetupCFRate] setting CFRate"<<endl;
  TString datapath=getenv("DATA_DIR");
  if(!IsExists(datapath+"/"+GetEra()+"/SMP/CFRate.root")){
    cout<<"[SMPAnalyzerCore::SetupCFRate] no CFRate.root"<<endl;
    return;
  }
  TFile f(datapath+"/"+GetEra()+"/SMP/CFRate.root");
  hcfrate_data=(TH2*)f.Get("cfdata");
  if(hcfrate_data){
    cout<<"[SMPAnalyzerCore::SetupCFRate] load hcfrate_data"<<endl;
    hcfrate_data->SetDirectory(0);
  }else cout<<"[SMPAnalyzerCore::SetupCFRate] no hcfrate_data"<<endl;
  hcfrate_mc=(TH2*)f.Get("cfmc");
  if(hcfrate_mc){
    cout<<"[SMPAnalyzerCore::SetupCFRate] load hcfrate_mc"<<endl;
    hcfrate_mc->SetDirectory(0);
  }else cout<<"[SMPAnalyzerCore::SetupCFRate] no hcfrate_mc"<<endl;
  hcfsf=(TH2*)f.Get("cfsf");
  if(hcfsf){
    cout<<"[SMPAnalyzerCore::SetupCFRate] load hcfsf"<<endl;
    hcfsf->SetDirectory(0);
  }else cout<<"[SMPAnalyzerCore::SetupCFRate] no hcfsf"<<endl;
  //hcfenergyscale=(TH2*)f.Get("cfenergyscale");
  //if(hcfenergyscale){
  //  cout<<"[SMPAnalyzerCore::SetupCFRate] load hcfenergyscale"<<endl;
  //  hcfenergyscale->SetDirectory(0);
  //}else cout<<"[SMPAnalyzerCore::SetupCFRate] no hcfenergyscale"<<endl;
  f.Close();
}
double SMPAnalyzerCore::GetCFSF(const Lepton* l,int sys){
  if(IsDATA) return 1.;
  if(!hcfsf) return 1.;
  if(!l) return 1.;
  if(l->LeptonFlavour()!=Lepton::ELECTRON) return 1.;
  return GetBinContentUser(hcfsf,l->Eta(),l->Pt(),sys);
}
double SMPAnalyzerCore::GetCFSF(const Parameter& p,int sys){
  if(IsDATA) return 1.;
  double sf=1.;
  if(p.lepton0&&!p.truth_lepton0.IsEmpty())
    if(p.lepton0->Charge()*p.truth_lepton0.Charge()<0) sf*=GetCFSF(p.lepton0,sys);
  if(p.lepton1&&!p.truth_lepton1.IsEmpty())
    if(p.lepton1->Charge()*p.truth_lepton1.Charge()<0) sf*=GetCFSF(p.lepton1,sys);
  return sf;
}
void SMPAnalyzerCore::DeleteCFRate(){
  if(hcfrate_data) delete hcfrate_data;
  if(hcfrate_mc) delete hcfrate_mc;
  if(hcfsf) delete hcfsf;
}
void SMPAnalyzerCore::SetupMuonTrackingSF(){
  jSetupMuonTrackingSF=true;
  cout<<"[SMPAnalyzerCore::SetupMuonTrackingSF] setup"<<endl;
  TString datapath=getenv("DATA_DIR");
  TString era=GetEra();
  if(IsExists(datapath+"/"+era+"/SMP/muonPOGtrackingSF.root")){
    TFile f(datapath+"/"+era+"/SMP/muonPOGtrackingSF.root");
    if(GetEra()=="2016preVFP"){
      fMuonTrackingSF=(TH1*)f.Get("muonPOGtrackingSF_preVFP");
      if(fMuonTrackingSF){
	cout<<"[SMPAnalyzerCore::SetupMuonTrackingSF] load muonPOGtrackingSF_preVFP"<<endl;
	fMuonTrackingSF->SetDirectory(NULL);
      }
    }else if(GetEra()=="2016postVFP"){
      fMuonTrackingSF=(TH1*)f.Get("muonPOGtrackingSF_postVFP");
      if(fMuonTrackingSF){
	cout<<"[SMPAnalyzerCore::SetupMuonTrackingSF] load muonPOGtrackingSF_postVFP"<<endl;
	fMuonTrackingSF->SetDirectory(NULL);
      }
    }
  }
}
double SMPAnalyzerCore::GetMuonTrackingSF(double eta,int sys){
  double sf=1.;
  if(IsDATA) return sf;
  if(!jSetupMuonTrackingSF) SetupMuonTrackingSF();
  if(fMuonTrackingSF){
    sf=GetBinContentUser(fMuonTrackingSF,eta,sys);
  }
  return sf;
}
void SMPAnalyzerCore::DeleteMuonTrackingSF(){
  if(fMuonTrackingSF){
    delete fMuonTrackingSF;
  }
  fMuonTrackingSF=NULL;
  jSetupMuonTrackingSF=false;
}
double SMPAnalyzerCore::GetDYWeakWeight(double mass){
  if(IsDATA) return 1.;
  if(!IsDYSample) return 1.;
  if(mass<55) return 0.988939;
  else if(mass<60) return 0.992556;
  else if(mass<65) return 0.996362;
  else if(mass<70) return 1.00086;
  else if(mass<75) return 1.00593;
  else if(mass<80) return 1.00989;
  else if(mass<85) return 1.01263;
  else if(mass<90) return 1.01373;
  else if(mass<95) return 1.01338;
  else if(mass<100) return 1.01242;
  else if(mass<110) return 1.01078;
  else if(mass<120) return 1.00839;
  else if(mass<130) return 1.00628;
  else if(mass<140) return 1.00461;
  else if(mass<150) return 1.0033;
  else if(mass<170) return 1.00201;
  else if(mass<200) return 0.999256;
  else if(mass<250) return 0.995825;
  else if(mass<300) return 0.992451;
  else if(mass<400) return 0.986289;
  else if(mass<500) return 0.979024;
  else if(mass<600) return 0.972292;
  else if(mass<700) return 0.967596;
  else if(mass<800) return 0.959725;
  else if(mass<1000) return 0.953025;
  else if(mass<1500) return 0.935142;
  else if(mass<2000) return 0.909548;
  else if(mass<3000) return 0.8895;
  else return 0.900657;
}
void SMPAnalyzerCore::SetupFakeRate(){
  cout<<"[SMPAnalyzerCore::SetupFakeRate] setup"<<endl;
  TString datapath=getenv("DATA_DIR");
  TString era=GetEra();
  if(IsExists(datapath+"/"+era+"/SMP/FakeRate.root")){
    TFile f(datapath+"/"+era+"/SMP/FakeRate.root");
    fFakeRate_electron=(TH2*)f.Get("ee"+GetEra());
    if(fFakeRate_electron){
      cout<<"[SMPAnalyzerCore::SetupFakeRate] load ee"+GetEra()<<endl;
      fFakeRate_electron->SetDirectory(NULL);
    }
    fFakeRate_muon=(TH2*)f.Get("mm"+GetEra());
    if(fFakeRate_muon){
      cout<<"[SMPAnalyzerCore::SetupFakeRate] load mm"+GetEra()<<endl;
      fFakeRate_muon->SetDirectory(NULL);
    }
  }
  if(IsExists(datapath+"/"+era+"/SMP/FakeTF.root")){
    TFile f(datapath+"/"+era+"/SMP/FakeTF.root");
    for(const auto& obj:*(f.GetListOfKeys())){
      TKey* key=(TKey*)obj;
      TString name=key->GetName();
      if(name.Contains(GetEra())){
	TH2* hist=(TH2*)f.Get(name);
	if(hist){
	  cout<<"[SMPAnalyzerCore::SetupFakeTF] load "+name<<endl;
	  hist->SetDirectory(NULL);
	  fFakeTF[name]=hist;
	}
      }
    }
  }
}
double SMPAnalyzerCore::GetFakeTF(Parameter& p,TString option,int sys){
  TH2* fFakeTF_l0=NULL;
  TH2* fFakeTF_l1=NULL;
  TString key="";
  if(option.Contains("lj")){
    if(p.lepton0->InheritsFrom("Electron")){
      key+="ej";
    }else if(p.lepton0->InheritsFrom("Muon")){
      key+="mj";
    }
  }else{
    if(p.lepton0->InheritsFrom("Electron")){
      key+="ee";
    }else if(p.lepton0->InheritsFrom("Muon")){
      key+="mm";
    }
  }
  key+=GetEra();
  if(option.Contains("cpt")) key+="_cpt";
  else if(option.Contains("mpt")) key+="_mpt";
  
  if(option.Contains("0bjet")) key+="_0bjet";
  else if(option.Contains("nbjet")) key+="_nbjet";

  if(!option.Contains("lj")) key+="_noZ";

  if(fFakeTF.find(key+"_l0tf")!=fFakeTF.end()) fFakeTF_l0=fFakeTF[key+"_l0tf"];
  if(fFakeTF.find(key+"_l1tf")!=fFakeTF.end()) fFakeTF_l1=fFakeTF[key+"_l1tf"];

  if(!fFakeTF_l0||!fFakeTF_l1) return 0.;

  double tf=1.;
  tf*=GetBinContentUser(fFakeTF_l0,fabs(p.lepton0->Eta()),p.lepton0->Pt(),sys);
  if(tf<0) tf=0.;
  tf*=GetBinContentUser(fFakeTF_l1,fabs(p.lepton1->Eta()),p.lepton1->Pt(),sys);
  if(tf<0) tf=0.;
  return tf;
}
double SMPAnalyzerCore::GetFakeRate(const Lepton* lep){
  if(!lep) return 0.;
  if(lep->LeptonFlavour()==Lepton::ELECTRON)
    return GetFakeRate(lep->LeptonFlavour(),lep->Eta(),lep->Pt());
  if(lep->LeptonFlavour()==Lepton::MUON)
    return GetFakeRate(lep->LeptonFlavour(),lep->Eta(),lep->Pt());
  return 0.;
  //return GetFakeRate(lep->LeptonFlavour(),lep->Eta(),lep->Pt()*(1+lep->RelIso()*TMath::Max(0.,TMath::Min(1.,(lep->Pt()-30)/30))));
}
double SMPAnalyzerCore::GetFakeRate(Lepton::Flavour flavour,double eta,double pt){
  TH2* fFakeRate=NULL;
  if(flavour==Lepton::ELECTRON) fFakeRate=fFakeRate_electron;
  else if(flavour==Lepton::MUON) fFakeRate=fFakeRate_muon;
  if(!fFakeRate) return 0.;
  eta=fabs(eta);
  double etamin=fFakeRate->GetXaxis()->GetBinLowEdge(1);
  double etamax=fFakeRate->GetXaxis()->GetBinUpEdge(fFakeRate->GetNbinsX());
  double ptmin=fFakeRate->GetYaxis()->GetBinLowEdge(1);
  double ptmax=fFakeRate->GetYaxis()->GetBinUpEdge(fFakeRate->GetNbinsY());
  if(eta<etamin) eta=etamin+1e-6;
  if(eta>etamax) eta=etamax-1e-6;
  if(pt<ptmin) pt=ptmin+1e-6;
  if(pt>ptmax) pt=ptmax-1e-6;
  eta=fFakeRate->GetXaxis()->GetBinCenter(fFakeRate->GetXaxis()->FindBin(eta));
  return fFakeRate->Interpolate(eta,pt);
}
void SMPAnalyzerCore::DeleteFakeRate(){
  if(fFakeRate_electron) delete fFakeRate_electron;
  if(fFakeRate_muon) delete fFakeRate_muon;

  for(auto& [key,hist]:fFakeTF){
    if(hist) delete hist;
  }
  fFakeTF.clear();
}

void SMPAnalyzerCore::PrintGens(const vector<Gen>& gens){
  cout<<"index\tpid\tmother\tstatus\tpropt\thard\n";
  for(int i=0;i<(int)gens.size();i++){
    gens[i].Print();
    cout<<gens.at(i).Index()<<"\t"<<gens.at(i).PID()<<"\t"<<gens.at(i).MotherIndex()<<"\t"<<gens.at(i).Status()<<"\t"<<gens.at(i).isPrompt()<<"\t"<<gens.at(i).isHardProcess()<<endl;
  }
}

double SMPAnalyzerCore::GetBinContentUser(TH1* hist,double valx,int sys){
  double xmin=hist->GetXaxis()->GetXmin();
  double xmax=hist->GetXaxis()->GetXmax();
  if(xmin>=0) valx=fabs(valx);
  if(valx<xmin) valx=xmin+0.001;
  if(valx>=xmax) valx=xmax-0.001;
  return hist->GetBinContent(hist->FindBin(valx))+sys*hist->GetBinError(hist->FindBin(valx));
}
double SMPAnalyzerCore::GetBinContentUser(TH2* hist,double valx,double valy,int sys){
  double xmin=hist->GetXaxis()->GetXmin();
  double xmax=hist->GetXaxis()->GetXmax();
  double ymin=hist->GetYaxis()->GetXmin();
  double ymax=hist->GetYaxis()->GetXmax();
  if(xmin>=0) valx=fabs(valx);
  if(valx<xmin) valx=xmin+0.001;
  if(valx>=xmax) valx=xmax-0.001;
  if(ymin>=0) valy=fabs(valy);
  if(valy<ymin) valy=ymin+0.001;
  if(valy>=ymax) valy=ymax-0.001;
  return hist->GetBinContent(hist->FindBin(valx,valy))+sys*hist->GetBinError(hist->FindBin(valx,valy));
}
double SMPAnalyzerCore::GetBinContentUser(TH3* hist,double valx,double valy,double valz,int sys){
  double xmin=hist->GetXaxis()->GetXmin();
  double xmax=hist->GetXaxis()->GetXmax();
  double ymin=hist->GetYaxis()->GetXmin();
  double ymax=hist->GetYaxis()->GetXmax();
  double zmin=hist->GetZaxis()->GetXmin();
  double zmax=hist->GetZaxis()->GetXmax();
  if(xmin>=0) valx=fabs(valx);
  if(valx<xmin) valx=xmin+0.001;
  if(valx>xmax) valx=xmax-0.001;
  if(ymin>=0) valy=fabs(valy);
  if(valy<ymin) valy=ymin+0.001;
  if(valy>ymax) valy=ymax-0.001;
  if(zmin>=0) valz=fabs(valz);
  if(valz<zmin) valz=zmin+0.001;
  if(valz>zmax) valz=zmax-0.001;
  return hist->GetBinContent(hist->FindBin(valx,valy,valz))+sys*hist->GetBinError(hist->FindBin(valx,valy,valz));
}
void SMPAnalyzerCore::GetAFBLHEParticles(const vector<LHE>& lhes,LHE& p0,LHE& p1,LHE& l0,LHE& l1,LHE& j0){
  if(!IsDYSample&&!MCSample.Contains("GamGamToLL")&&!MCSample.Contains("TTLL")){
    cout <<"[AFBAnalyzer::GetAFBLHEParticles] this is only for dilepton event"<<endl;
    exit(EXIT_FAILURE);
  }
  p0=LHE();
  p1=LHE();
  l0=LHE();
  l1=LHE();
  j0=LHE();
  if(!lhes.size()) return;
  for(int i=0;i<(int)lhes.size();i++){
    if(p0.ID()==0&&lhes[i].Status()==-1&&lhes[i].Eta()>0) p0=lhes[i];
    if(p1.ID()==0&&lhes[i].Status()==-1&&lhes[i].Eta()<0) p1=lhes[i];
    if(l0.ID()==0&&(abs(lhes[i].ID())==11||abs(lhes[i].ID())==13||abs(lhes[i].ID())==15)) l0=lhes[i];
    if(l0.ID()&&(abs(lhes[i].ID())==11||abs(lhes[i].ID())==13||abs(lhes[i].ID())==15)) l1=lhes[i];
    if(lhes[i].Status()==1)
      if(abs(lhes[i].ID())<=6||lhes[i].ID()==21)
	if(lhes[i].Pt()>j0.Pt()) j0=lhes[i];
  }
  if(p0.ID()==0||p1.ID()==0||l0.ID()==0||l1.ID()==0){
    cout <<"[AFBAnalyzer::GetLHEParticles] something is wrong"<<endl;
    exit(EXIT_FAILURE);
  }
  if(l0.Pt()<l1.Pt()){
    LHE temp=l0;
    l0=l1;
    l1=temp;
  }
}


void SMPAnalyzerCore::GetAFBGenParticles(const vector<Gen>& gens,Gen& parton0,Gen& parton1,Gen& l0,Gen& l1,int mode){
  //mode 0:bare 1:dressed01 2:dressed04 3:beforeFSR
  if(!IsDYSample&&!MCSample.Contains("GamGamToLL")&&!MCSample.Contains("TTLL")){
    cout <<"[SMPAnalyzerCore::GetAFBGenParticles] this is only for dilepton event"<<endl;
    exit(EXIT_FAILURE);
  }
  parton0=Gen();
  parton1=Gen();
  l0=Gen();
  l1=Gen();
  vector<const Gen*> leptons;
  vector<const Gen*> photons;
  int ngen=gens.size();
  for(int i=0;i<ngen;i++){
    if(!gens.at(i).isPrompt()) continue;
    int genpid=gens.at(i).PID();
    if(gens.at(i).isHardProcess()){
      if(abs(genpid)<7||genpid==21||genpid==22){
	if(parton0.IsEmpty()) parton0=gens[i];
	else if(parton1.IsEmpty()) parton1=gens[i];
      }
    }
    if(gens.at(i).Status()==1){
      if(abs(genpid)==11||abs(genpid)==13) leptons.push_back(&gens[i]);
      else if(gens.at(i).PID()==22) photons.push_back(&gens[i]);
    }
  }
  int nlepton=leptons.size();
  const double maxdr=0.4;
  for(int i=0;i<nlepton;i++){
    if(leptons[i]->PID()!=lhe_l0.ID()) continue;
    if(leptons[i]->DeltaR(lhe_l0)>maxdr) continue;
    if( fabs(leptons[i]->E()-lhe_l0.E()) < fabs(l0.E()-lhe_l0.E()) ){
      l0=*leptons[i];
    }
  }
  if(l0.PID()==0){
    for(int i=0;i<nlepton;i++){
      if(leptons[i]->PID()!=lhe_l0.ID()) continue;
      if(l0.PID()==0 || leptons[i]->DeltaR(lhe_l0)<l0.DeltaR(lhe_l0)){
	l0=*leptons[i];
      }
    }
  }
  for(int i=0;i<nlepton;i++){
    if(leptons[i]->PID()!=lhe_l1.ID()) continue;
    if(leptons[i]->DeltaR(lhe_l1)>maxdr) continue;
    if( fabs(leptons[i]->E()-lhe_l1.E()) < fabs(l1.E()-lhe_l1.E()) ){
      l1=*leptons[i];
    }
  }
  if(l1.PID()==0){
    for(int i=0;i<nlepton;i++){
      if(leptons[i]->PID()!=lhe_l1.ID()) continue;
      if(l1.PID()==0 || leptons[i]->DeltaR(lhe_l1)<l1.DeltaR(lhe_l1)){
	l1=*leptons[i];
      }
    }
  }
  if(l0.Pt()<l1.Pt()){
    Gen tmp=l0;
    l0=l1;
    l1=tmp;
  }
  if(mode>=3){
    if(nlepton>=4){
      for(int i=0;i<nlepton;i++){
	if(leptons[i]->Index()==l0.Index()||leptons[i]->Index()==l1.Index()) continue;
	for(int j=i+1;j<nlepton;j++){
	  if(leptons[j]->Index()==l0.Index()||leptons[j]->Index()==l1.Index()) continue;
	  if(!(leptons[i]->PID()+leptons[j]->PID()==0)) continue;
	  vector<int> history_i=TrackGenSelfHistory(*leptons[i],gens);
	  vector<int> history_j=TrackGenSelfHistory(*leptons[j],gens);
	  if(history_i.at(1)==history_j.at(1)){
	    photons.push_back(leptons[i]);
	    photons.push_back(leptons[j]);
	  }
	}
      }
    }	
    for(const auto& photon:photons){
      vector<int> history=TrackGenSelfHistory(*photon,gens);
      if(gens[history.at(1)].PID()==l0.PID()) l0+=*photon;
      else if(gens[history.at(1)].PID()==l1.PID()) l1+=*photon;
      else if(gens[history.at(1)].PID()==23){ // for minnlo+photos
	if(photon->DeltaR(l0)<photon->DeltaR(l1)) l0+=*photon;
	else l1+=*photon;
      }
    }    
  }else if(mode>=1){
    double delr=mode==1?0.1:0.4;
    for(const auto& photon:photons){
      if(l0.DeltaR(*photon)>delr&&l1.DeltaR(*photon)>delr) continue;
      if(l0.DeltaR(*photon)<l1.DeltaR(*photon)) l0+=*photon;
      else l1+=*photon;
    }
  }
}

Gen SMPAnalyzerCore::SMPGetGenMatchedLepton(const Lepton& lep,const std::vector<Gen>& gens,int mode){
  //0: default
  //1: dressed 0.1
  Gen gen_lepton=GetGenMatchedLepton(lep,gens);
  if(gen_lepton.IsEmpty()) return gen_lepton;
  if(mode==1){ // dressed 0.1 cone
    for(const auto& gen: gens){
      if(gen.Status()!=1) continue;
      if(gen.PID()!=22) continue;
      if(gen.DeltaR(gen_lepton)>0.1) continue;
      gen_lepton+=gen;
    }
  }
    
  return gen_lepton;
}

std::vector<Electron> SMPAnalyzerCore::SMPGetElectrons(TString id, double ptmin, double fetamax){
  std::vector<Electron> out;
  if(id=="passMediumID_SelQ"){
    std::vector<Electron> electrons = GetAllElectrons();
    for(unsigned int i=0; i<electrons.size(); i++){
      Electron this_electron= electrons.at(i);
      if(!( this_electron.Pt()>ptmin ))	continue;
      if(!( fabs(this_electron.scEta())<fetamax )) continue;
      if(!( this_electron.PassID("passMediumID") ))	continue;
      if(!electron_isGsfCtfScPixChargeConsistent->at(i)) continue;
      out.push_back(this_electron);
    }
  }else if(id=="passTightID_SelQ"){
    std::vector<Electron> electrons = GetAllElectrons();
    for(unsigned int i=0; i<electrons.size(); i++){
      Electron this_electron= electrons.at(i);
      if(!( this_electron.Pt()>ptmin ))	continue;
      if(!( fabs(this_electron.scEta())<fetamax )) continue;
      if(!( this_electron.PassID("passTightID") )) continue;
      if(!electron_isGsfCtfScPixChargeConsistent->at(i)) continue;
      out.push_back(this_electron);
    }
  }else if(id=="passMediumIDWithAntiIso"){
    vector<Electron> electrons = GetAllElectrons();
    for(unsigned int i=0; i<electrons.size(); i++){
      Electron el= electrons.at(i);
      if(!( el.Pt()>ptmin ))	continue;
      if(!( fabs(el.scEta())<fetamax )) continue;
      if( el.etaRegion()==Electron::GAP ) continue;
      if( fabs(el.scEta()) <= 1.479 ){
	if(! (el.Full5x5_sigmaIetaIeta() < 0.0106) ) continue;
	if(! (fabs(el.dEtaSeed()) < 0.0032) ) continue;
	if(! (fabs(el.dPhiIn()) < 0.0547) ) continue;
	if( (el.HoverE() < 0.046 + 1.16/el.scE() + 0.0324*el.Rho()/el.scE()) && (el.RelIso() < 0.0478+0.506/el.UncorrPt()) ) continue;
	if(! (fabs(el.InvEminusInvP()) < 0.184) ) continue;
	if(! (el.NMissingHits() <= 1) ) continue;
	if(! (el.PassConversionVeto()) ) continue;
      }else{
	if(! (el.Full5x5_sigmaIetaIeta() < 0.0387) ) continue;
	if(! (fabs(el.dEtaSeed()) < 0.00632) ) continue;
	if(! (fabs(el.dPhiIn()) <  0.0394 ) ) continue;
	if( (el.HoverE() < 0.0275 + 2.52/el.scE() + 0.183*el.Rho()/el.scE()) && (el.RelIso() < 0.0658+0.963/el.UncorrPt()) ) continue;
	if(! (fabs(el.InvEminusInvP()) < 0.0721) ) continue;
	if(! (el.NMissingHits() <= 1) ) continue;
	if(! (el.PassConversionVeto()) ) continue;
      }
      out.push_back(el);
    }
  }else if(id=="passAntiMediumID"){
    vector<Electron> electrons = GetAllElectrons();
    for(unsigned int i=0; i<electrons.size(); i++){
      Electron el= electrons.at(i);
      if(!( el.Pt()>ptmin ))	continue;
      if(!( fabs(el.scEta())<fetamax )) continue;
      if( el.etaRegion()==Electron::GAP ) continue;
      if( el.PassID("passMediumID") ) continue;
      out.push_back(el);
    }
  }else if(id=="passMediumIDSideBand"){
    vector<Electron> electrons = GetAllElectrons();
    for(unsigned int i=0; i<electrons.size(); i++){
      Electron el= electrons.at(i);
      if(!( el.Pt()>ptmin ))	continue;
      if(!( fabs(el.scEta())<fetamax )) continue;
      if( el.etaRegion()==Electron::GAP ) continue;
      if( el.PassID("passMediumID") ) continue;
      bool passtrigger=true;
      if(fabs(el.scEta())<=1.479){
	if(el.Full5x5_sigmaIetaIeta()>0.013) passtrigger=false;
	if(el.HoverE()>0.13*el.scE()) passtrigger=false;
	if(el.ecalPFClusterIso()>0.5*el.Et()) passtrigger=false;
	if(el.hcalPFClusterIso()>0.3*el.Et()) passtrigger=false;
	if(el.ecalPFClusterIso()>(0.5+0.29*el.Rho())*el.Et()) passtrigger=false;
	if(el.hcalPFClusterIso()>(0.3+0.2*el.Rho())*el.Et()) passtrigger=false;
	if(!el.PassConversionVeto()) passtrigger=false;
	if(el.dEtaSeed()>0.01) passtrigger=false;
	if(el.dPhiIn()>0.07) passtrigger=false;
	if(el.TrkIso()>0.2*el.Et()) passtrigger=false;						
      }else{
	if(el.Full5x5_sigmaIetaIeta()>0.035) passtrigger=false;
	if(el.HoverE()>0.13*el.scE()) passtrigger=false;
	if(el.ecalPFClusterIso()>(0.5+0.21*el.Rho())*el.Et()) passtrigger=false;
	if(el.hcalPFClusterIso()>(0.3+0.25*el.Rho())*el.Et()) passtrigger=false;
	if(!el.PassConversionVeto()) passtrigger=false;
	if(el.dEtaSeed()>0.015) passtrigger=false;
	if(el.dPhiIn()>0.1) passtrigger=false;
	if(el.TrkIso()>0.2*el.Et()) passtrigger=false;
      }
      if(!passtrigger) continue;
      out.push_back(el);
    }
  }else if(id=="passAntiLooseID"){
    vector<Electron> electrons = GetAllElectrons();
    for(unsigned int i=0; i<electrons.size(); i++){
      Electron el= electrons.at(i);
      if(!( el.Pt()>ptmin ))	continue;
      if(!( fabs(el.scEta())<fetamax )) continue;
      if( el.etaRegion()==Electron::GAP ) continue;
      if( el.PassID("passLooseID") ) continue;
      out.push_back(el);
    }
  }else out=GetElectrons(id,ptmin,fetamax);
  std::sort(out.begin(),out.end(),PtComparing);
  return out;
}    
std::vector<Muon> SMPAnalyzerCore::SMPGetMuons(TString id,double ptmin,double fetamax){
  vector<Muon> out;
  if(id=="POGTightWithLooseTrkIso"){
    vector<Muon> muons=GetMuons("POGTight",ptmin,fetamax);
    for(auto const& muon: muons){
      //if(muon.TrkIso()/muon.Pt()<0.1) out.push_back(muon);
      if(muon.PassSelector(Muon::Selector::TkIsoLoose)) out.push_back(muon);
    }
  }else if(id=="POGMediumWithLooseTrkIso"){
    vector<Muon> muons;
    if(DataEra=="2016preVFP") muons=GetMuons("POGMedium_hip",ptmin,fetamax);
    else muons=GetMuons("POGMedium",ptmin,fetamax);
    for(auto const& muon: muons){
      if(muon.PassSelector(Muon::Selector::TkIsoLoose)) out.push_back(muon);
    }
  }else if(id=="POGMediumWithLooseTrkIsoSideBand"){
    vector<Muon> muons;
    if(DataEra=="2016preVFP") muons=GetMuons("POGMedium_hip",ptmin,fetamax);
    else muons=GetMuons("POGMedium",ptmin,fetamax);
    for(auto const& muon: muons){
      if(muon.PassSelector(Muon::Selector::TkIsoLoose)) continue;
      if(muon.TrkIso()/muon.Pt()>0.4) continue;
      out.push_back(muon);
    }
  }else if(id=="POGMediumWithAntiLooseTrkIso"){
    vector<Muon> muons;
    if(DataEra=="2016preVFP") muons=GetMuons("POGMedium_hip",ptmin,fetamax);
    else muons=GetMuons("POGMedium",ptmin,fetamax);
    for(auto const& muon: muons){
      if(muon.PassSelector(Muon::Selector::TkIsoLoose)) continue;
      out.push_back(muon);
    }
  }else if(id=="NotMediumWithLooseTrkIso"){
    vector<Muon> muons=GetMuons("NOCUT",ptmin,fetamax);
    for(auto const& muon: muons){
      if(DataEra=="2016preVFP"){
	if(muon.PassID("POGMedium_hip")&&muon.PassSelector(Muon::Selector::TkIsoLoose)) continue;
      }else{
	if(muon.PassID("POGMedium")&&muon.PassSelector(Muon::Selector::TkIsoLoose)) continue;	
      }
      out.push_back(muon);
    }
  }else if(id=="POGMediumNoHipWithLooseTrkIso"){
    vector<Muon> muons=GetMuons("POGMedium_nohip",ptmin,fetamax);
    for(auto const& muon: muons){
      if(muon.PassSelector(Muon::Selector::TkIsoLoose)) out.push_back(muon);
    }
  }else if(id=="POGTightWithAntiIso"){
    vector<Muon> muons=GetMuons("POGTight",ptmin,fetamax);
    for(auto const& muon: muons){
      if(muon.RelIso()>0.3) out.push_back(muon);
    }
  }else if(id=="POGTightWithAntiMediumIso"){
    vector<Muon> muons=GetMuons("POGTight",ptmin,fetamax);
    for(auto const& muon: muons){
      if(muon.RelIso()>0.2) out.push_back(muon);
    }
  }else out=GetMuons(id,ptmin,fetamax);
  return out;
}
double SMPAnalyzerCore::MuonMomentumCorrection(const Muon& muon,int set,int member){
  double rc=1.;
  if(IsDATA){
    rc=roc->kScaleDT(muon.Charge(),muon.MiniAODPt(),muon.Eta(),muon.Phi(),set,member);
  }else{
    Gen gen=GetGenMatchedLepton(muon,gens);
    if(gen.IsEmpty()){
      gRandom->SetSeed((run<<15)+(lumi<<10)+(event<<5)+muon.Eta()*100);
      double u=gRandom->Rndm();
      rc=roc->kSmearMC(muon.Charge(),muon.MiniAODPt(),muon.Eta(),muon.Phi(),muon.TrackerLayers(),u,set,member);
    }else{
      rc=roc->kSpreadMC(muon.Charge(),muon.MiniAODPt(),muon.Eta(),muon.Phi(),gen.Pt(),set,member);
    }
  }
  return rc;
}

std::vector<Muon> SMPAnalyzerCore::MuonMomentumCorrection(const vector<Muon>& muons,int set,int member,bool sort){
  if(!roc) return std::vector<Muon>(muons);
  std::vector<Muon> out;
  for(auto muon:muons){
    if(set>=0){
      double rc=MuonMomentumCorrection(muon,set,member);
      muon.SetPtEtaPhiM(muon.MiniAODPt()*rc,muon.Eta(),muon.Phi(),muon.M());
    }else if(set==-2){
      double rc=MuonMomentumCorrection(muon,0,0);
      if(fRoccorResidual.find("muon_scale")!=fRoccorResidual.end() && fRoccorResidual["muon_scale"]){
	double eta=fabs(muon.Eta());
	if(eta>=2.4) eta=2.39;
	int ibin=fRoccorResidual["muon_scale"]->FindBin(eta);
	if(IsDATA){
	  double scale=1/(1+fRoccorResidual["muon_scale"]->GetBinContent(ibin));
	  rc*=scale;
	}
	gRandom->SetSeed((run<<15)+(lumi<<10)+(event<<5)+muon.Eta()*101);
	double u=gRandom->Rndm();
	double centralRes=fRoccorResidual["muon_centralRes"]->GetBinContent(ibin);
	double leftRes=fRoccorResidual["muon_leftRes"]->GetBinContent(ibin);
	double rightRes=fRoccorResidual["muon_rightRes"]->GetBinContent(ibin);
	if(u<fRoccorResidual["muon_leftFrac"]->GetBinContent(ibin) && (leftRes>0)^IsDATA ){
	  rc*=1-fabs(gRandom->Gaus(0,fabs(leftRes)));
	}else if(1-u<fRoccorResidual["muon_rightFrac"]->GetBinContent(ibin) && (rightRes>0)^IsDATA ){
	  rc*=1+fabs(gRandom->Gaus(0,fabs(rightRes)));
	}else if((centralRes>0)^IsDATA){
	  rc*=gRandom->Gaus(1,fabs(centralRes));
	}
      }
      muon.SetPtEtaPhiM(muon.MiniAODPt()*rc,muon.Eta(),muon.Phi(),muon.M());
    }
    out.push_back(muon);
  }
  if(sort) std::sort(out.begin(),out.end(),PtComparing);
  return out;
}
double SMPAnalyzerCore::ElectronEnergyCorrection(const Electron& electron,int set,int member){
  double rc=1.;
  //double rcerr=0.;
  double el_eta=electron.Eta();
  double el_phi=electron.Phi();
  if(IsDATA){
    rc=rocele->kScaleDT(electron.UncorrPt(),el_eta,el_phi,electron.R9(),run,set,member);
  }else{	
    Gen gen;
    if(!gen_l0_dressed.IsEmpty()&&gen_l0_dressed.DeltaR(electron)<0.1){
      gen=gen_l0_dressed;
    }else if(!gen_l1_dressed.IsEmpty()&&gen_l1_dressed.DeltaR(electron)<0.1){
      gen=gen_l1_dressed;
    }else{
      gen=SMPGetGenMatchedLepton(electron,gens,1);
    }
    gRandom->SetSeed((run<<15)+(lumi<<10)+(event<<5)+electron.Eta()*100);
    double u=gRandom->Rndm();
    if(!gen.IsEmpty()&&fabs(electron.Pt()/gen.Pt()-1.)<0.5){
      rc=rocele->kSpreadMC(electron.UncorrPt(),el_eta,el_phi,electron.R9(),u,gen.Pt(),set,member);
    }else{
      rc=rocele->kScaleMC(electron.UncorrPt(),el_eta,el_phi,electron.R9(),set,member);
    }
  }      
  if(TMath::IsNaN(rc)) rc=1.;
  return rc*electron.UncorrE()/electron.E();
}
std::vector<Electron> SMPAnalyzerCore::ElectronEnergyCorrection(const vector<Electron>& electrons,int set,int member,bool sort){
  if(!rocele) return std::vector<Electron>(electrons);
  std::vector<Electron> out;
  for(auto electron:electrons){
    if(set>=0){
      electron*=ElectronEnergyCorrection(electron,set,member);
    }else if(set==-1){ //no energe cor
      electron*=electron.UncorrE()/electron.E();
    }else if(set==-2){ //residual energe cor
      electron*=ElectronEnergyCorrection(electron,0,0);
      if(fRoccorResidual.find("electron_scale")!=fRoccorResidual.end() && fRoccorResidual["electron_scale"]){
	double eta=fabs(electron.Eta());
	if(eta>=2.4) eta=2.39;
	int ibin=fRoccorResidual["electron_scale"]->FindBin(eta);
	if(IsDATA){
	  double scale=1/(1+fRoccorResidual["electron_scale"]->GetBinContent(ibin));
	  electron*=scale;
	}
	gRandom->SetSeed((run<<15)+(lumi<<10)+(event<<5)+electron.Eta()*101);
	double u=gRandom->Rndm();
	double centralRes=fRoccorResidual["electron_centralRes"]->GetBinContent(ibin);
	double leftRes=fRoccorResidual["electron_leftRes"]->GetBinContent(ibin);
	double rightRes=fRoccorResidual["electron_rightRes"]->GetBinContent(ibin);
	if(u<fRoccorResidual["electron_leftFrac"]->GetBinContent(ibin) && (leftRes>0)^IsDATA ){
	  electron*=1-fabs(gRandom->Gaus(0,fabs(leftRes)));
	}else if(1-u<fRoccorResidual["electron_rightFrac"]->GetBinContent(ibin) && (rightRes>0)^IsDATA ){
	  electron*=1+fabs(gRandom->Gaus(0,fabs(rightRes)));
	}else if((centralRes>0)^IsDATA){
	  electron*=gRandom->Gaus(1,fabs(centralRes));
	}
      }
    }else{
      cout<<"[SMPAnalyzerCore::ElectronEnergyCorrection] wrong set "<<set<<endl;
      exit(ENODATA);
    }
    out.push_back(electron);
  }
  if(sort){
    std::sort(out.begin(),out.end(),PtComparing);
  }
  return out;
}
  
void SMPAnalyzerCore::FillCutflow(TString histname,TString label,double weight){
  TH1D* hist=NULL;
  auto it=maphist_TH1D.find(histname);
  if(it==maphist_TH1D.end()){
    hist=new TH1D(histname,"",1,0,1);
    hist->SetDirectory(NULL);
    maphist_TH1D[histname]=hist;    
  }else hist=it->second;

  hist->Fill(label,weight);
}
    
TString SMPAnalyzerCore::Replace(TString str,TRegexp reg,TString repl){
  int extent;
  int start=str.Index(reg,&extent);
  if(start>=0) return str.Replace(start,extent,repl);
  else return str;
}
map<TString,double> SMPAnalyzerCore::SelectWeights(map<TString,double> origin,vector<TString> keys){
  map<TString,double> out;
  for(const auto& key:keys){
    if(origin.find(key)!=origin.end()){
      out[key]=origin[key];
    }
  }
  return out;
}

SMPAnalyzerCore::Parameter::Parameter(){
}
SMPAnalyzerCore::Parameter::~Parameter(){
}
void SMPAnalyzerCore::Apply(Parameter& p,TString vsuf,unique_ptr<Variation>& v){
  p.vsuffix=vsuf;

  if(typeid(*v)==typeid(VariationWeight)){
    const VariationWeight* vv=static_cast<const VariationWeight*>(&(*v));
    p.weight=vv->weight;
  }else{
    p.weight=p.default_weight;
  }
  
  if(typeid(*v)==typeid(VariationMuonMomentum)){
    const VariationMuonMomentum* vv=static_cast<const VariationMuonMomentum*>(&(*v));
    if(p.muonmomentum_set!=vv->set||p.muonmomentum_mem!=vv->mem){
      p.muons=MuonMomentumCorrection(p.muons,vv->set,vv->mem,false);
      p.amuons=MuonMomentumCorrection(p.amuons,vv->set,vv->mem,false);
      p.SetLeptons();
      p.muonmomentum_set=vv->set;
      p.muonmomentum_mem=vv->mem;
    }
  }else{
    if(p.muonmomentum_set!=p.default_muonmomentum_set||p.muonmomentum_mem!=p.default_muonmomentum_mem){
      p.muons=MuonMomentumCorrection(p.muons,p.default_muonmomentum_set,p.default_muonmomentum_mem,false);
      p.amuons=MuonMomentumCorrection(p.amuons,p.default_muonmomentum_set,p.default_muonmomentum_mem,false);
      p.SetLeptons();
      p.muonmomentum_set=p.default_muonmomentum_set;
      p.muonmomentum_mem=p.default_muonmomentum_mem;
    }
  }

  if(typeid(*v)==typeid(VariationElectronEnergy)){
    const VariationElectronEnergy* vv=static_cast<const VariationElectronEnergy*>(&(*v));
    if(p.electronenergy_set!=vv->set||p.electronenergy_mem!=vv->mem){
      p.electrons=ElectronEnergyCorrection(p.electrons,vv->set,vv->mem,false);
      p.aelectrons=ElectronEnergyCorrection(p.aelectrons,vv->set,vv->mem,false);
      p.SetLeptons();
      p.electronenergy_set=vv->set;
      p.electronenergy_mem=vv->mem;
    }    
  }else{
    if(p.electronenergy_set!=p.default_electronenergy_set||p.electronenergy_mem!=p.default_electronenergy_mem){
      p.electrons=ElectronEnergyCorrection(p.electrons,p.default_electronenergy_set,p.default_electronenergy_mem,false);
      p.aelectrons=ElectronEnergyCorrection(p.aelectrons,p.default_electronenergy_set,p.default_electronenergy_mem,false);
      p.SetLeptons();
      p.electronenergy_set=p.default_electronenergy_set;
      p.electronenergy_mem=p.default_electronenergy_mem;
    }
  }

  if(typeid(*v)==typeid(VariationJES)){
    const VariationJES* vv=static_cast<const VariationJES*>(&(*v));
    if(p.JES_direction!=vv->direction){
      p.jets=ScaleJets(p.jets,vv->direction);
      p.bjets=ScaleJets(p.bjets,vv->direction);
      p.JES_direction=vv->direction;
    }    
  }else{
    if(p.JES_direction!=p.default_JES_direction){
      p.jets=ScaleJets(p.jets,p.default_JES_direction);
      p.bjets=ScaleJets(p.bjets,p.default_JES_direction);
      p.JES_direction=p.default_JES_direction;
    }
  }

  if(typeid(*v)==typeid(VariationJER)){
    const VariationJER* vv=static_cast<const VariationJER*>(&(*v));
    if(p.JER_direction!=vv->direction){
      p.jets=SmearJets(p.jets,vv->direction);
      p.bjets=SmearJets(p.bjets,vv->direction);
      p.JER_direction=vv->direction;
    }    
  }else{
    if(p.JER_direction!=p.default_JER_direction){
      p.jets=SmearJets(p.jets,p.default_JER_direction);
      p.bjets=SmearJets(p.bjets,p.default_JER_direction);
      p.JER_direction=p.default_JER_direction;
    }
  }
}
void SMPAnalyzerCore::Parameter::SetChannel(TString ch){
  vector<TString> availables={"el","ee","eE","Ee","EE","mu","mm","mM","Mm","MM","em","me","en","mn","ej","Ej","mj","Mj"};
  bool pass=false;
  for(const TString& avail:availables)
    if(ch==avail) pass=true;
  if(!pass){
    cout<<"[SMPAnalyzerCore::Parameter::SetChannel] not available channel "<<ch<<endl;
    exit(EXIT_FAILURE);
  }
  channel=ch;
  SetLeptons();
}
void SMPAnalyzerCore::Parameter::SetElectronKeys(TString elID,vector<TString> trig){
  k.electronRECOSF="Electron_RECO";
  k.electronIDSF=elID;
  k.electronIDSF2="";
  k.triggerSF=trig;
}
void SMPAnalyzerCore::Parameter::SetElectronKeys(TString elID,TString elID2,vector<TString> trig){
  k.electronRECOSF="Electron_RECO";
  k.electronIDSF=elID;
  k.electronIDSF2=elID2;
  k.triggerSF=trig;
}
void SMPAnalyzerCore::Parameter::SetMuonKeys(TString muID,TString muISO,vector<TString> trig){
  k.muonIDSF=muID;
  k.muonISOSF=muISO;
  k.triggerSF=trig;
}
void SMPAnalyzerCore::Parameter::SetLeptonPtCut(double l0pt,double l1pt){
  c.lepton0pt=l0pt;
  c.lepton1pt=l1pt;
}
void SMPAnalyzerCore::Parameter::SetLeptons(){
  leptons={};
  lepton0=NULL; 
  lepton1=NULL;
  truth_lepton0=Gen();
  truth_lepton1=Gen();
  if(channel=="") return;
  unsigned int ie=0,iae=0,im=0,iam=0;
  int nc=channel.Length();
  for(int ic=0;ic<nc;ic++){
    char c=channel[ic];
    if(c=='e'||c=='l'){
      if(electrons.size()>ie){
	leptons.push_back(&electrons.at(ie));
	ie++;
      }else leptons.push_back(NULL);
    }else if(c=='m'||c=='u'){
      if(muons.size()>im){
	leptons.push_back(&muons.at(im));
	im++;
      }else leptons.push_back(NULL);
    }else if(c=='E'){
      if(aelectrons.size()>iae){
	leptons.push_back(&aelectrons.at(iae));
	iae++;
      }else leptons.push_back(NULL);
    }else if(c=='M'){
      if(amuons.size()>iam){
	leptons.push_back(&amuons.at(iam));
	iam++;
      }else leptons.push_back(NULL);
    }else if (c=='n'){
      //neutrino
    }else if (c=='j'){
      //jet
    }
  }
  //should not sort for em or me channel
  //std::sort(leptons.begin(),leptons.end(),PtComparingPtr);
  if(leptons.size()>0) lepton0=leptons.at(0);
  if(leptons.size()>1) lepton1=leptons.at(1);
  if(lepton0) truth_lepton0=SMPGetGenMatchedLepton(*lepton0,gens);
  if(lepton1) truth_lepton1=SMPGetGenMatchedLepton(*lepton1,gens);
  if(lepton0&&lepton1){
    if(lepton0->Charge()*lepton1->Charge()>0){
      if(!hprefix.Contains("ss_")) hprefix+="ss_";
    }else{
      if(hprefix.Contains("ss_")) hprefix.ReplaceAll("ss_","");
    }
  }

}
void SMPAnalyzerCore::Parameter::SetGens(vector<Gen> gs){
  gens=gs;
  SetLeptons();
}
void SMPAnalyzerCore::Parameter::SetElectrons(vector<Electron> els){
  electrons=els;
  SetLeptons();
}    
void SMPAnalyzerCore::Parameter::SetMuons(vector<Muon> mus){
  muons=mus;
  SetLeptons();
}    
void SMPAnalyzerCore::Parameter::SetAElectrons(vector<Electron> els){
  aelectrons=els;
  SetLeptons();
}    
void SMPAnalyzerCore::Parameter::SetAMuons(vector<Muon> mus){
  amuons=mus;
  SetLeptons();
}
void SMPAnalyzerCore::AddVariationWeight(Variations& v,TString suffix,double weight){
  v[suffix]=std::make_unique<VariationWeight>(weight);
}
void SMPAnalyzerCore::AddVariationMuonMomentum(Variations& v,TString suffix,int set,int mem){
  v[suffix]=std::make_unique<VariationMuonMomentum>(set,mem);
}
void SMPAnalyzerCore::AddVariationElectronEnergy(Variations& v,TString suffix,int set,int mem){
  v[suffix]=std::make_unique<VariationElectronEnergy>(set,mem);
}
void SMPAnalyzerCore::AddVariationJES(Variations& v,TString suffix,int direction){
  v[suffix]=std::make_unique<VariationJES>(direction);
}
void SMPAnalyzerCore::AddVariationJER(Variations& v,TString suffix,int direction){
  v[suffix]=std::make_unique<VariationJER>(direction);
}

SMPAnalyzerCore::Parameter SMPAnalyzerCore::MakeParameter(TString channel,TString option){
  Parameter p;
  p.SetChannel(channel);
  p.option=option;
  p.SetGens(gens);
  p.hprefix="";
  p.w.lumiweight=reductionweight;
  p.w.PUweight=1;  p.w.PUweight_up=1;  p.w.PUweight_down=1;
  p.w.prefireweight=1;  p.w.prefireweight_up=1;  p.w.prefireweight_down=1;
  p.w.z0weight=1;
  p.w.zptweight=1;
  p.w.weakweight=1;
  p.w.topptweight=1;
  if(!IsDATA){
    p.w.lumiweight*=MCweight()*_event.GetTriggerLumi("Full");
    p.w.PUweight=mcCorr->GetPileUpWeight(nPileUp,0);
    p.w.PUweight_up=mcCorr->GetPileUpWeight(nPileUp,1);
    p.w.PUweight_down=mcCorr->GetPileUpWeight(nPileUp,-1);
    p.w.z0weight=GetZ0Weight(vertex_Z);
    p.w.prefireweight=L1PrefireReweight_Central;
    p.w.prefireweight_up=L1PrefireReweight_Up;
    p.w.prefireweight_down=L1PrefireReweight_Down;
    if(IsDYSample){
      if(abs(lhe_l0.ID())==11||abs(lhe_l0.ID())==13){
	TLorentzVector genZ=(gen_l0+gen_l1);
	p.w.zptweight=fZptCorrection->GetZptWeight(genZ.Pt(),genZ.Rapidity());
	p.w.zptweight_g=fZptCorrection->GetZptWeight(genZ.Pt());
	p.w.zptweight_gy=fZptCorrection->GetZptWeight(genZ.Pt(),genZ.Rapidity());
	p.w.zptweight_gym=fZptCorrection->GetZptWeight(genZ.Pt(),genZ.Rapidity(),genZ.M());
	p.w.weakweight=GetDYWeakWeight(genZ.M());
      }else p.hprefix+="tau_";
    }
    if(IsTTSample){
      p.w.topptweight=mcCorr->GetTopPtReweight(gens);
    }
  }

  p.prefix=p.channel+GetEraShort()+"/";
  
  if(p.channel=="mu"){
    p.k.muonTrackingSF="Muon_Tracking";
    p.k.muonRECOSF="Muon_RECO";
    p.SetMuonKeys("Muon_MediumID_trkIsoLoose","",{"IsoMu24_MediumID_trkIsoLoose"});
    p.SetMuons(MuonMomentumCorrection(SMPGetMuons("POGMediumWithLooseTrkIso",8.0,2.4),p.muonmomentum_set,p.muonmomentum_mem));
    p.SetLeptonPtCut(27,10);
    if(GetEraShort()=="2016a"){
      p.triggers={"HLT_IsoMu24_v","HLT_IsoTkMu24_v"};
    }else if(GetEraShort()=="2016b"){
      p.triggers={"HLT_IsoMu24_v","HLT_IsoTkMu24_v"};
    }else if(GetEraShort()=="2017"){
      p.triggers={"HLT_IsoMu24_v","HLT_IsoMu27_v"};
      p.k.triggerSF={"IsoMu24_MediumID_trkIsoLoose","IsoMu27_MediumID_trkIsoLoose"};
    }else if(GetEraShort()=="2018"){
      p.triggers={"HLT_IsoMu24_v"};
    }
  }else if(p.channel=="mm"){
    p.k.muonTrackingSF="Muon_Tracking";
    p.k.muonRECOSF="Muon_RECO";
    p.SetMuonKeys("Muon_MediumID_trkIsoLoose","",{"Mu17Leg1_MediumID_trkIsoLoose","Mu8Leg2_MediumID_trkIsoLoose"});
    p.SetMuons(MuonMomentumCorrection(SMPGetMuons("POGMediumWithLooseTrkIso",8.0,2.4),p.muonmomentum_set,p.muonmomentum_mem));
    p.SetLeptonPtCut(20,10);
    if(GetEraShort()=="2016a"){
      p.triggers={"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v","HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v","HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v","HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",};
    }else if(GetEraShort()=="2016b"){
      p.triggers={"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v","HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v","HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v","HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",};
      p.k.DZSF="DZ_MediumID_trkIsoLoose";
    }else if(GetEraShort()=="2017"){
      p.triggers={"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v"};
      p.k.DZSF="DZ_MediumID_trkIsoLoose";
    }else if(GetEraShort()=="2018"){
      p.triggers={"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"};
      p.k.DZSF="DZ_MediumID_trkIsoLoose";
    }
  }else if(p.option.Contains("SelQ")&&p.channel=="el"){
    p.prefix="selq/"+p.prefix;
    p.SetElectronKeys("Electron_MediumID","Electron_SelQ_MediumID",{"Ele27_SelQ_MediumID"});
    p.SetElectrons(ElectronEnergyCorrection(SMPGetElectrons("passMediumID_SelQ",8.0,2.5),p.electronenergy_set,p.electronenergy_mem));
    p.SetLeptonPtCut(30,10);
    if(GetEraShort()=="2016a"){
      p.triggers={"HLT_Ele27_WPTight_Gsf_v"};
    }else if(GetEraShort()=="2016b"){
      p.triggers={"HLT_Ele27_WPTight_Gsf_v"};
    }else if(GetEraShort()=="2017"){
      p.triggers={"HLT_Ele27_WPTight_Gsf_v","HLT_Ele32_WPTight_Gsf_v"};
      p.k.triggerSF={"Ele27_SelQ_MediumID","Ele32_SelQ_MediumID"};
    }else if(GetEraShort()=="2018"){
      p.triggers={"HLT_Ele28_WPTight_Gsf_v","HLT_Ele32_WPTight_Gsf_v"};
      p.k.triggerSF={"Ele28_SelQ_MediumID","Ele32_SelQ_MediumID"};
    }
  }else if(p.channel=="el"){
    p.SetElectronKeys("Electron_MediumID",{"Ele27_MediumID"});
    p.SetElectrons(ElectronEnergyCorrection(SMPGetElectrons("passMediumID",8.0,2.5),p.electronenergy_set,p.electronenergy_mem));
    p.SetLeptonPtCut(30,10);
    if(GetEraShort()=="2016a"){
      p.triggers={"HLT_Ele27_WPTight_Gsf_v"};
    }else if(GetEraShort()=="2016b"){
      p.triggers={"HLT_Ele27_WPTight_Gsf_v"};
    }else if(GetEraShort()=="2017"){
      p.triggers={"HLT_Ele27_WPTight_Gsf_v","HLT_Ele32_WPTight_Gsf_v"};
      p.k.triggerSF={"Ele27_MediumID","Ele32_MediumID"};
    }else if(GetEraShort()=="2018"){
      p.triggers={"HLT_Ele28_WPTight_Gsf_v","HLT_Ele32_WPTight_Gsf_v"};
      p.k.triggerSF={"Ele28_MediumID","Ele32_MediumID"};
    }
  }else if(p.channel=="ee"){
    p.SetElectronKeys("Electron_MediumID",{"Ele23Leg1_MediumID","Ele12Leg2_MediumID"});
    p.SetElectrons(ElectronEnergyCorrection(SMPGetElectrons("passMediumID",8.0,2.5),p.electronenergy_set,p.electronenergy_mem));
    p.SetLeptonPtCut(25,15);
    if(GetEraShort()=="2016a"){
      p.triggers={"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"};
      p.k.DZSF="DZ_MediumID";
    }else if(GetEraShort()=="2016b"){
      p.triggers={"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"};
      p.k.DZSF="DZ_MediumID";
    }else if(GetEraShort()=="2017") p.triggers={"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"};
    else if(GetEraShort()=="2018") p.triggers={"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"};
  }else if(p.channel=="me"){
    p.SetMuonKeys("Muon_MediumID_trkIsoLoose","",{"IsoMu24_MediumID_trkIsoLoose"});
    p.k.electronIDSF="Electron_MediumID";
    p.SetMuons(MuonMomentumCorrection(SMPGetMuons("POGMediumWithLooseTrkIso",8.0,2.4),p.muonmomentum_set,p.muonmomentum_mem));
    p.SetElectrons(ElectronEnergyCorrection(SMPGetElectrons("passMediumID",8.0,2.5),p.electronenergy_set,p.electronenergy_mem));
    p.SetLeptonPtCut(27,10);
    if(GetEraShort()=="2016a"){
      p.triggers={"HLT_IsoMu24_v","HLT_IsoTkMu24_v"};
    }else if(GetEraShort()=="2016b"){
      p.triggers={"HLT_IsoMu24_v","HLT_IsoTkMu24_v"};
    }else if(GetEraShort()=="2017"){
      p.triggers={"HLT_IsoMu24_v","HLT_IsoMu27_v"};
      p.k.triggerSF={"IsoMu24_MediumID_trkIsoLoose","IsoMu27_MediumID_trkIsoLoose"};
    }else if(GetEraShort()=="2018"){
      p.triggers={"HLT_IsoMu24_v"};
    }
  }else if(p.channel=="em"){
    p.SetElectronKeys("Electron_MediumID",{"Ele27_MediumID"});
    p.k.muonIDSF="Muon_MediumID_trkIsoLoose";
    p.SetElectrons(ElectronEnergyCorrection(SMPGetElectrons("passMediumID",8.0,2.5),p.electronenergy_set,p.electronenergy_mem));
    p.SetMuons(MuonMomentumCorrection(SMPGetMuons("POGMediumWithLooseTrkIso",8.0,2.4),p.muonmomentum_set,p.muonmomentum_mem));
    p.SetLeptonPtCut(30,10);
    if(GetEraShort()=="2016a"){
      p.triggers={"HLT_Ele27_WPTight_Gsf_v"};
    }else if(GetEraShort()=="2016b"){
      p.triggers={"HLT_Ele27_WPTight_Gsf_v"};
    }else if(GetEraShort()=="2017"){
      p.triggers={"HLT_Ele27_WPTight_Gsf_v","HLT_Ele32_WPTight_Gsf_v"};
      p.k.triggerSF={"Ele27_MediumID","Ele32_MediumID"};
    }else if(GetEraShort()=="2018"){
      p.triggers={"HLT_Ele28_WPTight_Gsf_v","HLT_Ele32_WPTight_Gsf_v"};
      p.k.triggerSF={"Ele28_MediumID","Ele32_MediumID"};
    }
  }else if(p.channel=="mM"||p.channel=="Mm"){
    p.k.muonIDSF="Muon_MediumID_trkIsoLoose";
    p.SetMuons(MuonMomentumCorrection(SMPGetMuons("POGMediumWithLooseTrkIso",8.0,2.4),p.muonmomentum_set,p.muonmomentum_mem));
    p.SetAMuons(MuonMomentumCorrection(SMPGetMuons("POGMediumWithLooseTrkIsoSideBand",8.0,2.4),p.muonmomentum_set,p.muonmomentum_mem));
    //p.SetAMuons(MuonMomentumCorrection(SMPGetMuons("NotMediumWithLooseTrkIso",8.0,2.4),p.muonmomentum_set,p.muonmomentum_mem));
    p.SetLeptonPtCut(20,10);
    //p.c.nmuonmax=1;
    p.option+=" triggermatching strictorder";
    if(GetEraShort()=="2016a"){
      p.triggers={"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v","HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v","HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v","HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",};
    }else if(GetEraShort()=="2016b"){
      p.triggers={"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v","HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v","HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v","HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",};
    }else if(GetEraShort()=="2017"){
      p.triggers={"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v"};
    }else if(GetEraShort()=="2018"){
      p.triggers={"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"};
    }
    /*
    p.SetMuonKeys("Muon_MediumID_trkIsoLoose","",{"IsoMu24_MediumID_trkIsoLoose"});
    p.SetMuons(MuonMomentumCorrection(SMPGetMuons("POGMediumWithLooseTrkIso",8.0,2.4),p.muonmomentum_set,p.muonmomentum_mem));
    p.SetAMuons(MuonMomentumCorrection(SMPGetMuons("POGMediumWithLooseTrkIsoSideBand",8.0,2.4),p.muonmomentum_set,p.muonmomentum_mem));
    //p.SetAMuons(MuonMomentumCorrection(SMPGetMuons("NotMediumWithLooseTrkIso",8.0,2.4),p.muonmomentum_set,p.muonmomentum_mem));
    p.SetLeptonPtCut(27,10);
    if(GetEraShort()=="2016a"){
      p.triggers={"HLT_IsoMu24_v","HLT_IsoTkMu24_v"};
    }else if(GetEraShort()=="2016b"){
      p.triggers={"HLT_IsoMu24_v","HLT_IsoTkMu24_v"};
    }else if(GetEraShort()=="2017"){
      p.triggers={"HLT_IsoMu24_v","HLT_IsoMu27_v"};
      p.k.triggerSF={"IsoMu24_MediumID_trkIsoLoose","IsoMu27_MediumID_trkIsoLoose"};
    }else if(GetEraShort()=="2018"){
      p.triggers={"HLT_IsoMu24_v"};
    }
    */
  }else if(p.channel=="MM"){
    p.k.muonIDSF="Muon_MediumID_trkIsoLoose";
    p.SetMuons(MuonMomentumCorrection(SMPGetMuons("POGMediumWithLooseTrkIso",8.0,2.4),p.muonmomentum_set,p.muonmomentum_mem));
    p.SetAMuons(MuonMomentumCorrection(SMPGetMuons("POGMediumWithLooseTrkIsoSideBand",8.0,2.4),p.muonmomentum_set,p.muonmomentum_mem));
    //p.SetAMuons(MuonMomentumCorrection(SMPGetMuons("NotMediumWithLooseTrkIso",8.0,2.4),p.muonmomentum_set,p.muonmomentum_mem));
    p.SetLeptonPtCut(20,10);
    //p.c.nmuonmax=0;
    p.option+=" triggermatching strictorder";
    if(GetEraShort()=="2016a"){
      p.triggers={"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v","HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v","HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v","HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",};
    }else if(GetEraShort()=="2016b"){
      p.triggers={"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v","HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v","HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v","HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",};
    }else if(GetEraShort()=="2017"){
      p.triggers={"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v"};
    }else if(GetEraShort()=="2018"){
      p.triggers={"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"};
    }
  }else if(p.channel=="eE"||p.channel=="Ee"){
    p.k.electronIDSF="Electron_MediumID";
    p.SetElectrons(ElectronEnergyCorrection(SMPGetElectrons("passMediumID",8.0,2.5),p.electronenergy_set,p.electronenergy_mem));
    p.SetAElectrons(ElectronEnergyCorrection(SMPGetElectrons("passMediumIDSideBand",8.0,2.5),p.electronenergy_set,p.electronenergy_mem));
    p.SetLeptonPtCut(25,15);
    //p.c.nelectronmax=1;
    p.option+=" triggermatching strictorder";
    if(GetEraShort()=="2016a") p.triggers={"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"};
    else if(GetEraShort()=="2016b") p.triggers={"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"};
    else if(GetEraShort()=="2017") p.triggers={"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"};
    else if(GetEraShort()=="2018") p.triggers={"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"};
    /*
    p.SetElectronKeys("Electron_MediumID",{"Ele27_MediumID"});
    p.SetElectrons(ElectronEnergyCorrection(SMPGetElectrons("passMediumID",8.0,2.5),p.electronenergy_set,p.electronenergy_mem));
    p.SetAElectrons(ElectronEnergyCorrection(SMPGetElectrons("passMediumIDSideBand",8.0,2.5),p.electronenergy_set,p.electronenergy_mem));
    p.SetLeptonPtCut(30,10);
    if(GetEraShort()=="2016a"){
      p.triggers={"HLT_Ele27_WPTight_Gsf_v"};
    }else if(GetEraShort()=="2016b"){
      p.triggers={"HLT_Ele27_WPTight_Gsf_v"};
    }else if(GetEraShort()=="2017"){
      p.triggers={"HLT_Ele27_WPTight_Gsf_v","HLT_Ele32_WPTight_Gsf_v"};
      p.k.triggerSF={"Ele27_MediumID","Ele32_MediumID"};
    }else if(GetEraShort()=="2018"){
      p.triggers={"HLT_Ele28_WPTight_Gsf_v","HLT_Ele32_WPTight_Gsf_v"};
      p.k.triggerSF={"Ele28_MediumID","Ele32_MediumID"};
    }
    */
  }else if(p.channel=="EE"){
    p.k.electronIDSF="Electron_MediumID";
    p.SetElectrons(ElectronEnergyCorrection(SMPGetElectrons("passMediumID",8.0,2.5),p.electronenergy_set,p.electronenergy_mem));
    p.SetAElectrons(ElectronEnergyCorrection(SMPGetElectrons("passMediumIDSideBand",8.0,2.5),p.electronenergy_set,p.electronenergy_mem));
    //p.c.nelectronmax=0;
    p.option+=" triggermatching strictorder";
    p.SetLeptonPtCut(25,15);
    if(GetEraShort()=="2016a") p.triggers={"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"};
    else if(GetEraShort()=="2016b") p.triggers={"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"};
    else if(GetEraShort()=="2017") p.triggers={"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"};
    else if(GetEraShort()=="2018") p.triggers={"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"};
  }else if(p.channel=="mn"){
    p.SetMuonKeys("Muon_MediumID_trkIsoLoose","",{"IsoMu24_MediumID_trkIsoLoose"});
    p.SetMuons(MuonMomentumCorrection(SMPGetMuons("POGMediumWithLooseTrkIso",8.0,2.4),p.muonmomentum_set,p.muonmomentum_mem));
    p.SetLeptonPtCut(27,-1);
    p.c.nleptonmin=1;
    if(GetEraShort()=="2016a"){
      p.triggers={"HLT_IsoMu24_v","HLT_IsoTkMu24_v"};
    }else if(GetEraShort()=="2016b"){
      p.triggers={"HLT_IsoMu24_v","HLT_IsoTkMu24_v"};
    }else if(GetEraShort()=="2017"){
      p.triggers={"HLT_IsoMu24_v","HLT_IsoMu27_v"};
      p.k.triggerSF={"IsoMu24_MediumID_trkIsoLoose","IsoMu27_MediumID_trkIsoLoose"};
    }else if(GetEraShort()=="2018"){
      p.triggers={"HLT_IsoMu24_v"};
    }
  }else if(p.channel=="en"){
    p.SetElectronKeys("Electron_MediumID",{"Ele27_MediumID"});
    p.SetElectrons(ElectronEnergyCorrection(SMPGetElectrons("passMediumID",8.0,2.5),p.electronenergy_set,p.electronenergy_mem));
    p.SetLeptonPtCut(30,-1);
    p.c.nleptonmin=1;
    if(GetEraShort()=="2016a"){
      p.triggers={"HLT_Ele27_WPTight_Gsf_v"};
    }else if(GetEraShort()=="2016b"){
      p.triggers={"HLT_Ele27_WPTight_Gsf_v"};
    }else if(GetEraShort()=="2017"){
      p.triggers={"HLT_Ele27_WPTight_Gsf_v","HLT_Ele32_WPTight_Gsf_v"};
      p.k.triggerSF={"Ele27_MediumID","Ele32_MediumID"};
    }else if(GetEraShort()=="2018"){
      p.triggers={"HLT_Ele28_WPTight_Gsf_v","HLT_Ele32_WPTight_Gsf_v"};
      p.k.triggerSF={"Ele28_MediumID","Ele32_MediumID"};
    }
  }else if(p.channel=="ej"||p.channel=="Ej"){
    p.k.electronIDSF="Electron_MediumID";
    p.SetElectrons(ElectronEnergyCorrection(SMPGetElectrons("passMediumID",8.0,2.5),p.electronenergy_set,p.electronenergy_mem));
    p.SetAElectrons(ElectronEnergyCorrection(SMPGetElectrons("passMediumIDSideBand",8.0,2.5),p.electronenergy_set,p.electronenergy_mem));
    p.SetLeptonPtCut(10,-1);
    p.c.nleptonmin=1;
    if(p.channel=="ej") p.c.nelectronmax=1;
    else if(p.channel=="Ej") p.c.nelectronmax=0;
    p.option+=" strictorder";
    if(GetEraShort()=="2016a") p.triggers={"HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v","HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v","HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v"};
    else if(GetEraShort()=="2016b") p.triggers={"HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v","HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v","HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v"};
    else if(GetEraShort()=="2017") p.triggers={"HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v","HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v","HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v"};
    else if(GetEraShort()=="2018") p.triggers={"HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v","HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v","HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v"};
  }else if(p.channel=="mj"||p.channel=="Mj"){
    p.k.muonIDSF="Muon_MediumID_trkIsoLoose";
    p.SetMuons(MuonMomentumCorrection(SMPGetMuons("POGMediumWithLooseTrkIso",8.0,2.4),p.muonmomentum_set,p.muonmomentum_mem));
    p.SetAMuons(MuonMomentumCorrection(SMPGetMuons("POGMediumWithLooseTrkIsoSideBand",8.0,2.4),p.muonmomentum_set,p.muonmomentum_mem));
    p.SetLeptonPtCut(10,-1);
    p.c.nleptonmin=1;
    if(p.channel=="mj") p.c.nmuonmax=1;
    else if(p.channel=="Mj") p.c.nmuonmax=0;    
    p.option+=" strictorder triggermatching";
    p.triggers={};
    if(GetEraShort()=="2016a"){
      if(!IsDATA||DataStream.Contains("DoubleMuon")){
	//p.triggers.push_back("HLT_Mu3_PFJet40_v");
	p.triggers.push_back("HLT_Mu8_TrkIsoVVL_v");
	p.triggers.push_back("HLT_Mu17_TrkIsoVVL_v");
      }
    }else if(GetEraShort()=="2016b"){
      if(!IsDATA||DataStream.Contains("DoubleMuon")){
	//p.triggers.push_back("HLT_Mu3_PFJet40_v");
	p.triggers.push_back("HLT_Mu8_TrkIsoVVL_v");
	p.triggers.push_back("HLT_Mu17_TrkIsoVVL_v");
      }
    }else if(GetEraShort()=="2017"){
      if(!IsDATA||DataStream.Contains("SingleMuon")){
	//p.triggers.push_back("HLT_Mu3_PFJet40_v");
	//p.triggers.push_back("HLT_Mu50_v");
      }
      if(!IsDATA||DataStream.Contains("DoubleMuon")){
	p.triggers.push_back("HLT_Mu8_TrkIsoVVL_v");
	p.triggers.push_back("HLT_Mu17_TrkIsoVVL_v");
      }
    }else if(GetEraShort()=="2018"){
      if(!IsDATA||DataStream.Contains("SingleMuon")){
	//p.triggers.push_back("HLT_Mu3_PFJet40_v");
	//p.triggers.push_back("HLT_Mu50_v");
      }
      if(!IsDATA||DataStream.Contains("DoubleMuon")){
	p.triggers.push_back("HLT_Mu8_TrkIsoVVL_v");
	p.triggers.push_back("HLT_Mu17_TrkIsoVVL_v");
      }
    }
  }

  p.jets.clear();
  if(p.option.Contains("jet_scale_up")){
    p.suffix+="_jet_scale_up";
    p.jets=SelectJets(ScaleJets(GetAllJets(),1),"tightLepVeto",p.c.jetpt*0.8,2.4);
  }else if(p.option.Contains("jet_scale_down")){
    p.suffix+="_jet_scale_down";
    p.jets=SelectJets(ScaleJets(GetAllJets(),-1),"tightLepVeto",p.c.jetpt*0.8,2.4);
  }else if(p.option.Contains("jet_smear_up")){
    p.suffix+="_jet_smear_up";
    p.jets=SelectJets(SmearJets(GetAllJets(),1),"tightLepVeto",p.c.jetpt*0.8,2.4);
  }else if(p.option.Contains("jet_smear_down")){
    p.suffix+="_jet_smear_down";
    p.jets=SelectJets(SmearJets(GetAllJets(),-1),"tightLepVeto",p.c.jetpt*0.8,2.4);
  }else{
    p.jets=SelectJets(GetAllJets(),"tightLepVeto",p.c.jetpt*0.8,2.4);
  }    
  std::sort(p.jets.begin(),p.jets.end(),PtComparing);
  JetTagging::Parameters jtp;
  if(p.option.Contains("DeepCSV::Medium")){
    jtp=JetTagging::Parameters(JetTagging::DeepCSV,JetTagging::Medium,JetTagging::incl,JetTagging::comb);
  }else if(p.option.Contains("DeepCSV")){
    jtp=JetTagging::Parameters(JetTagging::DeepCSV,JetTagging::Tight,JetTagging::incl,JetTagging::comb);
  }else if(p.option.Contains("DeepJet::Medium")){
    jtp=JetTagging::Parameters(JetTagging::DeepJet,JetTagging::Medium,JetTagging::incl,JetTagging::comb);
  }else if(p.option.Contains("DeepJet::Tight::mujets")){
    jtp=JetTagging::Parameters(JetTagging::DeepJet,JetTagging::Tight,JetTagging::incl,JetTagging::mujets);
  }else{
    jtp=JetTagging::Parameters(JetTagging::DeepJet,JetTagging::Tight,JetTagging::incl,JetTagging::comb);
  }
  p.bjets.clear();
  vector<Muon> allmuons=GetAllMuons();
  vector<Electron> allelectrons=GetAllElectrons();
  for(const auto& jet:p.jets){
    if(jet.GetTaggerResult(jtp.j_Tagger) < mcCorr->GetJetTaggingCutValue(jtp.j_Tagger, jtp.j_WP)) continue;
    if(!p.option.Contains("nobjetcleaning")){
      if(p.lepton0&&jet.DeltaR(*p.lepton0)<0.4) continue;
      if(p.lepton1&&jet.DeltaR(*p.lepton1)<0.4) continue;
    }
    Jet bjet=jet;
    double jetCharge = bjet.Charge();

    vector<Muon> bmuon;
    vector<Electron> belectron;

    for(int l=0,n=allmuons.size(); l<n; l++){
      if(allmuons.at(l).P()*sin(allmuons.at(l).Angle(bjet.Vect())) <0.6) continue; // original, 1GeV
      if(allmuons.at(l).TrkIso()/allmuons.at(l).Pt() <0.05) continue; // original, 0.1
      if(abs(allmuons.at(l).IP3D())/allmuons.at(l).IP3Derr() <2.) continue; // original, 2.5
      if(bjet.DeltaR(allmuons.at(l))<0.4) bmuon.push_back(allmuons.at(l));
    }

    //belectron Trial
    for(unsigned int l=0; l<allelectrons.size(); l++){
      if(allelectrons.at(l).P()*sin(allelectrons.at(l).Angle(bjet.Vect())) <0.6) continue;
      if(allelectrons.at(l).ecalPFClusterIso()/allelectrons.at(l).Pt() == 0.) continue;
      if(abs(allelectrons.at(l).IP3D())/allelectrons.at(l).IP3Derr() <2.0) continue;
      if(!allelectrons.at(l).IsGsfCtfScPixChargeConsistent()) continue;
      if(bjet.DeltaR(allelectrons.at(l))<0.4) belectron.push_back(allelectrons.at(l));
    }

    //The jet has soft muon inside, and its charge will determine the jet charge
    if(bmuon.size() > 0) jetCharge += 2 * bmuon.at(0).Charge();
    else if(belectron.size() > 0) jetCharge += 4 * belectron.at(0).Charge();
    bjet.userFloat["AFBCharge"]=jetCharge;
    p.bjets.push_back(bjet);
  }
  p.w.btagSF=1.;
  p.w.btagSF_hup=1.;
  p.w.btagSF_hdown=1.;
  p.w.btagSF_hcorr=1.;
  p.w.btagSF_huncorr=1.;
  p.w.btagSF_lup=1.;
  p.w.btagSF_ldown=1.;
  p.w.btagSF_lcorr=1.;
  p.w.btagSF_luncorr=1.;
  if(!IsDATA){
    p.w.btagSF=mcCorr->GetBTaggingReweight_1a(p.jets,jtp);
    p.w.btagSF_hup=mcCorr->GetBTaggingReweight_1a(p.jets,jtp,"SystUpHTag");
    p.w.btagSF_hdown=mcCorr->GetBTaggingReweight_1a(p.jets,jtp,"SystDownHTag");
    p.w.btagSF_hcorr=mcCorr->GetBTaggingReweight_1a(p.jets,jtp,"SystUpHTagCorr");
    p.w.btagSF_huncorr=mcCorr->GetBTaggingReweight_1a(p.jets,jtp,"SystUpHTagUnCorr");
    p.w.btagSF_lup=mcCorr->GetBTaggingReweight_1a(p.jets,jtp,"SystUpLTag");
    p.w.btagSF_ldown=mcCorr->GetBTaggingReweight_1a(p.jets,jtp,"SystDownLTag");
    p.w.btagSF_lcorr=mcCorr->GetBTaggingReweight_1a(p.jets,jtp,"SystUpLTagCorr");
    p.w.btagSF_luncorr=mcCorr->GetBTaggingReweight_1a(p.jets,jtp,"SystUpLTagUnCorr");
    //cout<<p.w.btagSF<<" "<<p.w.btagSF_hup<<" "<<p.w.btagSF_hcorr<<" "<<p.w.btagSF_huncorr<<endl;
    p.w.bchargeSF=GetBchargeSF(p);
    p.w.bchargeSF_s0m0=GetBchargeSF(p,0,0);
    p.w.bchargeSF_s0m1=GetBchargeSF(p,0,1);
  }
  return p;
}
vector<TString> SMPAnalyzerCore::Split(TString s,TString del){
  TObjArray* array=s.Tokenize(del);
  vector<TString> out;
  for(const auto& obj:*array){
    out.push_back(((TObjString*)obj)->String());
  }
  array->Delete();
  return out;
}
double SMPAnalyzerCore::GetPFMET_T1Smear() const {
  if(isnan(pfMET_Type1_pt)) return pfMET_Type1_pt;
  if(IsDATA) return pfMET_Type1_pt;
  TLorentzVector out(pfMET_Type1_pt*cos(pfMET_Type1_phi),pfMET_Type1_pt*sin(pfMET_Type1_phi),0,0);
  for(unsigned int i=0;i<jet_pt->size();i++){
    if(jet_neutralEmEnergyFraction->at(i)+jet_chargedEmEnergyFraction->at(i)>0.9) continue;
    if(fabs(jet_eta->at(i))>9.9) continue;
    TLorentzVector jet;
    jet.SetPtEtaPhiM(jet_pt->at(i), jet_eta->at(i), jet_phi->at(i), jet_m->at(i));
    jet*=(1-jet_muonEnergyFraction->at(i));
    if(jet.Pt()<15) continue;
    TLorentzVector jet_smear=jet*jet_smearedRes->at(i);
    out-=jet_smear-jet;
  }
  return out.Pt();
}
TString SMPAnalyzerCore::GetSkimName() const {
  TString skimname="";
  if(fChain->GetListOfFiles()->GetEntries()){
    TString filename=fChain->GetListOfFiles()->At(0)->GetTitle();
    if(filename.Contains("SkimTree_")){
      skimname=((TObjString*)TPRegexp("SkimTree_([^_/]*)").MatchS(filename)->At(1))->GetString();
    }
  }else{
    cout<<"[SMPAnalyzerCore::GetSkimName] no input file"<<endl;
    exit(EXIT_FAILURE);
  }
  return skimname;
}
void SMPAnalyzerCore::SetupL1PrefiringWeight(){
  TString datapath=(TString)getenv("DATA_DIR")+"/"+GetEra()+"/SMP/";
  TString L1PrefiringMaps="L1PrefiringMaps.root";
  TString L1MuonPrefiringParametriations="L1MuonPrefiringParametriations.root";
  if(IsExists(datapath+L1PrefiringMaps)&&IsExists(datapath+L1MuonPrefiringParametriations)){
    cout<<"[SMPAnalyzerCore::SetupL1PrefiringWeight] using file "+datapath+L1PrefiringMaps<<endl;
    cout<<"[SMPAnalyzerCore::SetupL1PrefiringWeight] using file "+datapath+L1MuonPrefiringParametriations<<endl;
  }else{
    cout<<"[SMPAnalyzerCore::SetupL1PrefiringWeight] no "+datapath+L1PrefiringMaps+" or "+datapath+L1MuonPrefiringParametriations<<endl;
    return;
  }
  DeleteL1PrefiringWeight();
  if(DataYear<2018){
    TString era=DataEra;
    if(era=="2017") era="2017BtoF";
    TFile f(datapath+L1PrefiringMaps);
    fL1Prefiring_photon=(TH2*)f.Get("L1prefiring_photonptvseta_UL"+era);
    fL1Prefiring_jet=(TH2*)f.Get("L1prefiring_jetptvseta_UL"+era);
    if(!fL1Prefiring_photon){
      cout<<"[SMPAnalyzerCore::SetupL1PrefiringWeight] no hist L1prefiring_photonptvseta_UL"+era<<endl;
      exit(ENODATA);
    }
    if(!fL1Prefiring_jet){
      cout<<"[SMPAnalyzerCore::SetupL1PrefiringWeight] no hist L1prefiring_jetptvseta_UL"+era<<endl;
      exit(ENODATA);
    }
    fL1Prefiring_photon->SetDirectory(0);
    fL1Prefiring_jet->SetDirectory(0);
  }
  {
    TString era=DataEra;
    if(era=="2017"||era=="2018") era="20172018";
    TFile f(datapath+L1MuonPrefiringParametriations);
    fL1Prefiring_muon[0]=(TF1*)f.Get("L1prefiring_muonparam_0.0To0.2_"+era);
    fL1Prefiring_muon[1]=(TF1*)f.Get("L1prefiring_muonparam_0.2To0.3_"+era);
    fL1Prefiring_muon[2]=(TF1*)f.Get("L1prefiring_muonparam_0.3To0.55_"+era);
    fL1Prefiring_muon[3]=(TF1*)f.Get("L1prefiring_muonparam_0.55To0.83_"+era);
    fL1Prefiring_muon[4]=(TF1*)f.Get("L1prefiring_muonparam_0.83To1.24_"+era);
    fL1Prefiring_muon[5]=(TF1*)f.Get("L1prefiring_muonparam_1.24To1.4_"+era);
    fL1Prefiring_muon[6]=(TF1*)f.Get("L1prefiring_muonparam_1.4To1.6_"+era);
    fL1Prefiring_muon[7]=(TF1*)f.Get("L1prefiring_muonparam_1.6To1.8_"+era);
    fL1Prefiring_muon[8]=(TF1*)f.Get("L1prefiring_muonparam_1.8To2.1_"+era);
    fL1Prefiring_muon[9]=(TF1*)f.Get("L1prefiring_muonparam_2.1To2.25_"+era);
    fL1Prefiring_muon[10]=(TF1*)f.Get("L1prefiring_muonparam_2.25To2.4_"+era);
    if(era.Contains("2016")){
      fL1Prefiring_muon[11]=(TF1*)f.Get("L1prefiring_muonparam_HotSpot_"+era);
    }
    for(int i=0;i<11+era.Contains("2016")?1:0;i++){
      if(!fL1Prefiring_muon[i]){
	cout<<"[SMPAnalyzerCore::SetupL1PrefiringWeight] no hist for fL1Prefiring_muon["<<i<<"]"<<endl;
	exit(ENODATA);
      }
    }
  }
  {
    TString infile=(TString)getenv("SKFlat_WD")+"/external/RocPFProb/prefiring_table_v1.txt";
    if(IsExists(infile)){
      cout<<"[SMPAnalyzerCore::SetupL1PrefiringWeight] using file "+infile<<endl;
      rocpfprob=new RocPFProb(infile.Data());
    }else{
      cout<<"[SMPAnalyzerCore::SetupL1PrefiringWeight] no "+infile<<endl;
    }
  }
  {
    TString infile=(TString)getenv("SKFlat_WD")+"/external/RocPFProb/fine_grain_map.root";
    if(IsExists(infile)){
      cout<<"[SMPAnalyzerCore::SetupL1PrefiringWeight] using file "+infile<<endl;
      TFile f(infile);
      fFGPP=(TH2*)f.Get("fgpp_map");
      fFGPM=(TH2*)f.Get("fgpm_map");
      if(fFGPP) fFGPP->SetDirectory(NULL);
      if(fFGPM) fFGPM->SetDirectory(NULL);
    }else{
      cout<<"[SMPAnalyzerCore::SetupL1PrefiringWeight] no "+infile<<endl;
    }
  }    
  return;
}
void SMPAnalyzerCore::DeleteL1PrefiringWeight(){
  if(fL1Prefiring_photon) delete fL1Prefiring_photon;
  if(fL1Prefiring_jet) delete fL1Prefiring_jet;
  for(int i=0;i<12;i++)
    if(fL1Prefiring_muon[i]) delete fL1Prefiring_muon[i];
  if(rocpfprob) delete rocpfprob;
  if(fFGPP) delete fFGPP;
  if(fFGPM) delete fFGPM;
  return;
}
double SMPAnalyzerCore::getPrefiringRateEcal(double eta, double pt, TH2* h_prefmap, int sys, int mode) const {
  double prefiringRateSystUncEcal_=0.2;
  //Check pt is not above map overflow
  int nbinsy = h_prefmap->GetNbinsY();
  double maxy = h_prefmap->GetYaxis()->GetBinLowEdge(nbinsy + 1);
  if (pt >= maxy)
    pt = maxy - 0.01;
  int thebin = h_prefmap->FindBin(eta, pt);

  double prefrate = h_prefmap->GetBinContent(thebin);
  double abseta=fabs(eta);
  if(mode==1){
    prefrate=h_prefmap->Interpolate(eta,pt);
  }else if(mode==2&&2.0<abseta&&abseta<2.5){
    int sign=eta>0?1:-1;
    double y1=h_prefmap->GetBinContent(h_prefmap->FindBin(sign*2.1, pt));
    double y2=h_prefmap->GetBinContent(h_prefmap->FindBin(sign*2.4, pt));
    double e12=h_prefmap->GetBinError(h_prefmap->FindBin(sign*2.1, pt)); e12=e12*e12;
    double e22=h_prefmap->GetBinError(h_prefmap->FindBin(sign*2.4, pt)); e22=e22*e22;
    double x12=0.015625;
    double x22=0.140625;
      
    double a=(y1*e22*x12+y2*e12*x22)/(x12*x12*e22+x22*x22*e12);
    prefrate=a*pow(abseta-2,2);
    //cout<<" "<<eta<<" "<<pt<<" "<<h_prefmap->GetBinContent(thebin)<<" "<<h_prefmap->Interpolate(eta,pt)<<" "<<prefrate<<endl;
  }
    
  double statuncty = h_prefmap->GetBinError(thebin);
  double systuncty = prefiringRateSystUncEcal_ * prefrate;

  if (sys == 1)
    prefrate = std::min(1., prefrate + sqrt(pow(statuncty, 2) + pow(systuncty, 2)));
  else if (sys == -1)
    prefrate = std::max(0., prefrate - sqrt(pow(statuncty, 2) + pow(systuncty, 2)));
  if (prefrate > 1.) {
    //edm::LogWarning("L1PrefireWeightProducer") << "Found a prefiring probability > 1. Setting to 1." << std::endl;
    return 1.;
  }
  return prefrate;
}
double SMPAnalyzerCore::getPrefiringRatePhoton(double eta, double pt, int sys, int mode) const {
  if(mode>2){
    int iera=0;
    if(DataEra=="2016preVFP") iera=1;
    else if(DataEra=="2016postVPF") iera=2;
    else if(DataEra=="2017") iera=3;
    if(mode==3) //RocPFProb
      return rocpfprob->getPrefireProb(iera,1,eta,pt);
    else if(mode==4){ //RocPFProb linear interpolation
      double abseta=fabs(eta);
      int sign=eta>0?1:-1;
      double x0,y0,x1,y1;
      if(abseta<2) return 0;
      else if(abseta<2.125){
	x0=2; y0=0;
	x1=2.125; y1=rocpfprob->getPrefireProb(iera,1,sign*x1,pt);
      }else if(abseta<2.375){
	x0=2.125; y0=rocpfprob->getPrefireProb(iera,1,sign*x0,pt);
	x1=2.375; y1=rocpfprob->getPrefireProb(iera,1,sign*x1,pt);
      }else if(abseta<2.625){
	x0=2.375; y0=rocpfprob->getPrefireProb(iera,1,sign*x0,pt);
	x1=2.625; y1=rocpfprob->getPrefireProb(iera,1,sign*x1,pt);
      }else{
	x0=2.625; y0=rocpfprob->getPrefireProb(iera,1,sign*x0,pt);
	x1=2.875; y1=rocpfprob->getPrefireProb(iera,1,sign*x1,pt);
      }
      return (y1-y0)/(x1-x0)*(abseta-x0)+y0;
    }else if(mode==5){ //RocPFProb linear interpolation only |eta|<2.5
      double abseta=fabs(eta);
      if(abseta>2.5)
	return rocpfprob->getPrefireProb(iera,1,eta,pt);
      int sign=eta>0?1:-1;
      double x0,y0,x1,y1;
      if(abseta<2) return 0;
      else if(abseta<2.125){
	x0=2; y0=0;
	x1=2.125; y1=rocpfprob->getPrefireProb(iera,1,sign*x1,pt);
      }else{
	x0=2.125; y0=rocpfprob->getPrefireProb(iera,1,sign*x0,pt);
	x1=2.375; y1=rocpfprob->getPrefireProb(iera,1,sign*x1,pt);
      }
      return (y1-y0)/(x1-x0)*(abseta-x0)+y0;
    }else if(mode==6){
      double abseta=fabs(eta);
      if(abseta<2) return 0;
      TH2* fg_map=NULL;
      if(eta>0) fg_map=fFGPP;
      else fg_map=fFGPM;
      int ibin=fg_map->FindBin(eta,pt);
      return fg_map->GetBinContent(ibin);
    }
  }
  return getPrefiringRateEcal(eta, pt, fL1Prefiring_photon, sys, mode);
}
double SMPAnalyzerCore::getPrefiringRateJet(double eta, double pt, int sys, int mode) const {
  if(mode>2){
    int iera=0;
    if(DataEra=="2016preVFP") iera=1;
    else if(DataEra=="2016postVPF") iera=2;
    else if(DataEra=="2017") iera=3;
    return rocpfprob->getPrefireProb(iera,2,eta,pt);
  }    
  return getPrefiringRateEcal(eta, pt, fL1Prefiring_jet, sys, mode);
}
double SMPAnalyzerCore::getPrefiringRateMuon(double eta, double phi, double pt, int sys) const {
  double prefiringRateSystUncMuon_=0.2;
  double prefrate;
  double statuncty;
  if ((DataYear==2016) && (eta > 1.24 && eta < 1.6) &&
      (phi > 2.44346 && phi < 2.79253)) {
    prefrate = fL1Prefiring_muon[11]->Eval(pt);
    statuncty = fL1Prefiring_muon[11]->GetParError(2);
  } else if (std::abs(eta) < 0.2) {
    prefrate = fL1Prefiring_muon[0]->Eval(pt);
    statuncty = fL1Prefiring_muon[0]->GetParError(2);
  } else if (std::abs(eta) < 0.3) {
    prefrate = fL1Prefiring_muon[1]->Eval(pt);
    statuncty = fL1Prefiring_muon[1]->GetParError(2);
  } else if (std::abs(eta) < 0.55) {
    prefrate = fL1Prefiring_muon[2]->Eval(pt);
    statuncty = fL1Prefiring_muon[2]->GetParError(2);
  } else if (std::abs(eta) < 0.83) {
    prefrate = fL1Prefiring_muon[3]->Eval(pt);
    statuncty = fL1Prefiring_muon[3]->GetParError(2);
  } else if (std::abs(eta) < 1.24) {
    prefrate = fL1Prefiring_muon[4]->Eval(pt);
    statuncty = fL1Prefiring_muon[4]->GetParError(2);
  } else if (std::abs(eta) < 1.4) {
    prefrate = fL1Prefiring_muon[5]->Eval(pt);
    statuncty = fL1Prefiring_muon[5]->GetParError(2);
  } else if (std::abs(eta) < 1.6) {
    prefrate = fL1Prefiring_muon[6]->Eval(pt);
    statuncty = fL1Prefiring_muon[6]->GetParError(2);
  } else if (std::abs(eta) < 1.8) {
    prefrate = fL1Prefiring_muon[7]->Eval(pt);
    statuncty = fL1Prefiring_muon[7]->GetParError(2);
  } else if (std::abs(eta) < 2.1) {
    prefrate = fL1Prefiring_muon[8]->Eval(pt);
    statuncty = fL1Prefiring_muon[8]->GetParError(2);
  } else if (std::abs(eta) < 2.25) {
    prefrate = fL1Prefiring_muon[9]->Eval(pt);
    statuncty = fL1Prefiring_muon[9]->GetParError(2);
  } else if (std::abs(eta) < 2.4) {
    prefrate = fL1Prefiring_muon[10]->Eval(pt);
    statuncty = fL1Prefiring_muon[10]->GetParError(2);
  } else {
    //LogDebug("L1PrefireWeightProducer") << "Muon outside of |eta| <= 2.4. Prefiring weight set to 0." << std::endl;
    return 0.;
  }
  double systuncty = prefiringRateSystUncMuon_ * prefrate;

  if (sys == 1)
    prefrate = std::min(1., prefrate + sqrt(pow(statuncty, 2) + pow(systuncty, 2)));
  else if (sys == -1)
    prefrate = std::max(0., prefrate - sqrt(pow(statuncty, 2) + pow(systuncty, 2)));
  //else if (fluctuation == upSyst)
  //  prefrate = std::min(1., prefrate + systuncty);
  //else if (fluctuation == downSyst)
  //  prefrate = std::max(0., prefrate - systuncty);
  //else if (fluctuation == upStat)
  //  prefrate = std::min(1., prefrate + statuncty);
  //else if (fluctuation == downStat)
  //  prefrate = std::max(0., prefrate - statuncty);

  if (prefrate > 1.) {
    //edm::LogWarning("L1PrefireWeightProducer") << "Found a prefiring probability > 1. Setting to 1." << std::endl;
    return 1.;
  }
  return prefrate;
}
double SMPAnalyzerCore::GetL1PrefiringWeight(int mode) const {
  double jetMaxMuonFraction_=0.5;

  //Photons
  vector<Photon> thePhotons = GetAllPhotons();

  //Jets
  vector<Jet> theJets = GetAllJets();

  //Muons
  vector<Muon> theMuons = GetAllMuons();

  //Probability for the event NOT to prefire, computed with the prefiring maps per object.
  //Up and down values correspond to the resulting value when shifting up/down all prefiring rates in prefiring maps.
  double nonPrefiringProba[3] = {1., 1., 1.};      //0: central, 1: up, 2: down
  double nonPrefiringProbaECAL[3] = {1., 1., 1.};  //0: central, 1: up, 2: down
  double nonPrefiringProbaMuon[7] = {
    1., 1., 1., 1., 1., 1., 1.};  //0: central, 1: up, 2: down, 3: up stat, 4: down stat, 5: up syst, 6: down syst

  for (const auto sys : {0, +1, -1}) {
    if(DataYear<2018){
      for (const auto& photon : thePhotons) {
	double pt_gam = photon.UncorrPt();
	double eta_gam = photon.Eta();
	if (pt_gam < 20.)
	  continue;
	if (fabs(eta_gam) < 2.)
	  continue;
	if (fabs(eta_gam) > 3.)
	  continue;
	double prefiringprob_gam = getPrefiringRatePhoton(eta_gam, pt_gam, sys, mode);
	nonPrefiringProbaECAL[sys] *= (1. - prefiringprob_gam);
      }
      
      //Now applying the prefiring maps to jets in the affected regions.
      for (const auto& jet : theJets) {
	double pt_jet = jet.Pt()/jet.JER();
	double eta_jet = jet.Eta();
	//double phi_jet = jet.Phi();
	if (pt_jet < 20.)
	  continue;
	if (fabs(eta_jet) < 2.)
	  continue;
	if (fabs(eta_jet) > 3.)
	  continue;
	if (jetMaxMuonFraction_ > 0 && jet.muonEnergyFraction() > jetMaxMuonFraction_)
	  continue;
	//Loop over photons to remove overlap
	double nonprefiringprobfromoverlappingphotons = 1.;
	bool foundOverlappingPhotons = false;
	for (const auto& photon : thePhotons) {
	  double pt_gam = photon.UncorrPt();
	  double eta_gam = photon.Eta();
	  //double phi_gam = photon.Phi();
	  if (pt_gam < 20.)
	  continue;
	  if (fabs(eta_gam) < 2.)
	    continue;
	  if (fabs(eta_gam) > 3.)
	    continue;
	  double dR = jet.DeltaR(photon);
	  if (dR > 0.4)
	    continue;
	  double prefiringprob_gam = getPrefiringRatePhoton(eta_gam, pt_gam, sys, mode);
	  nonprefiringprobfromoverlappingphotons *= (1. - prefiringprob_gam);
	  foundOverlappingPhotons = true;
	}
	//useEMpt =true if one wants to use maps parametrized vs Jet EM pt instead of pt.
	//if (useEMpt_) pt_jet *= (jet.neutralEmEnergyFraction() + jet.chargedEmEnergyFraction());
	double nonprefiringprobfromoverlappingjet = 1. - getPrefiringRateJet(eta_jet, pt_jet, sys);
	
	if (!foundOverlappingPhotons) {
	  nonPrefiringProbaECAL[sys] *= nonprefiringprobfromoverlappingjet;
	}
	//If overlapping photons have a non prefiring rate larger than the jet, then replace these weights by the jet one
	else if (nonprefiringprobfromoverlappingphotons > nonprefiringprobfromoverlappingjet) {
	  if (nonprefiringprobfromoverlappingphotons > 0.) {
            nonPrefiringProbaECAL[sys] *= nonprefiringprobfromoverlappingjet / nonprefiringprobfromoverlappingphotons;
	  } else {
	    nonPrefiringProbaECAL[sys] = 0.;
	  }
	}
	//Last case: if overlapping photons have a non prefiring rate smaller than the jet, don't consider the jet in the event weight, and do nothing.
      }
    }
    //Now calculate prefiring weights for muons
    for (const auto& muon : theMuons) {
      double pt = muon.MiniAODPt();
      double phi = muon.Phi();
      double eta = muon.Eta();
      // Remove crappy tracker muons which would not have prefired the L1 trigger
      if (pt < 5 || !muon.isPOGLoose())
	continue;
      double prefiringprob_mu = getPrefiringRateMuon(eta, phi, pt, sys);
      nonPrefiringProbaMuon[sys] *= (1. - prefiringprob_mu);
    }
  }
  // Calculate combined weight as product of the weight for individual objects
  for (const auto sys : {0, +1, -1}) {
    nonPrefiringProba[sys] = nonPrefiringProbaECAL[sys] * nonPrefiringProbaMuon[sys];
  }
  // Calculate statistical and systematic uncertainty separately in the muon case
  /*
  for (const auto fluct :
    {fluctuations::upSyst, fluctuations::downSyst, fluctuations::upStat, fluctuations::downStat}) {
    if (!missingInputMuon_ && doMuons_) {
      for (const auto& muon : theMuons) {
        double pt = muon.pt();
        double phi = muon.phi();
        double eta = muon.eta();
        // Remove crappy tracker muons which would not have prefired the L1 trigger
        if (pt < 5 || !muon.isLooseMuon())
          continue;
        double prefiringprob_mu = getPrefiringRateMuon(eta, phi, pt, fluct);
        nonPrefiringProbaMuon[fluct] *= (1. - prefiringprob_mu);
      }
    }
  }
  */
  //if(fabs(L1PrefireReweight_Central-nonPrefiringProba[0])<1e-6) cout<<"OK"<<endl;
  //else cout<<L1PrefireReweight_Central<<" "<<nonPrefiringProba[0]<<" "<<nonPrefiringProbaECAL[0]<<" "<<nonPrefiringProbaMuon[0]<<endl;
  return nonPrefiringProba[0];
}

bool SMPAnalyzerCore::PassSLT1(const Lepton* lep) const{
  if(!lep) return false;
  if(lep->LeptonFlavour()==Lepton::MUON){
    const Muon* muon=(Muon*)lep;
    if(DataYear==2016)
      return muon->PassFilterOR({"hltL3fL1sMu22L1f0Tkf24QL3trkIsoFiltered0p09","hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09"});
    else if(DataYear==2017||DataYear==2018)
      return muon->PassFilter("hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07");
    else{
      cout<<"[SMPAnalyzerCore::PassSLT1] unknown DataYear '"<<DataYear<<"'"<<endl;
      exit(EXIT_FAILURE);
    } 
  }else if(lep->LeptonFlavour()==Lepton::ELECTRON){
    cout<<"[SMPAnalyzerCore::PassSLT1] not yet implemented"<<endl;
    exit(EXIT_FAILURE);
  }else{
    cout<<"[SMPAnalyzerCore::PassSLT1] invalid lepton flavour"<<endl;
    exit(EXIT_FAILURE);
  }
  return false;
}
bool SMPAnalyzerCore::PassSLT2(const Lepton* lep) const{
  if(!lep) return false;
  if(lep->LeptonFlavour()==Lepton::MUON){
    const Muon* muon=(Muon*)lep;
    if(DataYear==2016)
      return false;
    else if(DataYear==2017)
      return muon->PassFilter("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07");
    else if(DataYear==2018)
      return false;
    else{
      cout<<"[SMPAnalyzerCore::PassSLT2] unknown DataYear '"<<DataYear<<"'"<<endl;
      exit(EXIT_FAILURE);
    } 
  }else if(lep->LeptonFlavour()==Lepton::ELECTRON){
    cout<<"[SMPAnalyzerCore::PassSLT2] not yet implemented"<<endl;
    exit(EXIT_FAILURE);
  }else{
    cout<<"[SMPAnalyzerCore::PassSLT2] invalid lepton flavour"<<endl;
    exit(EXIT_FAILURE);
  }
  return false;
}
bool SMPAnalyzerCore::PassDLT1(const Lepton* lep) const{
  if(!lep) return false;
  if(lep->LeptonFlavour()==Lepton::MUON){
    const Muon* muon=(Muon*)lep;
    if(DataYear==2016)
      return (muon->PassFilter("hltL3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17")&&muon->PassFilter("hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4"))||
	(muon->PassFilter("hltL3fL1sDoubleMu114L1f0L2f10L3Filtered17")&&muon->PassFilter("hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4"))||
	(muon->PassFilter("hltL3fL1sDoubleMu114TkFiltered17Q")&&muon->PassFilter("hltDiMuonTrk17Trk8RelTrkIsoFiltered0p4"));
    else if(DataYear==2017||DataYear==2018)
      return (muon->PassFilter("hltL3fL1DoubleMu155fFiltered17")&&muon->PassFilter("hltDiMuon178RelTrkIsoFiltered0p4"));
    else{
      cout<<"[SMPAnalyzerCore::PassDLT1] unknown DataYear '"<<DataYear<<"'"<<endl;
      exit(EXIT_FAILURE);
    } 
  }else if(lep->LeptonFlavour()==Lepton::ELECTRON){
    const Electron* electron=(Electron*)lep;
    return electron->PassFilter("hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter");
  }else{
    cout<<"[SMPAnalyzerCore::PassDLT1] invalid lepton flavour"<<endl;
    exit(EXIT_FAILURE);
  }
  return false;
}
bool SMPAnalyzerCore::PassDLT2(const Lepton* lep) const{
  if(!lep) return false;
  if(lep->LeptonFlavour()==Lepton::MUON){
    const Muon* muon=(Muon*)lep;
    if(DataYear==2016)
      return (muon->PassFilter("hltL3pfL1sDoubleMu114ORDoubleMu125L1f0L2pf0L3PreFiltered8")&&muon->PassFilter("hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4"))||
	(muon->PassFilter("hltDiMuonGlbFiltered17TrkFiltered8")&&muon->PassFilter("hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4"))||
	(muon->PassFilter("hltDiTkMuonTkFiltered17TkFiltered8")&&muon->PassFilter("hltDiMuonTrk17Trk8RelTrkIsoFiltered0p4"));
    else if(DataYear==2017||DataYear==2018)
      return (muon->PassFilter("hltL3fL1DoubleMu155fPreFiltered8")&&muon->PassFilter("hltDiMuon178RelTrkIsoFiltered0p4"));
    else{
      cout<<"[SMPAnalyzerCore::PassDLT2] unknown DataYear '"<<DataYear<<"'"<<endl;
      exit(EXIT_FAILURE);
    } 
  }else if(lep->LeptonFlavour()==Lepton::ELECTRON){
    const Electron* electron=(Electron*)lep;
    return electron->PassFilter("hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter");
  }else{
    cout<<"[SMPAnalyzerCore::PassDLT2] invalid lepton flavour"<<endl;
    exit(EXIT_FAILURE);
  }
  return false;
}
vector<vector<Weight>> SMPAnalyzerCore::Make2DWeights(const vector<int>& structure){
  vector<vector<Weight>> rt;
  for(const int& nmem:structure){
    rt.push_back(vector<Weight>(nmem,1.0));
  }
  return rt;
}
double SMPAnalyzerCore::GetBchargeSF(const Jet& bjet,int set,int mem) const {
  double sf=1.;
  double charge=bjet.GetUserFloat("AFBCharge");
  int flavour=bjet.GenHFHadronMatcherFlavour();
  int origin=bjet.GenHFHadronMatcherOrigin();
  if(flavour!=5) return sf;
  if(origin==-999) return sf;
  double data_m_accuracy=0.615651;
  double data_p_accuracy=0.629591;
  //data (BB SS method (mm)): 0.619894+-0.004496
  //data - (BB SS method (me)): 0.615651+-0.003372
  //data + (BB SS method (me)): 0.629591+-0.003385
  double sim_m_accuracy=0.632788;
  double sim_p_accuracy=0.645423;
  //sim (BB SS method): 0.639483+-0.000594
  //sim - (BB SS method (me)): 0.632788+-0.000486
  //sim + (BB SS method (me)): 0.645423+-0.000485
  if(set==0){
    if(mem==0){
      data_m_accuracy+=-0.001086;
      data_p_accuracy+=0.001081;
    }else if(mem==1){
      data_m_accuracy+=-0.003192;
      data_p_accuracy+=-0.003208;
    }
  }
  if(origin*charge<0){
    if(origin>0){
      sf=data_m_accuracy/sim_m_accuracy;
    }else{
      sf=data_p_accuracy/sim_p_accuracy;
    }      
  }else{
    if(origin>0){
      sf=(1-data_m_accuracy)/(1-sim_m_accuracy);
    }else{
      sf=(1-data_p_accuracy)/(1-sim_p_accuracy);
    }      
  }  
  return sf;
}
double SMPAnalyzerCore::GetBchargeSF(const Parameter& p,int set,int mem) const {
  double sf=1.;
  if(!p.bjets.size()) return sf;
  for(int i=0,n=p.bjets.size();i<n;i++){
    if(p.bjets.at(i).Pt()<p.c.jetpt) return sf;
    sf*=GetBchargeSF(p.bjets.at(i),set,mem);
  }
  return sf;
}
