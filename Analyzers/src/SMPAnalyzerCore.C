#include "SMPAnalyzerCore.h"

SMPAnalyzerCore::SMPAnalyzerCore(){}
SMPAnalyzerCore::~SMPAnalyzerCore(){
  if(roc) delete roc;
  if(rocele) delete rocele;
  for(std::map< TString, TH4D* >::iterator mapit = maphist_TH4D.begin(); mapit!=maphist_TH4D.end(); mapit++){
    delete mapit->second;
  }
  maphist_TH4D.clear();
  DeleteEfficiency();
  DeleteZptWeight();
  DeleteCFRate();
  DeleteFakeRate();
}

void SMPAnalyzerCore::initializeAnalyzer(){
  if(MaxEvent>0) reductionweight=1.*fChain->GetEntries()/MaxEvent;
  else reductionweight=1.;
  SetupEfficiency();
  SetupRoccoR();
  SetupPUJetWeight();
  SetupCFRate();
  SetupFakeRate();
  IsDYSample=false;
  if(MCSample.Contains("DYJets")||MCSample.Contains("ZToEE")||MCSample.Contains("ZToMuMu")||MCSample.Contains(TRegexp("DY[0-9]Jets"))) IsDYSample=true;
  if(IsDYSample) SetupZptWeight();

  vector<JetTagging::Parameters> jtps={JetTagging::Parameters(JetTagging::DeepJet,JetTagging::Tight,JetTagging::incl,JetTagging::comb),
                                       JetTagging::Parameters(JetTagging::DeepJet,JetTagging::Loose,JetTagging::incl,JetTagging::comb),
                                       JetTagging::Parameters(JetTagging::DeepCSV,JetTagging::Tight,JetTagging::incl,JetTagging::comb),
                                       JetTagging::Parameters(JetTagging::DeepCSV,JetTagging::Loose,JetTagging::incl,JetTagging::comb),

                                       JetTagging::Parameters(JetTagging::DeepJet_CvsB,JetTagging::Tight,JetTagging::incl,JetTagging::wcharm),
                                       JetTagging::Parameters(JetTagging::DeepJet_CvsB,JetTagging::Loose,JetTagging::incl,JetTagging::wcharm),
                                       JetTagging::Parameters(JetTagging::DeepJet_CvsL,JetTagging::Tight,JetTagging::incl,JetTagging::wcharm),
                                       JetTagging::Parameters(JetTagging::DeepJet_CvsL,JetTagging::Loose,JetTagging::incl,JetTagging::wcharm),
                                       JetTagging::Parameters(JetTagging::DeepCSV_CvsB,JetTagging::Tight,JetTagging::incl,JetTagging::wcharm),
                                       JetTagging::Parameters(JetTagging::DeepCSV_CvsB,JetTagging::Loose,JetTagging::incl,JetTagging::wcharm),
                                       JetTagging::Parameters(JetTagging::DeepCSV_CvsL,JetTagging::Tight,JetTagging::incl,JetTagging::wcharm),
                                       JetTagging::Parameters(JetTagging::DeepCSV_CvsL,JetTagging::Loose,JetTagging::incl,JetTagging::wcharm)};
  mcCorr->SetJetTaggingParameters(jtps);
}
void SMPAnalyzerCore::beginEvent(){
  _event=GetEvent();
  allmus=GetAllMuons();
  std::sort(allmus.begin(),allmus.end(),PtComparing);
  allels=GetAllElectrons();
  std::sort(allels.begin(),allels.end(),PtComparing);
  alljets=GetJets("tightLepVeto",20,2.4); // These part shows ERROR only in data when EfficiencyValidation
  std::sort(alljets.begin(),alljets.end(),PtComparing);

  if(!IsDATA){
    lhes=GetLHEs();
    gens=GetGens();
    //PrintLHEs(lhes);
    //PrintGens(gens);
    if(IsDYSample||MCSample.Contains("GamGamToLL")){
      GetDYLHEParticles(lhes,lhe_p0,lhe_p1,lhe_l0,lhe_l1,lhe_j0);
      GetDYGenParticles(gens,gen_p0,gen_p1,gen_l0,gen_l1,gen_j0,3);
      //GetDYGenParticles(gens,gen_p0,gen_p1,gen_l0_dressed,gen_l1_dressed,1);
      //GetDYGenParticles(gens,gen_p0,gen_p1,gen_l0_bare,gen_l1_bare,0);
    }
  }
}
void SMPAnalyzerCore::executeEventWithParameter(Parameter p){
  p.SetLeptons();
  if(p.channel.Length() > 2) p.SetJets();
  if(p.weightbit&NominalWeight) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"lumi",p.w.lumiweight);
  
  if(!_event.PassTrigger(p.triggers)) return;
  if(p.weightbit&NominalWeight) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"passTrig",p.w.lumiweight);

  if(p.weightbit&NominalWeight){
    FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"PU",p.w.lumiweight*p.w.PUweight);
    FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"prefire",p.w.lumiweight*p.w.PUweight*p.w.prefireweight);
    FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"zpt",p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight);
    FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"z0",p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight);
    FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"weak",p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight);
  }

  double eventweight=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.zptweight*p.w.weakweight;
  if(p.weightbit&NominalWeight) FillHist(p.prefix+p.hprefix+"nlepton"+p.suffix,p.muons.size()+p.electrons.size(),eventweight,10,0,10);

  /////////////////////// lepton selection ///////////////////////
  if(!SMPAnalyzerCore::PassSelection(p)) return;
  ///////////////// efficiency scale factors ///////////////////
  EvalIDSF(p);
  EvalTriggerSF(p);
  if(p.weightbit&NominalWeight){
    FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"RECOSF",eventweight*p.w.electronRECOSF);
    FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"IDSF",eventweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF);
    FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"ISOSF",eventweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF);
    FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"triggerSF",eventweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF);
    FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"CFSF",eventweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF);
  }

  /////////////////////// jet selection ///////////////////////
  if(p.channel.Length() > 2){ if(!PassSelection(p)) return;}
  ////// Fill histograms //////////
  FillHists(p);
}
void SMPAnalyzerCore::EvalIDSF(Parameter& p){
  p.doublemap["muontrackingSF"]=1.;
  if(!IsDATA){
    if(p.weightbit&EfficiencyWeight){
      p.w.electronRECOSF_sys=fEff->GetStructure(p.k.electronRECOSF);
      p.w.electronIDSF_sys=fEff->GetStructure(p.k.electronIDSF);
      p.w.muonIDSF_sys=fEff->GetStructure(p.k.muonIDSF);
      p.w.muonISOSF_sys=fEff->GetStructure(p.k.muonISOSF);
    }
    for(const auto& electron:p.electrons){
      p.w.electronRECOSF*=fEff->GetEfficiencySF(p.k.electronRECOSF,&electron,0,0);
      if(p.weightbit&EfficiencyWeight){
        int nset=p.w.electronRECOSF_sys.size();
        for(int s=0;s<nset;s++){
          int nmem=p.w.electronRECOSF_sys[s].size();
          for(int m=0;m<nmem;m++){
            p.w.electronRECOSF_sys[s][m]*=fEff->GetEfficiencySF(p.k.electronRECOSF,&electron,s,m);
          }
        }
      }
      p.w.electronIDSF*=fEff->GetEfficiencySF(p.k.electronIDSF,&electron,0,0);
      p.w.electronIDSF*=fEff->GetEfficiencySF(p.k.electronIDSF2,&electron,0,0);
      if(p.weightbit&EfficiencyWeight){
        int nset=p.w.electronIDSF_sys.size();
        for(int s=0;s<nset;s++){
          int nmem=p.w.electronIDSF_sys[s].size();
          for(int m=0;m<nmem;m++){
            p.w.electronIDSF_sys[s][m]*=fEff->GetEfficiencySF(p.k.electronIDSF,&electron,s,m);
            p.w.electronIDSF_sys[s][m]*=fEff->GetEfficiencySF(p.k.electronIDSF2,&electron,s,m);
          }
        }
      }
    }
    for(const auto& muon:p.muons){
      p.doublemap["muontrackingSF"]*=GetMuonTrackingSF(muon.Eta());
      p.w.muonIDSF*=fEff->GetEfficiencySF(p.k.muonIDSF,&muon,0,0);
      if(p.weightbit&EfficiencyWeight){
        int nset=p.w.muonIDSF_sys.size();
        for(int s=0;s<nset;s++){
          int nmem=p.w.muonIDSF_sys[s].size();
          for(int m=0;m<nmem;m++){
            p.w.muonIDSF_sys[s][m]*=fEff->GetEfficiencySF(p.k.muonIDSF,&muon,s,m);
          }
        }
      }
      p.w.muonISOSF*=fEff->GetEfficiencySF(p.k.muonISOSF,&muon,0,0);
      if(p.weightbit&EfficiencyWeight){
        int nset=p.w.muonISOSF_sys.size();
        for(int s=0;s<nset;s++){
          int nmem=p.w.muonISOSF_sys[s].size();
          for(int m=0;m<nmem;m++){
            p.w.muonISOSF_sys[s][m]*=fEff->GetEfficiencySF(p.k.muonISOSF,&muon,s,m);
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
    if(p.k.triggerSF.size()){
      if(p.weightbit&EfficiencyWeight){
        p.w.triggerSF_sys=fEff->GetStructure(p.k.triggerSF[0]);
      }
      if(p.k.triggerSF.at(0).Contains("Mu")) triggerables=MakeLeptonPointerVector(p.muons);
      else if(p.k.triggerSF.at(0).Contains("Ele")) triggerables=MakeLeptonPointerVector(p.electrons);
    }
    if(p.k.triggerSF.size()==1){
      p.w.triggerSF*=GetLeptonTriggerSF(p.k.triggerSF[0],triggerables,0,0);
      if(p.weightbit&EfficiencyWeight){
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
        p.w.triggerSF*=GetLeptonTriggerORSF(p.k.triggerSF[0],p.k.triggerSF[1],triggerables,0,0);
        if(p.weightbit&EfficiencyWeight){
          int nset=p.w.triggerSF_sys.size();
          for(int s=0;s<nset;s++){
            int nmem=p.w.triggerSF_sys[s].size();
            for(int m=0;m<nmem;m++){
              p.w.triggerSF_sys[s][m]*=GetLeptonTriggerORSF(p.k.triggerSF[0],p.k.triggerSF[1],triggerables,s,m);
            }
          }
        }
      }else{
        p.w.triggerSF*=GetDileptonTriggerSF(p.k.triggerSF[0],p.k.triggerSF[1],triggerables,0,0);
        if(p.weightbit&EfficiencyWeight){
          int nset=p.w.triggerSF_sys.size();
          for(int s=0;s<nset;s++){
            int nmem=p.w.triggerSF_sys[s].size();
            for(int m=0;m<nmem;m++){
              p.w.triggerSF_sys[s][m]*=GetDileptonTriggerSF(p.k.triggerSF[0],p.k.triggerSF[1],triggerables,s,m);
            }
          }
        }
      }
    }
  }
}
bool SMPAnalyzerCore::PassSelection(Parameter& p){
  double weight=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.zptweight*p.w.weakweight;

  if(!PassMETFilter()) return false;
  if(p.weightbit&NominalWeight) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"METfilter",weight);  

  if(p.c.nelectronmax>=0&&(int)p.electrons.size()>p.c.nelectronmax) return false;
  if(p.c.nmuonmax>=0&&(int)p.muons.size()>p.c.nmuonmax) return false;
    
  if(!p.lepton0||!p.lepton1) return false;
  if(p.weightbit&NominalWeight) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"Dilepton",weight);

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
  if(p.weightbit&NominalWeight) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"LepPtCut",weight);

  if(p.lepton0->Charge()*p.lepton1->Charge()>0) p.hprefix+="ss_";
  if(p.weightbit&NominalWeight) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"charge",weight);

  if((*p.lepton0+*p.lepton1).M()<52) return false;
  if(p.weightbit&NominalWeight) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"Mass52",weight);

  if(p.option.Contains("triggermatching")){
    if(p.triggers.size()){
      if(p.triggers.at(0).Contains(TRegexp("HLT_[Tk]*Mu17_TrkIsoVVL_[Tk]*Mu8_TrkIsoVVL"))){
	if(p.lepton0->LeptonFlavour()!=Lepton::MUON) return false;
	if(p.lepton1->LeptonFlavour()!=Lepton::MUON) return false;
	Muon *muon0=(Muon*)p.lepton0,*muon1=(Muon*)p.lepton1;
	if(GetEra()=="2016preVFP"||GetEra()=="2016postVFP"){
	  if(!muon0->PassFilterOR({"hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4","hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4","hltDiMuonTrk17Trk8RelTrkIsoFiltered0p4"})) return false;
	  if(!muon1->PassFilterOR({"hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4","hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4","hltDiMuonTrk17Trk8RelTrkIsoFiltered0p4"})) return false;
	  if(!muon0->PassFilterOR({"hltL3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17","hltL3fL1sDoubleMu114L1f0L2f10L3Filtered17","hltL3fL1sDoubleMu114TkFiltered17Q"})
	     &&!muon1->PassFilterOR({"hltL3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17","hltL3fL1sDoubleMu114L1f0L2f10L3Filtered17","hltL3fL1sDoubleMu114TkFiltered17Q"})) return false;
	  if(!muon0->PassFilterOR({"hltL3pfL1sDoubleMu114ORDoubleMu125L1f0L2pf0L3PreFiltered8","hltDiMuonGlbFiltered17TrkFiltered8","hltDiTkMuonTkFiltered17TkFiltered8"})
	     ||!muon1->PassFilterOR({"hltL3pfL1sDoubleMu114ORDoubleMu125L1f0L2pf0L3PreFiltered8","hltDiMuonGlbFiltered17TrkFiltered8","hltDiTkMuonTkFiltered17TkFiltered8"})) return false;
	}else if(GetEra()=="2017"||GetEra()=="2018"){
	  if(!muon0->PassFilter("hltDiMuon178RelTrkIsoFiltered0p4")) return false;
	  if(!muon1->PassFilter("hltDiMuon178RelTrkIsoFiltered0p4")) return false;
	  if(!muon0->PassFilter("hltL3fL1DoubleMu155fFiltered17")&&!muon1->PassFilter("hltL3fL1DoubleMu155fFiltered17")) return false;
	  if(!muon0->PassFilter("hltL3fL1DoubleMu155fPreFiltered8")||!muon1->PassFilter("hltL3fL1DoubleMu155fPreFiltered8")) return false;
	}else{
	  cout<<"[SMPAnalyzerCore::PassSelection] trigger matching for "<<GetEra()<<" "<<p.triggers.at(0)<<" is not implemented"<<endl;
	  exit(EXIT_FAILURE);
	}
      }else if(p.triggers.at(0).Contains("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL")){
	if(p.lepton0->LeptonFlavour()!=Lepton::ELECTRON) return false;
	if(p.lepton1->LeptonFlavour()!=Lepton::ELECTRON) return false;
	Electron *electron0=(Electron*)p.lepton0,*electron1=(Electron*)p.lepton1;
	if(!electron0->PassFilter("hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter")&&!electron1->PassFilter("hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter")) return false;
	if(!electron0->PassFilter("hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter")||!electron1->PassFilter("hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter")) return false;
      }else{
	cout<<"[SMPAnalyzerCore::PassSelection] trigger matching for "<<p.triggers.at(0)<<" is not implemented"<<endl;
	exit(EXIT_FAILURE);
      }
    }else{
      cout<<"[SMPAnalyzerCore::PassSelection] no trigger for trigger matching"<<endl;
      exit(EXIT_FAILURE);
    }
    if(p.weightbit&NominalWeight) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"TriggerMatching",weight);
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
    if(leptons.at(0)!=p.lepton0) return false;
    if(leptons.at(1)!=p.lepton1) return false;
    if(p.weightbit&NominalWeight) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"StrictOrder",weight);
  }
  return true;
}

TH4D* SMPAnalyzerCore::GetHist4D(TString histname){
  TH4D *h = NULL;
  std::map<TString, TH4D*>::iterator mapit = maphist_TH4D.find(histname);
  if(mapit != maphist_TH4D.end()) return mapit->second;
  return h;
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
void SMPAnalyzerCore::FillHists(Parameter& p){
  double weight=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF;
  TLorentzVector dilepton=(*p.lepton0)+(*p.lepton1);
  double dimass=dilepton.M();
  if(dimass>=60&&dimass<120){
    FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"m60to120",weight);
    FillHist(p.prefix+"m60to120/"+p.hprefix+"dimass"+p.suffix,dimass,weight,60,60,120);
    if(dimass>=80&&dimass<100){
      FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"m80to100",weight);
      FillHist(p.prefix+"m80to100/"+p.hprefix+"dimass"+p.suffix,dimass,weight,40,80,100);
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
double SMPAnalyzerCore::GetLeptonTriggerSF(TString triggerSF_key,const vector<Lepton*>& leps,int set,int mem){
  if(IsDATA) return 1;
  if(triggerSF_key=="") return 1;
  if(triggerSF_key=="Default") return 1;

  double data_eff=1.,sim_eff=1.;
  for(const auto& lep:leps){
    data_eff*=1-fEff->GetDataEfficiency(triggerSF_key,lep,set,mem);
    sim_eff*=1-fEff->GetSimEfficiency(triggerSF_key,lep,set,mem);
  }
  data_eff=1-data_eff;
  sim_eff=1-sim_eff;
  if(sim_eff==0) return 1.;
  else return data_eff/sim_eff;
}
double SMPAnalyzerCore::GetLeptonTriggerORSF(TString triggerSF_key0,TString triggerSF_key1,const vector<Lepton*>& leps,int set,int mem){
  if(IsDATA) return 1;

  double lumi0=1.; //trigger0 on
  double lumi1=0.; //only trigger1 on
  double lumi2=0.; //both off
  if(DataYear==2017&&triggerSF_key0.Contains("IsoMu24")&&triggerSF_key1.Contains("IsoMu27")){
    lumi0=37997.005;    lumi1=3480.873;    lumi2=0.;
  }else if(DataYear==2017&&triggerSF_key0.Contains("Ele27")&&triggerSF_key1.Contains("Ele32")){
    lumi0=31661.026;    lumi1=9522.208;    lumi2=295.;
  }else if(DataYear==2018&&triggerSF_key0.Contains("Ele28")&&triggerSF_key1.Contains("Ele32")){
    lumi0=23687.253;    lumi1=36140.626;   lumi2=0.;
  }else{
    cout<<"[SMPAnalyzerCore::LeptonTriggerOR_SF] not available combination "<<triggerSF_key0<<"||"<<triggerSF_key1<<" for "<<DataEra<<endl;
    exit(EXIT_FAILURE);
  }

  double data_eff_key0=1.,data_eff_key1=1.,sim_eff=1.;
  for(const auto& lep:leps){
    data_eff_key0*=1-fEff->GetDataEfficiency(triggerSF_key0,lep,set,mem);
    data_eff_key1*=1-fEff->GetDataEfficiency(triggerSF_key1,lep,set,mem);
    sim_eff*=1-fEff->GetSimEfficiency(triggerSF_key0,lep,set,mem);
  }
  data_eff_key0=1-data_eff_key0;
  data_eff_key1=1-data_eff_key1;
  sim_eff=1-sim_eff;
  if(sim_eff==0) return 1.;
  else return (lumi0*data_eff_key0+lumi1*data_eff_key1)/((lumi0+lumi1+lumi2)*sim_eff);
}
double SMPAnalyzerCore::GetDileptonTriggerSF(TString triggerSF_key0,TString triggerSF_key1,const vector<Lepton*>& leps,int set,int mem){
  if(IsDATA) return 1;
  if((triggerSF_key0==""||triggerSF_key0=="Default")&&(triggerSF_key1==""||triggerSF_key1=="Default")) return 1;
  int nlep=leps.size();
  if(nlep<2){
    cout<<"[SMPAnalyzerCore::DileptonTrigger_SF] nlep < 2. return 1."<<endl;
    return 1.;
  }
  double data_noleg1=1.,sim_noleg1=1.;
  vector<double> data_oneleg1_noleg2(nlep,1.);
  vector<double> sim_oneleg1_noleg2(nlep,1.);
  for(int i=0;i<nlep;i++){
    double data_eff_leg1=fEff->GetDataEfficiency(triggerSF_key0,leps.at(i),set,mem);
    double data_eff_leg2=fEff->GetDataEfficiency(triggerSF_key1,leps.at(i),set,mem);
    double sim_eff_leg1=fEff->GetSimEfficiency(triggerSF_key0,leps.at(i),set,mem);
    double sim_eff_leg2=fEff->GetSimEfficiency(triggerSF_key1,leps.at(i),set,mem);
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
  if(sim_eff==0) return 1.;
  else return data_eff/sim_eff;
}

// ZptWeight
void SMPAnalyzerCore::SetupZptWeight(){
  TString _MCSample=MCSample;
  if(MCSample.Contains("MiNNLO")) _MCSample="MiNNLO";
  TString zptpath=(TString)getenv("DATA_DIR")+"/"+GetEra()+"/SMP/ZptWeight_"+_MCSample+".root";
  if(IsExists(zptpath)){
    cout<<"[SMPAnalyzerCore::SetupZptWeight] using file "+zptpath<<endl;
  }else{
    cout<<"[SMPAnalyzerCore::SetupZptWeight] no "+zptpath<<endl;
    return;
  }
  DeleteZptWeight();
  TFile f(zptpath);
  fZptWeightG=(TF1*)f.Get("zptweight_g");
  if(!fZptWeightG){
    cout<<"[SMPAnalyzerCore::SetupZptWeight] no zptweight_g"<<endl;
    exit(ENODATA);
  }
  fZptWeightYaxis=(TAxis*)f.Get("yaxis");
  if(!fZptWeightYaxis){
    cout<<"[SMPAnalyzerCore::SetupZptWeight] no yaxis"<<endl;
    exit(ENODATA);
  }
  fZptWeightY.resize(fZptWeightYaxis->GetNbins()+2,NULL);
  for(int i=1;i<fZptWeightYaxis->GetNbins()+1;i++){
    fZptWeightY[i]=(TF1*)f.Get(Form("zptweight_y%d",i));
    if(!fZptWeightY[i]){
      cout<<"[SMPAnalyzerCore::SetupZptWeight] no zptweight_y"+TString(i)<<endl;
      exit(ENODATA);
    }
  }
  fZptWeightMaxis=(TAxis*)f.Get("maxis");
  if(!fZptWeightMaxis){
    cout<<"[SMPAnalyzerCore::SetupZptWeight] no maxis"<<endl;
    exit(ENODATA);
  }
  fZptWeightM.resize(fZptWeightMaxis->GetNbins()+2,NULL);
  for(int i=1;i<fZptWeightMaxis->GetNbins()+1;i++){
    fZptWeightM[i]=(TF1*)f.Get(Form("zptweight_m%d",i));
    if(!fZptWeightM[i]){
      cout<<"[SMPAnalyzerCore::SetupZptWeight] no zptweight_m"+TString(i)<<endl;
      exit(ENODATA);
    }
  }
}
double SMPAnalyzerCore::GetZptWeight(double mass,double rapidity,double pt,TString opt){
  if(!fZptWeightG) return 1.;
  if(mass==0) return 1.;
  if(isnan(rapidity)) return 1.;
  double m=mass;
  if(m<fZptWeightMaxis->GetXmin()) m=fZptWeightMaxis->GetXmin();
  if(m>=fZptWeightMaxis->GetXmax()) m=fZptWeightMaxis->GetXmax()-1e-6;
  double y=fabs(rapidity);
  if(y>=fZptWeightYaxis->GetXmax()) y=fZptWeightYaxis->GetXmax()-1e-6;
  if(pt<0) pt=0;
  if(pt>=650) pt=649.9;
  double sf=1.;

  opt.ToUpper();
  if(opt.Contains("G")) sf*=fZptWeightG->Eval(pt);

  if(opt.Contains("Y")){
    double ymin=fZptWeightYaxis->GetBinCenter(1);
    double ymax=fZptWeightYaxis->GetBinCenter(fZptWeightYaxis->GetNbins());
    int biny1,biny2;
    if(y<ymin){
      biny1=1;
      biny2=2;
    }else if(y>=ymax){
      biny1=fZptWeightYaxis->GetNbins()-1;
      biny2=fZptWeightYaxis->GetNbins();
    }else{
      int biny=fZptWeightYaxis->FindBin(y);
      if(y>=fZptWeightYaxis->GetBinCenter(biny)){
        biny1=biny;
        biny2=biny+1;
      }else{
        biny1=biny-1;
        biny2=biny;
      }
    }
    double y1=fZptWeightYaxis->GetBinCenter(biny1);
    double y2=fZptWeightYaxis->GetBinCenter(biny2);
    sf*=( (y2-y)*fZptWeightY[biny1]->Eval(pt) + (y-y1)*fZptWeightY[biny2]->Eval(pt) )/(y2-y1);
  }

  if(opt.Contains("M")){
    double mmin=fZptWeightMaxis->GetBinCenter(1);
    double mmax=fZptWeightMaxis->GetBinCenter(fZptWeightMaxis->GetNbins());
    int binm1,binm2;
    if(m<mmin){
      binm1=1;
      binm2=2;
    }else if(m>=mmax){
      binm1=fZptWeightMaxis->GetNbins()-1;
      binm2=fZptWeightMaxis->GetNbins();
    }else{
      int binm=fZptWeightMaxis->FindBin(m);
      if(m>=fZptWeightMaxis->GetBinCenter(binm)){
        binm1=binm;
        binm2=binm+1;
      }else{
        binm1=binm-1;
        binm2=binm;
      }
    }
    double m1=fZptWeightMaxis->GetBinCenter(binm1);
    double m2=fZptWeightMaxis->GetBinCenter(binm2);
    sf*=( (m2-m)*fZptWeightM[binm1]->Eval(pt) + (m-m1)*fZptWeightM[binm2]->Eval(pt) )/(m2-m1);
  }
  return sf;
}
void SMPAnalyzerCore::DeleteZptWeight(){
  if(fZptWeightG){
    delete fZptWeightG;
    fZptWeightG=NULL;
  }
  if(fZptWeightYaxis){
    delete fZptWeightYaxis;
    fZptWeightYaxis=NULL;
  }
  for(auto f:fZptWeightY){
    if(f) delete f;
  }
  fZptWeightY.clear();  
  if(fZptWeightMaxis){
    delete fZptWeightMaxis;
    fZptWeightMaxis=NULL;
  }
  for(auto f:fZptWeightM){
    if(f) delete f;
  }
  fZptWeightM.clear();
}

void SMPAnalyzerCore::SetupRoccoR(){
  cout<<"[SMPAnalyzerCore::SetupRoccoR] setting Rocheseter Correction"<<endl;
  TString datapath=getenv("DATA_DIR");
  TString rocpath=datapath+"/"+GetEra()+"/RoccoR/RoccoR"+GetEraShort()+"UL.txt";
  if(IsExists(rocpath)) roc=new RoccoR(rocpath.Data());
  else cout<<"[SMPAnalyzerCore::SetupRoccoR] no "+rocpath<<endl;
  TString erashort=GetEraShort();
  TString rocelepath=datapath+"/"+GetEra()+"/RoccoR/e_"+erashort(2,3)+"UL.txt";
  if(DataYear==2016) rocelepath=datapath+"/"+GetEra()+"/RoccoR/e_"+erashort(2,3)+"UL_1.txt";
  if(IsExists(rocelepath)) rocele=new Aepcor(rocelepath.Data());
  else cout<<"[SMPAnalyzerCore::SetupRoccoR] no "+rocelepath<<endl;  
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
}

void SMPAnalyzerCore::PrintGens(const vector<Gen>& gens){
  cout<<"index\tpid\tstatus\tmother\tHard\tPrompt\tpt\tpz\teta\tphi\tmass\n";
  for(int i=0;i<(int)gens.size();i++){
    gens[i].Print();
  }
}
void SMPAnalyzerCore::PrintLHEs(vector<LHE>& lhes){
  cout<<"index\tpid\tstatus\tpt\teta\tphi\tmass\n";
  for(int i=0;i<(int)lhes.size();i++){
    lhes[i].Print();
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

// This function used for DY+b analysis (in especially, qG or Gq collisions)
void SMPAnalyzerCore::GetDYLHEParticles(const vector<LHE>& lhes,LHE& p0,LHE& p1,LHE& l0,LHE& l1,LHE& j0){
  if(!IsDYSample&&!MCSample.Contains("GamGamToLL")){
    cout <<"[AFBAnalyzer::GetDYLHEParticles] this is for DY event"<<endl;
    exit(EXIT_FAILURE);
  }
  p0=LHE();
  p1=LHE();
  l0=LHE();
  l1=LHE();
  j0=LHE();

  if(!lhes.size()) return;
  bool IsqG = false;
  if(lhes[0].ID() !=lhes[1].ID() && max(lhes[0].ID(),lhes[1].ID()) == 21) IsqG = true;
  int bnum=0;

  for(int i=0;i<(int)lhes.size();i++){
    //cout<<lhes[i].Index()<<"\t"<<lhes[i].ID()<<"\t"<<lhes[i].Status()<<"\t"<<lhes[i].E()<<"\t"<<lhes[i].Px()<<"\t"<<lhes[i].Py()<<"\t"<<lhes[i].Pz()<<"\t"<<lhes[i].Eta()<<"\t"<<lhes[i].M()<<"\t"<<endl;
    if(p0.ID()==0&&lhes[i].Status()==-1&&lhes[i].Eta()>0) p0=lhes[i];
    if(p1.ID()==0&&lhes[i].Status()==-1&&lhes[i].Eta()<0) p1=lhes[i];
    if(l0.ID()==0&&(abs(lhes[i].ID())==11||abs(lhes[i].ID())==13||abs(lhes[i].ID())==15)) l0=lhes[i];
    if(l0.ID()&&lhes[i].ID()==-l0.ID()) l1=lhes[i];
    if(IsqG){ 
      if(j0.ID()==0 && lhes[i].ID()==min(lhes[0].ID(),lhes[1].ID()) && lhes[i].Status()==1) j0=lhes[i]; //Among status=1 lhes, the first,second are always leptons, the third is quark.
    }else{                                                                                              // (if gluon radiation only, then gluon)
      if(j0.ID()==0 && (abs(lhes[i].ID())<6 || lhes[i].ID()==21) && lhes[i].Status()==1) j0=lhes[i];
    }
    if(lhes[i].ID()==5&&lhes[i].Status()==1) bnum += 3;
    else if(lhes[i].ID()==-5&&lhes[i].Status()==1) bnum -= 2;
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
  /*
  if(bnum==1) p.fprefix += "bB_";
  else if(bnum==4) p.fprefix += "bbB_";
  else if(bnum==-1) p.fprefix += "BbB_";
  else if(bnum==3) p.fprefix += "b_";
  else if(bnum==-2) p.fprefix += "B_";
  else if(bnum==0) p.fprefix += "";
  else p.fprefix = p.fprefix+"b"+Form("%d",bnum)+"_";
  */
}
// This function used for DY+b analysis (in especially, qG or Gq collisions)
void SMPAnalyzerCore::GetDYGenParticles(const vector<Gen>& gens,Gen& parton0,Gen& parton1,Gen& l0,Gen& l1,Gen& j0,int mode){
  //mode 0:bare 1:dressed01 2:dressed04 3:beforeFSR
  if(!IsDYSample&&!MCSample.Contains("GamGamToLL")){
    cout <<"[SMPAnalyzerCore::GetDYGenParticles] this is for DY event"<<endl;
    exit(EXIT_FAILURE);
  }else{
    if(abs(lhe_l0.ID()) == 15) return; // ZtoTauTau is too complex. It sometimes include hadronic decay.
  }

  parton0=Gen();
  parton1=Gen();
  l0=Gen();
  l1=Gen();
  j0=Gen();
  vector<const Gen*> leptons;
  vector<const Gen*> photons;
  vector<const Gen*> jets;

  int ngen=gens.size();
  for(int i=0;i<ngen;i++){
    //gens.at(i).Print();
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
    if(gens.at(i).isHardProcess()&&(abs(genpid)<7||genpid==21)) jets.push_back(&gens[i]);
  }
  int nlepton=leptons.size();
  for(int i=0;i<nlepton;i++){
    for(int j=i+1;j<nlepton;j++){
      if(!(leptons[i]->PID()+leptons[j]->PID()==0)) continue;
      if((*leptons[i]+*leptons[j]).M()>(l0+l1).M()){
        if(leptons[i]->Pt()>leptons[j]->Pt()){
          l0=*leptons[i];
          l1=*leptons[j];
        }else{
          l0=*leptons[j];
          l1=*leptons[i];
        }
      }
    }
  }
  if(l0.PID()==0||l1.PID()==0){
    cout << "[AFBAnalyzer::GetGenParticles] something is wrong"<<endl;
    for(int i=0;i<ngen;i++){
      gens.at(i).Print();
      cout << "gen "<<i<<" : prompt? : "<<gens.at(i).isPrompt()<<endl;
    }
    cout << "l0 index, l1 index, l0l1mass : "<<l0.Index()<<", "<<l1.Index()<<", "<<(l0+l1).M()<<endl;
    exit(EXIT_FAILURE);
  }

  bool IsqG = false;
  if(parton0.PID() != parton1.PID() && max(parton0.PID(),parton1.PID()) == 21) IsqG = true;

  int njet=jets.size();
  for(int i=0;i<njet;i++){
    if(IsqG){
      if(jets[i]->PID()!=min(parton0.PID(),parton1.PID())) continue;
    }
    if((jets[i]->Pt()>j0.Pt())) j0=*jets[i];
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
  std::vector<Electron> electrons = allels;
  if(id=="passMediumID_SelQ"){
    for(unsigned int i=0; i<electrons.size(); i++){
      Electron this_electron= electrons.at(i);
      if(!( this_electron.Pt()>ptmin ))	continue;
      if(!( fabs(this_electron.scEta())<fetamax )) continue;
      if(!( this_electron.PassID("passMediumID") ))	continue;
      if(!electron_isGsfCtfScPixChargeConsistent->at(i)) continue;
      out.push_back(this_electron);
    }
  }else if(id=="passTightID_SelQ"){
    for(unsigned int i=0; i<electrons.size(); i++){
      Electron this_electron= electrons.at(i);
      if(!( this_electron.Pt()>ptmin ))	continue;
      if(!( fabs(this_electron.scEta())<fetamax )) continue;
      if(!( this_electron.PassID("passTightID") )) continue;
      if(!electron_isGsfCtfScPixChargeConsistent->at(i)) continue;
      out.push_back(this_electron);
    }
  }else if(id=="passMediumIDWithAntiIso"){
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
    for(unsigned int i=0; i<electrons.size(); i++){
      Electron el= electrons.at(i);
      if(!( el.Pt()>ptmin ))	continue;
      if(!( fabs(el.scEta())<fetamax )) continue;
      if( el.etaRegion()==Electron::GAP ) continue;
      if( el.PassID("passMediumID") ) continue;
      out.push_back(el);
    }
  }else if(id=="passAntiLooseID"){
    for(unsigned int i=0; i<electrons.size(); i++){
      Electron el= electrons.at(i);
      if(!( el.Pt()>ptmin ))	continue;
      if(!( fabs(el.scEta())<fetamax )) continue;
      if( el.etaRegion()==Electron::GAP ) continue;
      if( el.PassID("passLooseID") ) continue;
      out.push_back(el);
    }
  }else out=SelectElectrons(allels,id,ptmin,fetamax);
  std::sort(out.begin(),out.end(),PtComparing);
  return out;
}    
std::vector<Muon> SMPAnalyzerCore::SMPGetMuons(TString id,double ptmin,double fetamax){
  vector<Muon> out;
  if(id=="POGTightWithLooseTrkIso"){
    vector<Muon> muons=SelectMuons(allmus,"POGTight",ptmin,fetamax);
    for(auto const& muon: muons){
      //if(muon.TrkIso()/muon.Pt()<0.1) out.push_back(muon);
      if(muon.PassSelector(Muon::Selector::TkIsoLoose)) out.push_back(muon);
    }
  }else if(id=="POGMediumWithLooseTrkIso"){
    TString IDhip = GetEraShort()=="2016a"? "POGMedium" : "POGMedium_nohip";
    vector<Muon> muons=SelectMuons(allmus,IDhip,ptmin,fetamax);
    for(auto const& muon: muons){
      if(muon.PassSelector(Muon::Selector::TkIsoLoose)) out.push_back(muon);
    }
  }else if(id=="POGMediumWithAntiLooseTrkIso"){
    TString IDhip = GetEraShort()=="2016a"? "POGMedium": "POGMedium_nohip";
    vector<Muon> muons=SelectMuons(allmus,IDhip,ptmin,fetamax);
    for(auto const& muon: muons){
      if(muon.PassSelector(Muon::Selector::TkIsoLoose)) continue;
      out.push_back(muon);
    }
  }else if(id=="POGTightWithAntiIso"){
    vector<Muon> muons=SelectMuons(allmus,"POGTight",ptmin,fetamax);
    for(auto const& muon: muons){
      if(muon.RelIso()>0.3) out.push_back(muon);
    }
  }else if(id=="POGTightWithAntiMediumIso"){
    vector<Muon> muons=SelectMuons(allmus,"POGTight",ptmin,fetamax);
    for(auto const& muon: muons){
      if(muon.RelIso()>0.2) out.push_back(muon);
    }
  }else out=SelectMuons(allmus,id,ptmin,fetamax);
  return out;
}

std::vector<Muon> SMPAnalyzerCore::MuonMomentumCorrection(const vector<Muon>& muons,int sys,int set,int member){
  if(!roc) return std::vector<Muon>(muons);
  std::vector<Muon> out;
  for(auto muon:muons){
    double rc=1.;
    double rcerr=0.;
    if(set>=0){
      if(IsDATA){
        rc=roc->kScaleDT(muon.Charge(),muon.MiniAODPt(),muon.Eta(),muon.Phi(),set,member);
        rcerr=roc->kScaleDTerror(muon.Charge(),muon.MiniAODPt(),muon.Eta(),muon.Phi());
      }else{
	Gen gen=GetGenMatchedLepton(muon,gens);
	if(gen.IsEmpty()){
	  gRandom->SetSeed((run<<15)+(lumi<<10)+(event<<5)+muon.Eta()*100);
	  double u=gRandom->Rndm();
	  rc=roc->kSmearMC(muon.Charge(),muon.MiniAODPt(),muon.Eta(),muon.Phi(),muon.TrackerLayers(),u,set,member);
	  rcerr=roc->kSmearMCerror(muon.Charge(),muon.MiniAODPt(),muon.Eta(),muon.Phi(),muon.TrackerLayers(),u);
	}else{
	  rc=roc->kSpreadMC(muon.Charge(),muon.MiniAODPt(),muon.Eta(),muon.Phi(),gen.Pt(),set,member);
	  rcerr=roc->kSpreadMCerror(muon.Charge(),muon.MiniAODPt(),muon.Eta(),muon.Phi(),gen.Pt());
	}
      }      
    }
    muon.SetPtEtaPhiM(muon.MiniAODPt()*(rc+sys*rcerr),muon.Eta(),muon.Phi(),muon.M());
    out.push_back(muon);
  }
  std::sort(out.begin(),out.end(),PtComparing);
  return out;
}

std::vector<Electron> SMPAnalyzerCore::ElectronEnergyCorrection(const vector<Electron>& electrons,int set,int member){
  if(!rocele) return std::vector<Electron>(electrons);
  std::vector<Electron> out;
  for(auto electron:electrons){
    if(set>=0){
      double rc=1.;
      //double rcerr=0.;
      double el_eta=electron.scEta();
      double el_phi=electron.Phi();
      if(IsDATA){
	rc=rocele->kScaleDT(electron.UncorrPt(),el_eta,el_phi,electron.R9(),run,set,member);
      }else{	
	Gen gen=SMPGetGenMatchedLepton(electron,gens,1);
	gRandom->SetSeed((run<<15)+(lumi<<10)+(event<<5)+electron.Eta()*100);
	double u=gRandom->Rndm();
	if(!gen.IsEmpty()&&fabs(electron.Pt()/gen.Pt()-1.)<0.5){
	  rc=rocele->kSpreadMC(electron.UncorrPt(),el_eta,el_phi,electron.R9(),u,gen.Pt(),set,member);
	}else{
	  rc=rocele->kSmearMC(electron.UncorrPt(),el_eta,el_phi,electron.R9(),u,set,member);
	}
      }      
      //if(electron.Pt()>100) rc=1.;
      if(TMath::IsNaN(rc)) rc=1.;
      electron*=rc*electron.UncorrE()/electron.E();
    }else if(set==-1){ //no energe cor
      electron*=electron.UncorrE()/electron.E();
    }else{
      cout<<"[SMPAnalyzerCore::ElectronEnergyCorrection] wrong set "<<set<<endl;
      exit(ENODATA);
    }
    out.push_back(electron);
  }
  std::sort(out.begin(),out.end(),PtComparing);
  return out;
}

void SMPAnalyzerCore::FillCutflow(TString histname,TString label,double weight){
  TH1D* hist=NULL;
  auto it=maphist_TH1D.find(histname);
  if(it==maphist_TH1D.end()){
    hist=new TH1D(histname,"",1,0,1);
    hist->GetXaxis()->SetBinLabel(1,label);
    maphist_TH1D[histname]=hist;    
  }else hist=it->second;

  int nbin=hist->GetNbinsX();
  int ibin=0;
  for(int i=1;i<=nbin;i++){
    if(hist->GetXaxis()->GetBinLabel(i)==label){
      ibin=i;
    }
  }

  if(!ibin){
    hist->SetBins(nbin+1,0,nbin+1);
    ibin=nbin+1;
    hist->GetXaxis()->SetBinLabel(ibin,label);
  }
  hist->Fill(ibin-0.5,weight);
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

double SMPAnalyzerCore::GetBTaggingReweight_1a_2WP(const vector<Jet>& jets, JetTagging::Parameters jtpT, JetTagging::Parameters jtpL, string Syst){
  //Syst. usage ex.: "SystUpHTag"(all component variation for heavy flav(b,c).),
  //                 "SystUpHTagCorr"(variation of heavy flav(b,c) sf only for yearly correlated components)
  //change H->L for light flav., Up->Down for downward variation, Corr->UnCorr for yearly independent components

  if(IsDATA) return 1.;

  TString SystStr(Syst);
  double Prob_MC(1.), Prob_DATA(1.), SF(1.);
  bool Syst_HTag=false, Syst_LTag=false; int SystDir=0, CorrType=0;
  string SystKey;
  if(SystStr.Contains("Syst")){
    if     (SystStr.Contains("HTag")) Syst_HTag=true;
    else if(SystStr.Contains("LTag")) Syst_LTag=true;
    if     (SystStr.Contains("Up")  ) SystDir= 1;
    else if(SystStr.Contains("Down")) SystDir=-1;
    if     (SystStr.Contains("UnCorr")) CorrType=-1;
    else if(SystStr.Contains("Corr"))   CorrType= 1;
    if(SystDir==0){ cout<<"SystStr in not correct form"<<endl; exit(ENODATA); }
    if(!(Syst_HTag or Syst_LTag)){ cout<<"SystMode but no H/L mode assigned"<<endl; exit(ENODATA); }
  }

  for(unsigned int i=0; i<jets.size(); i++){
    int JetHadFlav = jets.at(i).hadronFlavour();
    bool ApplySyst=false;
    if     (Syst_HTag && (JetHadFlav==4 or JetHadFlav==5)){ ApplySyst=true; }
    else if(Syst_LTag && (JetHadFlav==0                 )){ ApplySyst=true; }

    if     (ApplySyst && CorrType==0) SystKey=SystDir>0? "up":"down";
    else if(ApplySyst && CorrType >0) SystKey=SystDir>0? "up_correlated":"down_correlated";
    else if(ApplySyst && CorrType <0) SystKey=SystDir>0? "up_uncorrelated":"down_uncorrelated";
    else                              SystKey="central";


    double this_MC_EffT = mcCorr->GetMCJetTagEff(jtpT.j_Tagger, jtpT.j_WP, jets.at(i).hadronFlavour(), jets.at(i).Pt(), jets.at(i).Eta());
    double this_MC_EffL = mcCorr->GetMCJetTagEff(jtpL.j_Tagger, jtpL.j_WP, jets.at(i).hadronFlavour(), jets.at(i).Pt(), jets.at(i).Eta());
    double this_SFT = mcCorr->GetJetTaggingSF(jtpT,
                                              jets.at(i).hadronFlavour(),
                                              jets.at(i).Pt(),
                                              jets.at(i).Eta(),
                                              jets.at(i).GetTaggerResult(jtpT.j_Tagger),
                                              SystKey );
    double this_SFL = mcCorr->GetJetTaggingSF(jtpL,
                                              jets.at(i).hadronFlavour(),
                                              jets.at(i).Pt(),
                                              jets.at(i).Eta(),
                                              jets.at(i).GetTaggerResult(jtpL.j_Tagger),
                                              SystKey );
    double this_DATA_EffT = this_MC_EffT*this_SFT;
    double this_DATA_EffL = this_MC_EffL*this_SFL;

    bool isTaggedT = jets.at(i).GetTaggerResult(jtpT.j_Tagger) > mcCorr->GetJetTaggingCutValue(jtpT.j_Tagger, jtpT.j_WP);
    bool isTaggedL = jets.at(i).GetTaggerResult(jtpL.j_Tagger) > mcCorr->GetJetTaggingCutValue(jtpL.j_Tagger, jtpL.j_WP);
    if(isTaggedT){
      Prob_MC *= this_MC_EffT;
      Prob_DATA *= this_DATA_EffT;
    }
    else if(isTaggedL){
      if(this_MC_EffL == this_MC_EffT) this_MC_EffL += 1E-10;
      Prob_MC *= this_MC_EffL - this_MC_EffT;
      Prob_DATA *= this_DATA_EffL - this_DATA_EffT;
    }
    else{
      Prob_MC *= 1.-this_MC_EffL;
      Prob_DATA *= 1.-this_DATA_EffL;
    }
  }

  if(Prob_MC>0. && Prob_DATA>0.) SF=Prob_DATA/Prob_MC;
  else SF=0.;

  return SF;
}

bool SMPAnalyzerCore::PUJetIDPass(Jet jet, TString ID){
  if(jet.Pt() >= 50) return true;

  if(DataEra=="2016preVFP" || DataEra=="2016postVFP"){
    if(ID=="Tight"){
      if(jet.Pt() >= 40){
        if(jet.PileupJetId() > 0.97) return true;
      }else if(jet.Pt() >= 30){
        if(jet.PileupJetId() > 0.94) return true;
      }else if(jet.Pt() >= 20){
        if(jet.PileupJetId() > 0.87) return true;
      }
    }
    else if(ID=="Medium"){
      if(jet.Pt() >= 40){
        if(jet.PileupJetId() > 0.93) return true;
      }else if(jet.Pt() >= 30){
        if(jet.PileupJetId() > 0.86) return true;
      }else if(jet.Pt() >= 20){
        if(jet.PileupJetId() > 0.62) return true;
      }
    }
    else{
      if(jet.Pt() >= 40){
        if(jet.PileupJetId() > -0.42) return true;
      }else if(jet.Pt() >= 30){
        if(jet.PileupJetId() > -0.71) return true;
      }else if(jet.Pt() >= 20){
        if(jet.PileupJetId() > -0.90) return true;
      }
    }
  }

  if(DataEra=="2017"){
    if(ID=="Tight"){
      if(jet.Pt() >= 40){
        if(jet.PileupJetId() > 0.98) return true;
      }else if(jet.Pt() >= 30){
        if(jet.PileupJetId() > 0.96) return true;
      }else if(jet.Pt() >= 20){
        if(jet.PileupJetId() > 0.90) return true;
      }
    }
    else if(ID=="Medium"){
      if(jet.Pt() >= 40){
        if(jet.PileupJetId() > 0.96) return true;
      }else if(jet.Pt() >= 30){
        if(jet.PileupJetId() > 0.90) return true;
      }else if(jet.Pt() >= 20){
        if(jet.PileupJetId() > 0.68) return true;
      }
    }
    else{
      if(jet.Pt() >= 40){
        if(jet.PileupJetId() > -0.19) return true;
      }else if(jet.Pt() >= 30){
        if(jet.PileupJetId() > -0.63) return true;
      }else if(jet.Pt() >= 20){
        if(jet.PileupJetId() > -0.88) return true;
      }
    }
  }

  if(DataEra=="2018"){
    if(ID=="Tight"){
      if(jet.Pt() >= 40){
        if(jet.PileupJetId() > 0.98) return true;
      }else if(jet.Pt() >= 30){
        if(jet.PileupJetId() > 0.96) return true;
      }else if(jet.Pt() >= 20){
        if(jet.PileupJetId() > 0.90) return true;
      }
    }
    else if(ID=="Medium"){
      if(jet.Pt() >= 40){
        if(jet.PileupJetId() > 0.96) return true;
      }else if(jet.Pt() >= 30){
        if(jet.PileupJetId() > 0.90) return true;
      }else if(jet.Pt() >= 20){
        if(jet.PileupJetId() > 0.68) return true;
      }
    }
    else{
      if(jet.Pt() >= 40){
        if(jet.PileupJetId() > -0.19) return true;
      }else if(jet.Pt() >= 30){
        if(jet.PileupJetId() > -0.63) return true;
      }else if(jet.Pt() >= 20){
        if(jet.PileupJetId() > -0.88) return true;
      }
    }
  }
  return false;
}

void SMPAnalyzerCore::SetupPUJetWeight(){
  TString datapath = getenv("DATA_DIR");
  TFile fPUID(datapath+"/"+GetEra()+"/ID/PUJet/PUID.root");
  vector<TString> IDs = {"T", "M", "L"};
  for(unsigned int i=0; i<IDs.size(); i++){
    cout<<"[SMPAnalyzerCore::SetupPUJetWeight] setting PUJetWeight with ID : "+IDs.at(i)<<endl;

    TString era = GetEra();
    if(era == "2016postVFP") era = "2016";
    else if(era == "2016preVFP") era = "2016APV";

    heff_data = (TH2F*)fPUID.Get("h2_eff_dataUL"+era+"_"+IDs.at(i));
    heff_mc   = (TH2F*)fPUID.Get("h2_eff_mcUL"+era+"_"+IDs.at(i));
    hmistag_data = (TH2F*)fPUID.Get("h2_mistag_dataUL"+era+"_"+IDs.at(i));
    hmistag_mc   = (TH2F*)fPUID.Get("h2_mistag_mcUL"+era+"_"+IDs.at(i));

    heff_data->SetDirectory(0);
    heff_mc->SetDirectory(0);
    hmistag_data->SetDirectory(0);
    hmistag_mc->SetDirectory(0);
  }

  fPUID.Close();
}

double SMPAnalyzerCore::GetPUJetWeight(const vector<Jet>& jets, TString ID, int sys){
  sys = 0;
  if(IsDATA) return 1.;

  vector<Gen> gens=GetGens();

  double Prob_MC(1.), Prob_DATA(1.);
  for(unsigned int i=0; i<jets.size(); i++){
    double jetpt = jets.at(i).Pt();
    double jeteta = jets.at(i).Eta();
    if(jets.at(i).Pt() < 20) cout<<"jet pt < 20GeV, something wrong"<<endl;;
    if(jets.at(i).Pt() > 50) continue;
    if(abs(jets.at(i).Eta()) > 2.5) continue;

    double this_DATA_eff = heff_data->GetBinContent(heff_data->FindBin(jetpt, jeteta));
    double this_MC_eff = heff_mc->GetBinContent(heff_mc->FindBin(jetpt, jeteta));
    double this_DATA_mistag = hmistag_data->GetBinContent(hmistag_data->FindBin(jetpt, jeteta));
    double this_MC_mistag = hmistag_mc->GetBinContent(hmistag_mc->FindBin(jetpt, jeteta));
    if(this_DATA_eff * this_MC_eff * this_DATA_mistag * this_MC_mistag == 0.) continue;

    bool isRealJet = false;
    isRealJet = (jets.at(i).GenHFHadronMatcherFlavour() >= 0.);
    bool isPassID = PUJetIDPass(jets.at(i), ID);

    if(isRealJet){
      if(isPassID){
        if(this_MC_eff == 0) this_MC_eff += 1E-4;
        Prob_DATA *= this_DATA_eff;
        Prob_MC *= this_MC_eff;
      }else{
        if(this_MC_eff == 1) this_MC_eff -= 1E-4;
        Prob_DATA *= 1.-this_DATA_eff;
        Prob_MC *= 1.-this_MC_eff;
      }
    }else{
      if(isPassID){
        if(this_MC_mistag == 0) this_MC_mistag += 1E-4;
        Prob_DATA *= this_DATA_mistag;
        Prob_MC *= this_MC_mistag;
      }else{
        if(this_MC_mistag == 1) this_MC_mistag -= 1E-4;
        Prob_DATA *= 1.-this_DATA_mistag;
        Prob_MC *= 1.-this_MC_mistag;
      }
    }
  }

  return Prob_DATA/Prob_MC;
}

bool SMPAnalyzerCore::isGenMatchedJet(const Jet& jet, const vector<Gen>& gens){
  
  if(IsDATA){
    cout<<"This is for MC, something wrong here" <<endl;
    return false;
  }

  //Before using Gen Jet, we should use anti-kt algorithm clustering gen particles
  //But now, I used parton as an axis of gen jet, so dR(jet, gen parton) will criteria for matching
  int NumNearGen = 0;
  double GenpTSum = 0.;
  for(unsigned int i=0; i<gens.size(); i++){
    //gens.at(i).Print();
    if(gens.at(i).Status() != 1) continue;
    if(abs(gens.at(i).PID()) ==12 || abs(gens.at(i).PID()) ==14 ||  abs(gens.at(i).PID()) ==16) continue;
    //if(abs(gens.at(i).PID()) >10 && abs(gens.at(i).PID()) <17) continue;
    //if(abs(gens.at(i).PID()) == 24 || gens.at(i).PID() == 22 || gens.at(i).PID() == 23 || gens.at(i).PID() == 25) continue;

    if(jet.DeltaR(gens.at(i)) <1.0) GenpTSum += gens.at(i).Pt();//return true;
    if(jet.DeltaR(gens.at(i)) <1.0) NumNearGen += 1;
  }
  cout<<"reco jet pT = "<<jet.Pt()<<", gen pT sum = "<<GenpTSum<<", and ratio = "<<GenpTSum/jet.Pt()<<" , #ofnearGen is "<<NumNearGen<<endl;
  if(GenpTSum > jet.Pt()*0.5 && GenpTSum < jet.Pt()*1.5) return true;
  return false;
}

//double SMPAnalyzerCore::bjetCharge(const Jet& jet, int mode){
double SMPAnalyzerCore::jetCharge(const Jet& jet, int mode, TString prefix, double eventweight, bool doFillHists){
  //In mode0, output is jet charge (Sum of pt weighted charge of tracks)
  double jetCharge = jet.Charge();
  if(mode == 0) return jetCharge;
  vector<Muon> bmuon;
  bmuon.clear();
  vector<Electron> belectron;
  belectron.clear();

  for(unsigned int l=0; l<allmus.size(); l++){
    if(allmus.at(l).P()*sin(allmus.at(l).Angle(jet.Vect())) <0.6) continue; // original, 1GeV
    if(allmus.at(l).TrkIso()/allmus.at(l).Pt() <0.05) continue; // original, 0.1
    if(abs(allmus.at(l).IP3D())/allmus.at(l).IP3Derr() <2.) continue; // original, 2.5
    if(jet.DeltaR(allmus.at(l))<0.4) bmuon.push_back(allmus.at(l));
  }

  //belectron Trial
  for(unsigned int l=0; l<allels.size(); l++){
    if(allels.at(l).P()*sin(allels.at(l).Angle(jet.Vect())) <0.6) continue;
    if(allels.at(l).ecalPFClusterIso()/allels.at(l).Pt() == 0.) continue;
    if(abs(allels.at(l).IP3D())/allels.at(l).IP3Derr() <2.0) continue;
    if(!allels.at(l).IsGsfCtfScPixChargeConsistent()) continue;
    if(jet.DeltaR(allels.at(l))<0.4) belectron.push_back(allels.at(l));
  }

  if(doFillHists){
    if(prefix!="") FillHist(prefix+"jetCharge_raw",jetCharge,eventweight,200, -2, 2);
  }
  //The jet has soft muon inside, and its charge will determine the jet charge
  if(bmuon.size() > 0) jetCharge += 2 * bmuon.at(0).Charge();
  else if(belectron.size() > 0) jetCharge += 4 * belectron.at(0).Charge();

  if(doFillHists){
    if(prefix!="") FillHist(prefix+"jetCharge_rawwide",jetCharge,eventweight,1000, -10, 10);
    if(prefix!="" && jetCharge > 0.) FillHist(prefix+"jetCharge_rawwide_plus",jetCharge,eventweight,500, 0, 10);
    if(prefix!="" && jetCharge > 0.1) FillHist(prefix+"jetCharge_rawwide_plus01",jetCharge,eventweight,500, 0, 10);
    if(prefix!="" && jetCharge > 0.2) FillHist(prefix+"jetCharge_rawwide_plus02",jetCharge,eventweight,500, 0, 10);
    if(prefix!="" && jetCharge > 0.3) FillHist(prefix+"jetCharge_rawwide_plus03",jetCharge,eventweight,500, 0, 10);
    if(prefix!="" && jetCharge > 0.4) FillHist(prefix+"jetCharge_rawwide_plus04",jetCharge,eventweight,500, 0, 10);
    if(prefix!="" && jetCharge < 0.) FillHist(prefix+"jetCharge_rawwide_minus",jetCharge,eventweight,500, -10, 0);
    if(prefix!="" && jetCharge < -0.1) FillHist(prefix+"jetCharge_rawwide_minus01",jetCharge,eventweight,500, -10, 0);
    if(prefix!="" && jetCharge < -0.2) FillHist(prefix+"jetCharge_rawwide_minus02",jetCharge,eventweight,500, -10, 0);
    if(prefix!="" && jetCharge < -0.3) FillHist(prefix+"jetCharge_rawwide_minus03",jetCharge,eventweight,500, -10, 0);
    if(prefix!="" && jetCharge < -0.4) FillHist(prefix+"jetCharge_rawwide_minus04",jetCharge,eventweight,500, -10, 0);
  }
  return jetCharge;
}

SMPAnalyzerCore::Parameter::Parameter(){
}
SMPAnalyzerCore::Parameter::~Parameter(){
}
void SMPAnalyzerCore::Parameter::SetChannel(TString ch){
  vector<TString> availables={"el","ee","eE","EE","mu","mm","mM","MM","em","me"};
  bool pass=false;
  for(const TString& avail:availables)
    if(ch(0,2)==avail) pass=true;
  if(!pass){
    cout<<"[SMPAnalyzerCore::Parameter::SetChannel] not available channel "<<ch<<endl;
    exit(EXIT_FAILURE);
  }
  channel=ch;
  SetLeptons();
  if(channel.Length()>2) SetJets();
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
void SMPAnalyzerCore::Parameter::SetJetPtCut(double j0pt,double j1pt){
  c.jet0pt=j0pt;
  c.jet1pt=j1pt;
}
void SMPAnalyzerCore::Parameter::SetLeptons(){
  leptons={};
  lepton0=NULL;
  lepton1=NULL;
  truth_lepton0=Gen();
  truth_lepton1=Gen();
  if(channel=="") return;
  unsigned int ie=0,iae=0,im=0,iam=0;
  for(int ic=0;ic<2;ic++){
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
    }
  }
  std::sort(leptons.begin(),leptons.end(),PtComparingPtr);
  if(leptons.size()>0) lepton0=leptons.at(0);
  if(leptons.size()>1) lepton1=leptons.at(1);
  if(lepton0) truth_lepton0=SMPGetGenMatchedLepton(*lepton0,gens);
  if(lepton1) truth_lepton1=SMPGetGenMatchedLepton(*lepton1,gens);
}
void SMPAnalyzerCore::Parameter::SetJets(){
  jets={};
  jet0=NULL;
  jet1=NULL;
  truth_jet0=Gen();
  truth_jet1=Gen();
  if(channel=="") return;
  unsigned int iB=0,iC=0,iL=0,iA=0;
  int nc=channel.Length();
  if(nc>2){
    for(int ic=2;ic<nc;ic++){
      char c=channel[ic];
      if(c=='b'||c=='B'){
        if(bjets.size()>iB){
          jets.push_back(&bjets.at(iB));
          iB++;
        }else jets.push_back(NULL);
      }if(c=='c'||c=='C'){
        if(cjets.size()>iC){
          jets.push_back(&cjets.at(iC));
          iC++;
        }else jets.push_back(NULL);
      }if(c=='l'||c=='L'||c=='j'){
        if(ljets.size()>iL){
          jets.push_back(&ljets.at(iL));
          iL++;
        }else jets.push_back(NULL);
      }if(c=='x'){
        if(ajets.size()>iA){
          jets.push_back(&ajets.at(iA));
          iA++;
        }else jets.push_back(NULL);
      }
    }
  }
  //std::sort(jets.begin(),jets.end(),PtComparingPtr);
  if(jets.size()>0) jet0=jets.at(0);
  if(jets.size()>1) jet1=jets.at(1);
  //if(jet0) truth_jet0=SMPGetGenMatchedJet(*jet0,gens);
  //if(jet1) truth_jet1=SMPGetGenMatchedJet(*jet1,gens);
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
void SMPAnalyzerCore::Parameter::SetBJets(vector<Jet> bs){
  bjets=bs;
  SetJets();
}
void SMPAnalyzerCore::Parameter::SetCJets(vector<Jet> cs){
  cjets=cs;
  SetJets();
}
void SMPAnalyzerCore::Parameter::SetLJets(vector<Jet> ls){
  ljets=ls;
  SetJets();
}
void SMPAnalyzerCore::Parameter::SetAJets(vector<Jet> as){
  ajets=as;
  SetJets();
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
  p.w.pujetSF=1;
  p.w.tagjetSF=1;
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
	p.w.zptweight=GetZptWeight(genZ.M(),genZ.Rapidity(),genZ.Pt());
	p.w.weakweight=GetDYWeakWeight(genZ.M());

	/*
        // Only qqbar collisions (LO DY)
        if(lhe_p0.ID()+lhe_p1.ID()==0) p.hprefix+="";
        // Only qG collisions (NLO DY)
        else if((abs(lhe_p0.ID())<=5&&lhe_p1.ID()==21) || (lhe_p0.ID()==21&&abs(lhe_p1.ID())<=5)){
          //if((abs(gen_p0.PID())<=5&&gen_p1.PID()==21)) p.hprefix+="qG_";
          //else p.hprefix+="Gq_";
          if(lhe_p0.ID()==5||lhe_p1.ID()==5) p.hprefix+="Dyb_";
          else if(lhe_p0.ID()==-5||lhe_p1.ID()==-5) p.hprefix+="Dybbar_";
          else if(lhe_p0.ID()==4||lhe_p1.ID()==4) p.hprefix+="Dyc_";
          else if(lhe_p0.ID()==-4||lhe_p1.ID()==-4) p.hprefix+="Dycbar_";
          else p.hprefix+="Dyudsg_";
        } // Only GG collisions (NNLO DY)
        else if(lhe_p0.ID()==21 && lhe_p1.ID()==21){
          Gen heavyparton = gens.at(0);
          int nheavyparton = 0;
          for(unsigned int i=0; i<gens.size(); i++){
            if(!gens.at(i).isHardProcess()) continue;
            if(abs(gens.at(i).PID())>=11 && abs(gens.at(i).PID())<=16) continue; // No Lepton
            if(gens.at(i).PID()==22 || gens.at(i).PID()==23) continue; // No gamma, Z
            if(gens.at(i).Pt() < 30 || abs(gens.at(i).Eta())>2.4) continue; // In the acceptance

            if(nheavyparton==0 && (abs(gens.at(i).PID())==4 || abs(gens.at(i).PID())==5)){
              heavyparton=gens.at(i);
              nheavyparton++;
              continue;
            }
            else if(nheavyparton>0 && (abs(gens.at(i).PID())==4 || abs(gens.at(i).PID())==5)){
              heavyparton=(heavyparton.Pt()>gens.at(i).Pt()?heavyparton:gens.at(i));
              break;
            }
          }

          if(nheavyparton>0 && heavyparton.PID()==5) p.hprefix+="Dyb_";//"Dyggb_";
          else if(nheavyparton>0 && heavyparton.PID()==-5) p.hprefix+="Dybbar_";//"Dyggbbar_";
          else if(nheavyparton>0 && heavyparton.PID()==4) p.hprefix+="Dyc_";//"Dyggc_";
          else if(nheavyparton>0 && heavyparton.PID()==-4) p.hprefix+="Dycbar_";//"Dyggcbar_";
          else p.hprefix+="";//"Dygg_";
        } // Only bq or cq collisions (NNLO DY)
        else if(abs(lhe_p0.ID())==4 || abs(lhe_p0.ID())==5 || abs(lhe_p1.ID())==4 || abs(lhe_p1.ID())==5){
          Gen heavyparton = gens.at(0);
          int nheavyparton = 0;
          for(unsigned int i=0; i<gens.size(); i++){
            if(!gens.at(i).isHardProcess()) continue;
            if(abs(gens.at(i).PID())>=11 && abs(gens.at(i).PID())<=16) continue; // No Lepton
            if(gens.at(i).PID()==22 || gens.at(i).PID()==23) continue; // No gamma, Z
            if(gens.at(i).Pt() < 30 || abs(gens.at(i).Eta())>2.4) continue; // In the acceptance

            if(nheavyparton==0 && (abs(gens.at(i).PID())==4 || abs(gens.at(i).PID())==5)){
              heavyparton=gens.at(i);
              nheavyparton++;
              continue;
            }
            else if(nheavyparton>0 && (abs(gens.at(i).PID())==4 || abs(gens.at(i).PID())==5)){
              heavyparton=(heavyparton.Pt()>gens.at(i).Pt()?heavyparton:gens.at(i));
              break;
            }
          }
          if(nheavyparton>0 && heavyparton.PID()==5) p.hprefix+="Dyb_";//"Dyqqb_";
          else if(nheavyparton>0 && heavyparton.PID()==-5) p.hprefix+="Dybbar_";//"Dyqqbbar_";
          else if(nheavyparton>0 && heavyparton.PID()==4) p.hprefix+="Dyc_";//"Dyqqc_";
          else if(nheavyparton>0 && heavyparton.PID()==-4) p.hprefix+="Dycbar_";//"Dyqqcbar_";
          else p.hprefix+="";//"Dyqq_";
        }
	  */
      }else p.hprefix+="tau_";
    }
  }

  p.prefix=p.channel+GetEraShort()+"/";
  if(p.channel(0,2)=="mu"){
    p.SetMuonKeys("Muon_MediumID_trkIsoLoose","",{"IsoMu24_MediumID_trkIsoLoose"});
    p.SetMuons(MuonMomentumCorrection(SMPGetMuons("POGMediumWithLooseTrkIso",0.0,2.4),0,0));
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
  }else if(p.channel(0,2)=="mm"){
    p.SetMuonKeys("Muon_MediumID_trkIsoLoose","",{"Mu17Leg1_MediumID_trkIsoLoose","Mu8Leg2_MediumID_trkIsoLoose"});
    p.SetMuons(MuonMomentumCorrection(SMPGetMuons("POGMediumWithLooseTrkIso",0.0,2.4),0,0));
    p.SetLeptonPtCut(20,10);
    if(GetEraShort()=="2016a"){
      p.triggers={"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v","HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v","HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v","HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",};
    }else if(GetEraShort()=="2016b"){
      p.triggers={"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v","HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v","HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v","HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",};
    }else if(GetEraShort()=="2017"){
      p.triggers={"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v"};
    }else if(GetEraShort()=="2018"){
      p.triggers={"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"};
    }
  }else if(p.option.Contains("SelQ")&&p.channel=="el"){
    p.prefix="selq/"+p.prefix;
    p.SetElectronKeys("Electron_MediumID","Electron_SelQ_MediumID",{"Ele27_SelQ_MediumID"});
    p.SetElectrons(ElectronEnergyCorrection(SMPGetElectrons("passMediumID_SelQ",0.0,2.5),0,0));
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
  }else if(p.channel(0,2)=="el"){
    p.SetElectronKeys("Electron_MediumID",{"Ele27_MediumID"});
    p.SetElectrons(ElectronEnergyCorrection(SMPGetElectrons("passMediumID",0.0,2.5),0,0));
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
  }else if(p.channel(0,2)=="ee"){
    p.SetElectronKeys("Electron_MediumID",{"Ele23Leg1_MediumID","Ele12Leg2_MediumID"});
    p.SetElectrons(ElectronEnergyCorrection(SMPGetElectrons("passMediumID",0.0,2.5),0,0));
    p.SetLeptonPtCut(25,15);
    if(GetEraShort()=="2016a") p.triggers={"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"};
    else if(GetEraShort()=="2016b") p.triggers={"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"};
    else if(GetEraShort()=="2017") p.triggers={"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"};
    else if(GetEraShort()=="2018") p.triggers={"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"};
  }else if(p.channel(0,2)=="me"){
    p.SetMuonKeys("Muon_MediumID_trkIsoLoose","",{"IsoMu24_MediumID_trkIsoLoose"});
    p.k.electronIDSF="Electron_MediumID";
    p.SetMuons(MuonMomentumCorrection(SMPGetMuons("POGMediumWithLooseTrkIso",0.0,2.4),0,0));
    p.SetElectrons(SMPGetElectrons("passMediumID",0.0,2.5));
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
  }else if(p.channel(0,2)=="em"){
    p.SetElectronKeys("Electron_MediumID",{"Ele27_MediumID"});
    p.k.muonIDSF="Muon_MediumID_trkIsoLoose";
    p.SetElectrons(SMPGetElectrons("passMediumID",0.0,2.5));
    p.SetMuons(MuonMomentumCorrection(SMPGetMuons("POGMediumWithLooseTrkIso",0.0,2.4),0,0));
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
  }else if(p.channel(0,2)=="mM"){
    p.k.muonIDSF="Muon_MediumID_trkIsoLoose";
    p.SetMuons(MuonMomentumCorrection(SMPGetMuons("POGMediumWithLooseTrkIso",0.0,2.4),0,0));
    p.SetAMuons(MuonMomentumCorrection(SMPGetMuons("POGMediumWithAntiLooseTrkIso",0.0,2.4),0,0));
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
    p.SetMuonKeys("Muon_MediumID_trkIsoLoose","",{"IsoMu24_MediumID_trkIsoLoose"});
    p.SetMuons(MuonMomentumCorrection(SMPGetMuons("POGMediumWithLooseTrkIso",0.0,2.4),0,0));
    p.SetAMuons(MuonMomentumCorrection(SMPGetMuons("POGMediumWithAntiLooseTrkIso",0.0,2.4),0,0));
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
  }else if(p.channel(0,2)=="MM"){
    p.k.muonIDSF="Muon_MediumID_trkIsoLoose";
    p.SetMuons(MuonMomentumCorrection(SMPGetMuons("POGMediumWithLooseTrkIso",0.0,2.4),0,0));
    p.SetAMuons(MuonMomentumCorrection(SMPGetMuons("POGMediumWithAntiLooseTrkIso",0.0,2.4),0,0));
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
  }else if(p.channel(0,2)=="eE"){
    p.k.electronIDSF="Electron_MediumID";
    p.SetElectrons(SMPGetElectrons("passMediumID",0.0,2.5));
    p.SetAElectrons(SMPGetElectrons("passAntiLooseID",0.0,2.5));
    p.SetLeptonPtCut(25,15);
    //p.c.nelectronmax=1;
    p.option+=" triggermatching strictorder";
    if(GetEraShort()=="2016a") p.triggers={"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"};
    else if(GetEraShort()=="2016b") p.triggers={"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"};
    else if(GetEraShort()=="2017") p.triggers={"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"};
    else if(GetEraShort()=="2018") p.triggers={"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"};
    p.SetElectronKeys("Electron_MediumID",{"Ele27_MediumID"});
    p.SetElectrons(SMPGetElectrons("passMediumID",0.0,2.5));
    p.SetAElectrons(SMPGetElectrons("passAntiLooseID",0.0,2.5));
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
  }else if(p.channel(0,2)=="EE"){
    p.k.electronIDSF="Electron_MediumID";
    p.SetElectrons(SMPGetElectrons("passMediumID",0.0,2.5));
    p.SetAElectrons(SMPGetElectrons("passAntiLooseID",0.0,2.5));
    //p.c.nelectronmax=0;
    p.option+=" triggermatching strictorder";
    p.SetLeptonPtCut(25,15);
    if(GetEraShort()=="2016a") p.triggers={"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"};
    else if(GetEraShort()=="2016b") p.triggers={"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"};
    else if(GetEraShort()=="2017") p.triggers={"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"};
    else if(GetEraShort()=="2018") p.triggers={"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"};
  }
  else{
    cout<<"[SMPAnalyzerCore::MakeParameter] not available setting "<<p.channel<<endl;
    exit(EXIT_FAILURE);
  }

  if(p.channel.Length()>2){
    vector<Jet> lepvetojets = {};
    realjets.clear();
    for(unsigned int l=0; l<alljets.size();l++){
      if(p.lepton0){if(p.lepton0->DeltaR(alljets.at(l)) <0.4) continue;}
      if(p.lepton1){if(p.lepton1->DeltaR(alljets.at(l)) <0.4) continue;}
      lepvetojets.push_back(alljets.at(l));
    }
    for(unsigned int l=0; l<lepvetojets.size();l++){
      if(!PUJetIDPass(alljets.at(l), "Loose")) continue;
      realjets.push_back(alljets.at(l));
    }

    p.w.pujetSF = GetPUJetWeight(lepvetojets, "Loose", 0);

    if(p.channel(2,2)=="bx"){
      p.SetJetPtCut(30,20);
      p.c.nbjetmax=1;

      JetTagging::Parameters DeepJet_Tight = JetTagging::Parameters(JetTagging::DeepJet,JetTagging::Tight,JetTagging::incl,JetTagging::comb);
      JetTagging::Parameters DeepJet_Loose = JetTagging::Parameters(JetTagging::DeepJet,JetTagging::Loose,JetTagging::incl,JetTagging::comb);
      p.bjets.clear();
      p.ajets.clear();

      for(const auto& jet:realjets){
        if(jet.Pt() > p.c.jet0pt && jet.GetTaggerResult(DeepJet_Tight.j_Tagger) > mcCorr->GetJetTaggingCutValue(DeepJet_Tight.j_Tagger, DeepJet_Tight.j_WP)) p.bjets.push_back(jet);
        else if(jet.Pt() > p.c.jet1pt && jet.GetTaggerResult(DeepJet_Loose.j_Tagger) > mcCorr->GetJetTaggingCutValue(DeepJet_Loose.j_Tagger, DeepJet_Loose.j_WP)) p.ajets.push_back(jet);
      }
      p.SetBJets(p.bjets);
      p.SetAJets(p.ajets);

      p.w.tagjetSF = GetBTaggingReweight_1a_2WP(realjets, DeepJet_Tight, DeepJet_Loose, "central");
      p.doublemap["btagSF_hup"] = GetBTaggingReweight_1a_2WP(realjets, DeepJet_Tight, DeepJet_Loose, "Syst HTag Up Corr");
      p.doublemap["btagSF_hdown"] = GetBTaggingReweight_1a_2WP(realjets, DeepJet_Tight, DeepJet_Loose, "Syst HTag Down Corr");
      p.doublemap["btagSF_lup"] = GetBTaggingReweight_1a_2WP(realjets, DeepJet_Tight, DeepJet_Loose, "Syst LTag Up Corr");
      p.doublemap["btagSF_ldown"] = GetBTaggingReweight_1a_2WP(realjets, DeepJet_Tight, DeepJet_Loose, "Syst LTag Down Corr");
      p.doublemap["btagSF_hup"+GetEra()] = GetBTaggingReweight_1a_2WP(realjets, DeepJet_Tight, DeepJet_Loose, "Syst HTag Up UnCorr");
      p.doublemap["btagSF_hdown"+GetEra()] = GetBTaggingReweight_1a_2WP(realjets, DeepJet_Tight, DeepJet_Loose, "Syst HTag Down UnCorr");
      p.doublemap["btagSF_lup"+GetEra()] = GetBTaggingReweight_1a_2WP(realjets, DeepJet_Tight, DeepJet_Loose, "Syst LTag Up UnCorr");
      p.doublemap["btagSF_ldown"+GetEra()] = GetBTaggingReweight_1a_2WP(realjets, DeepJet_Tight, DeepJet_Loose, "Syst LTag Down UnCorr");
    }else if(p.channel(2,2)=="Bx"){
      p.SetJetPtCut(30,20);
      p.c.nbjetmax=1;

      JetTagging::Parameters DeepCSV_Tight = JetTagging::Parameters(JetTagging::DeepCSV,JetTagging::Tight,JetTagging::incl,JetTagging::comb);
      JetTagging::Parameters DeepCSV_Loose = JetTagging::Parameters(JetTagging::DeepCSV,JetTagging::Loose,JetTagging::incl,JetTagging::comb);
      p.bjets.clear();
      p.ajets.clear();

      for(const auto& jet:realjets){
        if(jet.Pt() > p.c.jet0pt && jet.GetTaggerResult(DeepCSV_Tight.j_Tagger) > mcCorr->GetJetTaggingCutValue(DeepCSV_Tight.j_Tagger, DeepCSV_Tight.j_WP)) p.bjets.push_back(jet);
        else if(jet.Pt() > p.c.jet1pt && jet.GetTaggerResult(DeepCSV_Loose.j_Tagger) > mcCorr->GetJetTaggingCutValue(DeepCSV_Loose.j_Tagger, DeepCSV_Loose.j_WP)) p.ajets.push_back(jet);
      }
      p.SetBJets(p.bjets);
      p.SetAJets(p.ajets);

      p.w.tagjetSF = GetBTaggingReweight_1a_2WP(realjets, DeepCSV_Tight, DeepCSV_Loose, "central");
      p.doublemap["btagSF_hup"] = GetBTaggingReweight_1a_2WP(realjets, DeepCSV_Tight, DeepCSV_Loose, "Syst HTag Up Corr");
      p.doublemap["btagSF_hdown"] = GetBTaggingReweight_1a_2WP(realjets, DeepCSV_Tight, DeepCSV_Loose, "Syst HTag Down Corr");
      p.doublemap["btagSF_lup"] = GetBTaggingReweight_1a_2WP(realjets, DeepCSV_Tight, DeepCSV_Loose, "Syst LTag Up Corr");
      p.doublemap["btagSF_ldown"] = GetBTaggingReweight_1a_2WP(realjets, DeepCSV_Tight, DeepCSV_Loose, "Syst LTag Down Corr");
      p.doublemap["btagSF_hup"+GetEra()] = GetBTaggingReweight_1a_2WP(realjets, DeepCSV_Tight, DeepCSV_Loose, "Syst HTag Up UnCorr");
      p.doublemap["btagSF_hdown"+GetEra()] = GetBTaggingReweight_1a_2WP(realjets, DeepCSV_Tight, DeepCSV_Loose, "Syst HTag Down UnCorr");
      p.doublemap["btagSF_lup"+GetEra()] = GetBTaggingReweight_1a_2WP(realjets, DeepCSV_Tight, DeepCSV_Loose, "Syst LTag Up UnCorr");
      p.doublemap["btagSF_ldown"+GetEra()] = GetBTaggingReweight_1a_2WP(realjets, DeepCSV_Tight, DeepCSV_Loose, "Syst LTag Down UnCorr");
    }else if(p.channel(2,2)=="bb"){
      p.SetJetPtCut(30,30);
      p.c.nbjetmax=2;
      p.c.nbjetmin=2;

      JetTagging::Parameters DeepJet_Tight = JetTagging::Parameters(JetTagging::DeepJet,JetTagging::Tight,JetTagging::incl,JetTagging::comb);
      p.bjets.clear();

      for(const auto& jet:realjets){
        if(jet.Pt() > p.c.jet0pt && jet.GetTaggerResult(DeepJet_Tight.j_Tagger) > mcCorr->GetJetTaggingCutValue(DeepJet_Tight.j_Tagger, DeepJet_Tight.j_WP)) p.bjets.push_back(jet);
      }
      p.SetBJets(p.bjets);

      p.w.tagjetSF = mcCorr->GetBTaggingReweight_1a(realjets, DeepJet_Tight, "central");
      p.doublemap["btagSF_hup"] = mcCorr->GetBTaggingReweight_1a(realjets, DeepJet_Tight, "Syst HTag Up Corr");
      p.doublemap["btagSF_hdown"] = mcCorr->GetBTaggingReweight_1a(realjets, DeepJet_Tight, "Syst HTag Down Corr");
      p.doublemap["btagSF_lup"] = mcCorr->GetBTaggingReweight_1a(realjets, DeepJet_Tight, "Syst LTag Up Corr");
      p.doublemap["btagSF_ldown"] = mcCorr->GetBTaggingReweight_1a(realjets, DeepJet_Tight, "Syst LTag Down Corr");
      p.doublemap["btagSF_hup"+GetEra()] = mcCorr->GetBTaggingReweight_1a(realjets, DeepJet_Tight, "Syst HTag Up UnCorr");
      p.doublemap["btagSF_hdown"+GetEra()] = mcCorr->GetBTaggingReweight_1a(realjets, DeepJet_Tight, "Syst HTag Down UnCorr");
      p.doublemap["btagSF_lup"+GetEra()] = mcCorr->GetBTaggingReweight_1a(realjets, DeepJet_Tight, "Syst LTag Up UnCorr");
      p.doublemap["btagSF_ldown"+GetEra()] = mcCorr->GetBTaggingReweight_1a(realjets, DeepJet_Tight, "Syst LTag Down UnCorr");
    }else if(p.channel(2,2)=="BB"){
      p.SetJetPtCut(30,30);
      p.c.nbjetmax=2;
      p.c.nbjetmin=2;

      JetTagging::Parameters DeepCSV_Tight = JetTagging::Parameters(JetTagging::DeepCSV,JetTagging::Tight,JetTagging::incl,JetTagging::comb);
      p.bjets.clear();

      for(const auto& jet:realjets){
        if(jet.Pt() > p.c.jet0pt && jet.GetTaggerResult(DeepCSV_Tight.j_Tagger) > mcCorr->GetJetTaggingCutValue(DeepCSV_Tight.j_Tagger, DeepCSV_Tight.j_WP)) p.bjets.push_back(jet);
      }
      p.SetBJets(p.bjets);

      p.w.tagjetSF = mcCorr->GetBTaggingReweight_1a(realjets, DeepCSV_Tight, "central");
      p.doublemap["btagSF_hup"] = mcCorr->GetBTaggingReweight_1a(realjets, DeepCSV_Tight, "Syst HTag Up Corr");
      p.doublemap["btagSF_hdown"] = mcCorr->GetBTaggingReweight_1a(realjets, DeepCSV_Tight, "Syst HTag Down Corr");
      p.doublemap["btagSF_lup"] = mcCorr->GetBTaggingReweight_1a(realjets, DeepCSV_Tight, "Syst LTag Up Corr");
      p.doublemap["btagSF_ldown"] = mcCorr->GetBTaggingReweight_1a(realjets, DeepCSV_Tight, "Syst LTag Down Corr");
      p.doublemap["btagSF_hup"+GetEra()] = mcCorr->GetBTaggingReweight_1a(realjets, DeepCSV_Tight, "Syst HTag Up UnCorr");
      p.doublemap["btagSF_hdown"+GetEra()] = mcCorr->GetBTaggingReweight_1a(realjets, DeepCSV_Tight, "Syst HTag Down UnCorr");
      p.doublemap["btagSF_lup"+GetEra()] = mcCorr->GetBTaggingReweight_1a(realjets, DeepCSV_Tight, "Syst LTag Up UnCorr");
      p.doublemap["btagSF_ldown"+GetEra()] = mcCorr->GetBTaggingReweight_1a(realjets, DeepCSV_Tight, "Syst LTag Down UnCorr");
    }if(p.channel(2,2)=="cx"){
      p.SetJetPtCut(30,20);
      p.c.ncjetmax=1;

      JetTagging::Parameters DeepJet_Tight = JetTagging::Parameters(JetTagging::DeepJet,JetTagging::Tight,JetTagging::incl,JetTagging::comb);

      JetTagging::Parameters DeepJet_CvsB_Tight = JetTagging::Parameters(JetTagging::DeepJet_CvsB,JetTagging::Tight,JetTagging::incl,JetTagging::wcharm);
      JetTagging::Parameters DeepJet_CvsL_Tight = JetTagging::Parameters(JetTagging::DeepJet_CvsL,JetTagging::Tight,JetTagging::incl,JetTagging::wcharm);
      JetTagging::Parameters DeepJet_CvsB_Loose = JetTagging::Parameters(JetTagging::DeepJet_CvsB,JetTagging::Loose,JetTagging::incl,JetTagging::wcharm);
      JetTagging::Parameters DeepJet_CvsL_Loose = JetTagging::Parameters(JetTagging::DeepJet_CvsL,JetTagging::Loose,JetTagging::incl,JetTagging::wcharm);
      p.cjets.clear();
      p.ajets.clear();

      for(const auto& jet:realjets){
        if(jet.Pt() > p.c.jet0pt && jet.GetTaggerResult(DeepJet_CvsB_Tight.j_Tagger) > mcCorr->GetJetTaggingCutValue(DeepJet_CvsB_Tight.j_Tagger, DeepJet_CvsB_Tight.j_WP) && jet.GetTaggerResult(DeepJet_CvsL_Tight.j_Tagger) > mcCorr->GetJetTaggingCutValue(DeepJet_CvsL_Tight.j_Tagger, DeepJet_CvsL_Tight.j_WP)) p.cjets.push_back(jet);
        else if(jet.Pt() > p.c.jet1pt && jet.GetTaggerResult(DeepJet_CvsB_Loose.j_Tagger) > mcCorr->GetJetTaggingCutValue(DeepJet_CvsB_Loose.j_Tagger, DeepJet_CvsB_Loose.j_WP) && jet.GetTaggerResult(DeepJet_CvsL_Loose.j_Tagger) > mcCorr->GetJetTaggingCutValue(DeepJet_CvsL_Loose.j_Tagger, DeepJet_CvsL_Loose.j_WP)) p.ajets.push_back(jet);
      }
      for(const auto& jet:p.cjets){
	if(jet.Pt() > p.c.jet0pt && jet.GetTaggerResult(DeepJet_Tight.j_Tagger) > mcCorr->GetJetTaggingCutValue(DeepJet_Tight.j_Tagger, DeepJet_Tight.j_WP)) p.bjets.push_back(jet);
      }
      p.SetCJets(p.cjets);
      p.SetAJets(p.ajets);

      p.w.tagjetSF = 1.;
      p.doublemap["ctagSF_hup"] = 1.;
      p.doublemap["ctagSF_hdown"] = 1.;
      p.doublemap["ctagSF_lup"] = 1.;
      p.doublemap["ctagSF_ldown"] = 1.;
      p.doublemap["ctagSF_hup"+GetEra()] = 1.;
      p.doublemap["ctagSF_hdown"+GetEra()] = 1.;
      p.doublemap["ctagSF_lup"+GetEra()] = 1.;
      p.doublemap["ctagSF_ldown"+GetEra()] = 1.;
    }if(p.channel(2,2)=="lx"){
      p.SetJetPtCut(30,20);
      p.c.nljetmax=1;

      JetTagging::Parameters DeepJet_Tight = JetTagging::Parameters(JetTagging::DeepJet,JetTagging::Tight,JetTagging::incl,JetTagging::comb);
      JetTagging::Parameters DeepJet_Loose = JetTagging::Parameters(JetTagging::DeepJet,JetTagging::Loose,JetTagging::incl,JetTagging::comb);
      JetTagging::Parameters DeepJet_CvsB_Tight = JetTagging::Parameters(JetTagging::DeepJet_CvsB,JetTagging::Tight,JetTagging::incl,JetTagging::wcharm);
      JetTagging::Parameters DeepJet_CvsL_Tight = JetTagging::Parameters(JetTagging::DeepJet_CvsL,JetTagging::Tight,JetTagging::incl,JetTagging::wcharm);
      JetTagging::Parameters DeepJet_CvsB_Loose = JetTagging::Parameters(JetTagging::DeepJet_CvsB,JetTagging::Loose,JetTagging::incl,JetTagging::wcharm);
      JetTagging::Parameters DeepJet_CvsL_Loose = JetTagging::Parameters(JetTagging::DeepJet_CvsL,JetTagging::Loose,JetTagging::incl,JetTagging::wcharm);
      p.ljets.clear();
      p.ajets.clear();

      for(const auto& jet:realjets){
	if(jet.Pt() > p.c.jet0pt && jet.GetTaggerResult(DeepJet_Loose.j_Tagger) < mcCorr->GetJetTaggingCutValue(DeepJet_Loose.j_Tagger, DeepJet_Loose.j_WP) && jet.GetTaggerResult(DeepJet_CvsL_Loose.j_Tagger) < mcCorr->GetJetTaggingCutValue(DeepJet_CvsL_Loose.j_Tagger, DeepJet_CvsL_Loose.j_WP)) p.ljets.push_back(jet);
        else if(jet.Pt() > p.c.jet1pt && jet.GetTaggerResult(DeepJet_Tight.j_Tagger) < mcCorr->GetJetTaggingCutValue(DeepJet_Tight.j_Tagger, DeepJet_Tight.j_WP) && jet.GetTaggerResult(DeepJet_CvsB_Tight.j_Tagger) < mcCorr->GetJetTaggingCutValue(DeepJet_CvsB_Tight.j_Tagger, DeepJet_CvsB_Tight.j_WP)) p.ajets.push_back(jet);
      }
      for(const auto& jet:p.ljets){
        if(jet.Pt() > p.c.jet0pt && jet.GetTaggerResult(DeepJet_Tight.j_Tagger) > mcCorr->GetJetTaggingCutValue(DeepJet_Tight.j_Tagger, DeepJet_Tight.j_WP)) p.bjets.push_back(jet);
        if(jet.Pt() > p.c.jet0pt && jet.GetTaggerResult(DeepJet_CvsB_Tight.j_Tagger) > mcCorr->GetJetTaggingCutValue(DeepJet_CvsB_Tight.j_Tagger, DeepJet_CvsB_Tight.j_WP) && jet.GetTaggerResult(DeepJet_CvsL_Tight.j_Tagger) > mcCorr->GetJetTaggingCutValue(DeepJet_CvsL_Tight.j_Tagger, DeepJet_CvsL_Tight.j_WP)) p.cjets.push_back(jet);
      }
      p.SetLJets(p.ljets);
      p.SetAJets(p.ajets);

      p.w.tagjetSF = 1.;
      p.doublemap["ltagSF_hup"] = 1.;
      p.doublemap["ltagSF_hdown"] = 1.;
      p.doublemap["ltagSF_lup"] = 1.;
      p.doublemap["ltagSF_ldown"] = 1.;
      p.doublemap["ltagSF_hup"+GetEra()] = 1.;
      p.doublemap["ltagSF_hdown"+GetEra()] = 1.;
      p.doublemap["ltagSF_lup"+GetEra()] = 1.;
      p.doublemap["ltagSF_ldown"+GetEra()] = 1.;
    }if(p.channel(2,2)=="jx"){
      p.SetJetPtCut(30,20);
      p.c.nljetmax=1;

      JetTagging::Parameters DeepJet_Tight = JetTagging::Parameters(JetTagging::DeepJet,JetTagging::Tight,JetTagging::incl,JetTagging::comb);
      JetTagging::Parameters DeepJet_CvsB_Tight = JetTagging::Parameters(JetTagging::DeepJet_CvsB,JetTagging::Tight,JetTagging::incl,JetTagging::wcharm);
      JetTagging::Parameters DeepJet_CvsL_Tight = JetTagging::Parameters(JetTagging::DeepJet_CvsL,JetTagging::Tight,JetTagging::incl,JetTagging::wcharm);
      p.ljets.clear();
      p.ajets.clear();

      for(const auto& jet:realjets){
        if(jet.Pt() > p.c.jet0pt) p.ljets.push_back(jet);
        else if(jet.Pt() > p.c.jet1pt) p.ajets.push_back(jet);
      }
      for(const auto& jet:p.ljets){
        if(jet.Pt() > p.c.jet0pt && jet.GetTaggerResult(DeepJet_Tight.j_Tagger) > mcCorr->GetJetTaggingCutValue(DeepJet_Tight.j_Tagger, DeepJet_Tight.j_WP)) p.bjets.push_back(jet);
        if(jet.Pt() > p.c.jet0pt && jet.GetTaggerResult(DeepJet_CvsB_Tight.j_Tagger) > mcCorr->GetJetTaggingCutValue(DeepJet_CvsB_Tight.j_Tagger, DeepJet_CvsB_Tight.j_WP) && jet.GetTaggerResult(DeepJet_CvsL_Tight.j_Tagger) > mcCorr->GetJetTaggingCutValue(DeepJet_CvsL_Tight.j_Tagger, DeepJet_CvsL_Tight.j_WP)) p.cjets.push_back(jet);
      }
      p.SetLJets(p.ljets);
      p.SetAJets(p.ajets);

      p.w.tagjetSF = 1.;
      p.doublemap["jtagSF_hup"] = 1.;
      p.doublemap["jtagSF_hdown"] = 1.;
      p.doublemap["jtagSF_lup"] = 1.;
      p.doublemap["jtagSF_ldown"] = 1.;
      p.doublemap["jtagSF_hup"+GetEra()] = 1.;
      p.doublemap["jtagSF_hdown"+GetEra()] = 1.;
      p.doublemap["jtagSF_lup"+GetEra()] = 1.;
      p.doublemap["jtagSF_ldown"+GetEra()] = 1.;
    }
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
