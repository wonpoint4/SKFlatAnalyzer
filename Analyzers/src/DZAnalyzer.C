#include "DZAnalyzer.h"

DZAnalyzer::DZAnalyzer(){
}
DZAnalyzer::~DZAnalyzer(){
}
void DZAnalyzer::SetupDZSF(){
  cout<<"[DZAnalyzer::SetupDZSF] setup"<<endl;
  TString datapath=getenv("DATA_DIR");
  TString era=GetEra();
  TString erashort=GetEraShort();
  if(IsExists(datapath+"/"+era+"/ID/Muon/mm"+erashort+"_DZ.root")){
    TFile f(datapath+"/"+era+"/ID/Muon/mm"+erashort+"_DZ.root");
    fDZSF_muon=(TH2*)f.Get("sf");
    if(fDZSF_muon){
      cout<<"[DZAnalyzer::SetupDZSF] load muon "<<erashort<<endl;
      fDZSF_muon->SetDirectory(NULL);
    }
  }
  if(IsExists(datapath+"/"+era+"/ID/Electron/ee"+erashort+"_DZ.root")){
    TFile f(datapath+"/"+era+"/ID/Electron/ee"+erashort+"_DZ.root");
    fDZSF_electron=(TH2*)f.Get("sf");
    if(fDZSF_electron){
      cout<<"[DZAnalyzer::SetupDZSF] load electron "<<erashort<<endl;
      fDZSF_electron->SetDirectory(NULL);
    }
  }
  if(IsExists(datapath+"/"+era+"/ID/Muon/mm"+erashort+"_DZ_DZ.root")){
    TFile f(datapath+"/"+era+"/ID/Muon/mm"+erashort+"_DZ_DZ.root");
    fDZSF_muon_DZ=(TH1*)f.Get("sf");
    if(fDZSF_muon_DZ){
      cout<<"[DZAnalyzer::SetupDZSF] load muon "<<erashort<<endl;
      fDZSF_muon_DZ->SetDirectory(NULL);
    }
  }
  if(IsExists(datapath+"/"+era+"/ID/Electron/ee"+erashort+"_DZ_DZ.root")){
    TFile f(datapath+"/"+era+"/ID/Electron/ee"+erashort+"_DZ_DZ.root");
    fDZSF_electron_DZ=(TH1*)f.Get("sf");
    if(fDZSF_electron_DZ){
      cout<<"[DZAnalyzer::SetupDZSF] load electron "<<erashort<<endl;
      fDZSF_electron_DZ->SetDirectory(NULL);
    }
  }
}
double DZAnalyzer::GetDZSF(Lepton* lep){
  double sf=1.;
  if(IsDATA) return sf;
  double eta=lep->Eta();
  double pt=lep->Pt();
  if(eta>=2.5) eta=2.49;
  if(eta<-2.5) eta=-2.5;
  if(pt>=200) pt=199;
  if(lep->InheritsFrom("Electron")){
    if(fDZSF_electron){
      int bin=fDZSF_electron->FindBin(eta,pt);
      sf*=fDZSF_electron->GetBinContent(bin);
    }
  }
  if(lep->InheritsFrom("Muon")){
    if(fDZSF_muon){
      int bin=fDZSF_muon->FindBin(eta,pt);
      sf*=fDZSF_muon->GetBinContent(bin);
    }
  }
  return sf;
}
double DZAnalyzer::GetDZSF_DZ(Lepton* lep1,Lepton* lep2){
  double sf=1.;
  if(IsDATA) return sf;
  double dz=fabs(lep1->dZ()-lep2->dZ());
  if(lep1->InheritsFrom("Electron")){
    if(fDZSF_electron_DZ){
      int bin=fDZSF_electron_DZ->FindBin(dz);
      sf*=fDZSF_electron_DZ->GetBinContent(bin);
    }
  }
  if(lep1->InheritsFrom("Muon")){
    if(fDZSF_muon_DZ){
      int bin=fDZSF_muon_DZ->FindBin(dz);
      sf*=fDZSF_muon_DZ->GetBinContent(bin);
    }
  }
  return sf;
}
void DZAnalyzer::initializeAnalyzer(){
  SMPAnalyzerCore::initializeAnalyzer(); //setup zpt roc z0
  SetupDZSF();
  fChain->SetBranchStatus("pfMET_*",false);
  fChain->SetBranchStatus("pfMET_Type1_pt",true);
  fChain->SetBranchStatus("jet_*",false);
  fChain->SetBranchStatus("fatjet_*",false);
  fChain->SetBranchStatus("photon_*",false);
}
void DZAnalyzer::executeEvent(){
  //////// nominal channels //////////
  if(!IsDATA||DataStream.Contains("DoubleMuon")){
    Parameter p=MakeParameter("mm");
    p.triggers={"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v"};
    executeEventWithParameter(p);
    if(fabs(vertex_Z)<7){
      p.prefix+="z0/";
      executeEventWithParameter(p);
    }
  }
  if(!IsDATA||DataStream.Contains("DoubleEG")||DataStream.Contains("EGamma")){
    Parameter p=MakeParameter("ee");
    p.triggers={"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"};
    executeEventWithParameter(p);
    if(fabs(vertex_Z)<7){
      p.prefix+="z0/";
      executeEventWithParameter(p);
    }
  }
}
void DZAnalyzer::EvalWeights(Parameter& p){
  p.weightmap[""]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF;
  p.weightmap["_dzsf"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*GetDZSF(p.lepton0)*GetDZSF(p.lepton1);
  p.weightmap["_dzsf_dz"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.CFSF*GetDZSF_DZ(p.lepton0,p.lepton1);
}
void DZAnalyzer::FillHists(Parameter& p){
  TLorentzVector dilepton=(*p.lepton0)+(*p.lepton1);
  double dimass=dilepton.M();
  if(dimass>=76&&dimass<106){
    FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"m76to106",p.weightmap[""]);
    FillHistsDZ(p,"_den");
    FillHistsDZ(p,"_test_den");
    if(p.channel=="ee"){
      if(((Electron*)p.lepton0)->PassFilter("hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter")
	 && ((Electron*)p.lepton1)->PassFilter("hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter")){
	FillHistsDZ(p,"_num");
	FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"passDZ",p.weightmap[""]);
      }else{
	FillHistsDZ(p,"_fail");
      }	
    }else if(p.channel=="mm"){
      TString dzfilter;
      if(DataYear==2016) dzfilter="hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4DzFiltered0p2";
      else if(DataYear==2017) dzfilter="hltDiMuon178RelTrkIsoFiltered0p4DzFiltered0p2";
      else if(DataYear==2018) dzfilter="hltDiMuon178RelTrkIsoFiltered0p4DzFiltered0p2";
	      
      if(((Muon*)p.lepton0)->PassFilter(dzfilter) && ((Muon*)p.lepton1)->PassFilter(dzfilter)){
	FillHistsDZ(p,"_num");
	FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"passDZ",p.weightmap[""]);
      }else{
	FillHistsDZ(p,"_fail");
      }
    }
    if(_event.PassTrigger(Replace(p.triggers[0],"_v$","_DZ_v"))){
      FillHistsDZ(p,"_test_num");
    }else{
      FillHistsDZ(p,"_test_fail");
    }
  }
}

void DZAnalyzer::FillHistsDZ(Parameter& p,TString suffix){
  TString pre=p.prefix+p.hprefix;
  for(const auto& [wname,w]:p.weightmap){
    TString suf=p.suffix+suffix+wname;
    
    //for leptons
    for(int i=0;i<(int)p.leptons.size();i++){
      double pt=p.leptons.at(i)->Pt();
      double eta=p.leptons.at(i)->Eta();
      double dz=p.leptons.at(i)->dZ();
      double dze=p.leptons.at(i)->dZerr();
      TString charge=p.leptons.at(i)->Charge()>0?"p":"m";
      FillHist(Form("%sl%dpt%s",pre.Data(),i,suf.Data()),pt,w,fineptbinnum,fineptbin);
      FillHist(Form("%sl%deta%s",pre.Data(),i,suf.Data()),eta,w,fineetabinnum,fineetabin);
      FillHist(Form("%sl%ddz%s",pre.Data(),i,suf.Data()),fabs(dz),w,1000,0,1);
      FillHist(Form("%sl%ddze%s",pre.Data(),i,suf.Data()),dze,w,1000,0,1);
      FillHist(Form("%sl%ipz%s",pre.Data(),i,suf.Data()),fabs(dz)/dze,w,1000,0,100);
      
      FillHist(Form("%slpt%s",pre.Data(),suf.Data()),pt,w,fineptbinnum,fineptbin);
      FillHist(Form("%sleta%s",pre.Data(),suf.Data()),eta,w,fineetabinnum,fineetabin);
      FillHist(Form("%sldz%s",pre.Data(),suf.Data()),fabs(dz),w,1000,0,1);
      FillHist(Form("%sldze%s",pre.Data(),suf.Data()),dze,w,1000,0,1);
      FillHist(Form("%slipz%s",pre.Data(),suf.Data()),fabs(dz)/dze,w,1000,0,100);
      FillHist(Form("%setapt%s",pre.Data(),suf.Data()),p.leptons.at(i)->Eta(),p.leptons.at(i)->Pt(),w,etabinnum,etabin,ptbinnum,ptbin);
      FillHist(Form("%setaptfine%s",pre.Data(),suf.Data()),p.leptons.at(i)->Eta(),p.leptons.at(i)->Pt(),w,fineetabinnum,fineetabin,fineptbinnum,fineptbin);

    }

    double lldzmin=999;
    if(p.channel=="ee"){
      vector<Electron> electrons=GetAllElectrons();
      vector<Electron> leg1s,leg2s;
      for(int i=0,n=electrons.size();i<n;i++){
	if(electrons.at(i).PassFilter("hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter"))
	  leg1s.push_back(electrons.at(i));
	if(electrons.at(i).PassFilter("hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter"))
	  leg2s.push_back(electrons.at(i));
      }
      for(auto leg1:leg1s){
	for(auto leg2:leg2s){
	  if(leg1.DeltaR(leg2)<0.001) continue;
	  double temp_dz=leg2.dZ()-leg1.dZ();
	  if(fabs(lldzmin)>fabs(temp_dz)){
	    lldzmin=temp_dz;
	  }
	}
      }
    }else if(p.channel=="mm"){
      vector<Muon> muons=GetAllMuons();
      vector<Muon> leg1s,leg2s;
      for(int i=0,n=muons.size();i<n;i++){
	if(muons.at(i).PassFilter("hltL3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17")&&muons.at(i).PassFilter("hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4"))
	  leg1s.push_back(muons.at(i));
	if(muons.at(i).PassFilter("hltL3pfL1sDoubleMu114ORDoubleMu125L1f0L2pf0L3PreFiltered8")&&muons.at(i).PassFilter("hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4"))
	  leg2s.push_back(muons.at(i));
      }
      for(auto leg1:leg1s){
	for(auto leg2:leg2s){
	  if(leg1.DeltaR(leg2)<0.001) continue;
	  double temp_dz=leg2.dZ()-leg1.dZ();
	  if(fabs(lldzmin)>fabs(temp_dz)){
	    lldzmin=temp_dz;
	  }
	}
      }
    }
    FillHist(Form("%slldzmin%s",pre.Data(),suf.Data()),fabs(lldzmin),w,1000,0,1);
    FillHist(Form("%slldz%s",pre.Data(),suf.Data()),fabs(p.lepton1->dZ()-p.lepton0->dZ()),w,1000,0,1);
    Lepton *lepton0,*lepton1;
    if(p.lepton0->Charge()>0&&p.lepton1->Charge()<0){
      lepton0=p.lepton1;
      lepton1=p.lepton0;
    }else if(p.lepton0->Charge()<0&&p.lepton1->Charge()>0){
      lepton0=p.lepton0;
      lepton1=p.lepton1;
    }else if(p.lepton0->Phi()>p.lepton1->Phi()){
      lepton0=p.lepton1;
      lepton1=p.lepton0;
    }else{
      lepton0=p.lepton0;
      lepton1=p.lepton1;      
    }
    FillHist(Form("%sepep%s",pre.Data(),suf.Data()),lepton0->Eta(),lepton0->Pt(),lepton1->Eta(),lepton1->Pt(),w,etabinnum,etabin,ptbinnum,ptbin,etabinnum,etabin,ptbinnum,ptbin);
    FillHist(Form("%sepepfine%s",pre.Data(),suf.Data()),lepton0->Eta(),lepton0->Pt(),lepton1->Eta(),lepton1->Pt(),w,fineetabinnum,fineetabin,fineptbinnum,fineptbin,fineetabinnum,fineetabin,fineptbinnum,fineptbin);

    FillHist(Form("%sptpt%s",pre.Data(),suf.Data()),lepton0->Pt(),lepton1->Pt(),w,fineptbinnum,fineptbin,fineptbinnum,fineptbin);
    FillHist(Form("%setaeta%s",pre.Data(),suf.Data()),lepton0->Eta(),lepton1->Eta(),w,fineetabinnum,fineetabin,fineetabinnum,fineetabin);
      
    
    TLorentzVector dilepton=(*p.lepton0)+(*p.lepton1);
    double dimass=dilepton.M();
    double dipt=dilepton.Pt();
    double dirap=dilepton.Rapidity();
    FillHist(pre+"dimass"+suf,dimass,w,196,52,150);
    FillHist(pre+"dipt"+suf,dipt,w,400,0,400);
    FillHist(pre+"dirap"+suf,dirap,w,120,-3,3);
    FillHist(pre+"z0"+suf,vertex_Z,w,100,-20,20);
    FillHist(pre+"nPV"+suf,nPV,w,100,0,100);

    if(wname!="") continue;

  }
}

bool DZAnalyzer::PassSelection(Parameter& p){
  if(p.triggers.size()==2){
    if(GetEra()=="2017"&&p.triggers[0]=="HLT_Ele27_WPTight_Gsf_v"&&p.triggers[1]=="HLT_Ele32_WPTight_Gsf_v"){
      //if(!_event.PassTrigger("HLT_Ele27_WPTight_Gsf_v")) p.c.lepton0pt=35;
      if(p.lepton0&&p.lepton0->Pt()<35){
	p.triggers={"HLT_Ele27_WPTight_Gsf_v"};
	if(!_event.PassTrigger(p.triggers)) return false;
	p.k.triggerSF={"Ele27_MediumID"};
	if(!IsDATA) p.w.lumiweight*=31.72/41.54;
      }
    }
    if(GetEra()=="2018"&&p.triggers[0]=="HLT_Ele28_WPTight_Gsf_v"&&p.triggers[1]=="HLT_Ele32_WPTight_Gsf_v"){
      //if(!_event.PassTrigger("HLT_Ele28_WPTight_Gsf_v")) p.c.lepton0pt=35;
      if(p.lepton0&&p.lepton0->Pt()<35){
	p.triggers={"HLT_Ele28_WPTight_Gsf_v"};
	if(!_event.PassTrigger(p.triggers)) return false;
	p.k.triggerSF={"Ele28_MediumID"};
	if(!IsDATA) p.w.lumiweight*=23687.253/59827.879;
      }
    }
    if(GetEra()=="2017"&&p.triggers[0]=="HLT_IsoMu24_v"&&p.triggers[1]=="HLT_IsoMu27_v"){
      //if(!_event.PassTrigger("HLT_IsoMu24_v")) p.c.lepton0pt=30;
      if(p.lepton0&&p.lepton0->Pt()<30){
	p.triggers={"HLT_IsoMu24_v"};
	if(!_event.PassTrigger(p.triggers)) return false;
	p.k.triggerSF={"IsoMu24_MediumID_trkIsoLoose"};
	if(!IsDATA) p.w.lumiweight*=37997.005/41477.878;
      }
    }
  }
  return SMPAnalyzerCore::PassSelection(p);
}
