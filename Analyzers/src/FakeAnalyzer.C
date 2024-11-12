#include "FakeAnalyzer.h"

FakeAnalyzer::FakeAnalyzer(){
}
FakeAnalyzer::~FakeAnalyzer(){
}
void FakeAnalyzer::initializeAnalyzer(){
  SMPAnalyzerCore::initializeAnalyzer(); //setup zpt roc z0
}
void FakeAnalyzer::executeEvent(){
  if(!IsDATA||DataStream.Contains("DoubleMuon")){
    {Parameter p=MakeParameter("mm");
    p.SetMuons(MuonMomentumCorrection(SMPGetMuons("FakeAllID",8.0,2.4),0,0));
    p.option+=" triggermatching strictorder";
    executeEventWithParameter(p);}
  }
  if(!IsDATA||DataStream.Contains("DoubleEG")||DataStream.Contains("EGamma")){
    {Parameter p=MakeParameter("ee");
    p.SetElectrons(ElectronEnergyCorrection(SMPGetElectrons("FakeAllID",8.0,2.5),0,0));
    p.option+=" triggermatching strictorder";
    executeEventWithParameter(p);}
  }
}
int FakeAnalyzer::GetWP(const Lepton* lep){
  if(lep->LeptonFlavour()==Lepton::MUON){
    if(PassID(lep,"POGTightWithPFIsoVeryTight")&&PassID(lep,"POGMediumWithLooseTrkIso")) return 0;
    //if(PassID(lep,"POGMediumWithTightTrkIso")) return 0;
    else if(PassID(lep,"POGMediumWithLooseTrkIso")) return 1;
    else if(PassID(lep,"POGMediumWithTrkIsoVVL")) return 2;
    else if(PassID(lep,"FakeAllID")) return 3;
    else return 4;
  }else if(lep->LeptonFlavour()==Lepton::ELECTRON){
    const Electron* el=(const Electron*)lep;
    if(el->passTightID()) return 0;
    else if(el->passMediumID()) return 1;
    else if(el->passLooseID()) return 2;
    else if(el->passVetoID()) return 3;
    else if(PassID(el,"passMVAIso0")) return 4;
    else return 5;
  }
  return -1;
}
void FakeAnalyzer::FillHists(Parameter& p){
  TLorentzVector dilepton=(*p.lepton0)+(*p.lepton1);
  double dimass=dilepton.M();
  int l0wp=GetWP(p.lepton0);
  int l1wp=GetWP(p.lepton1);
  TString swp=Form("wp%d%d/",l0wp,l1wp);
  //if(dimass>52){
  if(p.leptons.at(0)->DeltaPhi(*p.leptons.at(1))>1.6){
    for(TString rb:{"","0bjet/","nbjet/"}){
      for(TString rz:{"","noZ/"}){
	for(TString rmet:{"","metcut/"}){
	  int nbjet=count_if(p.bjets.begin(),p.bjets.end(),[&p](auto& jet){return jet.Pt()>p.c.jetpt;});
	  if(rb=="0bjet/"&&nbjet>0) continue;
	  if(rb=="nbjet/"&&nbjet==0) continue;
	  if(rz=="noZ/"&&(dimass>76&&dimass<106)) continue;
	  if(rmet=="metcut/"&&pfMET_Type1_pt>50) continue;
	  TString region=rb+rz+rmet;
	  FillDileptonHists(p,swp+region);
	  if(p.channel=="mm"){
	    if(l0wp<2&&l1wp<2){	
	      FillDileptonHists(p,region);
	    }else if(l0wp==2&&l1wp==2){
	      map<TString,int> map_fakesuf={{"",0},{"_fakeTF_up",1},{"_fakeTF_down",-1}};
	      for(auto [fakesuf,sys]:map_fakesuf){
		double tf=GetFakeTF(p,"",sys);
		FillDileptonHists(p,region,"fake_",fakesuf,tf);
	      }
	      if(region.Contains("nbjet")){
		double tf=GetFakeTF(p,"nbjet",0);
		FillDileptonHists(p,region,"fake_","_fakeTF_nb",tf);
	      }
	    }
	  }else if(p.channel=="ee"){
	    if(l0wp<2&&l1wp<2){	
	      FillDileptonHists(p,region);
	    }else if(l0wp==5&&l1wp==5){
	      map<TString,int> map_fakesuf={{"",0},{"_fakeTF_up",1},{"_fakeTF_down",-1}};
	      for(auto [fakesuf,sys]:map_fakesuf){
		double tf=GetFakeTF(p,"",sys);
		FillDileptonHists(p,region,"fake_",fakesuf,tf);
	      }
	      if(region.Contains("nbjet")){
		double tf=GetFakeTF(p,"nbjet",0);
		FillDileptonHists(p,region,"fake_","_fakeTF_nb",tf);
	      }
	    }	    
	  }
	}
      }
    }
  } 
}
void FakeAnalyzer::FillDileptonHists(Parameter& p,TString region,TString hprefix,TString suffix,double w){
  double weight=p.weight*w;
  TString pre=p.prefix+region+hprefix+p.hprefix;
  TString suf=p.suffix+p.vsuffix+suffix;
  TLorentzVector dilepton=(*p.lepton0)+(*p.lepton1);
  
  FillHist(pre+"dimass"+suf,dilepton.M(),weight,150,0,150);
  FillHist(pre+"dimass_wide"+suf,dilepton.M(),weight,nmbin,mbins);
  FillHist(pre+"dipt"+suf,dilepton.Pt(),weight,100,0,100);
  FillHist(pre+"dirap"+suf,dilepton.Rapidity(),weight,60,-3,3);
  FillHist(pre+"met"+suf,pfMET_Type1_pt,weight,100,0,100);
  FillHist(pre+"drll"+suf,p.leptons.at(0)->DeltaR(*p.leptons.at(1)),weight,100,0,5);
  FillHist(pre+"dphill"+suf,p.leptons.at(0)->DeltaPhi(*p.leptons.at(1)),weight,80,-4,4);
  if(p.bjets.size()&&p.bjets.at(0).Pt()>p.c.jetpt){
    FillHist(pre+"b0pt"+suf,p.bjets.at(0).Pt(),weight,100,0,200);
    FillHist(pre+"b0eta"+suf,p.bjets.at(0).Eta(),weight,60,-3,3);
    FillHist(pre+"drbl"+suf,TMath::Min(p.bjets.at(0).DeltaR(*p.lepton0),p.bjets.at(0).DeltaR(*p.lepton1)),weight,60,-3,3);      
  }
  if(p.jets.size()&&p.jets.at(0).Pt()>p.c.jetpt){
    FillHist(pre+"j0pt"+suf,p.jets.at(0).Pt(),weight,100,0,200);
    FillHist(pre+"j0eta"+suf,p.jets.at(0).Eta(),weight,60,-3,3);
    FillHist(pre+"drjl"+suf,TMath::Min(p.jets.at(0).DeltaR(*p.lepton0),p.jets.at(0).DeltaR(*p.lepton1)),weight,60,-3,3);
  }      
  int nbjet=count_if(p.bjets.begin(),p.bjets.end(),[&p](auto& jet){return jet.Pt()>p.c.jetpt;});
  int njet=count_if(p.jets.begin(),p.jets.end(),[&p](auto& jet){return jet.Pt()>p.c.jetpt;});
  FillHist(pre+"bjets"+suf,nbjet,weight,10,0,10);
  FillHist(pre+"jets"+suf,njet,weight,10,0,10);
  
  for(int i=0,n=p.leptons.size();i<n;i++){
    double eta=fabs(p.leptons.at(i)->Eta());
    double pt=p.leptons.at(i)->Pt();
    double riso=p.leptons.at(i)->RelIso();
    double energy=p.leptons.at(i)->E();
    FillHist(pre+Form("l%detapt",i)+suf,eta,pt,weight,netabin,etabins,nptbin,ptbins);
    FillHist(pre+Form("l%detae",i)+suf,eta,energy,weight,netabin,etabins,nptbin,ptbins);
    FillHist(pre+Form("l%driso",i)+suf,riso,weight,100,0,1);
    if(p.leptons.at(i)->LeptonFlavour()==Lepton::Flavour::MUON){
      FillHist(pre+Form("l%drtkiso",i)+suf,((Muon*)p.leptons.at(i))->TrkIso()/pt,weight,100,0,1);
    }
    if(p.leptons.at(i)->LeptonFlavour()==Lepton::Flavour::ELECTRON){
      FillHist(pre+Form("l%dmvaiso",i)+suf,((Electron*)p.leptons.at(i))->MVAIso(),weight,1000,-1,1);
      FillHist(pre+Form("l%dmvanoiso",i)+suf,((Electron*)p.leptons.at(i))->MVANoIso(),weight,1000,-1,1);
    }
    FillHist(pre+"lpt"+suf,pt,weight,nptbin,ptbins);
    FillHist(pre+Form("l%dpt",i)+suf,pt,weight,nptbin,ptbins);
    FillHist(pre+"leta"+suf,eta,weight,netabin,etabins);
    FillHist(pre+Form("l%deta",i)+suf,eta,weight,netabin,etabins);
    double mt=sqrt(2*pt*pfMET_Type1_pt*(1-TMath::Cos(p.leptons.at(i)->Phi()-pfMET_Type1_phi)));
    FillHist(pre+Form("l%dmt",i)+suf,mt,weight,100,0,100);
  }
}
