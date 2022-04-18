#include "FakeAnalyzer.h"

FakeAnalyzer::FakeAnalyzer(){
}
FakeAnalyzer::~FakeAnalyzer(){
}
void FakeAnalyzer::initializeAnalyzer(){
  SMPAnalyzerCore::initializeAnalyzer(); //setup zpt roc z0
  fChain->SetBranchStatus("pfMET_*",false);
  fChain->SetBranchStatus("pfMET_Type1_pt",true);
  fChain->SetBranchStatus("jet_*",false);
  fChain->SetBranchStatus("fatjet_*",false);
  fChain->SetBranchStatus("photon_*",false);
}
void FakeAnalyzer::executeEvent(){
  if(!IsDATA||DataStream.Contains("DoubleMuon")){
    executeEventWithParameter(MakeParameter("mM"));
    executeEventWithParameter(MakeParameter("MM"));
  }
  if(!IsDATA||DataStream.Contains("DoubleEG")||DataStream.Contains("EGamma")){
    executeEventWithParameter(MakeParameter("eE"));
    executeEventWithParameter(MakeParameter("EE"));
  }
}

void FakeAnalyzer::FillHists(Parameter& p){
  for(auto [suf,weight]:p.weightmap){
    TLorentzVector dilepton=(*p.lepton0)+(*p.lepton1);
    double dimass=dilepton.M();
    if(dimass>52){
      if(p.channel.Contains(TRegexp("[eE]E"))){
	for(int i=0,n=p.aelectrons.size();i<n;i++){
	  double eta=fabs(p.aelectrons.at(i).Eta());
	  double pt=p.aelectrons.at(i).Pt();
	  double riso=p.aelectrons.at(i).RelIso();
	  double energy=p.aelectrons.at(i).E();
	  FillHist(p.prefix+p.hprefix+Form("al%detapt",i)+p.suffix+suf,eta,pt,weight,netabin,etabins,nptbin,ptbins);
	  FillHist(p.prefix+p.hprefix+Form("al%detajpt",i)+p.suffix+suf,eta,pt*(1+riso),weight,netabin,etabins,nptbin,ptbins);
	  FillHist(p.prefix+p.hprefix+Form("al%detampt",i)+p.suffix+suf,eta,pt*(1+riso*TMath::Max(0.,TMath::Min(1.,(pt-30)/30))),weight,netabin,etabins,nptbin,ptbins);
	  if(pt>15) FillHist(p.prefix+p.hprefix+Form("al%detae",i)+p.suffix+suf,eta,energy,weight,netabin,etabins,nptbin,ptbins);
	}
	for(int i=0,n=p.electrons.size();i<n;i++){
	  double eta=fabs(p.electrons.at(i).Eta());
	  double pt=p.electrons.at(i).Pt();
	  double riso=p.electrons.at(i).RelIso();
	  double energy=p.electrons.at(i).E();
	  FillHist(p.prefix+p.hprefix+Form("l%detapt",i)+p.suffix+suf,eta,pt,weight,netabin,etabins,nptbin,ptbins);
	  FillHist(p.prefix+p.hprefix+Form("l%detajpt",i)+p.suffix+suf,eta,pt*(1+riso),weight,netabin,etabins,nptbin,ptbins);
	  FillHist(p.prefix+p.hprefix+Form("l%detampt",i)+p.suffix+suf,eta,pt*(1+riso*TMath::Max(0.,TMath::Min(1.,(pt-30)/30))),weight,netabin,etabins,nptbin,ptbins);
	  if(pt>15) FillHist(p.prefix+p.hprefix+Form("l%detae",i)+p.suffix+suf,eta,energy,weight,netabin,etabins,nptbin,ptbins);
	}
      }else if(p.channel.Contains(TRegexp("[mM]M"))){
	for(int i=0,n=p.amuons.size();i<n;i++){
	  double eta=fabs(p.amuons.at(i).Eta());
	  double pt=p.amuons.at(i).Pt();
	  double riso=p.amuons.at(i).RelIso();
	  FillHist(p.prefix+p.hprefix+Form("al%detapt",i)+p.suffix+suf,eta,pt,weight,netabin,etabins,nptbin,ptbins);
	  FillHist(p.prefix+p.hprefix+Form("al%detajpt",i)+p.suffix+suf,eta,pt*(1+riso),weight,netabin,etabins,nptbin,ptbins);
	  FillHist(p.prefix+p.hprefix+Form("al%detampt",i)+p.suffix+suf,eta,pt*(1+riso*TMath::Max(0.,TMath::Min(1.,(pt-30)/30))),weight,netabin,etabins,nptbin,ptbins);
	}
	for(int i=0,n=p.muons.size();i<n;i++){
	  double eta=fabs(p.muons.at(i).Eta());
	  double pt=p.muons.at(i).Pt();
	  double riso=p.muons.at(i).RelIso();
	  FillHist(p.prefix+p.hprefix+Form("l%detapt",i)+p.suffix+suf,eta,pt,weight,netabin,etabins,nptbin,ptbins);
	  FillHist(p.prefix+p.hprefix+Form("l%detajpt",i)+p.suffix+suf,eta,pt*(1+riso),weight,netabin,etabins,nptbin,ptbins);
	  FillHist(p.prefix+p.hprefix+Form("l%detampt",i)+p.suffix+suf,eta,pt*(1+riso*TMath::Max(0.,TMath::Min(1.,(pt-30)/30))),weight,netabin,etabins,nptbin,ptbins);
	}
      }
      if(p.channel=="EE"){
	TString prefix=p.prefix;
	prefix.ReplaceAll("EE","ee");
	for(int i=0,n=p.aelectrons.size();i<n;i++){
	  for(int j=i+1;j<n;j++){
	    if(p.aelectrons.at(i).Pt()>p.c.lepton0pt&&p.aelectrons.at(j).Pt()>p.c.lepton1pt&&(p.aelectrons.at(i)+p.aelectrons.at(j)).M()>52){
	      TString hprefix=p.aelectrons.at(i).Charge()*p.aelectrons.at(j).Charge()>0?"ss_":"";
	      double fakerate=GetFakeRate(&p.aelectrons.at(i))*GetFakeRate(&p.aelectrons.at(j));
	      for(int k=j+1;k<n;k++) fakerate*=1+GetFakeRate(&p.aelectrons.at(k));
	      FillHist(prefix+hprefix+"fake_mass"+p.suffix+suf,(p.aelectrons.at(i)+p.aelectrons.at(j)).M(),weight*fakerate,98,52,150);
	      FillHist(prefix+hprefix+"fake_lpt"+p.suffix+suf,p.aelectrons.at(i).Pt(),weight*fakerate,100,0,100);
	      FillHist(prefix+hprefix+"fake_lpt"+p.suffix+suf,p.aelectrons.at(j).Pt(),weight*fakerate,100,0,100);
	      FillHist(prefix+hprefix+"fake_l0pt"+p.suffix+suf,p.aelectrons.at(i).Pt(),weight*fakerate,100,0,100);
	      FillHist(prefix+hprefix+"fake_l1pt"+p.suffix+suf,p.aelectrons.at(j).Pt(),weight*fakerate,100,0,100);
	      FillHist(prefix+hprefix+"fake_leta"+p.suffix+suf,p.aelectrons.at(i).Eta(),weight*fakerate,100,-5,5);
	      FillHist(prefix+hprefix+"fake_leta"+p.suffix+suf,p.aelectrons.at(j).Eta(),weight*fakerate,100,-5,5);
	      FillHist(prefix+hprefix+"fake_l0eta"+p.suffix+suf,p.aelectrons.at(i).Eta(),weight*fakerate,100,-5,5);
	      FillHist(prefix+hprefix+"fake_l1eta"+p.suffix+suf,p.aelectrons.at(j).Eta(),weight*fakerate,100,-5,5);
	    }
	  }
	}
      }else if(p.channel=="MM"){
	TString prefix=p.prefix;
	prefix.ReplaceAll("MM","mm");
	for(int i=0,n=p.amuons.size();i<n;i++){
	  for(int j=i+1;j<n;j++){
	    if(p.amuons.at(i).Pt()>p.c.lepton0pt&&p.amuons.at(j).Pt()>p.c.lepton1pt&&(p.amuons.at(i)+p.amuons.at(j)).M()>52){
	      double fakerate=GetFakeRate(&p.amuons.at(i))*GetFakeRate(&p.amuons.at(j));
	      TString hprefix=p.amuons.at(i).Charge()*p.amuons.at(j).Charge()>0?"ss_":"";
	      for(int k=j+1;k<n;k++) fakerate*=1+GetFakeRate(&p.amuons.at(k));
	      FillHist(prefix+hprefix+"fake_mass"+p.suffix+suf,(p.amuons.at(i)+p.amuons.at(j)).M(),weight*fakerate,98,52,150);
	      FillHist(prefix+hprefix+"fake_lpt"+p.suffix+suf,p.amuons.at(i).Pt(),weight*fakerate,100,0,100);
	      FillHist(prefix+hprefix+"fake_lpt"+p.suffix+suf,p.amuons.at(j).Pt(),weight*fakerate,100,0,100);
	      FillHist(prefix+hprefix+"fake_l0pt"+p.suffix+suf,p.amuons.at(i).Pt(),weight*fakerate,100,0,100);
	      FillHist(prefix+hprefix+"fake_l1pt"+p.suffix+suf,p.amuons.at(j).Pt(),weight*fakerate,100,0,100);
	      FillHist(prefix+hprefix+"fake_leta"+p.suffix+suf,p.amuons.at(i).Eta(),weight*fakerate,100,-5,5);
	      FillHist(prefix+hprefix+"fake_leta"+p.suffix+suf,p.amuons.at(j).Eta(),weight*fakerate,100,-5,5);
	      FillHist(prefix+hprefix+"fake_l0eta"+p.suffix+suf,p.amuons.at(i).Eta(),weight*fakerate,100,-5,5);
	      FillHist(prefix+hprefix+"fake_l1eta"+p.suffix+suf,p.amuons.at(j).Eta(),weight*fakerate,100,-5,5);
	    }
	  }
	}
      }
      if(dimass<76||dimass>106){
	if(p.channel.Contains(TRegexp("[eE]E"))){
	  for(int i=0,n=p.aelectrons.size();i<n;i++){
	    double eta=fabs(p.aelectrons.at(i).Eta());
	    double pt=p.aelectrons.at(i).Pt();
	    double riso=p.aelectrons.at(i).RelIso();
	    double energy=p.aelectrons.at(i).E();
	    FillHist(p.prefix+p.hprefix+Form("al%detapt_noZ",i)+p.suffix+suf,eta,pt,weight,netabin,etabins,nptbin,ptbins);
	    FillHist(p.prefix+p.hprefix+Form("al%detajpt_noZ",i)+p.suffix+suf,eta,pt*(1+riso),weight,netabin,etabins,nptbin,ptbins);
	    FillHist(p.prefix+p.hprefix+Form("al%detampt_noZ",i)+p.suffix+suf,eta,pt*(1+riso*TMath::Max(0.,TMath::Min(1.,(pt-30)/30))),weight,netabin,etabins,nptbin,ptbins);
	    if(pt>15) FillHist(p.prefix+p.hprefix+Form("al%detae_noZ",i)+p.suffix+suf,eta,energy,weight,netabin,etabins,nptbin,ptbins);
	  }
	  for(int i=0,n=p.electrons.size();i<n;i++){
	    double eta=fabs(p.electrons.at(i).Eta());
	    double pt=p.electrons.at(i).Pt();
	    double riso=p.electrons.at(i).RelIso();
	    double energy=p.electrons.at(i).E();
	    FillHist(p.prefix+p.hprefix+Form("l%detapt_noZ",i)+p.suffix+suf,eta,pt,weight,netabin,etabins,nptbin,ptbins);
	    FillHist(p.prefix+p.hprefix+Form("l%detajpt_noZ",i)+p.suffix+suf,eta,pt*(1+riso),weight,netabin,etabins,nptbin,ptbins);
	    FillHist(p.prefix+p.hprefix+Form("l%detampt_noZ",i)+p.suffix+suf,eta,pt*(1+riso*TMath::Max(0.,TMath::Min(1.,(pt-30)/30))),weight,netabin,etabins,nptbin,ptbins);
	    if(pt>15) FillHist(p.prefix+p.hprefix+Form("l%detae_noZ",i)+p.suffix+suf,eta,energy,weight,netabin,etabins,nptbin,ptbins);
	  }
	}else if(p.channel.Contains(TRegexp("[mM]M"))){
	  for(int i=0,n=p.amuons.size();i<n;i++){
	    double eta=fabs(p.amuons.at(i).Eta());
	    double pt=p.amuons.at(i).Pt();
	    double riso=p.amuons.at(i).RelIso();
	    FillHist(p.prefix+p.hprefix+Form("al%detapt_noZ",i)+p.suffix+suf,eta,pt,weight,netabin,etabins,nptbin,ptbins);
	    FillHist(p.prefix+p.hprefix+Form("al%detajpt_noZ",i)+p.suffix+suf,eta,pt*(1+riso),weight,netabin,etabins,nptbin,ptbins);
	    FillHist(p.prefix+p.hprefix+Form("al%detampt_noZ",i)+p.suffix+suf,eta,pt*(1+riso*TMath::Max(0.,TMath::Min(1.,(pt-30)/30))),weight,netabin,etabins,nptbin,ptbins);
	  }
	  for(int i=0,n=p.muons.size();i<n;i++){
	    double eta=fabs(p.muons.at(i).Eta());
	    double pt=p.muons.at(i).Pt();
	    double riso=p.muons.at(i).RelIso();
	    FillHist(p.prefix+p.hprefix+Form("l%detapt_noZ",i)+p.suffix+suf,eta,pt,weight,netabin,etabins,nptbin,ptbins);
	    FillHist(p.prefix+p.hprefix+Form("l%detajpt_noZ",i)+p.suffix+suf,eta,pt*(1+riso),weight,netabin,etabins,nptbin,ptbins);
	    FillHist(p.prefix+p.hprefix+Form("l%detampt_noZ",i)+p.suffix+suf,eta,pt*(1+riso*TMath::Max(0.,TMath::Min(1.,(pt-30)/30))),weight,netabin,etabins,nptbin,ptbins);
	  }
	}	
      }
    }
  }
}
