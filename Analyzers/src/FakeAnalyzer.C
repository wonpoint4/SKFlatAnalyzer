#include "FakeAnalyzer.h"

FakeAnalyzer::FakeAnalyzer(){
}
FakeAnalyzer::~FakeAnalyzer(){
}
void FakeAnalyzer::initializeAnalyzer(){
  SMPAnalyzerCore::initializeAnalyzer(); //setup zpt roc z0
  vector<JetTagging::Parameters> jtps={JetTagging::Parameters(JetTagging::DeepJet,JetTagging::Tight,JetTagging::incl,JetTagging::comb)};
  mcCorr->SetJetTaggingParameters(jtps);
}
//void FakeAnalyzer::executeEvent(){
//  if(GetSkimName()=="Dilepton") executeDileptonEvent();
//  else if(GetSkimName()=="HNFake") executeLeptonEvent();
//}
void FakeAnalyzer::executeEvent(){
  if(!IsDATA||DataStream.Contains("DoubleMuon")){
    {Parameter p=MakeParameter("mm");
    p.SetAMuons(MuonMomentumCorrection(SMPGetMuons("POGMediumWithAntiLooseTrkIso",0.0,2.4),0,0));
    p.option+=" triggermatching strictorder";
    executeEventWithParameter(p);}
      
    executeEventWithParameter(MakeParameter("mM"));
    executeEventWithParameter(MakeParameter("Mm"));
    executeEventWithParameter(MakeParameter("MM"));

    {Parameter p=MakeParameter("mM");
    p.SetAMuons(ToConePt(p.amuons));
    p.prefix+="cpt/";
    executeEventWithParameter(p);}

    {Parameter p=MakeParameter("Mm");
    p.SetAMuons(ToConePt(p.amuons));
    p.prefix+="cpt/";
    executeEventWithParameter(p);}

    {Parameter p=MakeParameter("MM");
    p.SetAMuons(ToConePt(p.amuons));
    p.prefix+="cpt/";
    executeEventWithParameter(p);}

    {Parameter p=MakeParameter("mM");
    p.SetAMuons(ToModifiedPt(p.amuons));
    p.prefix+="mpt/";
    executeEventWithParameter(p);}

    {Parameter p=MakeParameter("Mm");
    p.SetAMuons(ToModifiedPt(p.amuons));
    p.prefix+="mpt/";
    executeEventWithParameter(p);}

    {Parameter p=MakeParameter("MM");
    p.SetAMuons(ToModifiedPt(p.amuons));
    p.prefix+="mpt/";
    executeEventWithParameter(p);}

  }
  if(!IsDATA||DataStream.Contains("DoubleEG")||DataStream.Contains("EGamma")){
    {Parameter p=MakeParameter("ee");
    p.SetAElectrons(SMPGetElectrons("passAntiLooseID",0.0,2.5));
    p.option+=" triggermatching strictorder";
    UseSelectiveCharge(p);
    executeEventWithParameter(p);}
    {Parameter p=MakeParameter("eE");
    UseSelectiveCharge(p);
    executeEventWithParameter(p);}
    {Parameter p=MakeParameter("Ee");
    UseSelectiveCharge(p);
    executeEventWithParameter(p);}
    {Parameter p=MakeParameter("EE");
    UseSelectiveCharge(p);
    executeEventWithParameter(p);}

    {Parameter p=MakeParameter("eE");
    p.SetAElectrons(ToConePt(p.aelectrons));
    p.prefix+="cpt/";
    UseSelectiveCharge(p);
    executeEventWithParameter(p);}

    {Parameter p=MakeParameter("Ee");
    p.SetAElectrons(ToConePt(p.aelectrons));
    p.prefix+="cpt/";
    UseSelectiveCharge(p);
    executeEventWithParameter(p);}

    {Parameter p=MakeParameter("EE");
    p.SetAElectrons(ToConePt(p.aelectrons));
    p.prefix+="cpt/";
    UseSelectiveCharge(p);
    executeEventWithParameter(p);}

    {Parameter p=MakeParameter("eE");
    p.SetAElectrons(ToModifiedPt(p.aelectrons));
    p.prefix+="mpt/";
    UseSelectiveCharge(p);
    executeEventWithParameter(p);}

    {Parameter p=MakeParameter("Ee");
    p.SetAElectrons(ToModifiedPt(p.aelectrons));
    p.prefix+="mpt/";
    UseSelectiveCharge(p);
    executeEventWithParameter(p);}

    {Parameter p=MakeParameter("EE");
    p.SetAElectrons(ToModifiedPt(p.aelectrons));
    p.prefix+="mpt/";
    UseSelectiveCharge(p);
    executeEventWithParameter(p);}
  }
}

void FakeAnalyzer::FillHists(Parameter& p){
  int n_bjet=0;
  std::vector<Jet> jets=GetJets("tightLepVeto",40,2.4);
  std::sort(jets.begin(),jets.end(),PtComparing);
  JetTagging::Parameters jtp = JetTagging::Parameters(JetTagging::DeepJet,JetTagging::Tight,JetTagging::incl,JetTagging::comb);
  for(const auto& jet:jets)
    if(jet.GetTaggerResult(jtp.j_Tagger) > mcCorr->GetJetTaggingCutValue(jtp.j_Tagger, jtp.j_WP))
      n_bjet++;

  TLorentzVector dilepton=(*p.lepton0)+(*p.lepton1);
  double dimass=dilepton.M();
  if(dimass>52){
    if(n_bjet) FillFakeHists(p,"nbjet/");
    else FillFakeHists(p,"0bjet/");
    FillFakeHists(p,"");
    if(dimass<76||dimass>106){
      if(n_bjet) FillFakeHists(p,"nbjet/noZ/");
      else FillFakeHists(p,"0bjet/noZ/");
      FillFakeHists(p,"noZ/");
    }
  }    
}
void FakeAnalyzer::FillFakeHists(Parameter& p,TString region){
  for(auto [suf,weight]:p.weightmap){
    TLorentzVector dilepton=(*p.lepton0)+(*p.lepton1);
    FillHist(p.prefix+region+p.hprefix+"dimass"+p.suffix+suf,dilepton.M(),weight,98,52,150);
    FillHist(p.prefix+region+p.hprefix+"dimass_wide"+p.suffix+suf,dilepton.M(),weight,nmbin,mbins);
    FillHist(p.prefix+region+p.hprefix+"dipt"+p.suffix+suf,dilepton.Pt(),weight,100,0,100);
    FillHist(p.prefix+region+p.hprefix+"dirap"+p.suffix+suf,dilepton.Rapidity(),weight,60,-3,3);
    FillHist(p.prefix+region+p.hprefix+"lpt"+p.suffix+suf,p.lepton0->Pt(),weight,nptbin,ptbins);
    FillHist(p.prefix+region+p.hprefix+"lpt"+p.suffix+suf,p.lepton1->Pt(),weight,nptbin,ptbins);
    FillHist(p.prefix+region+p.hprefix+"l0pt"+p.suffix+suf,p.lepton0->Pt(),weight,nptbin,ptbins);
    FillHist(p.prefix+region+p.hprefix+"l1pt"+p.suffix+suf,p.lepton1->Pt(),weight,nptbin,ptbins);
    FillHist(p.prefix+region+p.hprefix+"leta"+p.suffix+suf,p.lepton0->Eta(),weight,netabin,etabins);
    FillHist(p.prefix+region+p.hprefix+"leta"+p.suffix+suf,p.lepton1->Eta(),weight,netabin,etabins);
    FillHist(p.prefix+region+p.hprefix+"l0eta"+p.suffix+suf,p.lepton0->Eta(),weight,netabin,etabins);
    FillHist(p.prefix+region+p.hprefix+"l1eta"+p.suffix+suf,p.lepton1->Eta(),weight,netabin,etabins);
    for(int i=0,n=p.leptons.size();i<n;i++){
      double eta=fabs(p.leptons.at(i)->Eta());
      double pt=p.leptons.at(i)->Pt();
      double riso=p.leptons.at(i)->RelIso();
      double energy=p.leptons.at(i)->E();
      FillHist(p.prefix+region+p.hprefix+Form("l%detapt",i)+p.suffix+suf,eta,pt,weight,netabin,etabins,nptbin,ptbins);
      FillHist(p.prefix+region+p.hprefix+Form("l%detae",i)+p.suffix+suf,eta,energy,weight,netabin,etabins,nptbin,ptbins);
      FillHist(p.prefix+region+p.hprefix+Form("l%driso",i)+p.suffix+suf,riso,weight,100,0,2);
    }
    if((p.channel=="EE"||p.channel=="MM")&&p.prefix.Contains("cpt/")){
      TString prefix=p.prefix+region;
      prefix.ReplaceAll("EE","ee");
      prefix.ReplaceAll("MM","mm");
      prefix.ReplaceAll("cpt/","");
      double tf=GetFakeTF(p,"cpt");
      FillHist(prefix+"fakecpt_"+p.hprefix+"dimass"+p.suffix+suf,(*p.lepton0+*p.lepton1).M(),weight*tf,98,52,150);
      FillHist(prefix+"fakecpt_"+p.hprefix+"dimass_wide"+p.suffix+suf,(*p.lepton0+*p.lepton1).M(),weight*tf,nmbin,mbins);
      FillHist(prefix+"fakecpt_"+p.hprefix+"dipt"+p.suffix+suf,(*p.lepton0+*p.lepton1).Pt(),weight*tf,100,0,100);
      FillHist(prefix+"fakecpt_"+p.hprefix+"dirap"+p.suffix+suf,(*p.lepton0+*p.lepton1).Rapidity(),weight*tf,60,-3,3);
      FillHist(prefix+"fakecpt_"+p.hprefix+"lpt"+p.suffix+suf,p.lepton0->Pt(),weight*tf,nptbin,ptbins);
      FillHist(prefix+"fakecpt_"+p.hprefix+"lpt"+p.suffix+suf,p.lepton1->Pt(),weight*tf,nptbin,ptbins);
      FillHist(prefix+"fakecpt_"+p.hprefix+"l0pt"+p.suffix+suf,p.lepton0->Pt(),weight*tf,nptbin,ptbins);
      FillHist(prefix+"fakecpt_"+p.hprefix+"l1pt"+p.suffix+suf,p.lepton1->Pt(),weight*tf,nptbin,ptbins);
      FillHist(prefix+"fakecpt_"+p.hprefix+"leta"+p.suffix+suf,p.lepton0->Eta(),weight*tf,netabin,etabins);
      FillHist(prefix+"fakecpt_"+p.hprefix+"leta"+p.suffix+suf,p.lepton1->Eta(),weight*tf,netabin,etabins);
      FillHist(prefix+"fakecpt_"+p.hprefix+"l0eta"+p.suffix+suf,p.lepton0->Eta(),weight*tf,netabin,etabins);
      FillHist(prefix+"fakecpt_"+p.hprefix+"l1eta"+p.suffix+suf,p.lepton1->Eta(),weight*tf,netabin,etabins);
    }else if((p.channel=="EE"||p.channel=="MM")&&p.prefix.Contains("mpt/")){
      TString prefix=p.prefix+region;
      prefix.ReplaceAll("EE","ee");
      prefix.ReplaceAll("MM","mm");
      prefix.ReplaceAll("mpt/","");
      double tf=GetFakeTF(p,"mpt");
      FillHist(prefix+"fakempt_"+p.hprefix+"dimass"+p.suffix+suf,(*p.lepton0+*p.lepton1).M(),weight*tf,98,52,150);
      FillHist(prefix+"fakempt_"+p.hprefix+"dimass_wide"+p.suffix+suf,(*p.lepton0+*p.lepton1).M(),weight*tf,nmbin,mbins);
      FillHist(prefix+"fakempt_"+p.hprefix+"dipt"+p.suffix+suf,(*p.lepton0+*p.lepton1).Pt(),weight*tf,100,0,100);
      FillHist(prefix+"fakempt_"+p.hprefix+"dirap"+p.suffix+suf,(*p.lepton0+*p.lepton1).Rapidity(),weight*tf,60,-3,3);
      FillHist(prefix+"fakempt_"+p.hprefix+"lpt"+p.suffix+suf,p.lepton0->Pt(),weight*tf,nptbin,ptbins);
      FillHist(prefix+"fakempt_"+p.hprefix+"lpt"+p.suffix+suf,p.lepton1->Pt(),weight*tf,nptbin,ptbins);
      FillHist(prefix+"fakempt_"+p.hprefix+"l0pt"+p.suffix+suf,p.lepton0->Pt(),weight*tf,nptbin,ptbins);
      FillHist(prefix+"fakempt_"+p.hprefix+"l1pt"+p.suffix+suf,p.lepton1->Pt(),weight*tf,nptbin,ptbins);
      FillHist(prefix+"fakempt_"+p.hprefix+"leta"+p.suffix+suf,p.lepton0->Eta(),weight*tf,netabin,etabins);
      FillHist(prefix+"fakempt_"+p.hprefix+"leta"+p.suffix+suf,p.lepton1->Eta(),weight*tf,netabin,etabins);
      FillHist(prefix+"fakempt_"+p.hprefix+"l0eta"+p.suffix+suf,p.lepton0->Eta(),weight*tf,netabin,etabins);
      FillHist(prefix+"fakempt_"+p.hprefix+"l1eta"+p.suffix+suf,p.lepton1->Eta(),weight*tf,netabin,etabins);
    }else if((p.channel=="EE"||p.channel=="MM")){
      TString prefix=p.prefix+region;
      prefix.ReplaceAll("EE","ee");
      prefix.ReplaceAll("MM","mm");
      double tf=GetFakeTF(p);
      FillHist(prefix+"fake_"+p.hprefix+"dimass"+p.suffix+suf,(*p.lepton0+*p.lepton1).M(),weight*tf,98,52,150);
      FillHist(prefix+"fake_"+p.hprefix+"dimass_wide"+p.suffix+suf,(*p.lepton0+*p.lepton1).M(),weight*tf,nmbin,mbins);
      FillHist(prefix+"fake_"+p.hprefix+"dipt"+p.suffix+suf,(*p.lepton0+*p.lepton1).Pt(),weight*tf,100,0,100);
      FillHist(prefix+"fake_"+p.hprefix+"dirap"+p.suffix+suf,(*p.lepton0+*p.lepton1).Rapidity(),weight*tf,60,-3,3);
      FillHist(prefix+"fake_"+p.hprefix+"lpt"+p.suffix+suf,p.lepton0->Pt(),weight*tf,nptbin,ptbins);
      FillHist(prefix+"fake_"+p.hprefix+"lpt"+p.suffix+suf,p.lepton1->Pt(),weight*tf,nptbin,ptbins);
      FillHist(prefix+"fake_"+p.hprefix+"l0pt"+p.suffix+suf,p.lepton0->Pt(),weight*tf,nptbin,ptbins);
      FillHist(prefix+"fake_"+p.hprefix+"l1pt"+p.suffix+suf,p.lepton1->Pt(),weight*tf,nptbin,ptbins);
      FillHist(prefix+"fake_"+p.hprefix+"leta"+p.suffix+suf,p.lepton0->Eta(),weight*tf,netabin,etabins);
      FillHist(prefix+"fake_"+p.hprefix+"leta"+p.suffix+suf,p.lepton1->Eta(),weight*tf,netabin,etabins);
      FillHist(prefix+"fake_"+p.hprefix+"l0eta"+p.suffix+suf,p.lepton0->Eta(),weight*tf,netabin,etabins);
      FillHist(prefix+"fake_"+p.hprefix+"l1eta"+p.suffix+suf,p.lepton1->Eta(),weight*tf,netabin,etabins);
    }
  }
}
vector<Muon> FakeAnalyzer::ToConePt(vector<Muon> muons){
  for(auto& muon:muons){
    muon*=1+muon.RelIso();
  }
  std::sort(muons.begin(),muons.end(),PtComparing);
  return muons;
}
vector<Electron> FakeAnalyzer::ToConePt(vector<Electron> electrons){
  for(auto& electron:electrons){
    electron*=1+electron.RelIso();
  }
  std::sort(electrons.begin(),electrons.end(),PtComparing);
  return electrons;
}
vector<Muon> FakeAnalyzer::ToModifiedPt(vector<Muon> muons){
  for(auto& muon:muons){
    muon*=1+muon.RelIso()*TMath::Max(0.,TMath::Min(1.,(muon.Pt()-30)/30));
  }
  std::sort(muons.begin(),muons.end(),PtComparing);
  return muons;
}
vector<Electron> FakeAnalyzer::ToModifiedPt(vector<Electron> electrons){
  for(auto& electron:electrons){
    electron*=1+electron.RelIso()*TMath::Max(0.,TMath::Min(1.,(electron.Pt()-30)/30));
  }
  std::sort(electrons.begin(),electrons.end(),PtComparing);
  return electrons;
}
void FakeAnalyzer::UseSelectiveCharge(Parameter& p){
  vector<Electron> electrons;
  vector<Electron> aelectrons;
  for(auto e:p.electrons){
    if(e.IsGsfCtfScPixChargeConsistent()) electrons.push_back(e);
  }
  for(auto e:p.aelectrons){
    if(e.IsGsfCtfScPixChargeConsistent()) aelectrons.push_back(e);
  }
  p.SetElectrons(electrons);
  p.SetAElectrons(aelectrons);
  p.k.electronIDSF2("Electron_SelQ_MediumID");
}
