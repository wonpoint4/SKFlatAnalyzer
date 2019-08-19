#include "SMPAnalyzerCore.h"

SMPAnalyzerCore::SMPAnalyzerCore(){
  hzpt_muon=NULL;
  hzpt_electron=NULL;
  hzpt_norm_muon=NULL;
  hzpt_norm_electron=NULL;
}

SMPAnalyzerCore::~SMPAnalyzerCore(){
  if(hzpt_muon) delete hzpt_muon;
  if(hzpt_norm_muon) delete hzpt_norm_muon;
  if(hzpt_electron) delete hzpt_electron;
  if(hzpt_norm_electron) delete hzpt_norm_electron;
  if(roc) delete roc;
  if(rocele) delete rocele;
}

void SMPAnalyzerCore::initializeAnalyzer(){
  SetupZPtWeight();
  SetupRoccoR();
  IsDYSample=false;
  if(MCSample.Contains("DYJets")||MCSample.Contains("ZToEE")||MCSample.Contains("ZToMuMu")||MCSample.Contains(TRegexp("DY[0-9]Jets"))) IsDYSample=true;
}

void SMPAnalyzerCore::FillDileptonHists(TString pre,TString suf,Particle *l0,Particle *l1,double w){
  TLorentzVector dilepton=*l0+*l1;
  double dimass=dilepton.M();
  double dipt=dilepton.Pt();
  double dirap=dilepton.Rapidity();
  FillHist(pre+"dimass"+suf,dimass,w,400,0,400);
  FillHist(pre+"dipt"+suf,dipt,w,400,0,400);
  FillHist(pre+"dirap"+suf,dirap,w,120,-6,6);
  vector<Particle*> leps;
  if(l0->Pt()>l1->Pt()) leps={l0,l1};
  else leps={l1,l0};
  for(int i=0;i<(int)leps.size();i++){
    FillHist(Form("%sl%dpt%s",pre.Data(),i,suf.Data()),leps.at(i)->Pt(),w,400,0,400);
    FillHist(Form("%sl%deta%s",pre.Data(),i,suf.Data()),leps.at(i)->Eta(),w,200,-5,5);
    FillHist(Form("%sl%dphi%s",pre.Data(),i,suf.Data()),leps.at(i)->Phi(),w,160,-4,4);
  }
  FillHist(pre+"lldelR"+suf,l0->DeltaR(*l1),w,70,0,7);  
  FillHist(pre+"lldelphi"+suf,l0->DeltaPhi(*l1),w,80,-4,4);
  FillHist(pre+"lldeleta"+suf,fabs(l0->Eta()-l1->Eta()),w,100,-5,5);
}
double SMPAnalyzerCore::Lepton_SF(TString histkey,const Lepton* lep,int sys){
  if(IsDATA) return 1.;
  if(histkey=="") return 1.;
  if(histkey=="Default") return 1.;
  double this_pt,this_eta;
  TH2* this_hist=NULL;
  if(histkey.Contains(TRegexp("_Q$"))){
    if(lep->Charge()>0) histkey+="Plus";
    else histkey+="Minus";
  }
  if(lep->LeptonFlavour()==Lepton::MUON){
    this_pt=((Muon*)lep)->MiniAODPt();
    this_eta=lep->Eta();
    this_hist=mcCorr->map_hist_Muon[histkey];
  }else if(lep->LeptonFlavour()==Lepton::ELECTRON){
    this_pt=lep->Pt();
    this_eta=((Electron*)lep)->scEta();
    this_hist=mcCorr->map_hist_Electron[histkey];
  }else{
    cout <<"[SMPAnalyzerCore::Lepton_SF] It is not lepton"<<endl;
    exit(EXIT_FAILURE);
  }
  if(!this_hist){
    cout <<"[SMPAnalyzerCore::Lepton_SF] no hist "<<histkey<<endl;
    exit(EXIT_FAILURE);
  }    
  double this_x,this_y;
  if(this_hist->GetXaxis()->GetXmax()>this_hist->GetYaxis()->GetXmax()){
    this_x=this_pt;
    this_y=this_eta;
  }else{
    this_x=this_eta;
    this_y=this_pt;
  }
  return GetBinContentUser(this_hist,this_x,this_y,sys);
}

double SMPAnalyzerCore::DileptonTrigger_SF(TString triggerSF_key0,TString triggerSF_key1,const vector<Lepton*>& leps,int sys){
  if(IsDATA) return 1;
  if(triggerSF_key0==""&&triggerSF_key1=="") return 1;
  if(leps.size()<2){
    cout<<"[SMPAnalyzerCore::Trigger_SF] only dilepton algorithm"<<endl;
    return 1;
  }else if(leps.size()>2){
    cout<<"[SMPAnalyzerCore::Trigger_SF] only dilepton algorithm"<<endl;
  }
  TString histkeys[2]={triggerSF_key0,triggerSF_key1};
  if(!(DataYear==2016&&leps[0]->LeptonFlavour()==Lepton::MUON)){
    double triggerSF=1.;
    for(int i=0;i<2;i++){
      triggerSF*=Lepton_SF("Trigger_SF_"+histkeys[i],leps.at(i),sys);
    }
    return triggerSF;
  }else{
    double lumi_periodB = 5.929001722;
    double lumi_periodC = 2.645968083;
    double lumi_periodD = 4.35344881;
    double lumi_periodE = 4.049732039;
    double lumi_periodF = 3.157020934;
    double lumi_periodG = 7.549615806;
    double lumi_periodH = 8.545039549 + 0.216782873;
    double total_lumi = (lumi_periodB+lumi_periodC+lumi_periodD+lumi_periodE+lumi_periodF+lumi_periodG+lumi_periodH);

    double WeightBtoF = (lumi_periodB+lumi_periodC+lumi_periodD+lumi_periodE+lumi_periodF)/total_lumi;
    double WeightGtoH = (lumi_periodG+lumi_periodH)/total_lumi;
    
    double triggerEff[2][2]={{1.,1.},{1.,1.}};
    TString sdata[2]={"DATA","MC"};
    TString speriod[2]={"BCDEF","GH"};
    TString scharge[2]={"",""};
    for(int il=0;il<2;il++){
      TString scharge="";
      if(histkeys[il].Contains(TRegexp("_Q$"))){
	if(leps.at(il)->Charge()>0) scharge="Plus";
	else scharge="Minus";
      }
      for(int id=0;id<2;id++){
	for(int ip=0;ip<2;ip++){
	  triggerEff[id][ip]*=Lepton_SF("Trigger_Eff_"+sdata[id]+"_"+histkeys[il]+scharge+"_"+speriod[ip],leps.at(il),(id?-1.:1.)*sys);
	}
      }
    }
    return (triggerEff[0][0]*WeightBtoF+triggerEff[0][1]*WeightGtoH)/(triggerEff[1][0]*WeightBtoF+triggerEff[1][1]*WeightGtoH);
  }
}
void SMPAnalyzerCore::SetupZPtWeight(){
  cout<<"[SMPAnalyzerCore::SetupZPtWeight] setting zptcor"<<endl;
  TString datapath=getenv("DATA_DIR");
  TFile fzpt(datapath+"/"+TString::Itoa(DataYear,10)+"/ZPt/ZPtWeight.root");
  TString sflavour[2]={"muon","electron"};
  TH2D **hzpt=NULL,**hzpt_norm=NULL;
  for(int ifl=0;ifl<2;ifl++){
    if(ifl==0){
      hzpt=&hzpt_muon;
      hzpt_norm=&hzpt_norm_muon;
    }else if(ifl==1){
      hzpt=&hzpt_electron;
      hzpt_norm=&hzpt_norm_electron;
    }
    for(int i=0;i<20;i++){
      TH2D* this_hzpt=(TH2D*)fzpt.Get(Form("%s%d_iter%d",sflavour[ifl].Data(),DataYear,i));
      if(this_hzpt){
	if(*hzpt){
	  (*hzpt)->Multiply(this_hzpt);
	  cout<<"[SMPAnalyzerCore::SetupZPtWeight] setting "<<sflavour[ifl]<<" zptcor iter"<<i<<endl;
	}else{
	  (*hzpt)=this_hzpt;
	  cout<<"[SMPAnalyzerCore::SetupZPtWeight] setting first "<<sflavour[ifl]<<" zptcor"<<i<<endl;
	}
      }else break;
    }
    if(*hzpt) (*hzpt)->SetDirectory(0);
    *hzpt_norm=(TH2D*)fzpt.Get(Form("%s%d_norm",sflavour[ifl].Data(),DataYear));
    if(*hzpt_norm){
      (*hzpt_norm)->SetDirectory(0);
      cout<<"[SMPAnalyzerCore::SetupZPtWeight] setting "<<sflavour[ifl]<<" zptcor norm"<<endl;
    }
  }
}
void SMPAnalyzerCore::SetupRoccoR(){
  cout<<"[SMPAnalyzerCore::SetupRoccoR] setting Rocheseter Correction"<<endl;
  TString datapath=getenv("DATA_DIR");
  TString textfile=datapath+"/"+TString::Itoa(DataYear,10)+"/RoccoR/RoccoR"+TString::Itoa(DataYear,10)+".txt";
  roc=new RoccoR(textfile.Data());
  rocele=new RocelecoR((datapath+"/"+TString::Itoa(DataYear,10)+"/RoccoR/RocelecoR"+TString::Itoa(DataYear,10)+"_new.txt").Data());
}
double SMPAnalyzerCore::GetZPtWeight(double zpt,double zrap,Lepton::Flavour flavour){
  double valzptcor=1.;
  double valzptcor_norm=1.;
  TH2D* hzpt=NULL;
  TH2D* hzpt_norm=NULL;
  if(flavour==Lepton::MUON){
    hzpt=hzpt_muon;
    hzpt_norm=hzpt_norm_muon;
  }else if(flavour==Lepton::ELECTRON){
    hzpt=hzpt_electron;
    hzpt_norm=hzpt_norm_electron;
  }
  if(hzpt) valzptcor*=GetBinContentUser(hzpt,zpt,zrap,0);
  if(hzpt_norm) valzptcor_norm*=GetBinContentUser(hzpt_norm,zpt,zrap,0);
  return valzptcor*valzptcor_norm;
}
void SMPAnalyzerCore::PrintGens(const vector<Gen>& gens){
  cout<<"index\tpid\tmother\tstatus\tpropt\thard\n";
  for(int i=0;i<(int)gens.size();i++){
    cout<<gens.at(i).Index()<<"\t"<<gens.at(i).PID()<<"\t"<<gens.at(i).MotherIndex()<<"\t"<<gens.at(i).Status()<<"\t"<<gens.at(i).isPrompt()<<"\t"<<gens.at(i).isHardProcess()<<endl;
  }
}

double SMPAnalyzerCore::GetBinContentUser(TH2* hist,double valx,double valy,int sys){
  double xmin=hist->GetXaxis()->GetXmin();
  double xmax=hist->GetXaxis()->GetXmax();
  double ymin=hist->GetYaxis()->GetXmin();
  double ymax=hist->GetYaxis()->GetXmax();
  if(xmin>=0) valx=fabs(valx);
  if(valx<xmin) valx=xmin+0.001;
  if(valx>xmax) valx=xmax-0.001;
  if(ymin>=0) valy=fabs(valy);
  if(valy<ymin) valy=ymin+0.001;
  if(valy>ymax) valy=ymax-0.001;
  return hist->GetBinContent(hist->FindBin(valx,valy))+sys*hist->GetBinError(hist->FindBin(valx,valy));
}

void SMPAnalyzerCore::GetGenIndex(const vector<Gen>& gens,int& parton0,int& parton1,int& hardl0,int& hardl1,int& l0,int& l1,vector<int>& photons){
  parton0=0;parton1=0;hardl0=0;hardl1=0;l0=0;l1=0;photons.clear();
  if(!IsDYSample){
    cout <<"[SMPAnalyzerCore::GetGenIndex] this is for DY event"<<endl;
    return;
  }
  for(int i=0;i<(int)gens.size();i++){
    if(!gens.at(i).isPrompt()) continue;
    int genpid=gens.at(i).PID();
    if(gens.at(i).isHardProcess()){
      if(abs(genpid)<7||genpid==21){
	if(!parton0) parton0=i;
	else if(!parton1) parton1=i;
      }else if(abs(genpid)==11||abs(genpid)==13||abs(genpid)==15){
	if(!hardl0) hardl0=i;
	else if(!hardl1&&genpid==-gens[hardl0].PID()) hardl1=i;
      }
    }
    if(gens.at(i).Status()==1){
      if(abs(genpid)==11||abs(genpid)==13){
	if(!l0) l0=i;
	else if(!l1&&genpid==-gens[l0].PID()) l1=i;
      }else if(gens.at(i).PID()==22){
	photons.push_back(i);
      }
    }
  }
  if(abs(gens[hardl0].PID())!=15){
    if(!hardl0||!hardl1||!l0||!l1){
      PrintGens(gens);
      cout <<"[SMPAnalyzerCore::GetGenIndex] something is wrong"<<endl;
      exit(EXIT_FAILURE);
    }
  }   
  return;
}

std::vector<Electron> SMPAnalyzerCore::SMPGetElectrons(TString id, double ptmin, double fetamax){
  std::vector<Electron> out;
  if(id=="passMediumID_selective"){
    std::vector<Electron> electrons = GetAllElectrons();
    for(unsigned int i=0; i<electrons.size(); i++){
      Electron this_electron= electrons.at(i);
      if(!( this_electron.Pt()>ptmin ))	continue;
      if(!( fabs(this_electron.scEta())<fetamax )) continue;
      if(!( this_electron.PassID("passMediumID") ))	continue;
      if(!electron_isGsfCtfScPixChargeConsistent->at(i)) continue;
      out.push_back(this_electron);
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
      if(muon.TrkIso()/muon.Pt()<0.1) out.push_back(muon);
    }
  }else out=GetMuons(id,ptmin,fetamax);
  std::sort(out.begin(),out.end(),PtComparing);
  return out;
}

std::vector<Muon> SMPAnalyzerCore::MuonMomentumCorrection(const vector<Muon>& muons,int set,int member){
  std::vector<Muon> out;
  for(auto muon:muons){
    double rc=1.;
    if(set>=0){
      if(IsDATA){
	rc=roc->kScaleDT(muon.Charge(),muon.MiniAODPt(),muon.Eta(),muon.Phi(),set,member);
      }else{
	rc=roc->kSmearMC(muon.Charge(),muon.MiniAODPt(),muon.Eta(),muon.Phi(),muon.TrackerLayers(),gRandom->Rndm(),set,member);
      }      
    }
    muon.SetPtEtaPhiM(muon.MiniAODPt()*rc,muon.Eta(),muon.Phi(),muon.M());
    out.push_back(muon);
  }
  return out;
}

std::vector<Electron> SMPAnalyzerCore::ElectronEnergyCorrection(const vector<Electron>& electrons,int set,int member){
  std::vector<Electron> out;
  for(auto electron:electrons){
    double rc=1.;
    double pt=electron.Pt();
    if(set>=0){
      if(IsDATA){
	pt=electron.Pt();
	rc=rocele->kScaleDT(electron.Charge(),pt,electron.Eta(),electron.Phi(),set,member);
      }else{
	pt=electron.UncorrPt();
	rc=rocele->kScaleMC(electron.Charge(),pt,electron.Eta(),electron.Phi(),set,member);
      }      
    }
    electron.SetPtEtaPhiM(pt*rc,electron.Eta(),electron.Phi(),electron.M());
    out.push_back(electron);
  }
  return out;
}
  
