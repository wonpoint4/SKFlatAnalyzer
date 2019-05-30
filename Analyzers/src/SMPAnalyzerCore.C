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
}

void SMPAnalyzerCore::initializeAnalyzer(){
  SetupZPtWeight();
  IsDYSample=false;
  if(MCSample.Contains("DYJets")||MCSample.Contains("ZToEE")||MCSample.Contains("ZToMuMu")) IsDYSample=true;
}

void SMPAnalyzerCore::FillGenHists(TString pre,TString suf,TLorentzVector genl0,TLorentzVector genl1,TLorentzVector genfsr,double w){
  TLorentzVector genZ=genl0+genl1+genfsr;
  FillHist(pre+"genZmass"+suf,genZ.M(),w,400,0,400);
  FillHist(pre+"genZpt"+suf,genZ.Pt(),w,400,0,400);
  FillHist(pre+"genZrap"+suf,genZ.Rapidity(),w,60,-6,6);
  if(genl0.Pt()<genl1.Pt()){
    TLorentzVector temp=genl0;
    genl0=genl1;
    genl1=temp;
  }
  FillHist(Form("%sgenl0pt%s",pre.Data(),suf.Data()),genl0.Pt(),w,200,0,200);
  FillHist(Form("%sgenl0eta%s",pre.Data(),suf.Data()),genl0.Eta(),w,200,-5,5);
  FillHist(Form("%sgenl1pt%s",pre.Data(),suf.Data()),genl1.Pt(),w,200,0,200);
  FillHist(Form("%sgenl1eta%s",pre.Data(),suf.Data()),genl1.Eta(),w,200,-5,5);
  FillHist(pre+"lldelR"+suf,genl0.DeltaR(genl1),w,70,0,7);  
  FillHist(pre+"lldelphi"+suf,genl0.DeltaPhi(genl1),w,80,-4,4);
}

void SMPAnalyzerCore::FillBasicHists(TString pre,TString suf,const vector<Lepton*>& leps,double w){
  Particle dilepton=((*leps.at(0))+(*leps.at(1)));
  FillHist(pre+"dimass"+suf,dilepton.M(),w,400,0,400);
  FillHist(pre+"dipt"+suf,dilepton.Pt(),w,400,0,400);
  FillHist(pre+"dirap"+suf,dilepton.Rapidity(),w,50,-5,5);
  for(int i=0;i<(int)leps.size();i++){
    FillHist(Form("%sl%dpt%s",pre.Data(),i,suf.Data()),leps.at(i)->Pt(),w,1000,0,1000);
    FillHist(Form("%sl%deta%s",pre.Data(),i,suf.Data()),leps.at(i)->Eta(),w,200,-5,5);
    FillHist(Form("%sl%driso%s",pre.Data(),i,suf.Data()),leps.at(i)->RelIso(),w,30,0,0.3);
  }
  FillHist(pre+"lldelR"+suf,leps.at(0)->DeltaR(*leps.at(1)),w,70,0,7);  
  FillHist(pre+"lldelphi"+suf,leps.at(0)->DeltaPhi(*leps.at(1)),w,80,-4,4);
  FillHist(pre+"nPV"+suf,nPV,w,60,0,60);
}
void SMPAnalyzerCore::FillSystematicHists(TString pre,TString suf,const vector<Lepton*>& leps,map<TString,double> map_systematic){
  for(auto iter=map_systematic.begin();iter!=map_systematic.end();iter++){
    FillBasicHists(pre,"_"+iter->first+suf,leps,iter->second);
  }
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
  map<TString,TH2F*>* map_hist_Lepton=NULL;
  if(leps[0]->LeptonFlavour()==Lepton::MUON){
    map_hist_Lepton=&mcCorr->map_hist_Muon;
  }else if(leps[0]->LeptonFlavour()==Lepton::ELECTRON){
    map_hist_Lepton=&mcCorr->map_hist_Electron;
  }else{
    cout <<"[SMPAnalyzerCore::Trigger_SF] Wrong flavour"<<endl;
    exit(EXIT_FAILURE);
  }    
      
  double this_pt[2]={},this_eta[2]={};
  TString this_charge[2];
  for(int i=0;i<2;i++){
    if(leps[i]->LeptonFlavour()==Lepton::MUON){
      this_pt[i]=((Muon*)leps.at(i))->MiniAODPt();
      this_eta[i]=leps.at(i)->Eta();
    }else if(leps[i]->LeptonFlavour()==Lepton::ELECTRON){
      this_pt[i]=leps.at(i)->Pt();
      this_eta[i]=((Electron*)leps.at(i))->scEta();
    }else{
      cout << "[SMPAnalyzerCore::Trigger_SF] It is not lepton"<<endl;
      exit(EXIT_FAILURE);
    }
    if(triggerSF_key0.Contains(TRegexp("_Q$"))&&triggerSF_key1.Contains(TRegexp("_Q$"))){
      if(leps.at(i)->Charge()>0) this_charge[i]="Plus";
      else this_charge[i]="Minus";
    }
  }
  if(!(DataYear==2016&&leps[0]->LeptonFlavour()==Lepton::MUON)){
    TH2F* this_hist[2]={NULL,NULL};
    double triggerSF=1.;
    this_hist[0]=(*map_hist_Lepton)["Trigger_SF_"+triggerSF_key0+this_charge[0]];
    this_hist[1]=(*map_hist_Lepton)["Trigger_SF_"+triggerSF_key1+this_charge[1]];
    if(!this_hist[0]||!this_hist[1]){
      cout << "[SMPAnalyzerCore::Trigger_SF] No Trigger_SF_"<<triggerSF_key0<<" or Trigger_SF_"<<triggerSF_key1<<endl;
      exit(EXIT_FAILURE);
    }
    for(int i=0;i<2;i++){
      triggerSF*=GetBinContentUser(this_hist[i],this_eta[i],this_pt[i],sys);
    }
    return triggerSF;
  }else{
    TH2F* this_hist[8]={};
    TString sdata[2]={"DATA","MC"};
    TString speriod[2]={"BCDEF","GH"};
    for(int id=0;id<2;id++){
      for(int ip=0;ip<2;ip++){
	this_hist[4*id+2*ip]=(*map_hist_Lepton)["Trigger_Eff_"+sdata[id]+"_"+triggerSF_key0+this_charge[0]+"_"+speriod[ip]];
	this_hist[4*id+2*ip+1]=(*map_hist_Lepton)["Trigger_Eff_"+sdata[id]+"_"+triggerSF_key1+this_charge[1]+"_"+speriod[ip]];
      }
    }
    if(!this_hist[0]||!this_hist[1]||!this_hist[2]||!this_hist[3]){
      cout << "[SMPAnalyzerCore::Trigger_SF] No "<<triggerSF_key0<<" or "<<triggerSF_key1<<endl;
      exit(EXIT_FAILURE);
    }
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
    
    double triggerEff[4]={1.,1.,1.,1.};
    for(int i=0;i<4;i++){
      for(int il=0;il<2;il++){
	triggerEff[i]*=GetBinContentUser(this_hist[2*i+il],this_eta[il],this_pt[il],(i<2?1.:-1.)*sys);
      }
    }
    return (triggerEff[0]*WeightBtoF+triggerEff[1]*WeightGtoH)/(triggerEff[2]*WeightBtoF+triggerEff[3]*WeightGtoH);
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
  if(id=="passMediumID_selective"){
    std::vector<Electron> electrons = GetAllElectrons();
    std::vector<Electron> out;
    for(unsigned int i=0; i<electrons.size(); i++){
      Electron this_electron= electrons.at(i);
      if(!( this_electron.Pt()>ptmin ))	continue;
      if(!( fabs(this_electron.scEta())<fetamax )) continue;
      if(!( this_electron.PassID("passMediumID") ))	continue;
      if(!electron_isGsfCtfScPixChargeConsistent->at(i)) continue;
      out.push_back(this_electron);
    }
    return out;
  }else return GetElectrons(id,ptmin,fetamax);
}    
std::vector<Muon> SMPAnalyzerCore::SMPGetMuons(TString id,double ptmin,double fetamax){
  if(id=="POGTightWithLooseTrkIso"){
    vector<Muon> muons=GetMuons("POGTight",ptmin,fetamax);
    vector<Muon> out;
    for(auto const& muon: muons){
      if(muon.TrkIso()/muon.Pt()<0.1) out.push_back(muon);
    }
    return out;
  }else return GetMuons(id,ptmin,fetamax);
}
