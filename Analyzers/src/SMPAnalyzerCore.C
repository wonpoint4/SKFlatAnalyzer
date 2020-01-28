#include "SMPAnalyzerCore.h"

SMPAnalyzerCore::SMPAnalyzerCore(){
  roc=NULL;
  rocele=NULL;
  hz0=NULL;
}

SMPAnalyzerCore::~SMPAnalyzerCore(){
  for(auto& [key,hist]:map_hist_zpt){
    if(hist) delete hist;
  }
  if(roc) delete roc;
  if(rocele) delete rocele;
  if(hz0) delete hz0;
}

void SMPAnalyzerCore::initializeAnalyzer(){
  SetupZptWeight();
  SetupRoccoR();
  SetupZ0Weight();
  IsDYSample=false;
  if(MCSample.Contains("DYJets")||MCSample.Contains("ZToEE")||MCSample.Contains("ZToMuMu")||MCSample.Contains(TRegexp("DY[0-9]Jets"))) IsDYSample=true;
}
void SMPAnalyzerCore::FillHist(TString histname, double value, map<TString,double> weights, int n_bin, double x_min, double x_max){
  for(const auto& [suffix,weight]:weights) FillHist(histname+suffix,value,weight,n_bin,x_min,x_max);
}
void SMPAnalyzerCore::FillHist(TString histname, double value, map<TString,double> weights, int n_bin, double *xbins){
  for(const auto& [suffix,weight]:weights) FillHist(histname+suffix,value,weight,n_bin,xbins);
}
void SMPAnalyzerCore::FillHist(TString histname, double value_x, double value_y, map<TString,double> weights, int n_binx, double x_min, double x_max, int n_biny, double y_min, double y_max){
  for(const auto& [suffix,weight]:weights) FillHist(histname+suffix,value_x,value_y,weight,n_binx,x_min,x_max,n_biny,y_min,y_max);
}
void SMPAnalyzerCore::FillHist(TString histname, double value_x, double value_y, map<TString,double> weights, int n_binx, double *xbins, int n_biny, double *ybins){
  for(const auto& [suffix,weight]:weights) FillHist(histname+suffix,value_x,value_y,weight,n_binx,xbins,n_biny,ybins);
}
void SMPAnalyzerCore::FillHist(TString histname, double value_x, double value_y, double value_z, map<TString,double> weights, int n_binx, double x_min, double x_max, int n_biny, double y_min, double y_max, int n_binz, double z_min, double z_max){
  for(const auto& [suffix,weight]:weights) FillHist(histname+suffix,value_x,value_y,value_z,weight,n_binx,x_min,x_max,n_biny,y_min,y_max,n_binz,z_min,z_max);
}
void SMPAnalyzerCore::FillHist(TString histname, double value_x, double value_y, double value_z, map<TString,double> weights, int n_binx, double *xbins, int n_biny, double *ybins, int n_binz, double *zbins){
  for(const auto& [suffix,weight]:weights) FillHist(histname+suffix,value_x,value_y,value_z,weight,n_binx,xbins,n_biny,ybins,n_binz,zbins);
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

double SMPAnalyzerCore::LeptonTrigger_SF(TString triggerSF_key,const vector<Lepton*>& leps,int sys){
  if(IsDATA) return 1;
  if(triggerSF_key=="") return 1;
  if(triggerSF_key=="Default") return 1;

  double triggerSF=1.,data_eff=1.,mc_eff=1.;
  for(const auto& lep:leps){
    data_eff*=1-Lepton_SF("Trigger_Eff_DATA_"+triggerSF_key,lep,sys);
    mc_eff*=1-Lepton_SF("Trigger_Eff_MC_"+triggerSF_key,lep,-sys);
  }
  data_eff=1-data_eff;
  mc_eff=1-mc_eff;
  return data_eff/mc_eff;
}
double SMPAnalyzerCore::DileptonTrigger_SF(TString triggerSF_key0,TString triggerSF_key1,const vector<Lepton*>& leps,int sys){
  if(IsDATA) return 1;
  if(triggerSF_key0==""&&triggerSF_key1=="") return 1;
  if(leps.size()!=2){
    cout<<"[SMPAnalyzerCore::Trigger_SF] only dilepton algorithm"<<endl;
    return 1;
  }
  TString histkeys[2]={triggerSF_key0,triggerSF_key1};
  if(!(DataYear==2016&&leps[0]->LeptonFlavour()==Lepton::MUON)){
    double eff[2][2][2]={}; //[data/mc][l0/l1][leg1/leg2]
    TString sdata[2]={"DATA","MC"};
    for(int id=0;id<2;id++){
      for(int ilep=0;ilep<2;ilep++){
	for(int ileg=0;ileg<2;ileg++){
	  TString scharge="";
	  if(histkeys[ileg].Contains(TRegexp("_Q$"))){
	    if(leps.at(ilep)->Charge()>0) scharge="Plus";
	    else scharge="Minus";
	  }
	  eff[id][ilep][ileg]=Lepton_SF("Trigger_Eff_"+sdata[id]+"_"+histkeys[ileg]+scharge,leps.at(ilep),(id?-1.:1.)*sys);
	}
      }
    }
    double eff_data=eff[0][0][1]*eff[0][1][1]-(eff[0][0][1]-eff[0][0][0])*(eff[0][1][1]-eff[0][1][0]);
    double eff_mc=eff[1][0][1]*eff[1][1][1]-(eff[1][0][1]-eff[1][0][0])*(eff[1][1][1]-eff[1][1][0]);
    return eff_data/eff_mc;
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
    
    double eff[2][2][2][2]={}; //[period][data/mc][l0/l1][leg1/leg2]
    TString speriod[2]={"BCDEF","GH"};
    TString sdata[2]={"DATA","MC"};
    for(int ip=0;ip<2;ip++){
      for(int id=0;id<2;id++){
	for(int ilep=0;ilep<2;ilep++){
	  for(int ileg=0;ileg<2;ileg++){
	    TString scharge="";
	    if(histkeys[ileg].Contains(TRegexp("_Q$"))){
	      if(leps.at(ilep)->Charge()>0) scharge="Plus";
	      else scharge="Minus";
	    }
	    eff[ip][id][ilep][ileg]=Lepton_SF("Trigger_Eff_"+sdata[id]+"_"+histkeys[ileg]+scharge+"_"+speriod[ip],leps.at(ilep),(id?-1.:1.)*sys);
	  }
	}
      }
    }
    double eff_data_BtoF=eff[0][0][0][1]*eff[0][0][1][1]-(eff[0][0][0][1]-eff[0][0][0][0])*(eff[0][0][1][1]-eff[0][0][1][0]);
    double eff_data_GtoH=eff[1][0][0][1]*eff[1][0][1][1]-(eff[1][0][0][1]-eff[1][0][0][0])*(eff[1][0][1][1]-eff[1][0][1][0]);
    double eff_mc_BtoF=eff[0][1][0][1]*eff[0][1][1][1]-(eff[0][1][0][1]-eff[0][1][0][0])*(eff[0][1][1][1]-eff[0][1][1][0]);
    double eff_mc_GtoH=eff[1][1][0][1]*eff[1][1][1][1]-(eff[1][1][0][1]-eff[1][1][0][0])*(eff[1][1][1][1]-eff[1][1][1][0]);
    return (eff_data_BtoF*WeightBtoF+eff_data_GtoH*WeightGtoH)/(eff_mc_BtoF*WeightBtoF+eff_mc_GtoH*WeightGtoH);
  }
}
void SMPAnalyzerCore::SetupZptWeight(){
  cout<<"[SMPAnalyzerCore::SetupZptWeight] setting zptcor"<<endl;
  TString datapath=getenv("DATA_DIR");
  ifstream file_check(datapath+"/"+TString::Itoa(DataYear,10)+"/Zpt/ZptWeight.root");
  bool isexist=file_check.is_open();
  file_check.close();
  if(!isexist){
    cout<<"[SMPAnalyzerCore::SetupZptWeight] no ZptWeight.root"<<endl;
    return;
  }
  TFile fzpt(datapath+"/"+TString::Itoa(DataYear,10)+"/Zpt/ZptWeight.root");
  for(const auto&& key:*(fzpt.GetListOfKeys())){
    TH2D* this_hist=(TH2D*)((TKey*)key)->ReadObj();
    TString histname=this_hist->GetName();
    TString striphistname=histname;
    int extent=0,start=striphistname.Index(TRegexp("_iter[0-9]*$"),&extent);
    if(start>=0) striphistname.Replace(start,extent,"");
    cout<<"[SMPAnalyzerCore::SetupZptWeight] setting "<<histname<<endl;
    if(histname.Contains(TRegexp("_iter0$"))){
      map_hist_zpt[striphistname]=this_hist;
      this_hist->SetDirectory(0);
    }else if(histname.Contains(TRegexp("_iter[0-9]*$"))){
      map_hist_zpt[striphistname]->Multiply(this_hist);
    }else if(histname.Contains(TRegexp("_norm$"))){
      map_hist_zpt[histname]=this_hist;
      this_hist->SetDirectory(0);
    }
  }
}
void SMPAnalyzerCore::SetupRoccoR(){
  cout<<"[SMPAnalyzerCore::SetupRoccoR] setting Rocheseter Correction"<<endl;
  TString datapath=getenv("DATA_DIR");
  roc=new RoccoR((datapath+"/"+TString::Itoa(DataYear,10)+"/RoccoR/RoccoR"+TString::Itoa(DataYear,10)+".txt").Data());
  rocele=new RocelecoR((datapath+"/"+TString::Itoa(DataYear,10)+"/RoccoR/RocelecoR"+TString::Itoa(DataYear,10)+"_new.txt").Data());
}
void SMPAnalyzerCore::SetupZ0Weight(){
  cout<<"[SMPAnalyzerCore::SetupZ0Weight] setting Z0Weight"<<endl;
  TString datapath=getenv("DATA_DIR");
  TFile fz0(datapath+"/"+TString::Itoa(DataYear,10)+"/Z0/Z0Weight.root");
  hz0=(TH1D*)fz0.Get("z0weight");
  if(hz0) hz0->SetDirectory(0);
  fz0.Close();
}
double SMPAnalyzerCore::GetZ0Weight(double valx){
  double xmin=hz0->GetXaxis()->GetXmin();
  double xmax=hz0->GetXaxis()->GetXmax();
  if(xmin>=0) valx=fabs(valx);
  if(valx<xmin) valx=xmin+0.001;
  if(valx>xmax) valx=xmax-0.001;
  return hz0->GetBinContent(hz0->FindBin(valx));
}
double SMPAnalyzerCore::GetZptWeight(double zpt,double zrap,Lepton::Flavour flavour){
  double valzptcor=1.;
  double valzptcor_norm=1.;
  zrap=fabs(zrap);
  TString sflavour=flavour==Lepton::MUON?"muon":"electron";
  TString MCName=MCSample;
  MCName=Replace(MCName,"DY[0-9]Jets","DYJets");
  TString hzptname=MCName+"_"+sflavour;
  TString hnormname=hzptname+"_norm";
  auto it=map_hist_zpt.find(hzptname);
  if(it!=map_hist_zpt.end()){
    valzptcor*=GetBinContentUser(map_hist_zpt[hzptname],zpt,zrap,0);
  }
  it=map_hist_zpt.find(hnormname);
  if(it!=map_hist_zpt.end()){
    valzptcor_norm*=GetBinContentUser(map_hist_zpt[hnormname],zpt,zrap,0);
  }
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
  if(id=="passMediumID_Selective"){
    std::vector<Electron> electrons = GetAllElectrons();
    for(unsigned int i=0; i<electrons.size(); i++){
      Electron this_electron= electrons.at(i);
      if(!( this_electron.Pt()>ptmin ))	continue;
      if(!( fabs(this_electron.scEta())<fetamax )) continue;
      if(!( this_electron.PassID("passMediumID") ))	continue;
      if(!electron_isGsfCtfScPixChargeConsistent->at(i)) continue;
      out.push_back(this_electron);
    }
  }else if(id=="passTightID_Selective"){
    std::vector<Electron> electrons = GetAllElectrons();
    for(unsigned int i=0; i<electrons.size(); i++){
      Electron this_electron= electrons.at(i);
      if(!( this_electron.Pt()>ptmin ))	continue;
      if(!( fabs(this_electron.scEta())<fetamax )) continue;
      if(!( this_electron.PassID("passTightID") )) continue;
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
  std::sort(out.begin(),out.end(),PtComparing);
  return out;
}

std::vector<Muon> SMPAnalyzerCore::MuonMomentumCorrection(const vector<Muon>& muons,int sys,int set,int member){
  std::vector<Muon> out;
  vector<Gen> gens;
  if(!IsData) gens=GetGens();
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
  std::vector<Electron> out;
  vector<Gen> gens;
  if(!IsData) gens=GetGens();
  for(auto electron:electrons){
    double rc=1.;
    double rcerr=0.;
    if(set>=0){
      if(IsDATA){
	rc=rocele->kScaleDT(electron.Charge(),electron.E(),electron.scEta(),electron.Phi(),set,member);
	//rcerr=rocele->kScaleDTerror(electron.Charge(),electron.E(),electron.scEta(),electron.Phi());
	electron.SetPtEtaPhiM(electron.Pt()*rc,electron.Eta(),electron.Phi(),electron.M());
      }else{
	Gen gen=GetGenMatchedLepton(electron,gens);
	if(gen.IsEmpty()){
	  rc=rocele->kScaleMC(electron.Charge(),electron.UncorrPt(),electron.scEta(),electron.Phi(),set,member);
	  //rcerr=rocele->kScaleMCerror(electron.Charge(),electron.UncorrPt(),electron.scEta(),electron.Phi(),set,member);	  
	}else{
	  rc=rocele->kSpreadMC(electron.Charge(),electron.UncorrPt(),electron.scEta(),electron.Phi(),gen.Pt(),set,member);
	  //rcerr=rocele->kSpreadMCerror(electron.Charge(),electron.UncorrPt(),electron.scEta(),electron.Phi(),gen.Pt(),set,member);
	}
	electron.SetPtEtaPhiM(electron.UncorrPt()*rc,electron.Eta(),electron.Phi(),electron.M());
      }      
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
