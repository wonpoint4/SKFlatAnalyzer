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
  for(std::map< TString, TH4D* >::iterator mapit = maphist_TH4D.begin(); mapit!=maphist_TH4D.end(); mapit++){
    delete mapit->second;
  }
  maphist_TH4D.clear();
}

void SMPAnalyzerCore::initializeAnalyzer(){
  SetupZptWeight();
  SetupRoccoR();
  SetupZ0Weight();
  IsDYSample=false;
  if(MCSample.Contains("DYJets")||MCSample.Contains("ZToEE")||MCSample.Contains("ZToMuMu")||MCSample.Contains(TRegexp("DY[0-9]Jets"))) IsDYSample=true;
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
                            Int_t n_binx, Double_t *xbins,
                            Int_t n_biny, Double_t *ybins,
                            Int_t n_binz, Double_t *zbins,
                            Int_t n_binu, Double_t *ubins){

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
                            Int_t n_binx, Double_t *xbins,
                            Int_t n_biny, Double_t *ybins,
                            Int_t n_binz, Double_t *zbins,
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
void SMPAnalyzerCore::FillHist(TString histname, double value_x, double value_y, double value_z, double value_u, map<TString,double> weights, int n_binx, double x_min, double x_max, int n_biny, double y_min, double y_max, int n_binz, double z_min, double z_max, int n_binu, double u_min, double u_max){
  for(const auto& [suffix,weight]:weights) FillHist(histname+suffix,value_x,value_y,value_z,value_u,weight,n_binx,x_min,x_max,n_biny,y_min,y_max,n_binz,z_min,z_max,n_binu,u_min,u_max);
}
void SMPAnalyzerCore::FillHist(TString histname, double value_x, double value_y, double value_z, double value_u, map<TString,double> weights, int n_binx, double *xbins, int n_biny, double *ybins, int n_binz, double *zbins, int n_binu, double *ubins){
  for(const auto& [suffix,weight]:weights) FillHist(histname+suffix,value_x,value_y,value_z,value_u,weight,n_binx,xbins,n_biny,ybins,n_binz,zbins,n_binu,ubins);
}
void SMPAnalyzerCore::FillHist(TString histname, double value_x, double value_y, double value_z, double value_u, map<TString,double> weights, int n_binx, double *xbins, int n_biny, double *ybins, int n_binz, double *zbins, int n_binu, double u_min, double u_max){
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
    if(!this_hist && DataYear==2016 && !histkey.Contains("_BCDEF$") && !histkey.Contains("_GH$")){
      double lumi_periodB = 5750.490644035;
      double lumi_periodC = 2572.903488748;
      double lumi_periodD = 4242.291556970;
      double lumi_periodE = 4025.228136967;
      double lumi_periodF = 3104.509131800;
      double lumi_periodG = 7575.824256098;
      double lumi_periodH = 8650.628380028;
      double total_lumi = (lumi_periodB+lumi_periodC+lumi_periodD+lumi_periodE+lumi_periodF+lumi_periodG+lumi_periodH);
      
      double WeightBtoF = (lumi_periodB+lumi_periodC+lumi_periodD+lumi_periodE+lumi_periodF)/total_lumi;
      double WeightGtoH = (lumi_periodG+lumi_periodH)/total_lumi;

      if(histkey.Contains("_SF_")){
	TString histkey_data=histkey;
	histkey_data.ReplaceAll("_SF_","_Eff_DATA_");
	TString histkey_mc=histkey;
	histkey_mc.ReplaceAll("_SF_","_Eff_MC_");
	return (WeightBtoF*Lepton_SF(histkey_data+"_BCDEF",lep,sys)+WeightGtoH*Lepton_SF(histkey_data+"_GH",lep,sys))/(WeightBtoF*Lepton_SF(histkey_mc+"_BCDEF",lep,-sys)+WeightGtoH*Lepton_SF(histkey_mc+"_GH",lep,-sys));
      }else if(histkey.Contains("_Eff_")){
	return WeightBtoF*Lepton_SF(histkey+"_BCDEF",lep,sys)+WeightGtoH*Lepton_SF(histkey+"_GH",lep,sys);
      }
    }
  }else if(lep->LeptonFlavour()==Lepton::ELECTRON){
    this_pt=((Electron*)lep)->UncorrPt();
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

  double data_eff=1.,mc_eff=1.;
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
  double eff[2][2][2]={}; //[data/mc][l0/l1][leg1/leg2]
  TString sdata[2]={"DATA","MC"};
  for(int id=0;id<2;id++)
    for(int ilep=0;ilep<2;ilep++)
      for(int ileg=0;ileg<2;ileg++)
	eff[id][ilep][ileg]=Lepton_SF("Trigger_Eff_"+sdata[id]+"_"+histkeys[ileg],leps.at(ilep),(id?-1.:1.)*sys);
  double eff_data=eff[0][0][1]*eff[0][1][1]-(eff[0][0][1]-eff[0][0][0])*(eff[0][1][1]-eff[0][1][0]);
  double eff_mc=eff[1][0][1]*eff[1][1][1]-(eff[1][0][1]-eff[1][0][0])*(eff[1][1][1]-eff[1][1][0]);
  return eff_data/eff_mc;
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
    gens[i].Print();
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

void SMPAnalyzerCore::GetDYLHEParticles(const vector<LHE>& lhes,LHE& l0,LHE& l1){
  if(!IsDYSample){
    cout <<"[AFBAnalyzer::GetDYLHEParticles] this is for DY event"<<endl;
    exit(EXIT_FAILURE);
  }
  for(int i=0;i<(int)lhes.size();i++){
    if(l0.ID()==0&&(abs(lhes[i].ID())==11||abs(lhes[i].ID())==13||abs(lhes[i].ID())==15)) l0=lhes[i];
    if(l0.ID()&&lhes[i].ID()==-l0.ID()) l1=lhes[i];
  }
  if(l0.ID()==0||l1.ID()==0){
    cout <<"[AFBAnalyzer::GetLHEParticles] something is wrong"<<endl;
    exit(EXIT_FAILURE);
  }
  if(l0.Pt()<l1.Pt()){
    LHE temp=l0;
    l0=l1;
    l1=temp;
  }
}


void SMPAnalyzerCore::GetDYGenParticles(const vector<Gen>& gens,Gen& parton0,Gen& parton1,Gen& l0,Gen& l1,bool dressed){
  if(!IsDYSample){
    cout <<"[SMPAnalyzerCore::GetDYGenParticles] this is for DY event"<<endl;
    exit(EXIT_FAILURE);
  }
  vector<const Gen*> leptons;
  vector<const Gen*> photons;
  int ngen=gens.size();
  for(int i=0;i<ngen;i++){
    if(!gens.at(i).isPrompt()) continue;
    int genpid=gens.at(i).PID();
    if(gens.at(i).isHardProcess()){
      if(abs(genpid)<7||genpid==21){
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
  if(dressed){
    int nphoton=photons.size();
    for(int i=0;i<nphoton;i++){
      if(l0.DeltaR(*photons[i])>0.4&&l1.DeltaR(*photons[i])>0.4) continue;
      if(l0.DeltaR(*photons[i])<l1.DeltaR(*photons[i])) l0+=*photons[i];
      else l1+=*photons[i];
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
  }else if(id=="POGMediumWithLooseTrkIso"){
    vector<Muon> muons=GetMuons("POGMedium",ptmin,fetamax);
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
    if(set>=0){
      double rc=1.;
      //double rcerr=0.;
      double el_eta=electron.scEta();
      double el_phi=electron.Phi();
      // FIXME: rocelecor don't have the corrections for 2.4<|eta|<2.5
      if(el_eta>=2.4) el_eta=2.39;
      if(el_eta<=-2.4) el_eta=-2.39;
      if(fabs(el_eta)<2.4){
	if(IsDATA){
	  rc=rocele->kScaleDT(electron.Charge(),electron.E(),el_eta,el_phi,set,member);
	  //rcerr=rocele->kScaleDTerror(electron.Charge(),electron.E(),el_eta,el_phi);
	  electron*=rc;
	}else{	
	  Gen gen=SMPGetGenMatchedLepton(electron,gens,1);
	  if(gen.IsEmpty()){
	    rc=rocele->kScaleMC(electron.Charge(),electron.UncorrPt(),el_eta,el_phi,set,member);
	    //rcerr=rocele->kScaleMCerror(electron.Charge(),electron.UncorrPt(),el_eta,el_phi,set,member);	  
	  }else{
	    rc=rocele->kSpreadMC(electron.Charge(),electron.UncorrPt(),el_eta,el_phi,gen.Pt(),set,member);
	    //rcerr=rocele->kSpreadMCerror(electron.Charge(),electron.UncorrPt(),el_eta,el_phi,gen.Pt(),set,member);
	  }
	  electron*=rc*electron.UncorrE()/electron.E();
	}      
      }
    }else if(set==-1){ //no energe cor
      electron*=electron.UncorrE()/electron.E();
    }else{
      cout<<"[SMPAnalyzerCore::ElectronEnergyCorrection] wrong set "<<set<<endl;
      exit(EXIT_FAILURE);
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
