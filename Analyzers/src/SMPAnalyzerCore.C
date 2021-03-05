#include "SMPAnalyzerCore.h"

SMPAnalyzerCore::SMPAnalyzerCore(){}
SMPAnalyzerCore::~SMPAnalyzerCore(){
  for(auto& iter:map_hist_zpt){
    if(iter.second) delete iter.second;
  }
  if(roc) delete roc;
  if(rocele) delete rocele;
  if(hz0_data) delete hz0_data;
  if(hz0_mc) delete hz0_mc;
  for(std::map< TString, TH4D* >::iterator mapit = maphist_TH4D.begin(); mapit!=maphist_TH4D.end(); mapit++){
    delete mapit->second;
  }
  maphist_TH4D.clear();
}

void SMPAnalyzerCore::initializeAnalyzer(){
  if(MaxEvent>0) reductionweight=1.*fChain->GetEntries()/MaxEvent;
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
bool SMPAnalyzerCore::IsExists(TString filepath){
  ifstream fcheck(filepath);
  return fcheck.good();
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
  }else if(histkey.Contains("_Q_")){
    if(lep->Charge()>0) histkey.ReplaceAll("_Q_","_QPlus_");
    else histkey.ReplaceAll("_Q_","_QMinus_");;
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
	double data_eff=WeightBtoF*Lepton_SF(histkey_data+"_BCDEF",lep,sys)+WeightGtoH*Lepton_SF(histkey_data+"_GH",lep,sys);
	double mc_eff=WeightBtoF*Lepton_SF(histkey_mc+"_BCDEF",lep,-sys)+WeightGtoH*Lepton_SF(histkey_mc+"_GH",lep,-sys);
	if(mc_eff==0) return 1;
	else return data_eff/mc_eff;
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
    exit(ENODATA);
  }    
  double this_x,this_y;
  if(this_hist->GetXaxis()->GetXmax()>this_hist->GetYaxis()->GetXmax()){
    if(histkey.Contains("_Eff_") && this_pt<this_hist->GetXaxis()->GetXmin()) return 0;
    this_x=this_pt;
    this_y=this_eta;
  }else{
    if(histkey.Contains("_Eff_") && this_pt<this_hist->GetYaxis()->GetXmin()) return 0;
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
  if(mc_eff==0) return 1.;
  else return data_eff/mc_eff;
}
double SMPAnalyzerCore::DileptonTrigger_SF(TString triggerSF_key0,TString triggerSF_key1,const vector<Lepton*>& leps,int sys){
  if(IsDATA) return 1;
  if((triggerSF_key0==""||triggerSF_key0=="Default")&&(triggerSF_key1==""||triggerSF_key1=="Default")) return 1;
  int nlep=leps.size();
  if(nlep<2){
    cout<<"[SMPAnalyzerCore::DileptonTrigger_SF] nlep < 2. return 1."<<endl;
    return 1.;
  }
  double data_noleg1=1.,mc_noleg1=1.;
  vector<double> data_oneleg1_noleg2(nlep,1.);
  vector<double> mc_oneleg1_noleg2(nlep,1.);
  for(int i=0;i<nlep;i++){
    double data_eff_leg1=Lepton_SF("Trigger_Eff_DATA_"+triggerSF_key0,leps.at(i),sys);
    double data_eff_leg2=Lepton_SF("Trigger_Eff_DATA_"+triggerSF_key1,leps.at(i),sys);
    double mc_eff_leg1=Lepton_SF("Trigger_Eff_MC_"+triggerSF_key0,leps.at(i),-sys);
    double mc_eff_leg2=Lepton_SF("Trigger_Eff_MC_"+triggerSF_key1,leps.at(i),-sys);
    data_noleg1*=(1-data_eff_leg1);    
    mc_noleg1*=(1-mc_eff_leg1);
    for(int j=0;j<nlep;j++){
      if(i==j){
	data_oneleg1_noleg2[j]*=data_eff_leg1;
	mc_oneleg1_noleg2[j]*=mc_eff_leg1;
      }else{
	data_oneleg1_noleg2[j]*=(1-data_eff_leg2);
	mc_oneleg1_noleg2[j]*=(1-mc_eff_leg2);
      }
    }
  }
  double data_eff=1.-data_noleg1;
  double mc_eff=1.-mc_noleg1;
  for(int i=0;i<nlep;i++){
    data_eff-=data_oneleg1_noleg2[i];
    mc_eff-=mc_oneleg1_noleg2[i];
  }
  if(mc_eff==0) return 1.;
  else return data_eff/mc_eff;
}
void SMPAnalyzerCore::SetupZptWeight(){
  /// TODO: currently, temp Zpt weight from 2017 single muon for all other channels
  cout<<"[SMPAnalyzerCore::SetupZptWeight] setting zptcor"<<endl;
  TString datapath=getenv("DATA_DIR");
  if(!IsExists(datapath+"/"+GetEra()+"/Zpt/ZptWeight.root")){
    cout<<"[SMPAnalyzerCore::SetupZptWeight] no ZptWeight.root"<<endl;
    return;
  }
  TFile fzpt(datapath+"/"+GetEra()+"/Zpt/ZptWeight.root");
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
  TString rocpath=datapath+"/"+GetEra()+"/RoccoR/RoccoR"+GetEraShort()+"UL.txt";
  if(IsExists(rocpath)) roc=new RoccoR(rocpath.Data());
  else cout<<"[SMPAnalyzerCore::SetupRoccoR] no "+rocpath<<endl;
  TString rocelepath=datapath+"/"+GetEra()+"/RoccoR/RocelecoR"+GetEraShort()+"_new.txt";
  if(IsExists(rocelepath)) rocele=new RocelecoR(rocelepath.Data());
  else cout<<"[SMPAnalyzerCore::SetupRoccoR] no "+rocelepath<<endl;  
}
void SMPAnalyzerCore::SetupZ0Weight(){
  cout<<"[SMPAnalyzerCore::SetupZ0Weight] setting Z0Weight"<<endl;
  TString datapath=getenv("DATA_DIR");
  if(!IsExists(datapath+"/"+GetEra()+"/Z0/Z0Weight.root")){
    cout<<"[SMPAnalyzerCore::SetupZ0Weight] no Z0Weight.root"<<endl;
    return;
  }
  TFile fz0(datapath+"/"+GetEra()+"/Z0/Z0Weight.root");
  hz0_data=(TF1*)fz0.Get("data_fit");
  hz0_mc=(TF1*)fz0.Get("mc_fit");
  fz0.Close();
}
double SMPAnalyzerCore::GetZ0Weight(double valx){
  return 1.; /// FIXME: no Z0 weight for UL
  double data_val = hz0_data->Eval(valx);
  double mc_val = hz0_mc->Eval(valx);
  double norm = hz0_mc->Integral(-100,100)/hz0_data->Integral(-100,100);
  return norm*data_val/mc_val;
}
double SMPAnalyzerCore::GetZptWeight(double zpt,double zrap,Lepton::Flavour flavour){
  double valzptcor=1.;
  double valzptcor_norm=1.;
  zrap=fabs(zrap);
  TString sflavour=flavour==Lepton::MUON?"muon":"electron";
  TString MCName=MCSample;
  if(MCName.Contains(TRegexp("^DY[0-9]Jets$"))) MCName="DYJets";
  if(MCName.Contains(TRegexp("^DYJets_Pt-[0-9]*To[0-9Inf]*$"))) MCName="DYJets";
  if(MCName.Contains(TRegexp("^DYJets_M-[0-9]*to[0-9Inf]*$"))) MCName="DYJets";
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
void SMPAnalyzerCore::GetEventWeights(){
  lumiweight=1;
  PUweight=1;
  PUweight_up=1;
  PUweight_down=1;
  prefireweight=1;
  prefireweight_up=1;
  prefireweight_down=1;
  zptweight=1;
  tauprefix="";
  z0weight=1;
  if(!IsDATA){
    lumiweight=weight_norm_1invpb*event.MCweight()*event.GetTriggerLumi("Full")*reductionweight;
    PUweight=mcCorr->GetPileUpWeight(nPileUp,0);
    PUweight_up=mcCorr->GetPileUpWeight(nPileUp,1);
    PUweight_down=mcCorr->GetPileUpWeight(nPileUp,-1);
    if(DataYear<2018){
      prefireweight=L1PrefireReweight_Central;
      prefireweight_up=L1PrefireReweight_Up;
      prefireweight_down=L1PrefireReweight_Down;
    }
    if(IsDYSample){
      vector<Gen> gens=GetGens();
      Gen parton0,parton1,l0,l1;
      GetDYGenParticles(gens,parton0,parton1,l0,l1,3);
      if(abs(l0.PID())==11||abs(l0.PID())==13){
	TLorentzVector genZ=(l0+l1);
	zptweight=GetZptWeight(genZ.Pt(),genZ.Rapidity(),abs(l0.PID())==13?Lepton::Flavour::MUON:Lepton::Flavour::ELECTRON);
      }else tauprefix="tau_";
    }
    z0weight=GetZ0Weight(vertex_Z);
  }else{
    lumiweight=reductionweight;
  }
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


void SMPAnalyzerCore::GetDYGenParticles(const vector<Gen>& gens,Gen& parton0,Gen& parton1,Gen& l0,Gen& l1,int mode){
  //mode 0:bare 1:dressed01 2:dressed04 3:beforeFSR
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
  if(mode>=3){
    if(nlepton>=4){
      for(int i=0;i<nlepton;i++){
	if(leptons[i]->Index()==l0.Index()||leptons[i]->Index()==l1.Index()) continue;
	for(int j=i+1;j<nlepton;j++){
	  if(leptons[j]->Index()==l0.Index()||leptons[j]->Index()==l1.Index()) continue;
	  if(!(leptons[i]->PID()+leptons[j]->PID()==0)) continue;
	  vector<int> history_i=TrackGenSelfHistory(*leptons[i],gens);
	  vector<int> history_j=TrackGenSelfHistory(*leptons[j],gens);
	  if(history_i.at(1)==history_j.at(1)) photons.push_back(&gens[history_i.at(1)]);
	}
      }
    }	
    for(const auto& photon:photons){
      vector<int> history=TrackGenSelfHistory(*photon,gens);
      if(gens[history.at(1)].PID()==l0.PID()) l0+=*photon;
      else if(gens[history.at(1)].PID()==l1.PID()) l1+=*photon;
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
      //if(muon.TrkIso()/muon.Pt()<0.1) out.push_back(muon);
      if(muon.PassSelector(Muon::Selector::TkIsoLoose)) out.push_back(muon);
    }
  }else if(id=="POGMediumWithLooseTrkIso"){
    vector<Muon> muons=GetMuons("POGMedium",ptmin,fetamax);
    for(auto const& muon: muons){
      //if(muon.TrkIso()/muon.Pt()<0.1) out.push_back(muon);
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
  std::sort(out.begin(),out.end(),PtComparing);
  return out;
}

std::vector<Muon> SMPAnalyzerCore::MuonMomentumCorrection(const vector<Muon>& muons,int sys,int set,int member){
  if(!roc) return std::vector<Muon>(muons);
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
  if(!rocele) return std::vector<Electron>(electrons);
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

SMPAnalyzerCore::Parameter::Parameter(){
  electronIDSF="ID_SF_MediumID_Q";
  muonIDSF="IDISO_SF_MediumID_trkIsoLoose_Q";
}
SMPAnalyzerCore::Parameter::~Parameter(){
}
SMPAnalyzerCore::Parameter::Parameter(TString pre,TString suf,TString elID,vector<TString> Trig,double l0ptcut,double l1ptcut,vector<Lepton*> leps_){
  prefix=pre;
  suffix=suf;
  electronIDSF=elID;
  triggerSF=Trig;
  if(l0ptcut>0) lep0ptcut=l0ptcut;
  if(l1ptcut>0) lep1ptcut=l1ptcut;
  if(leps_.size()) leps=leps_;
}
SMPAnalyzerCore::Parameter::Parameter(TString elID,vector<TString> Trig,double l0ptcut,double l1ptcut,vector<Lepton*> leps_){
  electronIDSF=elID;
  triggerSF=Trig;
  if(l0ptcut>0) lep0ptcut=l0ptcut;
  if(l1ptcut>0) lep1ptcut=l1ptcut;
  if(leps_.size()) leps=leps_;
}
SMPAnalyzerCore::Parameter::Parameter(TString pre,TString suf,TString muID,TString muISO,vector<TString> Trig,double l0ptcut,double l1ptcut,vector<Lepton*> leps_){
  prefix=pre;
  suffix=suf;
  muonIDSF=muID;
  muonISOSF=muISO;
  triggerSF=Trig;
  if(l0ptcut>0) lep0ptcut=l0ptcut;
  if(l1ptcut>0) lep1ptcut=l1ptcut;
  if(leps_.size()) leps=leps_;
}
SMPAnalyzerCore::Parameter::Parameter(TString muID,TString muISO,vector<TString> Trig,double l0ptcut,double l1ptcut,vector<Lepton*> leps_){
  muonIDSF=muID;
  muonISOSF=muISO;
  triggerSF=Trig;
  if(l0ptcut>0) lep0ptcut=l0ptcut;
  if(l1ptcut>0) lep1ptcut=l1ptcut;
  if(leps_.size()) leps=leps_;
}
SMPAnalyzerCore::Parameter SMPAnalyzerCore::Parameter::Clone(vector<Lepton*> leps_,int weightbit_) const {
  Parameter out=*this;
  out.leps=leps_;
  if(weightbit_>=0) out.weightbit=weightbit_;
  return out;
}
