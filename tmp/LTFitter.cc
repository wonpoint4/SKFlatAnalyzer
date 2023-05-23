TH2* response;
TH1* hreco;
map<TString,TH1*> hists;
const int nabin=LTAnalyzer::nptbin*LTAnalyzer::njetbin;
double As_nominal[nabin*2];
TString outfilename;
void NormalizeReco(TH1* hist){
  for(int i=0;i<nabin;i++){
    double istart=LTAnalyzer::ncostbin*LTAnalyzer::nphibin*i;
    double iend=LTAnalyzer::ncostbin*LTAnalyzer::nphibin*(i+1);
    double sum=hist->Integral(1+istart,iend);
    if(sum==0) continue;
    for(int j=0;j<LTAnalyzer::ncostbin*LTAnalyzer::nphibin;j++){
      double val=hist->GetBinContent(1+istart+j);
      double err=hist->GetBinError(1+istart+j);
      hist->SetBinContent(1+istart+j,val/sum);
      hist->SetBinError(1+istart+j,err/sum);
    }
  }
}
void Init(TString sim="mi",TString data="mg"){
  TH1::AddDirectory(0);
  TH1::SetDefaultSumw2(true);
  TFile* f;
  if(sim=="mi") f=TFile::Open("/data6/Users/hsseo/SKFlatOutput//Run2UltraLegacy_v3/LTAnalyzer/2018/LTAnalyzer_DYJetsToMuMu_MiNNLO.root");
  else if(sim=="mg") f=TFile::Open("/data6/Users/hsseo/SKFlatOutput//Run2UltraLegacy_v3/LTAnalyzer/2018/LTAnalyzer_DYJets_MG.root");
  response=(TH2*)f->Get("mm2018/response");
  f->Close();
  
  if(data=="mi"){
    f=TFile::Open("/data6/Users/hsseo/SKFlatOutput//Run2UltraLegacy_v3/LTAnalyzer/2018/LTAnalyzer_DYJetsToMuMu_MiNNLO.root");
    outfilename="mi.root";
  }else if(data=="mg"){
    f=TFile::Open("/data6/Users/hsseo/SKFlatOutput//Run2UltraLegacy_v3/LTAnalyzer/2018/LTAnalyzer_DYJets_MG.root");
    outfilename="mg.root";
  }
  hreco=(TH1*)f->Get("mm2018/reco");
  NormalizeReco(hreco);
  if(f->Get("mm2018/response")){
    TH2* true_response=(TH2*)f->Get("mm2018/response");
    TH1* true_gen=true_response->ProjectionX();
    TH1D* gen_j0_a0=new TH1D("gen_j0_a0","gen_j0_a0",LTAnalyzer::nptbin,LTAnalyzer::ptbins);
    TH1D* gen_j0_a2=new TH1D("gen_j0_a2","gen_j0_a2",LTAnalyzer::nptbin,LTAnalyzer::ptbins);
    TH1D* gen_j1_a0=new TH1D("gen_j1_a0","gen_j1_a0",LTAnalyzer::nptbin,LTAnalyzer::ptbins);
    TH1D* gen_j1_a2=new TH1D("gen_j1_a2","gen_j1_a2",LTAnalyzer::nptbin,LTAnalyzer::ptbins);
    TF2* func=new TF2("func","[0]*(1+x*x+0.5*[1]*(1-3*x*x)+0.5*[2]*(1-x*x)*cos(2*y))",0,1,0,TMath::Pi()/2);
    for(int i=0;i<nabin;i++){
      TH2D* hist=new TH2D("costphi","costphi",LTAnalyzer::ncostbin,0,1,LTAnalyzer::nphibin,0,TMath::Pi()/2);
      for(int icost=0;icost<LTAnalyzer::ncostbin;icost++){
	for(int iphi=0;iphi<LTAnalyzer::nphibin;iphi++){
	  hist->SetBinContent(1+icost,1+iphi,true_gen->GetBinContent(1+LTAnalyzer::nphibin*(LTAnalyzer::ncostbin*i+icost)+iphi));
	  hist->SetBinError(1+icost,1+iphi,true_gen->GetBinError(1+LTAnalyzer::nphibin*(LTAnalyzer::ncostbin*i+icost)+iphi));
	}
      }
      hist->Fit(func);
      if(i<LTAnalyzer::nptbin){
	gen_j0_a0->SetBinContent(i+1,func->GetParameter(1));
	gen_j0_a0->SetBinError(i+1,func->GetParError(1));
	gen_j0_a2->SetBinContent(i+1,func->GetParameter(2));
	gen_j0_a2->SetBinError(i+1,func->GetParError(2));
      }else{
	gen_j1_a0->SetBinContent(i%LTAnalyzer::nptbin+1,func->GetParameter(1));
	gen_j1_a0->SetBinError(i%LTAnalyzer::nptbin+1,func->GetParError(1));
	gen_j1_a2->SetBinContent(i%LTAnalyzer::nptbin+1,func->GetParameter(2));
	gen_j1_a2->SetBinError(i%LTAnalyzer::nptbin+1,func->GetParError(2));
      }
      delete hist;
    }
    hists["gen_j0_a0"]=gen_j0_a0;
    hists["gen_j0_a2"]=gen_j0_a2;
    hists["gen_j1_a0"]=gen_j1_a0;
    hists["gen_j1_a2"]=gen_j1_a2;
    delete true_response;
  }

  TH1* response_gen=response->ProjectionX();
  TF2* func=new TF2("func","[0]*(1+x*x+0.5*[1]*(1-3*x*x)+0.5*[2]*(1-x*x)*cos(2*y))",0,1,0,TMath::Pi()/2);
  for(int i=0;i<nabin;i++){
    TH2D* hist=new TH2D("costphi","costphi",LTAnalyzer::ncostbin,0,1,LTAnalyzer::nphibin,0,TMath::Pi()/2);
    for(int icost=0;icost<LTAnalyzer::ncostbin;icost++){
      for(int iphi=0;iphi<LTAnalyzer::nphibin;iphi++){
	hist->SetBinContent(1+icost,1+iphi,response_gen->GetBinContent(1+LTAnalyzer::nphibin*(LTAnalyzer::ncostbin*i+icost)+iphi));
	hist->SetBinError(1+icost,1+iphi,response_gen->GetBinError(1+LTAnalyzer::nphibin*(LTAnalyzer::ncostbin*i+icost)+iphi));
      }
    }
    hist->Fit(func);
    As_nominal[i]=func->GetParameter(1);
    As_nominal[nabin+i]=func->GetParameter(2);
    delete hist;
  }
  cout<<"Init"<<endl;
  for(int i=0;i<nabin;i++){
    cout<<i<<" "<<As_nominal[i]<<" "<<As_nominal[nabin+i]<<endl;
  }
}
double AngularDistribution(double cost,double phi,double A0,double A2){
  return 1+cost*cost+0.5*A0*(1-3*cost*cost)+0.5*A2*(1-cost*cost)*cos(2*phi);
}
double GetChi2(const double* As){
  const double* A0s=As;
  const double* A2s=As+nabin;
    
  TH2* response_reweighted=(TH2*)response->Clone();
  for(int i=0;i<2*LTAnalyzer::nptbin*LTAnalyzer::ncostbin*LTAnalyzer::nphibin;i++){
    int remain=i;
    int iphi=remain%LTAnalyzer::nphibin;
    remain=remain/LTAnalyzer::nphibin;
    double phi=TMath::Pi()/2*(1.0*iphi/LTAnalyzer::nphibin+0.5/LTAnalyzer::nphibin);
    int icost=remain%LTAnalyzer::ncostbin;
    remain=remain/LTAnalyzer::ncostbin;
    double cost=1.0*icost/LTAnalyzer::ncostbin+0.5/LTAnalyzer::ncostbin;
    int ipt=remain%LTAnalyzer::nptbin;
    remain=remain/LTAnalyzer::nptbin;
    int ijet=remain%LTAnalyzer::njetbin;
    double A0=A0s[LTAnalyzer::nptbin*ijet+ipt];
    double A2=A2s[LTAnalyzer::nptbin*ijet+ipt];
    double A0_nominal=As_nominal[LTAnalyzer::nptbin*ijet+ipt];
    double A2_nominal=As_nominal[nabin+LTAnalyzer::nptbin*ijet+ipt];
    double weight=AngularDistribution(cost,phi,A0,A2)/AngularDistribution(cost,phi,A0_nominal,A2_nominal);
    //cout<<"weight: "<<weight<<endl;
    for(int j=0;j<2*LTAnalyzer::nptbin*LTAnalyzer::ncostbin*LTAnalyzer::nphibin;j++){
      double val=response_reweighted->GetBinContent(i+1,j+1);
      double err=response_reweighted->GetBinError(i+1,j+1);
      response_reweighted->SetBinContent(i+1,j+1,val*weight);
      response_reweighted->SetBinError(i+1,j+1,err*weight);
    }
  }
  TH1* response_reco=response_reweighted->ProjectionY();
  NormalizeReco(response_reco);
  double chi2=0;
  for(int i=0;i<2*LTAnalyzer::nptbin*LTAnalyzer::ncostbin*LTAnalyzer::nphibin;i++){
    double val0=hreco->GetBinContent(i+1);
    double err0=hreco->GetBinError(i+1);
    double val1=response_reco->GetBinContent(i+1);
    double err1=response_reco->GetBinError(i+1);
    double denom=err0*err0+err1*err1;
    //cout<<i<<" "<<pow(val1-val0,2)/denom<<endl;
    if(denom==0) continue;
    chi2+=pow(val1-val0,2)/denom;
  }
  delete response_reweighted;
  cout<<"chi2:"<<chi2<<endl;
  return chi2;
}
void Test(){
  TH1* response_reco=response->ProjectionY();
  NormalizeReco(response_reco);
  TCanvas* c=new TCanvas;
  hreco->Draw("hist e");
  response_reco->SetLineColor(2);
  response_reco->Draw("same hist e");
  double chi2=0;
  for(int i=0;i<hreco->GetNcells();i++){
    cout<<i<<" "<<hreco->GetBinContent(i)<<" "<<hreco->GetBinError(i)<<" "<<response_reco->GetBinContent(i)<<" "<<response_reco->GetBinError(i)<<" "<<pow(hreco->GetBinContent(i)-response_reco->GetBinContent(i),2)/(pow(hreco->GetBinError(i),2)+pow(response_reco->GetBinError(i),2))<<endl;
    chi2+=pow(hreco->GetBinContent(i)-response_reco->GetBinContent(i),2)/(pow(hreco->GetBinError(i),2)+pow(response_reco->GetBinError(i),2));
  }
  cout<<"chi2: "<<chi2<<endl;
}
void minimize(){
  ROOT::Math::Minimizer* min=ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");
  ROOT::Math::Functor f(&GetChi2,2*nabin);
  min->SetFunction(f);
  for(int i=0;i<nabin;i++){
    min->SetVariable(i,Form("A0_%d",i),As_nominal[i],0.01);
    min->SetVariable(nabin+i,Form("A2_%d",i),As_nominal[nabin+i],0.01);
  }
  min->Minimize();
  const double *xs=min->X();
  const double *errs=min->Errors();
  TH1D* unfolded_j0_a0=new TH1D("unfolded_j0_a0","unfolded_j0_a0",LTAnalyzer::nptbin,LTAnalyzer::ptbins);
  TH1D* unfolded_j0_a2=new TH1D("unfolded_j0_a2","unfolded_j0_a2",LTAnalyzer::nptbin,LTAnalyzer::ptbins);
  TH1D* unfolded_j1_a0=new TH1D("unfolded_j1_a0","unfolded_j1_a0",LTAnalyzer::nptbin,LTAnalyzer::ptbins);
  TH1D* unfolded_j1_a2=new TH1D("unfolded_j1_a2","unfolded_j1_a2",LTAnalyzer::nptbin,LTAnalyzer::ptbins);
  for(int i=0;i<nabin;i++){
    if(i<LTAnalyzer::nptbin){
      unfolded_j0_a0->SetBinContent(i+1,xs[i]);
      unfolded_j0_a0->SetBinError(i+1,errs[i]);
      unfolded_j0_a2->SetBinContent(i+1,xs[nabin+i]);
      unfolded_j0_a2->SetBinError(i+1,errs[nabin+i]);
    }else{
      unfolded_j1_a0->SetBinContent(i%LTAnalyzer::nptbin+1,xs[i]);
      unfolded_j1_a0->SetBinError(i%LTAnalyzer::nptbin+1,errs[i]);
      unfolded_j1_a2->SetBinContent(i%LTAnalyzer::nptbin+1,xs[nabin+i]);
      unfolded_j1_a2->SetBinError(i%LTAnalyzer::nptbin+1,errs[nabin+i]);
    }
  }
  hists["unfolded_j0_a0"]=unfolded_j0_a0;
  hists["unfolded_j0_a2"]=unfolded_j0_a2;
  hists["unfolded_j1_a0"]=unfolded_j1_a0;
  hists["unfolded_j1_a2"]=unfolded_j1_a2;
  for(int i=0;i<nabin;i++){
    cout<<i<<" "<<xs[i]<<" "<<xs[nabin+i]<<endl;
  }
  TFile fout(outfilename,"recreate");
  for(auto [name,hist]:hists){
    hist->Write();
  }
}
