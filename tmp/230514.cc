void plot(){
  TH1::AddDirectory(0);
  TFile f("/data6/Users/hsseo/SKFlatOutput//Run2UltraLegacy_v3/LTAnalyzer/2018/LTAnalyzer_DYJetsToMuMu_MiNNLO.root");
  //TFile f("/data6/Users/hsseo/SKFlatOutput//Run2UltraLegacy_v3/LTAnalyzer/2018/LTAnalyzer_DYJets_MG.root");
  TH1D* ha0=new TH1D("a0","a0",LTAnalyzer::nptbin,LTAnalyzer::ptbins);
  TH1D* ha2=new TH1D("a2","a2",LTAnalyzer::nptbin,LTAnalyzer::ptbins);
  TF2* func=new TF2("func","[0]*(1+x*x+0.5*[1]*(1-3*x*x)+0.5*[2]*(1-x*x)*cos(2*y))",0,1,0,TMath::Pi()/2);
  for(int i=1;i<LTAnalyzer::nptbin+2;i++){
    TH2* hist=(TH2*)f.Get(Form("mm2018/gen_pt%d_costphi",i));
    hist->Fit(func);
    ha0->SetBinContent(i,func->GetParameter(1));
    ha0->SetBinError(i,func->GetParError(1));
    ha2->SetBinContent(i,func->GetParameter(2));
    ha2->SetBinError(i,func->GetParError(2));
  }
  TCanvas* c=new TCanvas("p1");
  ha0->Draw("hist e");
  ha0->SetStats(0);
  ha2->Draw("same hist e");
  ha2->SetLineColor(2);
  TLegend *leg=new TLegend(0.89,0.12,0.6,0.5);
  leg->AddEntry(ha0,"A0");
  leg->AddEntry(ha2,"A2");
  leg->Draw();
  ha0->SetTitle("");
  ha0->GetXaxis()->SetTitle("p_{T}(ll) [GeV]");
  ha0->GetYaxis()->SetTitle("A_{i}");
  
}
void plot2(){
  TH1::AddDirectory(0);
  TFile f("/data6/Users/hsseo/SKFlatOutput//Run2UltraLegacy_v3/LTAnalyzer/2018/LTAnalyzer_DYJetsToMuMu_MiNNLO.root");
  //TFile f("/data6/Users/hsseo/SKFlatOutput//Run2UltraLegacy_v3/LTAnalyzer/2018/LTAnalyzer_DYJets_MG.root");
  TH1D* ha0=new TH1D("a0","a0",LTAnalyzer::nptbin,LTAnalyzer::ptbins);
  TH1D* ha2=new TH1D("a2","a2",LTAnalyzer::nptbin,LTAnalyzer::ptbins);
  TF2* func=new TF2("func","[0]*(1+x*x+0.5*[1]*(1-3*x*x)+0.5*[2]*(1-x*x)*cos(2*y))",0,1,0,TMath::Pi()/2);
  TF1* funcx=new TF1("funcx","[0]*(1+x*x+0.5*[1]*(1-3*x*x))",0,1);
  //TF1* funcy=new TF1("funcy","[0]*(1+[1]/4*(1-2*sin(x)*sin(x))+[2]*sin(x))",0,TMath::Pi()/2);
  TF1* funcy=new TF1("funcy","[0]*(1+[1]/4*(1-2*sin(x)*sin(x)))",0,TMath::Pi()/2);
  for(int i=1;i<LTAnalyzer::nptbin+2;i++){
    TH2* hist=(TH2*)f.Get(Form("mm2018/gen_pt%d_costphi",i));
    TH1* histx=hist->ProjectionX(Form("x%d",i));
    TH1* histy=hist->ProjectionY(Form("y%d",i));
    new TCanvas;
    histx->Fit(funcx);
    new TCanvas;
    histy->Fit(funcy);
    ha0->SetBinContent(i,funcx->GetParameter(1));
    ha0->SetBinError(i,funcx->GetParError(1));
    ha2->SetBinContent(i,funcy->GetParameter(1));
    ha2->SetBinError(i,funcy->GetParError(1));
  }
  TCanvas* c=new TCanvas("p2");
  ha0->Draw("hist e");
  ha2->Draw("same hist e");
}
void plot3(){
  TH1::AddDirectory(0);
  //TFile f("/data6/Users/hsseo/SKFlatOutput//Run2UltraLegacy_v3/LTAnalyzer/2018/LTAnalyzer_DYJetsToMuMu_MiNNLO.root");
  TFile f("/data6/Users/hsseo/SKFlatOutput//Run2UltraLegacy_v3/LTAnalyzer/2018/LTAnalyzer_DYJets_MG.root");
  TH1D* hj0a0=new TH1D("j0a0","j0a0",LTAnalyzer::nptbin,LTAnalyzer::ptbins);
  TH1D* hj0a2=new TH1D("j0a2","j0a2",LTAnalyzer::nptbin,LTAnalyzer::ptbins);
  TH1D* hj1a0=new TH1D("j1a0","j1a0",LTAnalyzer::nptbin,LTAnalyzer::ptbins);
  TH1D* hj1a2=new TH1D("j1a2","j1a2",LTAnalyzer::nptbin,LTAnalyzer::ptbins);
  TF2* func=new TF2("func","[0]*(1+x*x+0.5*[1]*(1-3*x*x)+0.5*[2]*(1-x*x)*cos(2*y))",0,1,0,TMath::Pi()/2);
  for(int i=1;i<LTAnalyzer::nptbin+2;i++){
    TH2* histj0=(TH2*)f.Get(Form("mm2018/jet0/gen_pt%d_costphi",i));
    histj0->Fit(func);
    hj0a0->SetBinContent(i,func->GetParameter(1));
    hj0a0->SetBinError(i,func->GetParError(1));
    hj0a2->SetBinContent(i,func->GetParameter(2));
    hj0a2->SetBinError(i,func->GetParError(2));
    TH2* histj1=(TH2*)f.Get(Form("mm2018/jet1/gen_pt%d_costphi",i));
    histj1->Fit(func);
    hj1a0->SetBinContent(i,func->GetParameter(1));
    hj1a0->SetBinError(i,func->GetParError(1));
    hj1a2->SetBinContent(i,func->GetParameter(2));
    hj1a2->SetBinError(i,func->GetParError(2));
  }
  TCanvas* c=new TCanvas("p1");
  hj0a0->Draw("hist e");
  hj0a2->Draw("same hist e");
  hj1a0->Draw("same hist e");
  hj1a2->Draw("same hist e");
}
