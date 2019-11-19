vector<TString> yregion={"_y1.2to1.6","_y1.6to2.0","_y2.0to2.4"};

double getchi2(AFBPlotter& a,TString histname){
  TH1* h=a.GetHist("data",histname);
  h->GetXaxis()->SetRangeUser(60,120);
  double out=a.GetChi2(h);
  delete h;
  return out;
}
  
void chi2table(){
  DEBUG=0;
  vector<TString> regions={"m60to400","m60to400_y0.0to0.4","m60to400_y0.4to0.8","m60to400_y0.8to1.2","m60to400_y1.2to1.6","m60to400_y1.6to2.0","m60to400_y2.0to2.4","m60to400_pt50","m60to400_y0.0to0.4_pt50","m60to400_y0.4to0.8_pt50","m60to400_y0.8to1.2_pt50","m60to400_y1.2to1.6_pt50","m60to400_y1.6to2.0_pt50","m60to400_y2.0to2.4_pt50"};
  AFBPlotter a;
  for(int ic=0;ic<2;ic++){
    for(int iy=2016;iy<2019;iy++){
      a.Setup(ic,iy,0);
      for(const auto& region:regions){
	TString suffix=ic==0?"_trkiso":"_selective";
	cout<<a.schannel<<a.syear<<"\t"<<region<<"\t"<<getchi2(a,a.schannel+a.syear+"/"+region+"/AFB")<<"\t"<<getchi2(a,a.schannel+a.syear+"/"+region+"/weightedAFB")<<"\t"<<getchi2(a,a.schannel+a.syear+"/"+region+"/AFB"+suffix)<<"\t"<<getchi2(a,a.schannel+a.syear+"/"+region+"/weightedAFB"+suffix)<<endl;
      }
    }
  }
}

void plot_noeffiSF(){
  AFBPlotter p;
  p.Setup(1,2017,10);
  p.GetCanvas(Plotter::GetCompareAndRatio,"electron2017/m60to400/dimass","_noefficiencySF",64,2,60,120,"2:widey");
  p.GetCanvas(Plotter::GetCompareAndRatio,"electron2017/m60to400/l0mpt","_noefficiencySF",64,2,0,70,"2:widey");
  p.GetCanvas(Plotter::GetCompareAndRatio,"electron2017/m60to400/l1mpt","_noefficiencySF",64,2,0,70,"2:widey");
  p.GetCanvas(Plotter::GetCompareAndRatio,"electron2017/m60to400/l0ppt","_noefficiencySF",64,2,0,70,"2:widey");
  p.GetCanvas(Plotter::GetCompareAndRatio,"electron2017/m60to400/l1ppt","_noefficiencySF",64,2,0,70,"2:widey");
  p.GetCanvas(Plotter::GetCompareAndRatio,"electron2017/m60to400/l0meta","_noefficiencySF",64,0,0,0,"2:widey");
  p.GetCanvas(Plotter::GetCompareAndRatio,"electron2017/m60to400/l1meta","_noefficiencySF",64,0,0,0,"2:widey");
  p.GetCanvas(Plotter::GetCompareAndRatio,"electron2017/m60to400/l0peta","_noefficiencySF",64,0,0,0,"2:widey");
  p.GetCanvas(Plotter::GetCompareAndRatio,"electron2017/m60to400/l1peta","_noefficiencySF",64,0,0,0,"2:widey");

  p.Setup(0,2016,10);
  p.GetCanvas(Plotter::GetCompareAndRatio,"muon2016/m60to400/dimass","_noefficiencySF",64,2,60,120,"2:widey");
  p.GetCanvas(Plotter::GetCompareAndRatio,"muon2016/m60to400/l0mpt","_noefficiencySF",64,2,0,70,"2:widey");
  p.GetCanvas(Plotter::GetCompareAndRatio,"muon2016/m60to400/l1mpt","_noefficiencySF",64,2,0,70,"2:widey");
  p.GetCanvas(Plotter::GetCompareAndRatio,"muon2016/m60to400/l0ppt","_noefficiencySF",64,2,0,70,"2:widey");
  p.GetCanvas(Plotter::GetCompareAndRatio,"muon2016/m60to400/l1ppt","_noefficiencySF",64,2,0,70,"2:widey");
  p.GetCanvas(Plotter::GetCompareAndRatio,"muon2016/m60to400/l0meta","_noefficiencySF",64,0,0,0,"2:widey");
  p.GetCanvas(Plotter::GetCompareAndRatio,"muon2016/m60to400/l1meta","_noefficiencySF",64,0,0,0,"2:widey");
  p.GetCanvas(Plotter::GetCompareAndRatio,"muon2016/m60to400/l0peta","_noefficiencySF",64,0,0,0,"2:widey");
  p.GetCanvas(Plotter::GetCompareAndRatio,"muon2016/m60to400/l1peta","_noefficiencySF",64,0,0,0,"2:widey");

  p.Setup(0,2017,10);
  p.GetCanvas(Plotter::GetCompareAndRatio,"muon2017/m60to400/dimass","_noefficiencySF",64,2,60,120,"2:widey");
  p.GetCanvas(Plotter::GetCompareAndRatio,"muon2017/m60to400/l0mpt","_noefficiencySF",64,2,0,70,"2:widey");
  p.GetCanvas(Plotter::GetCompareAndRatio,"muon2017/m60to400/l1mpt","_noefficiencySF",64,2,0,70,"2:widey");
  p.GetCanvas(Plotter::GetCompareAndRatio,"muon2017/m60to400/l0ppt","_noefficiencySF",64,2,0,70,"2:widey");
  p.GetCanvas(Plotter::GetCompareAndRatio,"muon2017/m60to400/l1ppt","_noefficiencySF",64,2,0,70,"2:widey");
  p.GetCanvas(Plotter::GetCompareAndRatio,"muon2017/m60to400/l0meta","_noefficiencySF",64,0,0,0,"2:widey");
  p.GetCanvas(Plotter::GetCompareAndRatio,"muon2017/m60to400/l1meta","_noefficiencySF",64,0,0,0,"2:widey");
  p.GetCanvas(Plotter::GetCompareAndRatio,"muon2017/m60to400/l0peta","_noefficiencySF",64,0,0,0,"2:widey");
  p.GetCanvas(Plotter::GetCompareAndRatio,"muon2017/m60to400/l1peta","_noefficiencySF",64,0,0,0,"2:widey");

  for(auto* obj:*gROOT->GetListOfCanvases()){
    TCanvas* c=(TCanvas*)obj;
    gSystem->Exec("mkdir -p plot_noeffiSF");
    c->SaveAs(TString("plot_noeffiSF/")+c->GetName()+TString(".png"));
    //c->Delete();
  }
}

void plot_selective(){
  AFBPlotter p;
  for(int j=2016;j<2019;j++){
    p.Setup(1,j,11);
    p.GetCanvas(Plotter::GetCompareAndRatio,p.schannel+p.syear+"/m60to400/l0pt","_selective",64,2,0,70,"2:widey");
    p.GetCanvas(Plotter::GetCompareAndRatio,p.schannel+p.syear+"/m60to400/l1pt","_selective",64,2,0,70,"2:widey");
    p.GetCanvas(Plotter::GetCompareAndRatio,p.schannel+p.syear+"/m60to400/l0eta","_selective",64,0,0,0,"2:widey");
    p.GetCanvas(Plotter::GetCompareAndRatio,p.schannel+p.syear+"/m60to400/l1eta","_selective",64,0,0,0,"2:widey");

    p.GetCanvas(Plotter::GetCompareAndDiff,p.schannel+p.syear+"/m60to400_y1.2to1.6/AFB","_selective",64,0,60,120,"1:leftleg");
    p.GetCanvas(Plotter::GetCompareAndDiff,p.schannel+p.syear+"/m60to400_y1.2to1.6/weightedAFB","_selective",64,0,60,120,"1:leftleg");
    p.GetCanvas(Plotter::GetCompareAndDiff,p.schannel+p.syear+"/m60to400_y1.2to1.6/AFB","_selective",64,0,0,0,"1:leftleg");
    p.GetCanvas(Plotter::GetCompareAndDiff,p.schannel+p.syear+"/m60to400_y1.2to1.6/weightedAFB","_selective",64,0,0,0,"1:leftleg");
    p.GetCanvas(Plotter::GetCompareAndDiff,p.schannel+p.syear+"/m60to400_y1.2to1.6_pt50/AFB","_selective",64,0,0,0,"1:leftleg");
    p.GetCanvas(Plotter::GetCompareAndDiff,p.schannel+p.syear+"/m60to400_y1.2to1.6_pt50/weightedAFB","_selective",64,0,0,0,"1:leftleg");
  }
  for(auto* obj:*gROOT->GetListOfCanvases()){
    TCanvas* c=(TCanvas*)obj;
    gSystem->Exec("mkdir -p plot_selective");
    c->SaveAs(TString("plot_selective/")+c->GetName()+TString(".png"));
    //c->Delete();
  }
}
void plot_trkiso(){
  AFBPlotter p;
  for(int j=2016;j<2018;j++){
    p.Setup(0,j,12);
    p.GetCanvas(Plotter::GetCompareAndRatio,p.schannel+p.syear+"/m60to400/l0pt","_trkiso",64,2,0,70,"2:widey");
    p.GetCanvas(Plotter::GetCompareAndRatio,p.schannel+p.syear+"/m60to400/l1pt","_trkiso",64,2,0,70,"2:widey");
    p.GetCanvas(Plotter::GetCompareAndRatio,p.schannel+p.syear+"/m60to400/l0eta","_trkiso",64,0,0,0,"2:widey");
    p.GetCanvas(Plotter::GetCompareAndRatio,p.schannel+p.syear+"/m60to400/l1eta","_trkiso",64,0,0,0,"2:widey");

    for(const auto& y:yregion){
      p.GetCanvas(Plotter::GetCompareAndDiff,p.schannel+p.syear+"/m60to400"+y+"/AFB","_trkiso",64,0,60,120,"1:leftleg");
      p.GetCanvas(Plotter::GetCompareAndDiff,p.schannel+p.syear+"/m60to400"+y+"/weightedAFB","_trkiso",64,0,60,120,"1:leftleg");
      p.GetCanvas(Plotter::GetCompareAndDiff,p.schannel+p.syear+"/m60to400"+y+"/AFB","_trkiso",64,0,0,0,"1:leftleg");
      p.GetCanvas(Plotter::GetCompareAndDiff,p.schannel+p.syear+"/m60to400"+y+"/weightedAFB","_trkiso",64,0,0,0,"1:leftleg");
      p.GetCanvas(Plotter::GetCompareAndDiff,p.schannel+p.syear+"/m60to400"+y+"_pt50/AFB","_trkiso",64,0,0,0,"1:leftleg");
      p.GetCanvas(Plotter::GetCompareAndDiff,p.schannel+p.syear+"/m60to400"+y+"_pt50/weightedAFB","_trkiso",64,0,0,0,"1:leftleg");
    }
  }
  for(auto* obj:*gROOT->GetListOfCanvases()){
    TCanvas* c=(TCanvas*)obj;
    gSystem->Exec("mkdir -p plot_trkiso");
    c->SaveAs(TString("plot_trkiso/")+c->GetName()+TString(".png"));
    //c->Delete();
  }
}
