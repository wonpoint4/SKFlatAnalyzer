#include"AFBPlotter.cc"

void plot(){
  AFBPlotter aa("mi mi mi");
  aa.entries[0].title="reco";
  aa.entries[1].title="gen fiducial dressed";
  aa.entries[1].styles[0].linecolor=4;
  aa.entries[1].SetHistPrefix("genfid_");
  aa.entries[1].replace["$"]="_dressed";
  aa.entries[2].title="gen";
  aa.entries[2].styles[0].linecolor=6;
  aa.entries[2].SetHistPrefix("gen_");
  aa.DrawPlot("mm2018/0bjet/AFB(m)","logx type:1 ymin:-0.29 ymax:0.49");
  aa.DrawPlot("ee2018/0bjet/AFB(m)","logx type:1 ymin:-0.29 ymax:0.49");
  
  aa.Setup("mi mi");
  aa.entries[0].title="ee";
  aa.entries[0].replace["LL"]="ee";
  aa.entries[1].title="#mu#mu";
  aa.entries[1].styles[0].linecolor=4;
  aa.entries[1].replace["LL"]="mm";
  aa.DrawPlot("LL201[678][ab]?/0bjet/m[52,3000]/genfid_AFB_dressed","logx project:x");
  aa.DrawPlot("LL201[678][ab]?/0bjet/m[52,3000]/genfid_AFB_bare","logx project:x");
  aa.DrawPlot("LL201[678][ab]?/0bjet/m[52,3000]/gen_AFB","logx project:x");


  aa.Setup("mi mi");
  aa.entries[0].title="reco";
  aa.entries[0].replace["unfold_afbm"]="AFB";
  aa.entries[1].title="unfolded";
  aa.entries[1].styles[0].linecolor=4;
  aa.DrawPlot("ee2018/0bjet/m[52,3000]/unfold_afbm","logx project:x type:1");
  aa.DrawPlot("mm2018/0bjet/m[52,3000]/unfold_afbm","logx project:x type:1");

  aa.Setup("mi mi");
  aa.entries[0].title="POWHEG truth";
  aa.entries[0].replace["unfold_afbm"]="genfid_AFB_dressed";
  aa.entries[1].title="unfolded";
  aa.entries[1].SetTags("data");
  aa.entries[1].styles[0].linecolor=4;
  aa.DrawPlot("ee2018/0bjet/m[52,3000]/unfold_afbm","logx project:x type:1");
  aa.DrawPlot("mm2018/0bjet/m[52,3000]/unfold_afbm","logx project:x type:1");

  aa.Setup("amc amc");
  aa.entries[0].title="aMC@NLO truth";
  aa.entries[0].replace["unfold_afbm"]="genfid_AFB_dressed";
  aa.entries[1].title="unfolded by POWHEG";
  aa.entries[1].SetTags("data");
  aa.entries[1].styles[0].linecolor=4;
  aa.DrawPlot("ee201[678][ab]?/0bjet/m[52,3000]/unfold_afbm","logx project:x type:1");
  aa.DrawPlot("mm201[678][ab]?/0bjet/m[52,3000]/unfold_afbm","logx project:x type:1");
  
  aa.Setup("mi mi amc amc");
  aa.entries[0].title="POWHEG truth";
  aa.entries[0].replace["unfold_afbm"]="genfid_AFB_dressed";
  aa.entries[1].title="POWHEG unfolded";
  aa.entries[1].SetTags("data");
  aa.entries[1].styles[0].linecolor=4;
  aa.entries[2].title="aMC@NLO truth";
  aa.entries[2].replace["unfold_afbm"]="genfid_AFB_dressed";
  aa.entries[2].styles[0].linecolor=3;
  aa.entries[3].title="aMC@NLO unfolded by POWHEG";
  aa.entries[3].SetTags("data");
  aa.entries[3].styles[0].linecolor=6;
  aa.DrawPlot("ee201[678][ab]?/0bjet/m[52,3000]/unfold_afbm","logx project:x type:1");
  aa.DrawPlot("mm201[678][ab]?/0bjet/m[52,3000]/unfold_afbm","logx project:x type:1");


}
  
void result(){
  AFBPlotter aa("mi");
  aa.entries[0].SetTags("data");
  TH1* nominal=aa.GetHist(0,"ee201[678][ab]?/0bjet/m[52,3000]/unfold_afbm");
  vector<TH1*> sys;

  vector<TH1*> hists;
  //scale
  for(int i=0;i<9;i++){
    if(i==5||i==7) continue;
    hists.push_back(aa.GetHist(0,"ee201[678][ab]?/0bjet/m[52,3000]/unfold_afbm",Form("suffix:_scalevariation%d:dy",i)));
  }
  sys.push_back(aa.GetEnvelope(nominal,hists));
  hists.clear();
  
  //pdf
  for(int i=0;i<100;i++){
    hists.push_back(aa.GetHist(0,"ee201[678][ab]?/0bjet/m[52,3000]/unfold_afbm",Form("suffix:_pdf%d:dy",i)));
  }
  sys.push_back(aa.GetHessianError(nominal,hists));
  hists.clear();

  //zpt
  hists.push_back(aa.GetHist(0,"ee201[678][ab]?/0bjet/m[52,3000]/unfold_afbm","suffix:_nozptweight:dy"));
  sys.push_back(aa.GetEnvelope(nominal,hists));
  hists.clear();

  //eff stat
  for(int i=0;i<20;i++){
    hists.push_back(aa.GetHist(0,"ee201[678][ab]?/0bjet/m[52,3000]/unfold_afbm",Form("suffix:_efficiencySF_stat%d:dy",i)));
  }
  sys.push_back(aa.GetRMSError(nominal,hists));
  hists.clear();

  //eff syst FIXME
  for(int i=1;i<5;i++) hists.push_back(aa.GetHist(0,"ee201[678][ab]?/0bjet/m[52,3000]/unfold_afbm",Form("suffix:_electronRECOSF_s%d_m0:dy",i)));
  for(int i=1;i<8;i++) hists.push_back(aa.GetHist(0,"ee201[678][ab]?/0bjet/m[52,3000]/unfold_afbm",Form("suffix:_electronIDSF_s%d_m0:dy",i)));
  for(int i=1;i<5;i++) hists.push_back(aa.GetHist(0,"ee201[678][ab]?/0bjet/m[52,3000]/unfold_afbm",Form("suffix:_muonIDSF_s%d_m0:dy",i)));
  for(int i=1;i<5;i++) hists.push_back(aa.GetHist(0,"ee201[678][ab]?/0bjet/m[52,3000]/unfold_afbm",Form("suffix:_triggerSF_s%d_m0:dy",i)));
  sys.push_back(aa.GetHessianError(nominal,hists));
  hists.clear();
  
  for(int i=1;i<nominal->GetNbinsX()+1;i++){
    cout<<"bin "<<i<<" "<<nominal->GetBinLowEdge(i)<<"<m<"<<nominal->GetBinLowEdge(i+1)<<endl;
    cout<<nominal->GetBinContent(i)<<endl;
    cout<<nominal->GetBinError(i)<<endl;
    for(int j=0;j<sys.size();j++){
      cout<<sys[j]->GetBinError(i)<<endl;
    }
  }
}
