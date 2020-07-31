#include <iostream>
#include "TString.h"
#include "TFile.h"
using namespace std;
void Save(TString channelname){
  bool delete_canvas=true;
  EfficiencyPlotter aa;
  aa.SetupPlots("fig/efficiency/plot.dat");
  Verbosity=VERBOSITY::WARNING;
  bool ismuon=channelname.BeginsWith("m");
  TString id="MediumID_Q";
  TString slepton="e";
  if(ismuon){
    id="MediumID_trkIsoLoose_Q";
    slepton="#mu";
  }
  TString sdilepton=slepton+slepton;
  
  TString default_option;
  default_option="sysname:efficiencySF ytitle:Events";
  aa.SavePlot("default/mass_"+channelname,"histname:"+channelname+"/m80to100/dimass_"+id+" xtitle:'m("+sdilepton+") [GeV]' "+default_option,delete_canvas);
  aa.SavePlot("default/rapidity_"+channelname,"histname:"+channelname+"/m80to100/dirap_"+id+" xtitle:'y("+sdilepton+")' rebin:2 "+default_option,delete_canvas);
  aa.SavePlot("default/pt_"+channelname,"histname:"+channelname+"/m80to100/dipt_"+id+" xtitle:'p_{T}("+sdilepton+") [GeV]' rebin:2 xmax:100 "+default_option,delete_canvas);
  aa.SavePlot("default/l0pt_"+channelname,"histname:"+channelname+"/m80to100/l0pt_"+id+" xtitle:'leading "+slepton+" p_{T} [GeV]' xmin:20 xmax:80 "+default_option,delete_canvas);
  aa.SavePlot("default/l1pt_"+channelname,"histname:"+channelname+"/m80to100/l1pt_"+id+" xtitle:'sub-leading "+slepton+" p_{T} [GeV]' xmin:10 xmax:60 "+default_option,delete_canvas);
  aa.SavePlot("default/leta_"+channelname,"histname:"+channelname+"/m80to100/leta_"+id+" xtitle:'#eta("+slepton+")' xmin:-2.4 xmax:2.4 rebin:2 "+default_option,delete_canvas);
  if(!ismuon) aa.SavePlot("default/lsceta_"+channelname,"histname:"+channelname+"/m80to100/lsceta_"+id+" xtitle:'#eta_{SC}("+slepton+")' rebin:2 "+default_option,delete_canvas);
  
  default_option="sysname:efficiencySF norm ytitle:Nomalized";
  aa.SavePlot("norm/mass_"+channelname,"histname:"+channelname+"/m80to100/dimass_"+id+" xtitle:'m("+sdilepton+") [GeV]' "+default_option,delete_canvas);
  aa.SavePlot("norm/rapidity_"+channelname,"histname:"+channelname+"/m80to100/dirap_"+id+" xtitle:'y("+sdilepton+")' rebin:2 "+default_option,delete_canvas);
  aa.SavePlot("norm/pt_"+channelname,"histname:"+channelname+"/m80to100/dipt_"+id+" xtitle:'p_{T}("+sdilepton+") [GeV]' rebin:2 xmax:100 "+default_option,delete_canvas);
  aa.SavePlot("norm/l0pt_"+channelname,"histname:"+channelname+"/m80to100/l0pt_"+id+" xtitle:'leading "+slepton+" p_{T} [GeV]' xmin:20 xmax:80 "+default_option,delete_canvas);
  aa.SavePlot("norm/l1pt_"+channelname,"histname:"+channelname+"/m80to100/l1pt_"+id+" xtitle:'sub-leading "+slepton+" p_{T} [GeV]' xmin:10 xmax:60 "+default_option,delete_canvas);
  aa.SavePlot("norm/leta_"+channelname,"histname:"+channelname+"/m80to100/leta_"+id+" xtitle:'#eta("+slepton+")' xmin:-2.4 xmax:2.4 rebin:2 "+default_option,delete_canvas);
  if(!ismuon) aa.SavePlot("norm/lsceta_"+channelname,"histname:"+channelname+"/m80to100/lsceta_"+id+" xtitle:'#eta_{SC}("+slepton+")' rebin:2 "+default_option,delete_canvas);
  
  aa.SetupEntries("data sim sim_noSF");

  default_option="sysname:efficiencySF ytitle:Events base:1";
  aa.SavePlot("nosf/mass_"+channelname,"histname:"+channelname+"/m80to100/dimass_"+id+" xtitle:'m("+sdilepton+") [GeV]' "+default_option,delete_canvas);
  aa.SavePlot("nosf/rapidity_"+channelname,"histname:"+channelname+"/m80to100/dirap_"+id+" xtitle:'y("+sdilepton+")' rebin:2 "+default_option,delete_canvas);
  aa.SavePlot("nosf/pt_"+channelname,"histname:"+channelname+"/m80to100/dipt_"+id+" xtitle:'p_{T}("+sdilepton+") [GeV]' rebin:2 xmax:100 "+default_option,delete_canvas);
  aa.SavePlot("nosf/l0pt_"+channelname,"histname:"+channelname+"/m80to100/l0pt_"+id+" xtitle:'leading "+slepton+" p_{T} [GeV]' xmin:20 xmax:80 "+default_option,delete_canvas);
  aa.SavePlot("nosf/l1pt_"+channelname,"histname:"+channelname+"/m80to100/l1pt_"+id+" xtitle:'sub-leading "+slepton+" p_{T} [GeV]' xmin:10 xmax:60 "+default_option,delete_canvas);
  aa.SavePlot("nosf/leta_"+channelname,"histname:"+channelname+"/m80to100/leta_"+id+" xtitle:'#eta("+slepton+")' xmin:-2.4 xmax:2.4 rebin:2 BMleg "+default_option,delete_canvas);
  if(!ismuon) aa.SavePlot("nosf/lsceta_"+channelname,"histname:"+channelname+"/m80to100/lsceta_"+id+" xtitle:'#eta_{SC}("+slepton+")' rebin:2 BMleg "+default_option,delete_canvas);
  
  default_option="sysname:efficiencySF norm ytitle:Nomalized base:1";
  aa.SavePlot("nosf_norm/mass_"+channelname,"histname:"+channelname+"/m80to100/dimass_"+id+" xtitle:'m("+sdilepton+") [GeV]' "+default_option,delete_canvas);
  aa.SavePlot("nosf_norm/rapidity_"+channelname,"histname:"+channelname+"/m80to100/dirap_"+id+" xtitle:'y("+sdilepton+")' rebin:2 "+default_option,delete_canvas);
  aa.SavePlot("nosf_norm/pt_"+channelname,"histname:"+channelname+"/m80to100/dipt_"+id+" xtitle:'p_{T}("+sdilepton+") [GeV]' rebin:2 xmax:100 "+default_option,delete_canvas);
  aa.SavePlot("nosf_norm/l0pt_"+channelname,"histname:"+channelname+"/m80to100/l0pt_"+id+" xtitle:'leading "+slepton+" p_{T} [GeV]' xmin:20 xmax:80 "+default_option,delete_canvas);
  aa.SavePlot("nosf_norm/l1pt_"+channelname,"histname:"+channelname+"/m80to100/l1pt_"+id+" xtitle:'sub-leading "+slepton+" p_{T} [GeV]' xmin:10 xmax:60 "+default_option,delete_canvas);
  aa.SavePlot("nosf_norm/leta_"+channelname,"histname:"+channelname+"/m80to100/leta_"+id+" xtitle:'#eta("+slepton+")' xmin:-2.4 xmax:2.4 rebin:2 BMleg "+default_option,delete_canvas);
  if(!ismuon) aa.SavePlot("nosf_norm/lsceta_"+channelname,"histname:"+channelname+"/m80to100/lsceta_"+id+" xtitle:'#eta_{SC}("+slepton+")' rebin:2 BMleg "+default_option,delete_canvas);
}
void SaveAll(){
  TString channelnames[]={"ee2016","ee2017","ee2018","mm2016","mm2017","mm2018","el2016","el2017","el2018","mu2016","mu2017","mu2018"};
  for(auto channelname:channelnames) Save(channelname);
}
void make_muon_link(std::ostream& out=std::cout){
  TString prefix="https://wjun.web.cern.ch/wjun/Muon/AFB_Efficiency/";
  vector<pair<TString,TString>> pairs={
    {"Medium trkIsoLoose","2016BF_IDISO"},
    {"Single muon","2016BF_IsoMu24"},
    {"Double muon Leg1 (Mu17)","2016BF_Mu17"},
    {"Double muon Leg2 (Mu8)","2016BF_Mu8"},
    {"Medium trkIsoLoose","2016GH_IDISO"},
    {"Single muon","2016GH_IsoMu24"},
    {"Double muon Leg1 (Mu17)","2016GH_Mu17"},
    {"Double muon Leg2 (Mu8)","2016GH_Mu8"},
    {"Medium trkIsoLoose","2017_IDISO"},
    {"Single muon","2017_IsoMu27"},
    {"Double muon Leg1 (Mu17)","2017_Mu17"},
    {"Double muon Leg2 (Mu8)","2017_Mu8"},
    {"Medium trkIsoLoose","2018_IDISO"},
    {"Single muon","2018_IsoMu27"},
    {"Double muon Leg1 (Mu17)","2018_Mu17"},
    {"Double muon Leg2 (Mu8)","2018_Mu8"},
  };
  cout<<"<table>"<<endl;
  for(auto [title,path]:pairs){
    cout<<"<tr>"<<endl;
    cout<<"<td>"+title+"</td>"<<endl;
    cout<<"<td>+: <a href=\""+prefix+path+"_+/Plots\">plots</a> <a href=\""+prefix+path+"_+/result.root\">root</a> <br>-: <a href=\""+prefix+path+"_-/Plots\">plots</a> <a href=\""+prefix+path+"_-/result.root\">root</a> </td>"<<endl;
    cout<<"</tr>"<<endl;
  }
  cout<<"</table>"<<endl;
}
   


