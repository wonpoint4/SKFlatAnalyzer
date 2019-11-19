#include "Plotter/EfficiencyPlotter.cc"
void SaveAll(){
  vector<tuple<int,int,TString>> availables={make_tuple(0,2016,"POGTight_PFIsoTight"),
                                             make_tuple(0,2017,"POGTight_PFIsoTight"),
                                             make_tuple(1,2016,"MediumID"),
                                             make_tuple(1,2016,"MediumID_selective_Q"),
                                             make_tuple(1,2017,"MediumID"),
                                             make_tuple(1,2017,"MediumID_Q"),
                                             make_tuple(1,2017,"TightID_Q"),};
  DEBUG=2;
  for(const auto& [channel,year,mode]:availables){
    EfficiencyPlotter p;
    p.Setup(channel,year,mode,true);
    p.RemovePlots(".*");
    p.AddPlots("/l[01mp]*pt$","xmax:100 rebin:2");
    //p.SetPlotsOption("/l[01mp]*pt$","norm","_norm");
    p.AddPlots("/l[01mp]*eta$","rebin:2");
    //p.SetPlotsOption("/l[01mp]*eta$","norm","_norm");
    p.AddPlots("/dimass$","xmin:80 xmax:100");
    //p.SetPlotsOption("/dimass$","norm","_norm");
    p.AddPlots("/costhetaCS$","");
    //p.SetPlotsOption("/costhetaCS$","norm","_norm");
    p.AddPlots("/dipt$","xmax:100 rebin:2");
    //p.SetPlotsOption("/dipt$","norm","_norm");
    p.AddPlots("/dirap$","");
    //p.SetPlotsOption("/dirap$","norm","_norm");
    p.SetPlotsOption(".*","norm","_norm");
    p.PrintPlots();
    p.SaveCanvases(".*");
    for(const auto& [key,name]:p.plots){
      if(key.Contains(TRegexp("l[01]?mpt"))&&!key.Contains("norm")){
	p.AddPlot(Replace(key,"mpt","[mp]pt"),Form("type:%d base:0",Plot::Type::DoubleRatio));
      }
      if(key.Contains(TRegexp("l[01]?meta"))&&!key.Contains("norm")){
	p.AddPlot(Replace(key,"meta","[mp]eta"),Form("type:%d base:0",Plot::Type::DoubleRatio));
      }
    }
    p.entries.pop_back();
    p.SaveCanvases("\\[mp\\]");
  }
}

    
    
