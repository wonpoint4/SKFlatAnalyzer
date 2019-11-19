#include <vector>

void Init(int channel,int year,TString mode){
  EfficiencyPlotter aa;
  aa.Setup(channel,year,mode,1);
  aa.AddPlots("m60to120/lpt","xmax:80 norm sysname:efficiencySF");
  aa.AddPlots("m60to120/l0pt","xmax:80 norm sysname:efficiencySF");
  aa.AddPlots("m60to120/l1pt","xmax:80 norm sysname:efficiencySF");

  aa.AddPlots("m60to120/lppt","xmax:80 norm sysname:efficiencySF");
  aa.AddPlots("m60to120/l0ppt","xmax:80 norm sysname:efficiencySF");
  aa.AddPlots("m60to120/l1ppt","xmax:80 norm sysname:efficiencySF");

  aa.AddPlots("m60to120/lmpt","xmax:80 norm sysname:efficiencySF");
  aa.AddPlots("m60to120/l0mpt","xmax:80 norm sysname:efficiencySF");
  aa.AddPlots("m60to120/l1mpt","xmax:80 norm sysname:efficiencySF");

  aa.AddPlots("m60to120/leta","norm sysname:efficiencySF");
  aa.AddPlots("m60to120/l0eta","norm sysname:efficiencySF");
  aa.AddPlots("m60to120/l1eta","norm sysname:efficiencySF");

  aa.AddPlots("m60to120/lpeta","norm sysname:efficiencySF");
  aa.AddPlots("m60to120/l0peta","norm sysname:efficiencySF");
  aa.AddPlots("m60to120/l1peta","norm sysname:efficiencySF");

  aa.AddPlots("m60to120/lmeta","norm sysname:efficiencySF");
  aa.AddPlots("m60to120/l0meta","norm sysname:efficiencySF");
  aa.AddPlots("m60to120/l1meta","norm sysname:efficiencySF");
  
  aa.AddPlots("m60to120/dimass","norm xmin:60 xmax:120 sysname:efficiencySF");
  aa.AddPlots("m60to120/dipt","norm xmax:100 sysname:efficiencySF");

  aa.AddPlots("m60to120/.*costhetaCS","norm sysname:efficiencySF");
  
  aa.RemovePlots("POGTight");
  aa.RemovePlots("MediumID");  
  aa.RemovePlots("TightID");  
}

void SaveAll(bool init=false){
  std::vector<std::tuple<int,int,TString>> targets;
  //  targets.push_back(make_tuple(0,2016,"POGTight_PFIsoTight"));
  targets.push_back(make_tuple(0,2017,"POGTight_PFIsoTight"));
  targets.push_back(make_tuple(0,2017,"POGTight_TrkIsoLoose_Q"));
  targets.push_back(make_tuple(1,2016,"MediumID"));
  targets.push_back(make_tuple(1,2016,"MediumID_selective_Q"));
  targets.push_back(make_tuple(1,2017,"MediumID"));
  targets.push_back(make_tuple(1,2017,"MediumID_Q"));
  targets.push_back(make_tuple(1,2017,"MediumID_selective_Q"));
  targets.push_back(make_tuple(1,2017,"TightID_Q"));
  
  for(const auto& [channel,year,mode]:targets){
    if(init) Init(channel,year,mode);
    EfficiencyPlotter aa;
    aa.Setup(channel,year,mode,1);
    aa.SaveCanvases("m60to120");
  }
}
