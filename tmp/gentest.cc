{
  AFBPlotter aa;
  aa.AddFile("ulmi","/data6/Users/hsseo/SKFlatOutput/Run2UltraLegacy_v1/AFBAnalyzer/2016preVFP/AFBAnalyzer_DYJetsToMuMu_MiNNLO.root");
  aa.AddFile("ulamc","/data6/Users/hsseo/SKFlatOutput/Run2UltraLegacy_v1/AFBAnalyzer/2018/AFBAnalyzer_DYJets.root");
  aa.AddFile("ulmg","/data6/Users/hsseo/SKFlatOutput/Run2UltraLegacy_v1/AFBAnalyzer/2018/AFBAnalyzer_DYJets_MG.root");
  aa.AddFile("preamc","/data6/Users/hsseo/SKFlatOutput/Run2Legacy_v4/AFBAnalyzer/2018/AFBAnalyzer_DYJets.root");
  aa.AddFile("premg","/data6/Users/hsseo/SKFlatOutput/Run2Legacy_v4/AFBAnalyzer/2018/AFBAnalyzer_DYJets_MG.root");


  aa.AddFile("ulamc16a","/data6/Users/hsseo/SKFlatOutput/Run2UltraLegacy_v1/AFBAnalyzer/2016preVFP/AFBAnalyzer_DYJets.root");
  aa.AddFile("ulamc16aS20","/data6/Users/hsseo/SKFlatOutput/Run2UltraLegacy_v1/AFBAnalyzer/2016preVFP/AFBAnalyzer_DYJets_Summer20.root");
  aa.AddFile("ulmi16a","/data6/Users/hsseo/SKFlatOutput/Run2UltraLegacy_v1/AFBAnalyzer/2016preVFP/AFBAnalyzer_DYJetsToMuMu_MiNNLO.root");
  aa.AddFile("ulmi16b","/data6/Users/hsseo/SKFlatOutput/Run2UltraLegacy_v1/AFBAnalyzer/2016postVFP/AFBAnalyzer_DYJetsToMuMu_MiNNLO.root");
  aa.AddFile("ulmg16a","/data6/Users/hsseo/SKFlatOutput/Run2UltraLegacy_v1/AFBAnalyzer/2016preVFP/AFBAnalyzer_DYJets_MG.root");
  aa.AddFile("ulmg16b","/data6/Users/hsseo/SKFlatOutput/Run2UltraLegacy_v1/AFBAnalyzer/2016postVFP/AFBAnalyzer_DYJets_MG.root");
  aa.AddFile("ulmg17","/data6/Users/hsseo/SKFlatOutput/Run2UltraLegacy_v1/AFBAnalyzer/2017/AFBAnalyzer_DYJets_MG.root");
  aa.AddFile("ulmg18","/data6/Users/hsseo/SKFlatOutput/Run2UltraLegacy_v1/AFBAnalyzer/2018/AFBAnalyzer_DYJets_MG.root");
  
  aa.Setup("ulmi16a ulmi16a ulmi16b ulmi16b");
  aa.entries[0].replace["gen_"]="lhe_";
  aa.entries[2].replace["gen_"]="lhe_";
  //aa.entries[2].replace["_noweight"]="_dressed_noweight";
  //aa.entries[3].replace["_noweight"]="_bare_noweight";
  //aa.Setup("ulamc preamc");
  //aa.Setup("ulmi ulamc ulmg preamc premg");
  
  //vector<TString> hprefixes={"gen_"};
  vector<TString> hprefixes={"lhe_","gen_","genfid_"};
  vector<TString> types={"","_bare"};
  //vector<TString> types={"","_dressed","_bare"};
  //vector<TString> suffixes={""};
  vector<TString> suffixes={"","_noweight"};
  
  for(auto hprefix:hprefixes){
    for(auto type:types){
      for(auto suffix:suffixes){
	//aa.DrawPlot("mm201[678ab]+/m[80,100]/"+hprefix+"dipt"+type+suffix,"norm widthweight xmax:200");
	//aa.DrawPlot("mm201[678ab]+/m[80,100]/y[0,1.2]/"+hprefix+"dipt"+type+suffix,"norm widthweight xmax:200 absy");
	//aa.DrawPlot("mm201[678ab]+/m[80,100]/y[1.2,2.4]/"+hprefix+"dipt"+type+suffix,"norm widthweight xmax:200 absy");
	aa.DrawPlot("mm201[678ab]+/m[60,400]/"+hprefix+"dimass"+type+suffix,"norm widthweight");
	//aa.DrawPlot("mm201[678ab]+/m[80,100]/"+hprefix+"lpt"+type+suffix,"norm widthweight xmax:100");
	//aa.DrawPlot("mm201[678ab]+/m[60,400]/"+hprefix+"lpt"+type+suffix,"norm widthweight xmax:100");
      }
    }
  }
}
