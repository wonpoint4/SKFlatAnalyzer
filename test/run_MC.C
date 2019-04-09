R__LOAD_LIBRARY(libPhysics.so)
R__LOAD_LIBRARY(libTree.so)
R__LOAD_LIBRARY(libHist.so)
R__LOAD_LIBRARY(libDataFormats.so)
R__LOAD_LIBRARY(libAnalyzerTools.so)
R__LOAD_LIBRARY(libAnalyzers.so)
R__LOAD_LIBRARY(/cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/lhapdf/6.2.1-fmblme/lib/libLHAPDF.so)

void run_MC(){

  //==== Anaylzer class name
  //==== change "ExampleRun" if you want to run another one
  ExampleRun m;

  m.SetTreeName("recoTree/SKFlat");

  //==== Not really important
  m.MCSample = "DYJets";
  m.IsDATA = false;

  //==== Xsec and weight information
  m.xsec = 6225.42;
  m.sumW = 123584524;

  //==== Data year
  m.DataYear = 2017;

  //==== If you want to activate usre flags, put them here as a vector<TString>

  m.Userflags = {

  };

  //==== Activate it if you want to skip first "NSkip" events
  //m.NSkipEvent = 10;

  //==== Activate it if you want to print [SKFlatNtuple::Loop RUNNING] XX / YY (ZZ %) every "LogEvery" events
  //==== Default if 1000
  //m.LogEvery = 1;

  //==== # of events you want to run
  //==== ex) if m.NSkipEvent = 2000 and m.MaxEvent = 1000,
  //====     it will skip first 2000 events and then run 1000 events 
  m.MaxEvent = 1000;

  //==== Add file here
  m.AddFile("root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/jskim/SKFlat/2017/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/SKFlat_Run2Legacy_v1/190101_204350/0000/SKFlatNtuple_2017_MC_1.root");
  m.AddFile("root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/jskim/SKFlat/2017/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/SKFlat_Run2Legacy_v1/190101_204350/0000/SKFlatNtuple_2017_MC_10.root");
  m.AddFile("root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/jskim/SKFlat/2017/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/SKFlat_Run2Legacy_v1/190101_204350/0000/SKFlatNtuple_2017_MC_100.root");

  //==== Output path
  m.SetOutfilePath("hists.root");

  m.Init();
  m.initializeAnalyzerTools();
  m.initializeAnalyzer();
  m.SwitchToTempDir();
  m.Loop();

  m.WriteHist();

}
