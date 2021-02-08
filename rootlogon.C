{
  //gSystem->Load("/cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/lhapdf/6.2.1-fmblme/lib/libLHAPDF.so");
  gROOT->ProcessLine(".L Plotter/AFBPlotter.cc");
  gROOT->ProcessLine(".L Plotter/AFBPlotter_QCD.cc");
  gROOT->ProcessLine(".L Plotter/EfficiencyPlotter.cc");
}
