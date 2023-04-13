#include"AFBPlotter.cc"
#include"EfficiencyPlotter.cc"

void plot1(){
  EfficiencyPlotter aa;
  aa.plotdir="fig/230330napp";
  aa.SavePlotCondor("mm201[678][ab]?/m80to100/lpt","xmax:100 xtitle:'muon p_{T} [GeV]' suffix:_noefficiencySF:sim preliminary save:eff_mmRun2_lpt_noefficiencySF");
  aa.SavePlotCondor("mm201[678][ab]?/m80to100/lpt","xmax:100 xtitle:'muon p_{T} [GeV]' sysname:efficiencySF preliminary save:eff_mmRun2_lpt");
  aa.SavePlotCondor("ee201[678][ab]?/m80to100/lpt","xmax:100 xtitle:'electron p_{T} [GeV]' suffix:_noefficiencySF:sim preliminary save:eff_eeRun2_lpt_noefficiencySF");
  aa.SavePlotCondor("ee201[678][ab]?/m80to100/lpt","xmax:100 xtitle:'electron p_{T} [GeV]' sysname:efficiencySF preliminary save:eff_eeRun2_lpt");
  aa.SavePlotCondor("ee2017/m52to150/dimass","xtitle:'m(ee) [GeV]' preliminary save:ee2017_dimass");
  aa.SavePlotCondor("ee2017/m52to150/dimass","xtitle:'m(ee) [GeV]' preliminary save:ee2017_dimass_noroccor suffix:_noroccor");
  aa.SavePlotCondor("mm2017/m52to150/dimass","xtitle:'m(#mu#mu) [GeV]' preliminary save:mm2017_dimass");
  aa.SavePlotCondor("mm2017/m52to150/dimass","xtitle:'m(#mu#mu) [GeV]' preliminary save:mm2017_dimass_noroccor suffix:_noroccor");
}
void plot2(){
  AFBPlotter aa;
  aa.plotdir="fig/230330napp";
  aa.SavePlotCondor("mm201[678][ab]?/[0n]bjet/m[52,3000]/lpt","xmax:100 xtitle:'muon p_{T} [GeV]' sysname:totalsys preliminary save:mmRun2_lpt 1:logy widthweight");
  aa.SavePlotCondor("mm201[678][ab]?/[0n]bjet/m[52,3000]/lpt","xmax:100 xtitle:'muon p_{T} [GeV]' suffix:_nozptweight:dy preliminary save:mmRun2_lpt_nozptweight 1:logy widthweight");
  aa.SavePlotCondor("mm201[678][ab]?/[0n]bjet/m[52,3000]/dipt","xmax:650 xtitle:'muon p_{T} [GeV]' sysname:totalsys preliminary save:mmRun2_dipt 1:logy logx widthweight");
  aa.SavePlotCondor("mm201[678][ab]?/[0n]bjet/m[52,3000]/dipt","xmax:650 xtitle:'muon p_{T} [GeV]' suffix:_nozptweight:dy preliminary save:mmRun2_dipt_nozptweight 1:logy logx widthweight");
}
void plot(){
  plot1();
  plot2();
}
