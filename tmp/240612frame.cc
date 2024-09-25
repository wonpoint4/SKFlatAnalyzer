#include"AFBSystPlotter.cc"

void plot_DY_CS(){
  AFBSystPlotter aa("mi mi mi mi");
  aa.entries[0].title="Parton-level";
  aa.entries[0].SetHistPrefix("correct_");
  aa.entries[0].styles[0].linecolor=1;
  aa.entries[0].styles[0].linewidth=2;
  
  aa.entries[1].title="Direction guess using rapidity";
  aa.entries[1].SetHistPrefix("gen_");
  aa.entries[1].styles[0].linecolor=2;
  aa.entries[1].styles[0].linewidth=2;

  aa.entries[2].title="Acceptance cut";
  aa.entries[2].SetHistPrefix("genfid_");
  aa.entries[2].styles[0].linecolor=4;
  aa.entries[2].styles[0].linewidth=2;

  aa.entries[3].title="Detector-level";
  aa.entries[3].styles[0].linecolor=kGreen+1;
  aa.entries[3].styles[0].linewidth=2;
  aa.entries[3].replace["CS"]="";

  aa.SavePlot("mm201[678][ab]?/0bjet/dimassCS","AFB logx type:1 save:test.png BRleg xtitle:'m(ll) [GeV]' ytitle:A_{FB}");

}

void frame2(){
  AFBSystPlotter aa("mi mi mi mi");
  aa.entries[0].title="Correct direction";
  aa.entries[0].SetHistPrefix("gen_");
  aa.entries[0].styles[0].linecolor=1;
  aa.entries[0].styles[0].linewidth=2;
  
  aa.entries[1].title="Rapidity method";
  aa.entries[1].SetHistPrefix("gen_");
  aa.entries[1].styles[0].linecolor=2;
  aa.entries[1].styles[0].linewidth=2;

  aa.entries[2].title="Acceptance";
  aa.entries[2].SetHistPrefix("genfid_");
  aa.entries[2].styles[0].linecolor=4;
  aa.entries[2].styles[0].linewidth=2;

  aa.entries[3].title="Detector-level";
  aa.entries[3].styles[0].linecolor=kGreen+1;
  aa.entries[3].styles[0].linewidth=2;
  aa.entries[3].replace["Recoil"]="";

  aa.SavePlot("mm201[678][ab]?/nbjet/dimassRecoil","AFB logx type:1 save:test.png BRleg xtitle:'m(ll) [GeV]' ytitle:A_{FB} xmax:200 rebin:2");

}

void plot_DY_CSRecoil(){
  AFBSystPlotter aa("mi mi");
  aa.entries[0].title="CS frame";
  aa.entries[0].SetSuffix("CS");
  aa.entries[0].styles[0].linecolor=1;
  aa.entries[0].styles[0].linewidth=2;
  //aa.entries[0].replace["genfid_"]="correct_";
  
  aa.entries[1].title="Recoil frame";
  //aa.entries[1].SetSuffix("Recoil");
  aa.entries[1].replace["genfid_"]="";
  aa.entries[1].styles[0].linecolor=2;
  aa.entries[1].styles[0].linewidth=2;

  aa.SavePlot("[em][em]201[678][ab]?/nbjet/genfid_dimass","AFB logx type:1 save:test.png BRleg xtitle:'m(ll) [GeV]' ytitle:A_{FB} xmax:200 rebin:5");

}

void plot_TT_CSRecoil(){
  AFBSystPlotter aa("ttll ttll");
  aa.entries[0].title="CS frame";
  aa.entries[0].SetSuffix("CS");
  aa.entries[0].styles[0].linecolor=1;
  aa.entries[0].styles[0].linewidth=2;
  //aa.entries[0].replace["genfid_"]="correct_";
  
  aa.entries[1].title="Recoil frame";
  //aa.entries[1].SetSuffix("Recoil");
  aa.entries[1].replace["genfid_"]="";
  aa.entries[1].styles[0].linecolor=2;
  aa.entries[1].styles[0].linewidth=2;

  aa.SavePlot("[em][em]201[678][ab]?/nbjet/genfid_dimass","AFB logx type:1 save:test.png BRleg xtitle:'m(ll) [GeV]' ytitle:A_{FB} xmax:1000");

}

void frame(){
  //plot_DY_CS();
  plot_DY_CSRecoil();
  //plot_TT_CSRecoil();
}
