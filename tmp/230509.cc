#include"AFBPlotter.cc"

void plot1(){
  AFBPlotter aa("mi+bx_mi bx_mi mi");
  aa.plotdir="fig/230509";
  aa.entries[0].title="DY";
  aa.entries[0].styles[0].linecolor=1;
  aa.entries[0].styles[0].linewidth=2;
  aa.entries[0].styles[0].drawoption="hist e";
  aa.entries[1].title="DY b-init";
  aa.entries[1].styles[0].linecolor=2;
  aa.entries[2].title="DY no b-init";
  aa.entries[2].styles[0].linecolor=4;
  aa.SavePlot("mm201[678][ab]?/nbjet/m[52,500]/AFB(m)","logx rebin:4 ymin:-0.049 ymax:0.199 type:1 text:0.2,0.8,'CS frame' text:0.2,0.73,'N_{bjet}>0' save:dy_cs");
  aa.SavePlot("mm201[678][ab]?/nbjet/m[52,500]/AFBRecoil(m)","logx rebin:4 ymin:-0.049 ymax:0.199 type:1 text:0.2,0.8,'Recoil frame' text:0.2,0.73,'N_{bjet}>0' save:dy_recoil");
}
void plot2(){
  AFBPlotter aa("ttll bbllpll bbllplr bbllprl bbllprr");
  aa.plotdir="fig/230509";
  aa.entries[0].styles[0].linewidth=2;
  aa.SavePlot("[em][em]2018/nbjet/m[280,3000]/AFB(m)","logx ymin:-0.19 ymax:0.49 type:1 text:0.2,0.8,'CS frame' text:0.2,0.73,'N_{bjet}>0' save:bbll_cs rebin:{280,400,800,1000,3000}");
  aa.SavePlot("[em][em]2018/nbjet/m[280,3000]/AFBRecoil(m)","logx ymin:-0.19 ymax:0.49 type:1 text:0.2,0.8,'Recoil frame' text:0.2,0.73,'N_{bjet}>0' save:bbll_recoil rebin:{280,400,800,1000,3000}");
}  
void plot3(){
  AFBPlotter aa("ttll ttll");
  aa.plotdir="fig/230509";
  aa.entries[0].title+=" (CS frame)";
  aa.entries[1].title+=" (Recoil frame)";
  aa.entries[1].styles[0].linecolor=kBlue;
  aa.entries[1].AddTag("recoil");
  aa.SavePlot("mm2018/nbjet/m[52,3000]/AFB(m)","logx ymin:-0.09 ymax:0.59 type:1 replace:costhetaCS->costhetaRecoil:recoil text:0.2,0.73,'N_{bjet}>0' save:ttll");
}  
void plot(){
  plot1();
  plot2();
  plot3();
}
