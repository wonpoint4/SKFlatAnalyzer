#include "Plotter/AFBPlotter.cc"

void page2(){
  AFBPlotter cf("appamm","CF");
  TLatex latex;
  latex.SetNDC();
  TCanvas *c=cf.DrawPlot("electron2016/m[60,120]/y[0.0,2.4]/pt[0,100]/dimass","widthweight widey nodata");
  c->cd(1);
  latex.DrawLatex(0.18,0.8,"#it{aMC@NLO} + #it{PYTHIA8}");
  //cf.DrawPlot("electron2017/m[60,120]/y[0.0,2.4]/pt[0,100]/dimass","widthweight widey nodata");
  //cf.DrawPlot("electron2018/m[60,120]/y[0.0,2.4]/pt[0,100]/dimass","widthweight widey nodata");
  AFBPlotter wr("wr");
  c=wr.DrawPlot("muon2016/m[400,5000]/y[0.0,2.4]/pt[100,1000]/dipt","widthweight norm nodata type:1");
  c=wr.DrawPlot("muon2016/m[400,5000]/y[0.0,2.4]/pt[100,1000]/AFB(pt)","nodata type:1");
}
void page5(){
  AFBPlotter aa("dajbs");
  TCanvas *c;
  c=aa.DrawPlot("muon2016/m[60,2000]/y[0.0,2.4]/pt[0,1000]/dimass","widthweight 1:logy norm logx widey");
  c=aa.DrawPlot("muon2016/m[60,2000]/y[0.0,2.4]/pt[0,1000]/AFB(m)","type:1 nodata bottomleg logx"); 
}
void page6(){
  AFBPlotter aa("dajb");
  TCanvas *c;
  c=aa.DrawPlot("electron2016/m[60,2000]/y[0.0,2.4]/pt[0,1000]/dimass","widthweight 1:logy norm logx widey");
  c=aa.DrawPlot("electron2016/m[60,2000]/y[0.0,2.4]/pt[0,1000]/AFB(m)","type:1 nodata bottomleg logx");
}
void page7(){
  AFBPlotter aa("dajbs");
  TCanvas *c;
  c=aa.DrawPlot("muon2016/m[60,2000]/y[0.0,2.4]/pt[0,1000]/dirap","widthweight 1:logy norm widey");
  c=aa.DrawPlot("muon2016/m[60,2000]/y[0.0,2.4]/pt[0,1000]/AFB(y)","type:1 nodata bottomleg"); 
}
void page8(){
  AFBPlotter aa("dajb");
  TCanvas *c;
  c=aa.DrawPlot("electron2016/m[60,2000]/y[0.0,2.4]/pt[0,1000]/dirap","widthweight 1:logy norm widey");
  c=aa.DrawPlot("electron2016/m[60,2000]/y[0.0,2.4]/pt[0,1000]/AFB(y)","type:1 nodata bottomleg"); 
}
void page9(){
  AFBPlotter aa("dajbs");
  TCanvas *c;
  c=aa.DrawPlot("muon2016/m[60,2000]/y[0.0,2.4]/pt[0,1000]/dipt","widthweight 1:logy norm widey");
  c=aa.DrawPlot("muon2016/m[60,2000]/y[0.0,2.4]/pt[0,1000]/AFB(pt)","type:1 nodata bottomleg"); 
}
void page10(){
  AFBPlotter aa("dajb");
  TCanvas *c;
  c=aa.DrawPlot("electron2016/m[60,2000]/y[0.0,2.4]/pt[0,1000]/AFB(m)","type:1 nodata bottomleg logx"); 
  c=aa.DrawPlot("electron2016/m[60,2000]/y[0.0,2.4]/pt[0,1000]/dipt","widthweight 1:logy norm widey");
}  

void page5_2(){
  TCanvas *c;
  AFBPlotter dajbs("dajbs");
  c=dajbs.DrawPlot("muon2018/m[60,2000]/y[0.0,2.4]/pt[0,1000]/dimass","widthweight 1:logy norm logx widey");
  /*
  AFBPlotter dajb("dajb");
  c=dajb.DrawPlot("electron2018/m[60,2000]/y[0.0,2.4]/pt[0,1000]/dimass","widthweight 1:logy norm logx widey");

  AFBPlotter dabs("dabs");
  c=dabs.DrawPlot("muon2017/m[60,2000]/y[0.0,2.4]/pt[0,1000]/dimass","widthweight 1:logy norm logx widey");

  AFBPlotter dab("dab");
  c=dab.DrawPlot("electron2017/m[60,2000]/y[0.0,2.4]/pt[0,1000]/dimass","widthweight 1:logy norm logx widey");
  */
}
