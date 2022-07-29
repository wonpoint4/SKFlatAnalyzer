#include"AFBPlotter.cc"
void plot_qcdss(){
  TString common="";
  common+=" preliminary";
  AFBPlotter mmss("1.8*ss_mi 2.1*ss_mi 1.5*ss_mi");
  mmss.entries[0].title="QCD (SS method)";
  mmss.entries[0].styles[0].linewidth=2;
  mmss.entries[1].title="up";
  mmss.entries[1].styles[0].markerstyle=0;
  mmss.entries[1].styles[0].drawoption="HIST";
  mmss.entries[2].title="down";
  mmss.entries[2].styles[0].markerstyle=0;
  mmss.entries[2].styles[0].drawoption="HIST";
  mmss.DrawPlot("mm201[678][ab]?/[0n]bjet/m[52,3000]/dimass","logx widthweight 2:widewidey "+common);
  AFBPlotter eess("ss_mi ss_mi ss_mi");
  eess.entries[0].title="QCD (SS method)";
  eess.entries[0].styles[0].linewidth=2;
  eess.entries[1].title="up";
  eess.entries[1].subs[1].weight=-0.97;
  eess.entries[1].styles[0].markerstyle=0;
  eess.entries[1].styles[0].drawoption="HIST";
  eess.entries[2].title="down";
  eess.entries[2].subs[1].weight=-1.03;
  eess.entries[2].styles[0].markerstyle=0;
  eess.entries[2].styles[0].drawoption="HIST";
  eess.DrawPlot("ee201[678][ab]?/[0n]bjet/m[52,3000]/dimass","logx widthweight 2:widewidey "+common);
}  
void plot_qcd(){
  TString common=" base:3";
  common+=" preliminary";
  AFBPlotter mm("1.8*ss_mi 2.1*ss_mi 1.5*ss_mi fake_mi");
  mm.entries[0].title="QCD (SS method)";
  mm.entries[0].styles[0].linewidth=2;
  mm.entries[1].title="up";
  mm.entries[1].styles[0].markerstyle=0;
  mm.entries[1].styles[0].drawoption="HIST";
  mm.entries[2].title="down";
  mm.entries[2].styles[0].markerstyle=0;
  mm.entries[2].styles[0].drawoption="HIST";
  mm.entries[3].title="QCD (Fake-rate method)";
  mm.entries[3].styles[0].linewidth=2;

  mm.DrawPlot("mm201[678][ab]?/[0n]bjet/m[52,3000]/dimass","logx widthweight 2:widewidey"+common);
  mm.DrawPlot("mm201[678][ab]?/0bjet/m[52,3000]/dimass","logx widthweight 2:widewidey"+common);
  mm.DrawPlot("mm201[678][ab]?/nbjet/m[52,3000]/dimass","logx widthweight 2:widewidey"+common);

  mm.DrawPlot("mm201[678][ab]?/[0n]bjet/m[52,3000]/lpt","logx widthweight 2:widewidey xtitle:'p_{T}(#mu) [GeV]'"+common);
  mm.DrawPlot("mm201[678][ab]?/0bjet/m[52,3000]/lpt","logx widthweight 2:widewidey xtitle:'p_{T}(#mu) [GeV]'"+common);
  mm.DrawPlot("mm201[678][ab]?/nbjet/m[52,3000]/lpt","logx widthweight 2:widewidey xtitle:'p_{T}(#mu) [GeV]'"+common);

  mm.DrawPlot("mm201[678][ab]?/[0n]bjet/m[52,3000]/leta","widthweight 2:widewidey xtitle:'#eta_{#mu}'"+common);
  mm.DrawPlot("mm201[678][ab]?/0bjet/m[52,3000]/leta","widthweight 2:widewidey xtitle:'#eta_{#mu}'"+common);
  mm.DrawPlot("mm201[678][ab]?/nbjet/m[52,3000]/leta","widthweight 2:widewidey xtitle:'#eta_{#mu}'"+common);
  
  AFBPlotter ee("ss_mi ss_mi ss_mi fake_mi");
  ee.entries[0].title="QCD (SS method)";
  ee.entries[0].styles[0].linewidth=2;
  ee.entries[1].title="up";
  ee.entries[1].subs[1].weight=-0.97;
  ee.entries[1].styles[0].markerstyle=0;
  ee.entries[1].styles[0].drawoption="HIST";
  ee.entries[2].title="down";
  ee.entries[2].subs[1].weight=-1.03;
  ee.entries[2].styles[0].markerstyle=0;
  ee.entries[2].styles[0].drawoption="HIST";
  ee.entries[3].title="QCD (Fake-rate method)";
  ee.entries[3].styles[0].linewidth=2;
  
  ee.DrawPlot("ee201[678][ab]?/[0n]bjet/m[52,3000]/dimass","logx widthweight 2:widewidey"+common);
  ee.DrawPlot("ee201[678][ab]?/0bjet/m[52,3000]/dimass","logx widthweight 2:widewidey"+common);
  ee.DrawPlot("ee201[678][ab]?/nbjet/m[52,3000]/dimass","logx widthweight 2:widewidey"+common);

  ee.DrawPlot("ee201[678][ab]?/[0n]bjet/m[52,3000]/lpt","logx widthweight 2:widewidey xtitle:'p_{T}(e) [GeV]'"+common);
  ee.DrawPlot("ee201[678][ab]?/0bjet/m[52,3000]/lpt","logx widthweight 2:widewidey xtitle:'p_{T}(e) [GeV]'"+common);
  ee.DrawPlot("ee201[678][ab]?/nbjet/m[52,3000]/lpt","logx widthweight 2:widewidey xtitle:'p_{T}(e) [GeV]'"+common);

  ee.DrawPlot("ee201[678][ab]?/[0n]bjet/m[52,3000]/leta","widthweight 2:widewidey xtitle:'#eta_{e}'"+common);
  ee.DrawPlot("ee201[678][ab]?/0bjet/m[52,3000]/leta","widthweight 2:widewidey xtitle:'#eta_{e}'"+common);
  ee.DrawPlot("ee201[678][ab]?/nbjet/m[52,3000]/leta","widthweight 2:widewidey xtitle:'#eta_{e}'"+common);

}
void plot_fake(){
  TString common=" 2:widewidey xmax:500 xmin:10";
  //TString common="";
  AFBPlotter fake("data ^mi+tau_mi+wjets+vv+tttw+aa","FakeAnalyzer");
  fake.entries[1].styles[0].drawoption="hist e";
  fake.entries[1].styles[0].fillstyle=0;
  fake.DrawPlot("EE201[678][ab]?/al[01]etapt","project:y xtitle:'p_{T}(e) [GeV]' logx widthweight base:0"+common);
  fake.DrawPlot("eE201[678][ab]?/l0etapt","project:y xtitle:'p_{T}(e) [GeV]' logx widthweight base:0"+common);
  fake.DrawPlot("MM201[678][ab]?/al[01]etapt","project:y xtitle:'p_{T}(#mu) [GeV]' logx widthweight base:0"+common);
  fake.DrawPlot("mM201[678][ab]?/l0etapt","project:y xtitle:'p_{T}(#mu) [GeV]' logx widthweight base:0"+common);

  fake.DrawPlot("EE201[678][ab]?/al[01]etapt_noZ","project:y xtitle:'p_{T}(e) [GeV]' logx widthweight base:0"+common);
  fake.DrawPlot("eE201[678][ab]?/l0etapt_noZ","project:y xtitle:'p_{T}(e) [GeV]' logx widthweight base:0"+common);
  fake.DrawPlot("MM201[678][ab]?/al[01]etapt_noZ","project:y xtitle:'p_{T}(#mu) [GeV]' logx widthweight base:0"+common);
  fake.DrawPlot("mM201[678][ab]?/l0etapt_noZ","project:y xtitle:'p_{T}(#mu) [GeV]' logx widthweight base:0"+common);

  fake.DrawPlot("EE201[678][ab]?/ss_al[01]etapt_noZ","project:y xtitle:'p_{T}(e) [GeV]' logx widthweight base:0"+common);
  fake.DrawPlot("eE201[678][ab]?/ss_l0etapt_noZ","project:y xtitle:'p_{T}(e) [GeV]' logx widthweight base:0"+common);
  fake.DrawPlot("MM201[678][ab]?/ss_al[01]etapt_noZ","project:y xtitle:'p_{T}(#mu) [GeV]' logx widthweight base:0"+common);
  fake.DrawPlot("mM201[678][ab]?/ss_l0etapt_noZ","project:y xtitle:'p_{T}(#mu) [GeV]' logx widthweight base:0"+common);

  AFBPlotter fake_mi("fake_mi");
  
  

}

void plot(){
  plot_qcd();
  plot_fake();
}
