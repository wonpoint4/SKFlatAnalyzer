// first result of bjet region
void PlotAll(){
  AFBPlotter miel("miel");
  /*
  miel.DrawPlot("ee201[678][ab]?/m[52,3000]/nbjet","1:logy preliminary");

  miel.DrawPlot("ee201[678][ab]?/0bjet/m[52,3000]/dimass","norm rebin:2 widthweight logx 1:logy preliminary");
  miel.DrawPlot("ee201[678][ab]?/0bjet/m[52,3000]/AFB(m)","rebin:2 logx preliminary");
  miel.DrawPlot("ee201[678][ab]?/0bjet/m[52,80]/AFB(pt)","rebin:2 logx preliminary xmax:650");
  miel.DrawPlot("ee201[678][ab]?/0bjet/m[80,100]/AFB(pt)","rebin:2 logx preliminary xmax:650");
  miel.DrawPlot("ee201[678][ab]?/0bjet/m[100,3000]/AFB(pt)","rebin:2 logx preliminary xmax:650");
  */
  miel.DrawPlot("ee201[678][ab]?/nbjet/m[52,3000]/dimass","norm rebin:2 widthweight logx 1:logy preliminary");
  miel.DrawPlot("ee201[678][ab]?/nbjet/m[52,1500]/AFB(m)","rebin:2 logx preliminary");
  miel.DrawPlot("ee201[678][ab]?/nbjet/m[52,80]/AFB(pt)","rebin:{2,10,20,40,60,100,650} logx preliminary xmax:650");
  miel.DrawPlot("ee201[678][ab]?/nbjet/m[80,100]/AFB(pt)","rebin:{2,10,20,40,60,100,650} logx preliminary xmax:650");
  miel.DrawPlot("ee201[678][ab]?/nbjet/m[100,3000]/AFB(pt)","rebin:{2,10,20,40,60,100,650} logx preliminary xmax:650");

  AFBPlotter mimu("mimu");

  mimu.DrawPlot("mm201[678][ab]?/nbjet/m[52,3000]/dimass","norm rebin:2 widthweight logx 1:logy preliminary");
  mimu.DrawPlot("mm201[678][ab]?/nbjet/m[52,1500]/AFB(m)","rebin:2 logx preliminary");
  mimu.DrawPlot("mm201[678][ab]?/nbjet/m[52,80]/AFB(pt)","rebin:{2,10,20,40,60,100,650} logx preliminary xmax:650");
  mimu.DrawPlot("mm201[678][ab]?/nbjet/m[80,100]/AFB(pt)","rebin:{2,10,20,40,60,100,650} logx preliminary xmax:650");
  mimu.DrawPlot("mm201[678][ab]?/nbjet/m[100,3000]/AFB(pt)","rebin:{2,10,20,40,60,100,650} logx preliminary xmax:650");

}


  
  
