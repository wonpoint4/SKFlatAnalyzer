#include <thread>

void Func(TString setup,TString histname,TString option){
  TString savename=setup+"_"+histname;
  option.ReplaceAll("'","'\"'\"'");
  TPRegexp("[][?()]").Substitute(savename,"","g");
  TPRegexp("[/,]").Substitute(savename,"_","g");
  savename="kps/"+savename;
  if(option.Contains("nodefault")){
    option.ReplaceAll("nodefault","");
  }else{
    option+=" preliminary sysname:totalsys";
    if(!histname.Contains("AFB")) option+=" sysdetail 2:sysleg";
    else option+=" noleg sysleg type:1 nolumi";
  }
  system("export ROOT_HIST=0; echo 'gROOT->ProcessLine(\".L Plotter/AFBPlotter.cc\");\n AFBPlotter mimu(\""+setup+"\");\n mimu.SavePlot(\""+savename+"\",\"histname:"+histname+" "+option+"\"); \n'|root -l -b");
  //system("export ROOT_HIST=0; echo 'gROOT->ProcessLine(\".L Plotter/AFBPlotter.cc\");\n AFBPlotter mimu(\""+setup+"\",\"AFBAnalyzer_backup\");\n mimu.SavePlot(\""+savename+"\",\"histname:"+histname+" "+option+"\"); \n'|root -l -b");
}
thread Run(TString setup,TString histname,TString option){
  int interval=1;
  while(gSystem->GetFromPipe(Form("mpstat %d 1|tail -n1|awk '{print $NF}'",interval)).Atof()<30){
    interval=10;
  }
  return thread(Func,setup,histname,option);
}

void PlotAll(){
  vector<thread> threads;

  threads.push_back(Run("mimu","mm201[678][ab]?/[0n]bjet/m[52,3000]/nbjet","1:logy 2:widey xtitle:'Number of bjets'"));
  threads.push_back(Run("mimu","mm201[678][ab]?/0bjet/m[52,3000]/dimass","logx rebin:{52,56,60,65,70,74,77,80,82,84,86,88,89,90,91,92,93,94,96,98,100,103,106,110,115,120,130,140,150,175,200,240,280,340,400,500,600,700,800,1000,3000} norm widthweight 1:logy"));
  threads.push_back(Run("mimu","mm201[678][ab]?/0bjet/m[52,3000]/dirap","norm widthweight xmin:-2.4 xmax:2.4 noleg"));
  threads.push_back(Run("mimu","mm201[678][ab]?/0bjet/m[52,3000]/dipt","norm widthweight logx xmax:650"));
  threads.push_back(Run("mimu","mm201[678][ab]?/0bjet/m[52,3000]/lpt","norm widthweight logx xmin:10 xmax:500 xtitle:'p_{T}(#mu) [GeV]' 1:logy"));
  threads.push_back(Run("mimu","mm201[678][ab]?/0bjet/m[52,3000]/leta","norm widthweight xtitle:#eta(#mu) noleg"));
  threads.push_back(Run("mimu_nodata","mm201[678][ab]?/0bjet/m[52,3000]/AFB(m)","logx rebin:{52,56,60,65,70,74,77,80,82,84,86,88,89,90,91,92,93,94,96,98,100,103,106,110,115,120,130,140,150,175,200,240,280,340,400,500,600,700,800,1000,3000} TLleg ymin:-0.19 ymax:0.39 blind:76,106"));
  threads.push_back(Run("mimu_nodata","mm201[678][ab]?/0bjet/m[52,80]/AFB(pt)","logx rebin:{0,2,4,6,8,10,12,14,16,20,32,60,120,650} xmax:650 ymax:0.098 ymin:-0.098"));
  threads.push_back(Run("mimu_nodata","mm201[678][ab]?/0bjet/m[80,100]/AFB(pt)","logx rebin:{0,2,4,6,8,10,12,14,16,20,32,60,120,650} xmax:650 ymax:0.098 ymin:-0.098"));
  threads.push_back(Run("mimu_nodata","mm201[678][ab]?/0bjet/m[100,3000]/AFB(pt)","logx rebin:{0,2,4,6,8,10,12,14,16,20,32,60,120,650} xmax:650 ymax:0.098 ymin:-0.098 BRleg"));
  threads.push_back(Run("mimu","mm201[678][ab]?/nbjet/m[52,3000]/dimass","logx rebin:{52,56,60,65,70,74,77,80,82,84,86,88,89,90,91,92,93,94,96,98,100,103,106,110,115,120,130,140,150,175,200,240,280,340,400,500,600,700,800,1000,3000} norm widthweight 2:widey"));
  threads.push_back(Run("mimu","mm201[678][ab]?/nbjet/m[52,3000]/dirap","norm widthweight xmin:-2.4 xmax:2.4 noleg"));
  threads.push_back(Run("mimu","mm201[678][ab]?/nbjet/m[52,3000]/dipt","norm widthweight logx xmax:650 2:widey"));
  threads.push_back(Run("mimu","mm201[678][ab]?/nbjet/m[52,3000]/lpt","norm widthweight logx xmin:10 xmax:500 xtitle:'p_{T}(#mu) [GeV]'"));
  threads.push_back(Run("mimu","mm201[678][ab]?/nbjet/m[52,3000]/leta","norm widthweight xtitle:#eta(#mu) noleg"));
  threads.push_back(Run("mimu_nodata","mm201[678][ab]?/nbjet/m[52,3000]/AFB(m)","logx rebin:{52,60,70,77,82,86,89,91,93,96,100,106,115,130,150,200,280,400,600,3000} TLleg ymin:-0.099 ymax:0.099 blind:76,106"));
  threads.push_back(Run("mimu_nodata","mm201[678][ab]?/nbjet/m[52,80]/AFB(pt)","logx rebin:{0,2,8,14,20,32,60,120,650} xmax:650 ymax:0.098 ymin:-0.098"));
  threads.push_back(Run("mimu_nodata","mm201[678][ab]?/nbjet/m[80,100]/AFB(pt)","logx rebin:{0,2,8,14,20,32,60,120,650} xmax:650 ymax:0.098 ymin:-0.098"));
  threads.push_back(Run("mimu_nodata","mm201[678][ab]?/nbjet/m[100,3000]/AFB(pt)","logx rebin:{0,2,8,14,20,32,60,120,650} xmax:650 ymax:0.098 ymin:-0.098"));

  threads.push_back(Run("miel","ee201[678][ab]?/[0n]bjet/m[52,3000]/nbjet","1:logy 2:widey xtitle:'Number of bjets'"));
  threads.push_back(Run("miel","ee201[678][ab]?/0bjet/m[52,3000]/dimass","logx rebin:{52,56,60,65,70,74,77,80,82,84,86,88,89,90,91,92,93,94,96,98,100,103,106,110,115,120,130,140,150,175,200,240,280,340,400,500,600,700,800,1000,3000} norm widthweight 1:logy"));
  threads.push_back(Run("miel","ee201[678][ab]?/0bjet/m[52,3000]/dirap","norm widthweight xmin:-2.4 xmax:2.4 noleg"));
  threads.push_back(Run("miel","ee201[678][ab]?/0bjet/m[52,3000]/dipt","norm widthweight logx xmax:650"));
  threads.push_back(Run("miel","ee201[678][ab]?/0bjet/m[52,3000]/lpt","norm widthweight logx xmin:10 xmax:500 xtitle:'p_{T}(e) [GeV]' 1:logy"));
  threads.push_back(Run("miel","ee201[678][ab]?/0bjet/m[52,3000]/leta","norm widthweight xmin:-2.4 xmax:2.4 xtitle:#eta(e) noleg"));
  threads.push_back(Run("miel_nodata","ee201[678][ab]?/0bjet/m[52,3000]/AFB(m)","logx rebin:{52,56,60,65,70,74,77,80,82,84,86,88,89,90,91,92,93,94,96,98,100,103,106,110,115,120,130,140,150,175,200,240,280,340,400,500,600,700,800,1000,3000} TLleg ymin:-0.19 ymax:0.39 blind:76,106"));
  threads.push_back(Run("miel_nodata","ee201[678][ab]?/0bjet/m[52,80]/AFB(pt)","logx rebin:{0,2,4,6,8,10,12,14,16,20,32,60,120,650} xmax:650 ymax:0.098 ymin:-0.098"));
  threads.push_back(Run("miel_nodata","ee201[678][ab]?/0bjet/m[80,100]/AFB(pt)","logx rebin:{0,2,4,6,8,10,12,14,16,20,32,60,120,650} xmax:650 ymax:0.098 ymin:-0.098"));
  threads.push_back(Run("miel_nodata","ee201[678][ab]?/0bjet/m[100,3000]/AFB(pt)","logx rebin:{0,2,4,6,8,10,12,14,16,20,32,60,120,650} xmax:650 ymax:0.098 ymin:-0.098 BRleg"));
  threads.push_back(Run("miel","ee201[678][ab]?/nbjet/m[52,3000]/dimass","logx rebin:{52,56,60,65,70,74,77,80,82,84,86,88,89,90,91,92,93,94,96,98,100,103,106,110,115,120,130,140,150,175,200,240,280,340,400,500,600,700,800,1000,3000} norm widthweight 2:widey"));
  threads.push_back(Run("miel","ee201[678][ab]?/nbjet/m[52,3000]/dirap","norm widthweight xmin:-2.4 xmax:2.4 noleg"));
  threads.push_back(Run("miel","ee201[678][ab]?/nbjet/m[52,3000]/dipt","norm widthweight logx xmax:650 2:widey"));
  threads.push_back(Run("miel","ee201[678][ab]?/nbjet/m[52,3000]/lpt","norm widthweight logx xmin:10 xmax:500 xtitle:'p_{T}(e) [GeV]'"));
  threads.push_back(Run("miel","ee201[678][ab]?/nbjet/m[52,3000]/leta","norm widthweight xmin:-2.4 xmax:2.4 xtitle:#eta(e) noleg"));
  threads.push_back(Run("miel_nodata","ee201[678][ab]?/nbjet/m[52,3000]/AFB(m)","logx rebin:{52,60,70,77,82,86,89,91,93,96,100,106,115,130,150,200,280,400,600,3000} TLleg ymin:-0.099 ymax:0.099 blind:76,106"));
  threads.push_back(Run("miel_nodata","ee201[678][ab]?/nbjet/m[52,80]/AFB(pt)","logx rebin:{0,2,8,14,20,32,60,120,650} xmax:650 ymax:0.098 ymin:-0.098"));
  threads.push_back(Run("miel_nodata","ee201[678][ab]?/nbjet/m[80,100]/AFB(pt)","logx rebin:{0,2,8,14,20,32,60,120,650} xmax:650 ymax:0.098 ymin:-0.098"));
  threads.push_back(Run("miel_nodata","ee201[678][ab]?/nbjet/m[100,3000]/AFB(pt)","logx rebin:{0,2,8,14,20,32,60,120,650} xmax:650 ymax:0.098 ymin:-0.098"));

  threads.push_back(Run("amcmu","mm201[678][ab]?/0bjet/m[52,3000]/dimass","logx rebin:{52,56,60,65,70,74,77,80,82,84,86,88,89,90,91,92,93,94,96,98,100,103,106,110,115,120,130,140,150,175,200,240,280,340,400,500,600,700,800,1000,3000} norm widthweight 1:logy"));
  threads.push_back(Run("amcmu","mm201[678][ab]?/0bjet/m[52,3000]/lpt","norm widthweight logx xmin:10 xmax:500 xtitle:'p_{T}(#mu) [GeV]' 1:logy"));
  threads.push_back(Run("mimu","mm2016a/0bjet/m[80,100]/leta","widthweight xtitle:#eta(#mu) preliminary sysname:prefireweight 2:sysleg 1:BMleg nodefault"));
  threads.push_back(Run("mimu","mm2016a/0bjet/m[80,100]/leta(u)","suffix:_noprefireweight:sim widthweight xtitle:#eta(#mu) preliminary 1:BMleg nodefault"));

  for(auto& t:threads) t.join();

}
