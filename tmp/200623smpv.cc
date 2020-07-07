{
  AFBPlotter aa;
  aa.plotdir="smpv/low";
  aa.SavePlot("ee2016_dimass","histname:ee2016/costhetaCS project:x Xmin:60 Xmax:120 Ymin:-2.4 Ymax:2.4 widthweight");
  aa.SavePlot("ee2017_dimass","histname:ee2017/costhetaCS project:x Xmin:60 Xmax:120 Ymin:-2.4 Ymax:2.4 widthweight");
  aa.SavePlot("ee2018_dimass","histname:ee2018/costhetaCS project:x Xmin:60 Xmax:120 Ymin:-2.4 Ymax:2.4 widthweight");
  aa.SavePlot("mm2016_dimass","histname:mm2016/costhetaCS project:x Xmin:60 Xmax:120 Ymin:-2.4 Ymax:2.4 widthweight");
  aa.SavePlot("mm2017_dimass","histname:mm2017/costhetaCS project:x Xmin:60 Xmax:120 Ymin:-2.4 Ymax:2.4 widthweight");
  aa.SavePlot("mm2018_dimass","histname:mm2018/costhetaCS project:x Xmin:60 Xmax:120 Ymin:-2.4 Ymax:2.4 widthweight");

  aa.SavePlot("ee2016_dirap","histname:ee2016/costhetaCS project:y Xmin:60 Xmax:120 Ymin:-2.4 Ymax:2.4 BMleg");
  aa.SavePlot("ee2017_dirap","histname:ee2017/costhetaCS project:y Xmin:60 Xmax:120 Ymin:-2.4 Ymax:2.4 BMleg");
  aa.SavePlot("ee2018_dirap","histname:ee2018/costhetaCS project:y Xmin:60 Xmax:120 Ymin:-2.4 Ymax:2.4 BMleg");
  aa.SavePlot("mm2016_dirap","histname:mm2016/costhetaCS project:y Xmin:60 Xmax:120 Ymin:-2.4 Ymax:2.4 BMleg");
  aa.SavePlot("mm2017_dirap","histname:mm2017/costhetaCS project:y Xmin:60 Xmax:120 Ymin:-2.4 Ymax:2.4 BMleg");
  aa.SavePlot("mm2018_dirap","histname:mm2018/costhetaCS project:y Xmin:60 Xmax:120 Ymin:-2.4 Ymax:2.4 BMleg");

  aa.SavePlot("ee2016_dipt","histname:ee2016/costhetaCS project:z Xmin:60 Xmax:120 Ymin:-2.4 Ymax:2.4 Zmax:650 widthweight 1:logy");
  aa.SavePlot("ee2017_dipt","histname:ee2017/costhetaCS project:z Xmin:60 Xmax:120 Ymin:-2.4 Ymax:2.4 Zmax:650 widthweight 1:logy");
  aa.SavePlot("ee2018_dipt","histname:ee2018/costhetaCS project:z Xmin:60 Xmax:120 Ymin:-2.4 Ymax:2.4 Zmax:650 widthweight 1:logy");
  aa.SavePlot("mm2016_dipt","histname:mm2016/costhetaCS project:z Xmin:60 Xmax:120 Ymin:-2.4 Ymax:2.4 Zmax:650 widthweight 1:logy");
  aa.SavePlot("mm2017_dipt","histname:mm2017/costhetaCS project:z Xmin:60 Xmax:120 Ymin:-2.4 Ymax:2.4 Zmax:650 widthweight 1:logy");
  aa.SavePlot("mm2018_dipt","histname:mm2018/costhetaCS project:z Xmin:60 Xmax:120 Ymin:-2.4 Ymax:2.4 Zmax:650 widthweight 1:logy");

  aa.SavePlot("ee2016_cost","histname:ee2016/costhetaCS project:u Xmin:60 Xmax:120 Ymin:-2.4 Ymax:2.4 widthweight BMleg");
  aa.SavePlot("ee2017_cost","histname:ee2017/costhetaCS project:u Xmin:60 Xmax:120 Ymin:-2.4 Ymax:2.4 widthweight BMleg");
  aa.SavePlot("ee2018_cost","histname:ee2018/costhetaCS project:u Xmin:60 Xmax:120 Ymin:-2.4 Ymax:2.4 widthweight BMleg");
  aa.SavePlot("mm2016_cost","histname:mm2016/costhetaCS project:u Xmin:60 Xmax:120 Ymin:-2.4 Ymax:2.4 widthweight BMleg");
  aa.SavePlot("mm2017_cost","histname:mm2017/costhetaCS project:u Xmin:60 Xmax:120 Ymin:-2.4 Ymax:2.4 widthweight BMleg");
  aa.SavePlot("mm2018_cost","histname:mm2018/costhetaCS project:u Xmin:60 Xmax:120 Ymin:-2.4 Ymax:2.4 widthweight BMleg");

  AFBPlotter bb("data ^amcM+tau_amcM+vv+wjets+tt+tw");
  bb.plotdir="smpv/high";
  bb.SavePlot("ee2016_dimass","histname:ee2016/costhetaCS project:x Xmin:120 Xmax:3000 Ymin:-2.4 Ymax:2.4 widthweight 2:widey 1:logy");
  bb.SavePlot("ee2017_dimass","histname:ee2017/costhetaCS project:x Xmin:120 Xmax:3000 Ymin:-2.4 Ymax:2.4 widthweight 2:widey 1:logy");
  bb.SavePlot("ee2018_dimass","histname:ee2018/costhetaCS project:x Xmin:120 Xmax:3000 Ymin:-2.4 Ymax:2.4 widthweight 2:widey 1:logy");
  bb.SavePlot("mm2016_dimass","histname:mm2016/costhetaCS project:x Xmin:120 Xmax:3000 Ymin:-2.4 Ymax:2.4 widthweight 2:widey 1:logy");
  bb.SavePlot("mm2017_dimass","histname:mm2017/costhetaCS project:x Xmin:120 Xmax:3000 Ymin:-2.4 Ymax:2.4 widthweight 2:widey 1:logy");
  bb.SavePlot("mm2018_dimass","histname:mm2018/costhetaCS project:x Xmin:120 Xmax:3000 Ymin:-2.4 Ymax:2.4 widthweight 2:widey 1:logy");

  bb.SavePlot("ee2016_dirap","histname:ee2016/costhetaCS project:y Xmin:120 Xmax:3000 Ymin:-2.4 Ymax:2.4 2:widey");
  bb.SavePlot("ee2017_dirap","histname:ee2017/costhetaCS project:y Xmin:120 Xmax:3000 Ymin:-2.4 Ymax:2.4 2:widey");
  bb.SavePlot("ee2018_dirap","histname:ee2018/costhetaCS project:y Xmin:120 Xmax:3000 Ymin:-2.4 Ymax:2.4 2:widey");
  bb.SavePlot("mm2016_dirap","histname:mm2016/costhetaCS project:y Xmin:120 Xmax:3000 Ymin:-2.4 Ymax:2.4 2:widey");
  bb.SavePlot("mm2017_dirap","histname:mm2017/costhetaCS project:y Xmin:120 Xmax:3000 Ymin:-2.4 Ymax:2.4 2:widey");
  bb.SavePlot("mm2018_dirap","histname:mm2018/costhetaCS project:y Xmin:120 Xmax:3000 Ymin:-2.4 Ymax:2.4 2:widey");

  bb.SavePlot("ee2016_dipt","histname:ee2016/costhetaCS project:z Xmin:120 Xmax:3000 Ymin:-2.4 Ymax:2.4 Zmax:650 widthweight 2:widey 1:logy");
  bb.SavePlot("ee2017_dipt","histname:ee2017/costhetaCS project:z Xmin:120 Xmax:3000 Ymin:-2.4 Ymax:2.4 Zmax:650 widthweight 2:widey 1:logy");
  bb.SavePlot("ee2018_dipt","histname:ee2018/costhetaCS project:z Xmin:120 Xmax:3000 Ymin:-2.4 Ymax:2.4 Zmax:650 widthweight 2:widey 1:logy");
  bb.SavePlot("mm2016_dipt","histname:mm2016/costhetaCS project:z Xmin:120 Xmax:3000 Ymin:-2.4 Ymax:2.4 Zmax:650 widthweight 2:widey 1:logy");
  bb.SavePlot("mm2017_dipt","histname:mm2017/costhetaCS project:z Xmin:120 Xmax:3000 Ymin:-2.4 Ymax:2.4 Zmax:650 widthweight 2:widey 1:logy");
  bb.SavePlot("mm2018_dipt","histname:mm2018/costhetaCS project:z Xmin:120 Xmax:3000 Ymin:-2.4 Ymax:2.4 Zmax:650 widthweight 2:widey 1:logy");
  
  bb.SavePlot("ee2016_cost","histname:ee2016/costhetaCS project:u Xmin:120 Xmax:3000 Ymin:-2.4 Ymax:2.4 2:widey BMleg");
  bb.SavePlot("ee2017_cost","histname:ee2017/costhetaCS project:u Xmin:120 Xmax:3000 Ymin:-2.4 Ymax:2.4 2:widey BMleg");
  bb.SavePlot("ee2018_cost","histname:ee2018/costhetaCS project:u Xmin:120 Xmax:3000 Ymin:-2.4 Ymax:2.4 2:widey BMleg");
  bb.SavePlot("mm2016_cost","histname:mm2016/costhetaCS project:u Xmin:120 Xmax:3000 Ymin:-2.4 Ymax:2.4 2:widey BMleg");
  bb.SavePlot("mm2017_cost","histname:mm2017/costhetaCS project:u Xmin:120 Xmax:3000 Ymin:-2.4 Ymax:2.4 2:widey BMleg");
  bb.SavePlot("mm2018_cost","histname:mm2018/costhetaCS project:u Xmin:120 Xmax:3000 Ymin:-2.4 Ymax:2.4 2:widey BMleg");
}
