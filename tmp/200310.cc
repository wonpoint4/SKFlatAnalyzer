void genlevel(){
  AFBPlotter aa("gen_amc");
  aa.DrawPlot("mm");
}
void ee(){
  AFBPlotter ee("ee2018 ee2017 ee2016");
  ee.SavePlot("ee201[6-8]/m[60,120]/y[-2.4,2.4]/pt[0,650]/dimass","norm widthweight");
  ee.SavePlot("ee201[6-8]/m[60,120]/y[-2.4,2.4]/pt[4,650]/dipt","norm logx widthweight");
  ee.SavePlot("ee201[6-8]/m[60,120]/y[-2.4,2.4]/pt[0,650]/dirap","norm 1:BMleg");
  ee.SavePlot("ee201[6-8]/m[60,120]/y[-2.4,2.4]/pt[0,650]/AFB(m)","1:BRleg");
  ee.SavePlot("ee201[6-8]/m[60,120]/y[-2.4,2.4]/pt[0,650]/AFB(y)","1:BRleg");
  ee.SavePlot("ee201[6-8]/m[60,120]/y[-2.4,2.4]/pt[4,650]/AFB(pt)","logx 1:BLleg rebin:2");

  ee.SavePlot("ee201[6-8]/m[120,3000]/y[-2.4,2.4]/pt[0,650]/dimass","norm logx widthweight");
  ee.SavePlot("ee201[6-8]/m[120,3000]/y[-2.4,2.4]/pt[4,650]/dipt","norm logx widthweight");
  ee.SavePlot("ee201[6-8]/m[120,3000]/y[-2.4,2.4]/pt[0,650]/dirap","norm 1:BMleg");
  ee.SavePlot("ee201[6-8]/m[120,3000]/y[-2.4,2.4]/pt[0,650]/AFB(m)","1:TLleg logx");
  ee.SavePlot("ee201[6-8]/m[120,3000]/y[-2.4,2.4]/pt[0,650]/AFB(y)","1:TMleg");
  ee.SavePlot("ee201[6-8]/m[120,3000]/y[-2.4,2.4]/pt[4,650]/AFB(pt)","logx 1:BLleg rebin:2");
}
void mm(){
  AFBPlotter mm("mm2018 mm2017 mm2016");
  mm.SavePlot("mm201[6-8]/m[60,120]/y[-2.4,2.4]/pt[0,650]/dimass","norm widthweight");
  mm.SavePlot("mm201[6-8]/m[60,120]/y[-2.4,2.4]/pt[4,650]/dipt","norm logx widthweight");
  mm.SavePlot("mm201[6-8]/m[60,120]/y[-2.4,2.4]/pt[0,650]/dirap","norm 1:BMleg");
  mm.SavePlot("mm201[6-8]/m[60,120]/y[-2.4,2.4]/pt[0,650]/AFB(m)","1:BRleg");
  mm.SavePlot("mm201[6-8]/m[60,120]/y[-2.4,2.4]/pt[0,650]/AFB(y)","1:TMleg");
  mm.SavePlot("mm201[6-8]/m[60,120]/y[-2.4,2.4]/pt[4,650]/AFB(pt)","logx 1:TLleg rebin:2");

  mm.SavePlot("mm201[6-8]/m[120,3000]/y[-2.4,2.4]/pt[0,650]/dimass","norm logx widthweight");
  mm.SavePlot("mm201[6-8]/m[120,3000]/y[-2.4,2.4]/pt[4,650]/dipt","norm logx widthweight");
  mm.SavePlot("mm201[6-8]/m[120,3000]/y[-2.4,2.4]/pt[0,650]/dirap","norm 1:BMleg");
  mm.SavePlot("mm201[6-8]/m[120,3000]/y[-2.4,2.4]/pt[0,650]/AFB(m)","1:TLleg logx");
  mm.SavePlot("mm201[6-8]/m[120,3000]/y[-2.4,2.4]/pt[0,650]/AFB(y)","1:TMleg");
  mm.SavePlot("mm201[6-8]/m[120,3000]/y[-2.4,2.4]/pt[4,650]/AFB(pt)","logx 1:BLleg rebin:2");
}

void mmee(){
  AFBPlotter mmee("mm ee");
  mmee.SavePlot("[em][em]201[6-8]/m[60,120]/y[-2.4,2.4]/pt[0,650]/dimass","norm widthweight");
  mmee.SavePlot("[em][em]201[6-8]/m[60,120]/y[-2.4,2.4]/pt[4,650]/dipt","norm logx widthweight");
  mmee.SavePlot("[em][em]201[6-8]/m[60,120]/y[-2.4,2.4]/pt[0,650]/dirap","norm 1:BMleg");
  mmee.SavePlot("[em][em]201[6-8]/m[60,120]/y[-2.4,2.4]/pt[0,650]/AFB(m)","1:BRleg");
  mmee.SavePlot("[em][em]201[6-8]/m[60,120]/y[-2.4,2.4]/pt[0,650]/AFB(y)","1:TMleg");
  mmee.SavePlot("[em][em]201[6-8]/m[60,120]/y[-2.4,2.4]/pt[4,650]/AFB(pt)","logx 1:BLleg rebin:2");

  mmee.SavePlot("[em][em]201[6-8]/m[120,3000]/y[-2.4,2.4]/pt[0,650]/dimass","norm logx widthweight 2:widey");
  mmee.SavePlot("[em][em]201[6-8]/m[120,3000]/y[-2.4,2.4]/pt[4,650]/dipt","norm logx widthweight");
  mmee.SavePlot("[em][em]201[6-8]/m[120,3000]/y[-2.4,2.4]/pt[0,650]/dirap","norm 1:BMleg");
  mmee.SavePlot("[em][em]201[6-8]/m[120,3000]/y[-2.4,2.4]/pt[0,650]/AFB(m)","1:TLleg logx");
  mmee.SavePlot("[em][em]201[6-8]/m[120,3000]/y[-2.4,2.4]/pt[0,650]/AFB(y)","1:TMleg");
  mmee.SavePlot("[em][em]201[6-8]/m[120,3000]/y[-2.4,2.4]/pt[4,650]/AFB(pt)","logx 1:BLleg rebin:2");
}

void compare(){
  AFBPlotter dm("data *amc+bg");
  dm.SavePlot("ee201[6-8]/m[60,120]/y[-2.4,2.4]/pt[0,650]/dimass","norm widthweight");
  dm.SavePlot("ee201[6-8]/m[60,120]/y[-2.4,2.4]/pt[4,650]/dipt","norm logx widthweight");
  dm.SavePlot("ee201[6-8]/m[60,120]/y[-2.4,2.4]/pt[0,650]/dirap","norm 1:BMleg");
  dm.SavePlot("ee201[6-8]/m[60,120]/y[-2.4,2.4]/pt[0,650]/AFB(m)","1:BRleg");
  dm.SavePlot("ee201[6-8]/m[60,120]/y[-2.4,2.4]/pt[0,650]/AFB(y)","1:TMleg");
  dm.SavePlot("ee201[6-8]/m[60,120]/y[-2.4,2.4]/pt[4,650]/AFB(pt)","logx 1:BLleg rebin:2");

  dm.SavePlot("mm201[6-8]/m[60,120]/y[-2.4,2.4]/pt[0,650]/dimass","norm widthweight");
  dm.SavePlot("mm201[6-8]/m[60,120]/y[-2.4,2.4]/pt[4,650]/dipt","norm logx widthweight");
  dm.SavePlot("mm201[6-8]/m[60,120]/y[-2.4,2.4]/pt[0,650]/dirap","norm 1:BMleg");
  dm.SavePlot("mm201[6-8]/m[60,120]/y[-2.4,2.4]/pt[0,650]/AFB(m)","1:BRleg");
  dm.SavePlot("mm201[6-8]/m[60,120]/y[-2.4,2.4]/pt[0,650]/AFB(y)","1:TMleg");
  dm.SavePlot("mm201[6-8]/m[60,120]/y[-2.4,2.4]/pt[4,650]/AFB(pt)","logx 1:BLleg rebin:2");

  dm.Setup("data *amcM+bg");
  dm.SavePlot("ee201[6-8]/m[120,3000]/y[-2.4,2.4]/pt[0,650]/dimass","norm logx widthweight");
  dm.SavePlot("ee201[6-8]/m[120,3000]/y[-2.4,2.4]/pt[4,650]/dipt","norm logx widthweight");
  dm.SavePlot("ee201[6-8]/m[120,3000]/y[-2.4,2.4]/pt[0,650]/dirap","norm 1:BMleg");
  dm.SavePlot("ee201[6-8]/m[120,3000]/y[-2.4,2.4]/pt[0,650]/AFB(m)","1:TLleg logx");
  dm.SavePlot("ee201[6-8]/m[120,3000]/y[-2.4,2.4]/pt[0,650]/AFB(y)","1:TMleg");
  dm.SavePlot("ee201[6-8]/m[120,3000]/y[-2.4,2.4]/pt[4,650]/AFB(pt)","logx 1:BLleg rebin:2");

  dm.SavePlot("mm201[6-8]/m[120,3000]/y[-2.4,2.4]/pt[0,650]/dimass","norm logx widthweight");
  dm.SavePlot("mm201[6-8]/m[120,3000]/y[-2.4,2.4]/pt[4,650]/dipt","norm logx widthweight");
  dm.SavePlot("mm201[6-8]/m[120,3000]/y[-2.4,2.4]/pt[0,650]/dirap","norm 1:BMleg");
  dm.SavePlot("mm201[6-8]/m[120,3000]/y[-2.4,2.4]/pt[0,650]/AFB(m)","1:TLleg logx");
  dm.SavePlot("mm201[6-8]/m[120,3000]/y[-2.4,2.4]/pt[0,650]/AFB(y)","1:TMleg");
  dm.SavePlot("mm201[6-8]/m[120,3000]/y[-2.4,2.4]/pt[4,650]/AFB(pt)","logx 1:BLleg rebin:2"); 
}
