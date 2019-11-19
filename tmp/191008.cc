#include"Plotter/AFBPlotter.cc"
void SaveAll(TString mode="dab",TString flag=""){
  DEBUG=1;
  AFBPlotter p(mode,flag);
  p.RemovePlots(".*");
  vector<TString> syears={"2016/","2017/","2018/"};
  if(mode=="dajb") syears={"2016/","2018/"};
  vector<TString> schannels={"muon","electron"};
  if(mode=="m"||mode=="e"){
    syears={""};
    schannels={""};
  }
  for(auto syear:syears){
    for(auto schannel:schannels){      
      if(flag==""){
	p.SavePlot(schannel+syear+"m[80,100]/y[0.0,2.4]/pt[0,100]/dimass","widthweight norm");
	p.SavePlot(schannel+syear+"m[80,100]/y[0.0,2.4]/pt[0,100]/dipt","widthweight norm");
	p.SavePlot(schannel+syear+"m[80,100]/y[0.0,2.4]/pt[0,100]/dirap","widthweight norm");
	p.SavePlot(schannel+syear+"m[80,100]/y[0.0,2.4]/pt[0,100]/costhetaCS","widthweight norm");
	p.SavePlot(schannel+syear+"m[80,100]/y[0.0,2.4]/pt[0,100]/AFB(m)","type:5");
	p.SavePlot(schannel+syear+"m[80,100]/y[0.0,2.4]/pt[0,100]/AFB(y)","type:5");
	p.SavePlot(schannel+syear+"m[80,100]/y[0.0,2.4]/pt[0,100]/AFB(pt)","type:5");

	p.SavePlot(schannel+syear+"m[60,120]/y[0.0,2.4]/pt[0,100]/dimass","widthweight norm");
	p.SavePlot(schannel+syear+"m[60,120]/y[0.0,2.4]/pt[0,100]/dipt","widthweight norm");
	p.SavePlot(schannel+syear+"m[60,120]/y[0.0,2.4]/pt[0,100]/dirap","widthweight norm");
	p.SavePlot(schannel+syear+"m[60,120]/y[0.0,2.4]/pt[0,100]/costhetaCS","widthweight norm");
	p.SavePlot(schannel+syear+"m[60,120]/y[0.0,2.4]/pt[0,100]/AFB(m)","type:5");
	p.SavePlot(schannel+syear+"m[60,120]/y[0.0,2.4]/pt[0,100]/AFB(y)","type:5");
	p.SavePlot(schannel+syear+"m[60,120]/y[0.0,2.4]/pt[0,100]/AFB(pt)","type:5");
	p.SavePlot(schannel+syear+"m[60,120]/y[0.0,2.4]/pt[0,1000]/AFB(pt)","type:5");
	
	/*
	p.SavePlot(schannel+syear+"m[120,2000]/y[0.0,2.4]/pt[0,100]/dimass","widthweight norm");
	p.SavePlot(schannel+syear+"m[120,2000]/y[0.0,2.4]/pt[0,100]/dipt","widthweight norm");
	p.SavePlot(schannel+syear+"m[120,2000]/y[0.0,2.4]/pt[0,100]/dirap","widthweight norm");
	p.SavePlot(schannel+syear+"m[120,2000]/y[0.0,2.4]/pt[0,100]/costhetaCS","widthweight norm");
	p.SavePlot(schannel+syear+"m[120,2000]/y[0.0,2.4]/pt[0,100]/AFB(m)","type:5");
	p.SavePlot(schannel+syear+"m[120,2000]/y[0.0,2.4]/pt[0,100]/AFB(y)","type:5");
	p.SavePlot(schannel+syear+"m[120,2000]/y[0.0,2.4]/pt[0,100]/AFB(pt)","type:5");
	*/
	p.SavePlot(schannel+syear+"m[120,1000]/y[1.2,2.4]/pt[0,1000]/dimass","widthweight norm");
	p.SavePlot(schannel+syear+"m[120,1000]/y[1.2,2.4]/pt[0,1000]/dipt","widthweight norm");
	p.SavePlot(schannel+syear+"m[120,1000]/y[1.2,2.4]/pt[0,1000]/dirap","widthweight norm");
	p.SavePlot(schannel+syear+"m[120,1000]/y[1.2,2.4]/pt[0,1000]/costhetaCS","widthweight norm");
	p.SavePlot(schannel+syear+"m[120,1000]/y[1.2,2.4]/pt[0,1000]/AFB(m)","type:5");
	p.SavePlot(schannel+syear+"m[120,1000]/y[1.2,2.4]/pt[0,1000]/AFB(y)","type:5");
	p.SavePlot(schannel+syear+"m[120,1000]/y[1.2,2.4]/pt[0,1000]/AFB(pt)","type:5");
      }else if(flag=="LEP"){
	p.SavePlot(schannel+syear+"m[80,100]/y[0.0,2.4]/pt[0,100]/lpt","widthweight norm");
	p.SavePlot(schannel+syear+"m[80,100]/y[0.0,2.4]/pt[0,100]/leta","widthweight norm");
	
	p.SavePlot(schannel+syear+"m[60,120]/y[0.0,2.4]/pt[0,100]/lpt","widthweight norm");
	p.SavePlot(schannel+syear+"m[120,1000]/y[0.0,2.4]/pt[0,100]/lpt","widthweight norm");
	p.SavePlot(schannel+syear+"m[120,1000]/y[1.2,2.4]/pt[0,1000]/lpt","widthweight norm");
	
      }else if(flag=="GEN"){
	p.SavePlot(schannel+syear+"m[80,100]/y[0.0,2.4]/pt[0,100]/dimass","widthweight norm");	
      }else if(flag=="ETC"){
	p.SavePlot(schannel+syear+"m[80,100]/y[0.0,2.4]/pt[0,100]/z0","2:widewidey norm");	
      }
      if(mode=="dabs"){
	p.DrawPlot(schannel+syear+"m[60,800]/y[0.0,2.4]/pt[0,1000]/dimass","norm widthweight 1:logy");
      }
    }
  }
}      
  
