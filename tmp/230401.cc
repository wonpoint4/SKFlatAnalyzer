#include"EfficiencyPlotter.cc"

void plot(){
  EfficiencyPlotter aa;
  aa.plotdir="fig/230401_v15";
  for(TString channel:{"ee","el","mm","mu"}){
    //for(TString channel:{"mm","mu"}){
    for(TString era:{"2016a","2016b","2017","2018","Run2"}){
      //for(TString era:{"2017"}){
      TString era_=era=="Run2"?"201[678][ab]?":era;
      TString lepton=channel[0]=='e'?"electron":"muon";
      TString common=" 1:ytitle:Events 2:ytitle:'data/Pred.' 2:xtitle:'"+lepton+" #eta' norm 2:ymin:0.92 2:ymax:1.08 rebin:2 chi2 ";
      if(channel.BeginsWith('e')) common+=" xmin:-2.5 xmax:2.5 ";
      else common+=" xmin:-2.4 xmax:2.4 ";
      aa.SavePlotCondor(channel+era_+"/m80to100/leta",common+"sysname:effAN sysdetail 2:sysleg save:"+channel+era+"_leta");
      aa.SavePlotCondor("v14/"+channel+era_+"/m80to100/leta",common+"sysname:effAN sysdetail 2:sysleg save:"+channel+era+"_leta_v14");
      if(channel[0]=='m'){
	//aa.SavePlotCondor(channel+era_+"/m80to100/leta",common+"suffix:_nomuonTrackingSF:sim save:"+channel+era+"_leta_nomuonTrackingSF");
	//aa.SavePlotCondor(channel+era_+"/m80to100/leta",common+"suffix:_nomuonRECOSF:sim save:"+channel+era+"_leta_nomuonRECOSF");
	//aa.SavePlotCondor("v11/"+channel+era_+"/m80to100/leta",common+"sysname:effAN sysdetail  2:sysleg 2:ymin:0.92 2:ymax:1.08 rebin:2 save:"+channel+era+"_leta_v11");
	//aa.SavePlotCondor("v11/"+channel+era_+"/m80to100/leta",common+"suffix:_noDZSF:sim save:"+channel+era+"_leta_noDZSF_v11");
	//aa.SavePlotCondor("v13/"+channel+era_+"/m80to100/leta",common+"sysname:effAN sysdetail 2:sysleg save:"+channel+era+"_leta_v13");
	//aa.SavePlotCondor("v13/"+channel+era_+"/m80to100/leta",common+"suffix:_noDZSF:sim save:"+channel+era+"_leta_noDZSF_v13");
	//aa.SavePlotCondor("v13/"+channel+era_+"/m80to100/leta",common+"suffix:_nomuonTrackingSF:sim save:"+channel+era+"_leta_nomuonTrackingSF_v13");
	//aa.SavePlotCondor("v13/"+channel+era_+"/m80to100/leta",common+"suffix:_nomuonRECOSF:sim save:"+channel+era+"_leta_nomuonRECOSF_v13");
      }else{
	//aa.SavePlotCondor(channel+era_+"/m80to100/lsceta",common+"sysname:effAN sysdetail 2:sysleg save:"+channel+era+"_lsceta 2:xtitle:'electron #eta_{SC}'");
	//aa.SavePlotCondor(channel+era_+"/m80to100/leta",common+"sysname:effAN sysdetail 2:sysleg save:"+channel+era+"_leta");
	//aa.SavePlotCondor("newreco/"+channel+era_+"/m80to100/leta",common+"sysname:effAN sysdetail 2:sysleg save:"+channel+era+"_leta_newreco");
	//aa.SavePlotCondor("newid/"+channel+era_+"/m80to100/leta",common+"sysname:effAN sysdetail 2:sysleg save:"+channel+era+"_leta_newid");
	//aa.SavePlotCondor("newrecoid/"+channel+era_+"/m80to100/leta",common+"sysname:effAN sysdetail 2:sysleg save:"+channel+era+"_leta_newrecoid");
      }
      common=" xmin:-2.4 xmax:2.4 1:ytitle:Events 2:ytitle:'data/Pred.' 2:xtitle:'y(ll)' norm 2:ymin:0.92 2:ymax:1.08 rebin:2 chi2 ";
      aa.SavePlotCondor(channel+era_+"/m80to100/dirap",common+"sysname:effAN sysdetail 2:sysleg save:"+channel+era+"_dirap");
      aa.SavePlotCondor("v14/"+channel+era_+"/m80to100/dirap",common+"sysname:effAN sysdetail 2:sysleg save:"+channel+era+"_dirap_v14");
    }
  }
}
