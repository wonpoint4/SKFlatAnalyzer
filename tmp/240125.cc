#include"Chi2Plotter.cc"

void save_compare(TString suffix=""){
  Chi2Plotter aa;
  aa.plotdir="fig/240125/";
  for(TString channel:{"el","mu","ee","mm"}){
    for(TString era:{"2016a","2016b","2017","2018","Run2"}){
      TString era_=era;
      if(era=="Run2") era_="201[678][ab]?";
      TString common=" chi2 pvalue movesys ";
      if(suffix=="_residual"){
	if(channel[0]=='e'){
	  common+=" norm suffix:_electronIDSF_s18m0:sim ";
	}else{
	  common+=" norm suffix:_muonIDSF_s17m0:sim ";
	}	  
      }else if(suffix=="_roccor_residual"){
	if(channel[0]=='e'){
	  common+=" norm suffix:_electronenergy_residual ";
	}else{
	  common+=" norm suffix:_muonmomentum_residual ";
	}	  
      }else{
	common+=" sysname:totalsys"+suffix+" ";
      }
      if(suffix.Contains("_table")){
	common+=" chi2detail ";
      }

      aa.SavePlotCondor(channel+era_+"/ym"," save:"+channel+era+"_ym"+suffix+".png "+common);
      //aa.SavePlotCondor(channel+era_+"/ym"," xmin:0 xmax:24 save:"+channel+era+"_ym_low"+suffix+".png "+common);
      //aa.SavePlotCondor(channel+era_+"/ym"," xmin:24 xmax:42 save:"+channel+era+"_ym_middle"+suffix+".png "+common);
      //aa.SavePlotCondor(channel+era_+"/ym"," xmin:42 xmax:66 save:"+channel+era+"_ym_high"+suffix+".png "+common);
      // aa.SavePlotCondor(channel+era_+"/ym"," project:x save:"+channel+era+"_y"+suffix+".png "+common);
      // aa.SavePlotCondor(channel+era_+"/ym"," widthweight project:y save:"+channel+era+"_m"+suffix+".png "+common);
      // aa.SavePlotCondor(channel+era_+"/l0etapt"," save:"+channel+era+"_l0etapt"+suffix+".png "+common);
      // aa.SavePlotCondor(channel+era_+"/l1etapt"," save:"+channel+era+"_l1etapt"+suffix+".png "+common);
      // aa.SavePlotCondor(channel+era_+"/l0etapt"," project:x save:"+channel+era+"_l0eta"+suffix+".png "+common);
      // aa.SavePlotCondor(channel+era_+"/l0etapt"," logx widthweight project:y save:"+channel+era+"_l0pt"+suffix+".png "+common);
    }
  }
}

void save(){
  save_compare("");
  save_compare("_noresidual");
  save_compare("_roccor_noresidual");
  //save_compare("_noresidual_correlated");
  //save_compare("_correlated");
  save_compare("_residual");
  save_compare("_roccor_residual");
}
