#include"plot.cc"
int MODE=0;
bool SSELECTRON=false;
////////////////////////////// Setup function for AFBAnalyzer/////////////////////////////
void SetupSamples(int channel,int year,TString skim);
void SetupSystematics(int channel,int year);
void SetupDirectories(int channel,int year);
TString Setup(int channel,int year,TString skim="SkimTree_SMP"){
  samples.clear();
  systematics.clear();
  directories.clear();

  SetupSamples(channel,year,skim);
  SetupSystematics(channel,year);
  SetupDirectories(channel,year);

  cout<<"[Setup] nsample: "<<samples.size()<<endl;
  cout<<"[Setup] nsys: "<<systematics.size()<<endl;
  cout<<"[Setup] ndirectories: "<<directories.size()<<endl;

  TString schannel=GetStringChannel((Channel)channel);
  TString syear=Form("%d",year);
  return schannel+syear;
}

void SetupSamples(int channel,int year,TString skim){
  TString syear=Form("%d",year);
  TString analyzer="AFBAnalyzer";
  if(skim!="") skim="_"+skim;
  cout<<"[SetupSamples] "<<analyzer<<" "<<GetStringChannel((Channel)channel)<<year<<endl;
  TString filedir=TString(getenv("SKFlatOutputDir"))+getenv("SKFlatV")+"/"+analyzer+"/";
  if(year==2017){
    if(channel==Channel::MUON){
      AddSample("data",SampleType::DATA,EColor::kBlack,filedir+"2017/DATA/"+analyzer+skim+"_DoubleMuon_B.root",filedir+"2017/DATA/"+analyzer+skim+"_DoubleMuon_C.root",filedir+"2017/DATA/"+analyzer+skim+"_DoubleMuon_D.root",filedir+"2017/DATA/"+analyzer+skim+"_DoubleMuon_E.root",filedir+"2017/DATA/"+analyzer+skim+"_DoubleMuon_F.root");
    }else if(channel==Channel::ELECTRON){
      AddSample("data",SampleType::DATA,EColor::kBlack,filedir+"2017/DATA/"+analyzer+skim+"_DoubleEG_B.root",filedir+"2017/DATA/"+analyzer+skim+"_DoubleEG_C.root",filedir+"2017/DATA/"+analyzer+skim+"_DoubleEG_D.root",filedir+"2017/DATA/"+analyzer+skim+"_DoubleEG_E.root",filedir+"2017/DATA/"+analyzer+skim+"_DoubleEG_F.root");
    }else return "";
  }else if(year==2016){
    if(channel==Channel::MUON){
      AddSample("data",SampleType::DATA,EColor::kBlack,filedir+"2016/DATA/"+analyzer+skim+"_DoubleMuon_B_ver2.root",filedir+"2016/DATA/"+analyzer+skim+"_DoubleMuon_C.root",filedir+"2016/DATA/"+analyzer+skim+"_DoubleMuon_D.root",filedir+"2016/DATA/"+analyzer+skim+"_DoubleMuon_E.root",filedir+"2016/DATA/"+analyzer+skim+"_DoubleMuon_F.root",filedir+"2016/DATA/"+analyzer+skim+"_DoubleMuon_G.root",filedir+"2016/DATA/"+analyzer+skim+"_DoubleMuon_H.root");
    }else if(channel==Channel::ELECTRON){
      AddSample("data",SampleType::DATA,EColor::kBlack,filedir+"2016/DATA/"+analyzer+skim+"_DoubleEG_B_ver2.root",filedir+"2016/DATA/"+analyzer+skim+"_DoubleEG_C.root",filedir+"2016/DATA/"+analyzer+skim+"_DoubleEG_D.root",filedir+"2016/DATA/"+analyzer+skim+"_DoubleEG_E.root",filedir+"2016/DATA/"+analyzer+skim+"_DoubleEG_F.root",filedir+"2016/DATA/"+analyzer+skim+"_DoubleEG_G.root",filedir+"2016/DATA/"+analyzer+skim+"_DoubleEG_H.root");
    }else return "";
  }else return "";
  TString dytitle=(channel==Channel::MUON)?"#gamma*/Z#rightarrow#mu#mu":"#gamma*/Z#rightarrowee";
  if(MODE==0) AddSample(dytitle,SampleType::SIGNAL,(EColor)(EColor::kRed),filedir+syear+"/"+analyzer+skim+"_DYJets.root");
  else if(MODE==1){
    AddSampleWithPrefixAndWeight(dytitle+"(qqbar->DY+0parton)",SampleType::BG,EColor::kRed,filedir+syear+"/"+analyzer+skim+"_DYJets.root","qqbar_0j_",1,filedir+syear+"/"+analyzer+skim+"_DYJets.root","qbarq_0j_",1);
    AddSampleWithPrefixAndWeight(dytitle+"(qqbar->DY+Npartons)",SampleType::SIGNAL,EColor::kOrange,filedir+syear+"/"+analyzer+skim+"_DYJets.root","qqbar_nj_",1,filedir+syear+"/"+analyzer+skim+"_DYJets.root","qbarq_nj_",1);
    AddSampleWithPrefixAndWeight(dytitle+"(gq)",SampleType::SIGNAL,EColor::kGreen,filedir+syear+"/"+analyzer+skim+"_DYJets.root","gq_",1,filedir+syear+"/"+analyzer+skim+"_DYJets.root","qg_",1);
    AddSampleWithPrefixAndWeight(dytitle+"(gqbar)",SampleType::SIGNAL,EColor::kBlue,filedir+syear+"/"+analyzer+skim+"_DYJets.root","gqbar_",1,filedir+syear+"/"+analyzer+skim+"_DYJets.root","qbarg_",1);
    AddSampleWithPrefixAndWeight(dytitle+"(etc)",SampleType::BG,EColor::kMagenta,filedir+syear+"/"+analyzer+skim+"_DYJets.root","qq_",1,filedir+syear+"/"+analyzer+skim+"_DYJets.root","qbarqbar_",1,filedir+syear+"/"+analyzer+skim+"_DYJets.root","gg_",1);
  }else if(MODE==2){
    AddSampleWithPrefixAndWeight(dytitle+"(qg)",SampleType::SIGNAL,EColor::kOrange,filedir+syear+"/"+analyzer+skim+"_DYJets.root","qg_",1);
    AddSampleWithPrefixAndWeight(dytitle+"(gq)",SampleType::SIGNAL,EColor::kRed,filedir+syear+"/"+analyzer+skim+"_DYJets.root","gq_",1);
    AddSampleWithPrefixAndWeight(dytitle+"(gqbar)",SampleType::SIGNAL,EColor::kBlue,filedir+syear+"/"+analyzer+skim+"_DYJets.root","gqbar_",1);
    AddSampleWithPrefixAndWeight(dytitle+"(qbarg)",SampleType::SIGNAL,EColor::kYellow,filedir+syear+"/"+analyzer+skim+"_DYJets.root","qbarg_",1);
    AddSampleWithPrefixAndWeight(dytitle+"(qqbar->DY+Npartons)",SampleType::SIGNAL,EColor::kGreen,filedir+syear+"/"+analyzer+skim+"_DYJets.root","qqbar_nj_",1);
    AddSampleWithPrefixAndWeight(dytitle+"(qbarq->DY+Npartons)",SampleType::SIGNAL,EColor::kRed,filedir+syear+"/"+analyzer+skim+"_DYJets.root","qbarq_nj_",1);
    AddSampleWithPrefixAndWeight(dytitle+"(qqbar->DY+0parton)",SampleType::BG,EColor::kMagenta,filedir+syear+"/"+analyzer+skim+"_DYJets.root","qqbar_0j_",1);
    AddSampleWithPrefixAndWeight(dytitle+"(qbarq->DY+0parton)",SampleType::BG,EColor::kRed,filedir+syear+"/"+analyzer+skim+"_DYJets.root","qbarq_0j_",1);
    AddSampleWithPrefixAndWeight(dytitle+"(etc)",SampleType::BG,EColor::kMagenta,filedir+syear+"/"+analyzer+skim+"_DYJets.root","qq_",1,filedir+syear+"/"+analyzer+skim+"_DYJets.root","qbarqbar_",1,filedir+syear+"/"+analyzer+skim+"_DYJets.root","gg_",1);
  }
  
  AddSampleWithPrefixAndWeight("#gamma*/Z#rightarrow#tau#tau",SampleType::BG,EColor::kGreen,filedir+syear+"/"+analyzer+skim+"_DYJets.root","tau_",1);

  if(channel==Channel::MUON||SSELECTRON){
    if(year==2017){
      AddSampleWithPrefixAndWeight("QCD dijet",SampleType::BG,EColor::kCyan,samples[0].files[0],"ss_",1,samples[0].files[0],"ss_",1,samples[0].files[1],"ss_",1,samples[0].files[2],"ss_",1,samples[0].files[3],"ss_",1,samples[0].files[4],"ss_",1,filedir+syear+"/"+analyzer+skim+"_DYJets.root","ss_",-1,filedir+syear+"/"+analyzer+skim+"_DYJets.root","ss_tau_",-1,filedir+syear+"/"+analyzer+skim+"_WW_pythia.root","ss_",-1,filedir+syear+"/"+analyzer+skim+"_WZ_pythia.root","ss_",-1,filedir+syear+"/"+analyzer+skim+"_ZZ_pythia.root","ss_",-1,filedir+syear+"/"+analyzer+skim+"_WJets_MG.root","ss_",-1,year==2017?filedir+"2017/"+analyzer+skim+"_TTLL_powheg.root":filedir+"2016/"+analyzer+skim+"_TT_powheg.root","ss_",-1);
    }else{
      AddSampleWithPrefixAndWeight("QCD dijet",SampleType::BG,EColor::kCyan,samples[0].files[0],"ss_",1,samples[0].files[0],"ss_",1,samples[0].files[1],"ss_",1,samples[0].files[2],"ss_",1,samples[0].files[3],"ss_",1,samples[0].files[4],"ss_",1,samples[0].files[5],"ss_",1,samples[0].files[6],"ss_",1,filedir+syear+"/"+analyzer+skim+"_DYJets.root","ss_",-1,filedir+syear+"/"+analyzer+skim+"_DYJets.root","ss_tau_",-1,filedir+syear+"/"+analyzer+skim+"_WW_pythia.root","ss_",-1,filedir+syear+"/"+analyzer+skim+"_WZ_pythia.root","ss_",-1,filedir+syear+"/"+analyzer+skim+"_ZZ_pythia.root","ss_",-1,filedir+syear+"/"+analyzer+skim+"_WJets_MG.root","ss_",-1,year==2017?filedir+"2017/"+analyzer+skim+"_TTLL_powheg.root":filedir+"2016/"+analyzer+skim+"_TT_powheg.root","ss_",-1);
    }
  }

  AddSample("Diboson",SampleType::BG,EColor::kBlue,filedir+syear+"/"+analyzer+skim+"_WW_pythia.root",filedir+syear+"/"+analyzer+skim+"_WZ_pythia.root",filedir+syear+"/"+analyzer+skim+"_ZZ_pythia.root");
  AddSample("W",SampleType::BG,EColor::kYellow,filedir+syear+"/"+analyzer+skim+"_WJets_MG.root");
  AddSample("t#bar{t}",SampleType::BG,EColor::kMagenta,year==2017?filedir+"2017/"+analyzer+skim+"_TTLL_powheg.root":filedir+"2016/"+analyzer+skim+"_TT_powheg.root");
}
void SetupSystematics(int channel,int year){
  cout<<"[SetupSystematics]"<<endl;
  if(channel==Channel::ELECTRON) AddSystematic("RECOSF",SystematicType::ENVELOPE,"_RECOSF_up _RECOSF_down",false,true,true);
  AddSystematic("IDSF",SystematicType::ENVELOPE,"_IDSF_up _IDSF_down",false,true,true);
  if(channel==Channel::MUON) AddSystematic("ISOSF",SystematicType::ENVELOPE,"_ISOSF_up _ISOSF_down",false,true,true);
  AddSystematic("triggerSF",SystematicType::ENVELOPE,"_triggerSF_up _triggerSF_down",false,true,true);
  AddSystematic("PUreweight",SystematicType::ENVELOPE,"_PUreweight_up _PUreweight_down",false,true,true);
  AddSystematic("prefireweight",SystematicType::ENVELOPE,"_prefireweight_up _prefireweight_down",false,true,true);
  AddSystematic("scale",SystematicType::ENVELOPE,"_scale_up _scale_down",false,true,true);
  if(channel==Channel::ELECTRON) AddSystematic("smear",SystematicType::ENVELOPE,"_smear_up _smear_down",false,true,true);
  AddSystematic("alphaS",SystematicType::ENVELOPE,"_alphaS_up _alphaS_down",false,true,false);
  AddSystematic("scalevariation",SystematicType::ENVELOPE,"_scalevariation0 _scalevariation1 _scalevariation2 _scalevariation3 _scalevariation5 _scalevariation6 _scalevariation7",false,true,false);

  vector<TString> prefixes;
  for(int i=0;i<100;i++) prefixes.push_back(Form("_pdf%d",i));
  if(year==2017) AddSystematic("pdf",SystematicType::HESSIAN,prefixes,false,true,false);
  else if(year==2016) AddSystematic("pdf",SystematicType::GAUSSIAN,prefixes,false,true,false);
  else cout<<"###WARNING### [SetupSystematics] wrong year"<<endl;

  if(channel==Channel::ELECTRON) AddSystematic("noRECOSF",SystematicType::ENVELOPE,"_noRECOSF",false,true,true);
  AddSystematic("noIDSF",SystematicType::ENVELOPE,"_noIDSF",false,true,true);
  if(channel==Channel::MUON) AddSystematic("noISOSF",SystematicType::ENVELOPE,"_noISOSF",false,true,true);
  AddSystematic("notriggerSF",SystematicType::ENVELOPE,"_notriggerSF",false,true,true);
  AddSystematic("noPUreweight",SystematicType::ENVELOPE,"_noPUreweight",false,true,true);
  AddSystematic("noprefireweight",SystematicType::ENVELOPE,"_noprefireweight",false,true,true);
  AddSystematic("nozptcor",SystematicType::ENVELOPE,"_nozptcor",false,true,false);
  AddSystematic("noefficiencySF",SystematicType::ENVELOPE,"_noefficiencySF",false,true,true);
  if(channel==Channel::ELECTRON) AddSystematic("IDSF_POG",SystematicType::ENVELOPE,"_IDSF_POG",false,true,true);
  if(channel==Channel::ELECTRON) AddSystematic("selective",SystematicType::ENVELOPE,"_selective",true,true,true);
  AddSystematic("efficiencySF",SystematicType::MULTI,"RECOSF IDSF ISOSF triggerSF",false,false,false);
  AddSystematic("totalsys",SystematicType::MULTI,"RECOSF IDSF ISOSF triggerSF PUreweight prefireweight scale smear alphaS scalevariation pdf nozptcor",false,false,false);
}
void SetupDirectories(int channel,int year){
  cout<<"[SetupDirectoriesAFBAnalyzer] AFBAnalyzer "<<GetStringChannel((Channel)channel)<<year<<endl;
  const int ybinnum=6;
  double yrange[ybinnum+1]={0.0,0.4,0.8,1.2,1.6,2.0,2.4};
  const int mbinnum=12;
  double massrange[mbinnum+1]={60,70,78,84,87,89,91,93,95,98,104,112,120};
  vector<TString> region;
  region.push_back(Form("m%.0fto%.0f",massrange[0],massrange[mbinnum]));
  region.push_back(Form("m%.0fto%.0f_pt50",massrange[0],massrange[mbinnum]));
  for(int iy=0;iy<ybinnum;iy++){
    region.push_back(Form("m%.0fto%.0f_y%.1fto%.1f",massrange[0],massrange[mbinnum],yrange[iy],yrange[iy+1]));
    region.push_back(Form("m%.0fto%.0f_y%.1fto%.1f_pt50",massrange[0],massrange[mbinnum],yrange[iy],yrange[iy+1]));
    for(int im=0;im<mbinnum;im++){
      //region.push_back(Form("OS_m%.0fto%.0f_y%.1fto%.1f",massrange[im],massrange[im+1],yrange[iy],yrange[iy+1]));
    }
  }
  for(int i=0;i<region.size();i++){
    Directory directory;
    directory.name=GetStringChannel((Channel)channel)+Form("%d",year)+"/"+region[i]+"/";
    cout<<directory.name<<" ";
    AddPlot(directory,"dimass",2,60,120);
    AddPlot(directory,"dipt",4,0,200);
    AddPlot(directory,"dirap",0,0,0);
    AddPlot(directory,"l0pt",2,0,100);
    AddPlot(directory,"l0eta",0,0,0);
    AddPlot(directory,"l0riso",0,0,0);
    AddPlot(directory,"l1pt",2,0,100);
    AddPlot(directory,"l1eta",0,0,0);
    AddPlot(directory,"l1riso",0,0,0);
    AddPlot(directory,"lldelR",0,0,0);
    AddPlot(directory,"lldelphi",0,0,0);
    AddPlot(directory,"nPV",0,0,0,"widey");
    AddPlot(directory,"costhetaCS",0,0,0,"1:noleg");
    AddPlot(directory,"abscosthetaCS",0,0,0,"1:noleg");
    if(region[i].Contains("m60to120")) AddPlot(directory,"AFB",0,0,0,"AFB");
    directories.push_back(directory);
  }
  cout<<endl;
}
void SaveAFB(TString outputdir="AFBAnalyzer_plot"){
  int oldlevel=gErrorIgnoreLevel;
  gErrorIgnoreLevel=kWarning;
  int smax=systematics.size();
  TCanvas* c=NULL;
  TString AFB[]={"AFB","AFB_B","weightedAFB","weightedAFB_B"};
  TString dirname=directories[0].name(0,directories[0].name.Index('/'))+"/AFB/";
  std::cout<<"mkdir -p "+outputdir+"/"+dirname<<endl;
  system("mkdir -p "+outputdir+"/"+dirname);
  for(int j=0;j<sizeof(AFB)/sizeof(TString);j++){
    c=GetCompareAFBAll("_y...to.../",AFB[j],0,AFB[j]);
    c->SaveAs(outputdir+"/"+dirname+AFB[j]+".png");
    delete c;
  }
  dirname=directories[0].name(0,directories[0].name.Index('/'))+"/summary_pt50/";
  std::cout<<"mkdir -p "+outputdir+"/"+dirname<<endl;
  system("mkdir -p "+outputdir+"/"+dirname);
  for(int j=0;j<sizeof(AFB)/sizeof(TString);j++){
    c=GetCompareAFBAll("_y...to..._pt50/",AFB[j],0,AFB[j]);
    c->SaveAs(outputdir+"/"+dirname+AFB[j]+".png");
    delete c;
  }
  gErrorIgnoreLevel=oldlevel;
}
