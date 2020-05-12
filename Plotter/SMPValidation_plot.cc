#inlcude"plot.cc"
////////////////////////////// Setup function for SMPValidation/////////////////////////////
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
  TString analyzer="SMPValidation";
  if(skim!="") skim="_"+skim;
  cout<<"[SetupSamples] "<<analyzer<<" "<<GetStringChannel((Channel)channel)<<year<<endl;
  TString filedir=TString(getenv("SKFlatOutputDir"))+getenv("SKFlatV")+"/"+analyzer+"/";
  if(year==2017){
    if(channel==Channel::MUON){
      AddSample("data",SampleType::DATA,EColor::kBlack,filedir+"2017/DATA/"+analyzer+skim+"_DoubleMuon_B.root",filedir+"2017/DATA/"+analyzer+skim+"_DoubleMuon_C.root",filedir+"2017/DATA/"+analyzer+skim+"_DoubleMuon_D.root",filedir+"2017/DATA/"+analyzer+skim+"_DoubleMuon_E.root",filedir+"2017/DATA/"+analyzer+skim+"_DoubleMuon_F.root");
      AddSample("#gamma*/Z#rightarrow#mu#mu",SampleType::SIGNAL,EColor::kRed,filedir+"2017/"+analyzer+skim+"_DYJets.root");
    }else if(channel==Channel::ELECTRON){
      AddSample("data",SampleType::DATA,EColor::kBlack,filedir+"2017/DATA/"+analyzer+skim+"_DoubleEG_B.root",filedir+"2017/DATA/"+analyzer+skim+"_DoubleEG_C.root",filedir+"2017/DATA/"+analyzer+skim+"_DoubleEG_D.root",filedir+"2017/DATA/"+analyzer+skim+"_DoubleEG_E.root",filedir+"2017/DATA/"+analyzer+skim+"_DoubleEG_F.root");
      AddSample("#gamma*/Z#rightarrowee",SampleType::SIGNAL,EColor::kRed,filedir+"2017/"+analyzer+skim+"_DYJets.root");
    }else return "";
  }else if(year==2016){
    if(channel==Channel::MUON){
      AddSample("data",SampleType::DATA,EColor::kBlack,filedir+"2016/DATA/"+analyzer+skim+"_DoubleMuon_B_ver2.root",filedir+"2016/DATA/"+analyzer+skim+"_DoubleMuon_C.root",filedir+"2016/DATA/"+analyzer+skim+"_DoubleMuon_D.root",filedir+"2016/DATA/"+analyzer+skim+"_DoubleMuon_E.root",filedir+"2016/DATA/"+analyzer+skim+"_DoubleMuon_F.root",filedir+"2016/DATA/"+analyzer+skim+"_DoubleMuon_G.root",filedir+"2016/DATA/"+analyzer+skim+"_DoubleMuon_H.root");
      AddSample("#gamma*/Z#rightarrow#mu#mu",SampleType::SIGNAL,(EColor)(EColor::kRed),filedir+"2016/"+analyzer+skim+"_DYJets.root");
    }else if(channel==Channel::ELECTRON){
      AddSample("data",SampleType::DATA,EColor::kBlack,filedir+"2016/DATA/"+analyzer+skim+"_DoubleEG_B_ver2.root",filedir+"2016/DATA/"+analyzer+skim+"_DoubleEG_C.root",filedir+"2016/DATA/"+analyzer+skim+"_DoubleEG_D.root",filedir+"2016/DATA/"+analyzer+skim+"_DoubleEG_E.root",filedir+"2016/DATA/"+analyzer+skim+"_DoubleEG_F.root",filedir+"2016/DATA/"+analyzer+skim+"_DoubleEG_G.root",filedir+"2016/DATA/"+analyzer+skim+"_DoubleEG_H.root");
      AddSample("#gamma*/Z#rightarrowee",SampleType::SIGNAL,EColor::kRed,filedir+"2016/"+analyzer+skim+"_DYJets.root");
    }else return "";
  }else return "";
  AddSample("#gamma*/Z#rightarrow#tau#tau",SampleType::BG,EColor::kGreen,filedir+syear+"/"+analyzer+skim+"_DYJets.root");
  samples.back().prefix="tau_";
  AddSample("Diboson",SampleType::BG,EColor::kBlue,filedir+syear+"/"+analyzer+skim+"_WW_pythia.root",filedir+syear+"/"+analyzer+skim+"_WZ_pythia.root",filedir+syear+"/"+analyzer+skim+"_ZZ_pythia.root");
  AddSample("W",SampleType::BG,EColor::kYellow,filedir+syear+"/"+analyzer+skim+"_WJets_MG.root");
  AddSample("t#bar{t}",SampleType::BG,EColor::kMagenta,year==2017?filedir+"2017/"+analyzer+skim+"_TTLL_powheg.root":filedir+"2016/"+analyzer+skim+"_TT_powheg.root");
}
void SetupSystematics(int channel,int year){
  cout<<"[SetupSystematics]"<<endl;
  AddSystematic("RECOSF",SystematicType::ENVELOPE,"_RECOSF_up _RECOSF_down",false,true,true);
  AddSystematic("IDSF",SystematicType::ENVELOPE,"_IDSF_up _IDSF_down",false,true,true);
  AddSystematic("ISOSF",SystematicType::ENVELOPE,"_ISOSF_up _ISOSF_down",false,true,true);
  AddSystematic("triggerSF",SystematicType::ENVELOPE,"_triggerSF_up _triggerSF_down",false,true,true);
  AddSystematic("PUreweight",SystematicType::ENVELOPE,"_PUreweight_up _PUreweight_down",false,true,true);
  AddSystematic("prefireweight",SystematicType::ENVELOPE,"_prefireweight_up _prefireweight_down",false,true,true);
  AddSystematic("alphaS",SystematicType::ENVELOPE,"_alphaS_up _alphaS_down",false,true,false);
  AddSystematic("scalevariation",SystematicType::ENVELOPE,"_scalevariation0 _scalevariation1 _scalevariation2 _scalevariation3 _scalevariation5 _scalevariation6 _scalevariation7",false,true,false);

  vector<TString> prefixes;
  for(int i=0;i<100;i++) prefixes.push_back(Form("_pdf%d",i));
  if(year==2017) AddSystematic("pdf",SystematicType::HESSIAN,prefixes,false,true,false);
  else if(year==2016) AddSystematic("pdf",SystematicType::GAUSSIAN,prefixes,false,true,false);
  else cout<<"###WARNING### [SetupSystematics] wrong year"<<endl;

  AddSystematic("noRECOSF",SystematicType::ENVELOPE,"_noRECOSF",false,true,true);
  AddSystematic("noIDSF",SystematicType::ENVELOPE,"_noIDSF",false,true,true);
  AddSystematic("noISOSF",SystematicType::ENVELOPE,"_noISOSF",false,true,true);
  AddSystematic("notriggerSF",SystematicType::ENVELOPE,"_notriggerSF",false,true,true);
  AddSystematic("noPUreweight",SystematicType::ENVELOPE,"_noPUreweight",false,true,true);
  AddSystematic("noprefireweight",SystematicType::ENVELOPE,"_noprefireweight",false,true,true);
  AddSystematic("nozptcor",SystematicType::ENVELOPE,"_nozptcor",false,true,false);
  AddSystematic("noefficiencySF",SystematicType::ENVELOPE,"_noefficiencySF",false,true,true);
  AddSystematic("IDSF_POG",SystematicType::ENVELOPE,"_IDSF_POG",false,true,true);
  AddSystematic("efficiencySF",SystematicType::MULTI,"RECOSF IDSF ISOSF triggerSF",false,false,false);
  AddSystematic("totalsys",SystematicType::MULTI,"RECOSF IDSF ISOSF triggerSF PUreweight prefireweight alphaS scalevariation pdf nozptcor",false,false,false);
}
void SetupDirectories(int channel,int year){
  cout<<"[SetupDirectoriesSMPValidation] SMPValidation "<<GetStringChannel((Channel)channel)<<year<<endl;
  TString region[]={"OS","OS_Z","OS_Z_y0.0to0.4","OS_Z_y0.4to0.8","OS_Z_y0.8to1.2","OS_Z_y1.2to1.6","OS_Z_y1.6to2.0","OS_Z_y2.0to2.4","SS"};
  for(int i=0;i<sizeof(region)/sizeof(TString);i++){
    //for(int i=0;i<1;i++){
    Directory directory;
    directory.name=GetStringChannel((Channel)channel)+Form("%d",year)+"/"+region[i]+"/";
    cout<<directory.name<<" ";
    AddPlot(directory,"dimass",0,80,100);
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
    directories.push_back(directory);
  }
  cout<<endl;
  directories[0].plots[0].rebin=4;
  directories[0].plots[0].xmin=0;
  directories[0].plots[0].xmax=0;
}
