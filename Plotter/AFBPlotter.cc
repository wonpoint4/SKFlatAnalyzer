#include"Plotter.cc"

class AFBPlotter:public Plotter{
public:
  void SetupSamples();
  void SetupSystematics();
  void SetupPlots();
  void SavePlotsCondor(int channel,int year,int njob);
  int Setup(int channel_,int year_,int mode_=0);
  int mode;
  int channel,year;
  TString schannel,syear;
  AFBPlotter();
};
AFBPlotter::AFBPlotter(){
  TString filedir=TString(getenv("SKFlatOutputDir"))+getenv("SKFlatV")+"/AFBAnalyzer/";
  TString channels[]={"DoubleMuon","DoubleEG"};
  vector<tuple<TString,TString>> periods={make_tuple("2016","B_ver2"),make_tuple("2016","C"),make_tuple("2016","D"),make_tuple("2016","E"),make_tuple("2016","F"),make_tuple("2016","G"),make_tuple("2016","H"),make_tuple("2017","B"),make_tuple("2017","C"),make_tuple("2017","D"),make_tuple("2017","E"),make_tuple("2017","F")};
  for(int ic=0;ic<sizeof(channel)/sizeof(TString);ic++){
    for(unsigned int ip=0;ip<periods.size();ip++){
      TString path=filedir+get<0>(periods[ip])+"/DATA/AFBAnalyzer_SkimTree_SMP_"+channels[ic]+"_";
      TString samplefragkey=channels[ic]+get<0>(periods[ip])+get<1>(periods[ip]);
      samplefrags[samplefragkey]=SampleFrag(samplefragkey,SampleFrag::Type::DATA,kBlack,make_tuple(path+get<1>(periods[ip])+".root","",1.));
    }
  }

  TString syears[]={"2016","2017"};
  for(int i=0;i<sizeof(syears)/sizeof(TString);i++){
    samplefrags["muon"+syears[i]]=SampleFrag("muon"+syears[i],SampleFrag::Type::DATA,kBlack);
    samplefrags["muon"+syears[i]].Add(TRegexp("DoubleMuon"+syears[i]+"?.*"),1);
    samplefrags["electron"+syears[i]]=SampleFrag("electron"+syears[i],SampleFrag::Type::DATA,kBlack);
    samplefrags["electron"+syears[i]].Add(TRegexp("DoubleElectron"+syears[i]+"?.*"),1);

    samplefrags["dy"+syears[i]]=SampleFrag("dy"+syears[i],SampleFrag::Type::SIGNAL,kRed,make_tuple(filedir+syears[i]+"/AFBAnalyzer_SkimTree_SMP_DYJets.root","",1.));
    
    samplefrags["dytt"+syears[i]]=SampleFrag("#gamma*/Z#rightarrow#tau#tau",SampleFrag::Type::BG,EColor::kGreen,make_tuple(filedir+syears[i]+"/AFBAnalyzer_SkimTree_SMP_DYJets.root","tau_",1.));
					       
    TString path=filedir+syear+"/AFBAnalyzer_SkimTree_SMP_";
    samplefrags["ww"+syears[i]]=SampleFrag("ww"+syears[i],SampleFrag::Type::BG,kBlue,make_tuple(path+"WW_pythia.root","",1.));
    samplefrags["wz"+syears[i]]=SampleFrag("wz"+syears[i],SampleFrag::Type::BG,kGreen,make_tuple(path+"WZ_pythia.root","",1.));
    samplefrags["zz"+syears[i]]=SampleFrag("zz"+syears[i],SampleFrag::Type::BG,kYellow,make_tuple(path+"ZZ_pythia.root","",1.));
    samplefrags["vv"+syears[i]]=SampleFrag("Diboson",SampleFrag::Type::BG,kBlue,make_tuple("ww"+syears[i],1.),make_tuple("wz"+syears[i],1.),make_tuple("zz"+syears[i],1.));

    samplefrags["wjets"+syears[i]]=SampleFrag("W",SampleFrag::Type::BG,kYellow,make_tuple(path+"WJets_pythia.root","",1.));
  }
  samplefrags["tt2016"]=SampleFrag("t#bar{t}",SampleFrag::Type::BG,kMagenta,make_tuple(filedir+"2016/AFBAnalyzer_SkimTree_SMP_TT_powheg.root","",1.));
  samplefrags["tt2017"]=SampleFrag("t#bar{t}",SampleFrag::Type::BG,kMagenta,make_tuple(filedir+"2017/AFBAnalyzer_SkimTree_SMP_TTLL_powheg.root","",1.));
}
int AFBPlotter::Setup(int channel_,int year_,int mode_){
  samples.clear();
  systematics.clear();
  plots.clear();
  channel=channel_;
  year=year_;
  mode=mode_;
  schannel=GetStringChannel((Channel)channel);
  syear=Form("%d",year);

  SetupSamples();
  SetupSystematics();
  //SetupPlots();

  cout<<"[Setup] nsample: "<<samples.size()<<endl;
  cout<<"[Setup] nsys: "<<systematics.size()<<endl;
  cout<<"[Setup] nplot: "<<plots.size()<<endl;

  return 1;
}
void AFBPlotter::SetupSamples(){
  cout<<"[AFBPlotter::SetupSamples]"<<endl;
  samples["data"]=Sample("data",Sample::Type::DATA,kBlack,make_tuple<schannel+syear,1.>);
  TString dytitle=(channel==Channel::MUON)?"#gamma*/Z#rightarrow#mu#mu":"#gamma*/Z#rightarrowee";
  samplefrags["dy"+syear].title=dytitle;
  samples["sim"]=Sample("sim",SampleFrag::Type::STACK,kBlue,make_tuple("dy"+syear,1.),make_tuple("dytt"+syear,1.),make_tuple("vv"+syear,1.),make_tuple("wjets"+syear,1.),make_tuple("tt"+syear,1.));

  /*
  if(MODE==0) 
    AddSample(dytitle,SampleFrag::Type::SIGNAL,(EColor)(EColor::kRed),filedir+"/AFBAnalyzer_SkimTree_SMP_DYJets.root");
  else if(MODE==1){
    AddSampleWithPrefixAndWeight(dytitle+"(qqbar->DY+0parton)",SampleFrag::Type::BG,EColor::kRed,filedir+syear+"/"+analyzer+skim+"_DYJets.root","qqbar_0j_",1,filedir+syear+"/"+analyzer+skim+"_DYJets.root","qbarq_0j_",1);
    AddSampleWithPrefixAndWeight(dytitle+"(qqbar->DY+Npartons)",SampleFrag::Type::SIGNAL,EColor::kOrange,filedir+syear+"/"+analyzer+skim+"_DYJets.root","qqbar_nj_",1,filedir+syear+"/"+analyzer+skim+"_DYJets.root","qbarq_nj_",1);
    AddSampleWithPrefixAndWeight(dytitle+"(gq)",SampleFrag::Type::SIGNAL,EColor::kGreen,filedir+syear+"/"+analyzer+skim+"_DYJets.root","gq_",1,filedir+syear+"/"+analyzer+skim+"_DYJets.root","qg_",1);
    AddSampleWithPrefixAndWeight(dytitle+"(gqbar)",SampleFrag::Type::SIGNAL,EColor::kBlue,filedir+syear+"/"+analyzer+skim+"_DYJets.root","gqbar_",1,filedir+syear+"/"+analyzer+skim+"_DYJets.root","qbarg_",1);
    AddSampleWithPrefixAndWeight(dytitle+"(etc)",SampleFrag::Type::BG,EColor::kMagenta,filedir+syear+"/"+analyzer+skim+"_DYJets.root","qq_",1,filedir+syear+"/"+analyzer+skim+"_DYJets.root","qbarqbar_",1,filedir+syear+"/"+analyzer+skim+"_DYJets.root","gg_",1);
  }else if(MODE==2){
    AddSampleWithPrefixAndWeight(dytitle+"(qg)",SampleFrag::Type::SIGNAL,EColor::kOrange,filedir+syear+"/"+analyzer+skim+"_DYJets.root","qg_",1);
    AddSampleWithPrefixAndWeight(dytitle+"(gq)",SampleFrag::Type::SIGNAL,EColor::kRed,filedir+syear+"/"+analyzer+skim+"_DYJets.root","gq_",1);
    AddSampleWithPrefixAndWeight(dytitle+"(gqbar)",SampleFrag::Type::SIGNAL,EColor::kBlue,filedir+syear+"/"+analyzer+skim+"_DYJets.root","gqbar_",1);
    AddSampleWithPrefixAndWeight(dytitle+"(qbarg)",SampleFrag::Type::SIGNAL,EColor::kYellow,filedir+syear+"/"+analyzer+skim+"_DYJets.root","qbarg_",1);
    AddSampleWithPrefixAndWeight(dytitle+"(qqbar->DY+Npartons)",SampleFrag::Type::SIGNAL,EColor::kGreen,filedir+syear+"/"+analyzer+skim+"_DYJets.root","qqbar_nj_",1);
    AddSampleWithPrefixAndWeight(dytitle+"(qbarq->DY+Npartons)",SampleFrag::Type::SIGNAL,EColor::kRed,filedir+syear+"/"+analyzer+skim+"_DYJets.root","qbarq_nj_",1);
    AddSampleWithPrefixAndWeight(dytitle+"(qqbar->DY+0parton)",SampleFrag::Type::BG,EColor::kMagenta,filedir+syear+"/"+analyzer+skim+"_DYJets.root","qqbar_0j_",1);
    AddSampleWithPrefixAndWeight(dytitle+"(qbarq->DY+0parton)",SampleFrag::Type::BG,EColor::kRed,filedir+syear+"/"+analyzer+skim+"_DYJets.root","qbarq_0j_",1);
    AddSampleWithPrefixAndWeight(dytitle+"(etc)",SampleFrag::Type::BG,EColor::kMagenta,filedir+syear+"/"+analyzer+skim+"_DYJets.root","qq_",1,filedir+syear+"/"+analyzer+skim+"_DYJets.root","qbarqbar_",1,filedir+syear+"/"+analyzer+skim+"_DYJets.root","gg_",1);
  }
  if(year==2017){
    AddSampleWithPrefixAndWeight("QCD dijet",SampleFrag::Type::BG,EColor::kCyan,samples[0].files[0],"ss_",1,samples[0].files[0],"ss_",1,samples[0].files[1],"ss_",1,samples[0].files[2],"ss_",1,samples[0].files[3],"ss_",1,samples[0].files[4],"ss_",1,filedir+syear+"/"+analyzer+skim+"_DYJets.root","ss_",-1,filedir+syear+"/"+analyzer+skim+"_DYJets.root","ss_tau_",-1,filedir+syear+"/"+analyzer+skim+"_WW_pythia.root","ss_",-1,filedir+syear+"/"+analyzer+skim+"_WZ_pythia.root","ss_",-1,filedir+syear+"/"+analyzer+skim+"_ZZ_pythia.root","ss_",-1,filedir+syear+"/"+analyzer+skim+"_WJets_MG.root","ss_",-1,year==2017?filedir+"2017/"+analyzer+skim+"_TTLL_powheg.root":filedir+"2016/"+analyzer+skim+"_TT_powheg.root","ss_",-1);
  }else{
    AddSampleWithPrefixAndWeight("QCD dijet",SampleFrag::Type::BG,EColor::kCyan,samples[0].files[0],"ss_",1,samples[0].files[0],"ss_",1,samples[0].files[1],"ss_",1,samples[0].files[2],"ss_",1,samples[0].files[3],"ss_",1,samples[0].files[4],"ss_",1,samples[0].files[5],"ss_",1,samples[0].files[6],"ss_",1,filedir+syear+"/"+analyzer+skim+"_DYJets.root","ss_",-1,filedir+syear+"/"+analyzer+skim+"_DYJets.root","ss_tau_",-1,filedir+syear+"/"+analyzer+skim+"_WW_pythia.root","ss_",-1,filedir+syear+"/"+analyzer+skim+"_WZ_pythia.root","ss_",-1,filedir+syear+"/"+analyzer+skim+"_ZZ_pythia.root","ss_",-1,filedir+syear+"/"+analyzer+skim+"_WJets_MG.root","ss_",-1,year==2017?filedir+"2017/"+analyzer+skim+"_TTLL_powheg.root":filedir+"2016/"+analyzer+skim+"_TT_powheg.root","ss_",-1);
  }
*/
}
void AFBPlotter::SetupPlots(){
  if(channel==0){
    set<TString> excludes={"^/electron..../","/qqbar","/qbarq","/gq","/qg","/gqbar","/qbarg","/qq","/gg"};
    AddPlotsAuto(excludes);
  }else if(channel==1){
    set<TString> excludes={"^/muon..../","/qqbar","/qbarq","/gq","/qg","/gqbar","/qbarg","/qq","/gg"};
    AddPlotsAuto(excludes);
  }
}
void AFBPlotter::SavePlotsCondor(int channel,int year,int njob){
  system("mkdir -p condor");
  for(int i=0;i<njob;i++){
    system(Form("cd condor;echo 'echo $SKFlat_WD;cd $SKFlat_WD;source setup.sh;echo \".L Plotter/AFBAnalyzer_plot.cc \nSetup(%d,%d) \nSavePlots(\\\"AFBAnalyzer_plot\\\",%d,%d)\n.q\"|root -b'|condor_qsub --cwd -V",channel,year,njob,i));
  }
} 
void AFBPlotter::SetupSystematics(int channel,int year){
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
  AddSystematic("scalevariation",SystematicType::ENVELOPE,"_scalevariation0 _scalevariation1 _scalevariation2 _scalevariation3 _scalevariation4 _scalevariation6 _scalevariation8",false,true,false);

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
/*
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
    c=GetCompareAFBAll("_y...to.../.*"+AFB[j],0,AFB[j]);
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

*/
