#include"Plotter.cc"

class AFBPlotter:public Plotter{
public:
  void SetupSamples();
  void SetupSystematics();
  //void SetupPlots();
  //void SavePlotsCondor(int channel,int year,int njob);
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
  for(int ic=0;ic<sizeof(channels)/sizeof(TString);ic++){
    for(unsigned int ip=0;ip<periods.size();ip++){
      TString path=filedir+get<0>(periods[ip])+"/DATA/AFBAnalyzer_SkimTree_SMP_"+channels[ic]+"_";
      TString samplefragkey=channels[ic]+get<0>(periods[ip])+get<1>(periods[ip]);
      samplefrags[samplefragkey]=MakeSampleFrag(samplefragkey,SampleFrag::Type::DATA,kBlack,make_tuple(path+get<1>(periods[ip])+".root","",1.));
    }
  }

  TString syears[]={"2016","2017"};
  for(int i=0;i<sizeof(syears)/sizeof(TString);i++){
    samplefrags["muon"+syears[i]]=MakeSampleFrag("muon"+syears[i],SampleFrag::Type::DATA,kBlack);
    samplefrags["muon"+syears[i]].Add(TRegexp("DoubleMuon"+syears[i]+"[A-Z]"),1);
    samplefrags["electron"+syears[i]]=MakeSampleFrag("electron"+syears[i],SampleFrag::Type::DATA,kBlack);
    samplefrags["electron"+syears[i]].Add(TRegexp("DoubleElectron"+syears[i]+"[A-Z]"),1);

    samplefrags["dy"+syears[i]]=MakeSampleFrag("dy"+syears[i],SampleFrag::Type::SIGNAL,kRed,make_tuple(filedir+syears[i]+"/AFBAnalyzer_SkimTree_SMP_DYJets.root","",1.));
    samplefrags["powheg"+syears[i]+"muon"]=MakeSampleFrag("powheg"+syears[i]+"muon",SampleFrag::Type::SIGNAL,kBlue,make_tuple(filedir+syears[i]+"/AFBAnalyzer_ZToMuMu_M_50_120.root","",1.),make_tuple(filedir+syears[i]+"/AFBAnalyzer_ZToMuMu_M_120_200.root","",1.),make_tuple(filedir+syears[i]+"/AFBAnalyzer_ZToMuMu_M_200_400.root","",1.));
    samplefrags["powheg"+syears[i]+"electron"]=MakeSampleFrag("powheg"+syears[i]+"electron",SampleFrag::Type::SIGNAL,kBlue,make_tuple(filedir+syears[i]+"/AFBAnalyzer_ZToEE_M_50_120.root","",1.),make_tuple(filedir+syears[i]+"/AFBAnalyzer_ZToEE_M_120_200.root","",1.),make_tuple(filedir+syears[i]+"/AFBAnalyzer_ZToEE_M_200_400.root","",1.));

    samplefrags["powheg"+syears[i]+"muon"]=MakeSampleFrag("powheg"+syears[i]+"muon",SampleFrag::Type::SIGNAL,kBlue,make_tuple(filedir+syears[i]+"/AFBAnalyzer_ZToMuMu_M_50_120.root","",1.),make_tuple(filedir+syears[i]+"/AFBAnalyzer_ZToMuMu_M_120_200.root","",1.),make_tuple(filedir+syears[i]+"/AFBAnalyzer_ZToMuMu_M_200_400.root","",1.));
    samplefrags["powheg"+syears[i]+"electron"]=MakeSampleFrag("powheg"+syears[i]+"electron",SampleFrag::Type::SIGNAL,kBlue,make_tuple(filedir+syears[i]+"/AFBAnalyzer_ZToEE_M_50_120.root","",1.),make_tuple(filedir+syears[i]+"/AFBAnalyzer_ZToEE_M_120_200.root","",1.),make_tuple(filedir+syears[i]+"/AFBAnalyzer_ZToEE_M_200_400.root","",1.));

    samplefrags["dytt"+syears[i]]=MakeSampleFrag("#gamma*/Z#rightarrow#tau#tau",SampleFrag::Type::BG,EColor::kGreen,make_tuple(filedir+syears[i]+"/AFBAnalyzer_SkimTree_SMP_DYJets.root","tau_",1.));
					       
    TString path=filedir+syears[i]+"/AFBAnalyzer_SkimTree_SMP_";
    samplefrags["ww"+syears[i]]=MakeSampleFrag("ww"+syears[i],SampleFrag::Type::BG,kBlue,make_tuple(path+"WW_pythia.root","",1.));
    samplefrags["wz"+syears[i]]=MakeSampleFrag("wz"+syears[i],SampleFrag::Type::BG,kGreen,make_tuple(path+"WZ_pythia.root","",1.));
    samplefrags["zz"+syears[i]]=MakeSampleFrag("zz"+syears[i],SampleFrag::Type::BG,kYellow,make_tuple(path+"ZZ_pythia.root","",1.));
    samplefrags["vv"+syears[i]]=MakeSampleFrag("Diboson",SampleFrag::Type::BG,kBlue,make_tuple("ww"+syears[i],1.),make_tuple("wz"+syears[i],1.),make_tuple("zz"+syears[i],1.));

    samplefrags["wjets"+syears[i]]=MakeSampleFrag("W",SampleFrag::Type::BG,kYellow,make_tuple(path+"WJets_MG.root","",1.));
  }
  samplefrags["tt2016"]=MakeSampleFrag("t#bar{t}",SampleFrag::Type::BG,kMagenta,make_tuple(filedir+"2016/AFBAnalyzer_SkimTree_SMP_TT_powheg.root","",1.));
  samplefrags["tt2017"]=MakeSampleFrag("t#bar{t}",SampleFrag::Type::BG,kMagenta,make_tuple(filedir+"2017/AFBAnalyzer_SkimTree_SMP_TTLL_powheg.root","",1.));
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
  if(mode==0){
    samples["data"]=MakeSample("data",Sample::Type::SUM,kBlack,make_tuple(schannel+syear,1.));
    samples["data"].markerstyle=20;
    samples["data"].markersize=0.7;
    TString dytitle=(channel==Channel::MUON)?"#gamma*/Z#rightarrow#mu#mu":"#gamma*/Z#rightarrowee";
    samplefrags["dy"+syear].title=dytitle;
    samples["sim"]=MakeSample("sim",Sample::Type::STACK,kBlue,make_tuple("dy"+syear,1.),make_tuple("dytt"+syear,1.),make_tuple("vv"+syear,1.),make_tuple("wjets"+syear,1.),make_tuple("tt"+syear,1.));
  }else if(mode==1){
    samples["data"]=MakeSample("data",Sample::Type::SUM,kBlack,make_tuple(schannel+syear,1.),make_tuple("dytt"+syear,-1.),make_tuple("vv"+syear,-1.),make_tuple("wjets"+syear,-1.),make_tuple("tt"+syear,-1.));
    samples["data"].markerstyle=20;
    samples["data"].markersize=0.7;
    TString dytitle=(channel==Channel::MUON)?"#gamma*/Z#rightarrow#mu#mu":"#gamma*/Z#rightarrowee";
    samples["dy"]=MakeSample(dytitle,Sample::Type::SUM,kRed,make_tuple("dy"+syear,1.));
  }else if(mode==2){
    samples["data"]=MakeSample("data",Sample::Type::SUM,kBlack,make_tuple(schannel+syear,1.),make_tuple("dytt"+syear,-1.),make_tuple("vv"+syear,-1.),make_tuple("wjets"+syear,-1.),make_tuple("tt"+syear,-1.));
    samples["data"].markerstyle=20;
    samples["data"].markersize=0.7;
    samples["amc"]=MakeSample("aMC@NLO",Sample::Type::SUM,kRed,make_tuple("dy"+syear,1.));
    samples["powheg"]=MakeSample("powheg",Sample::Type::SUM,kBlue,make_tuple("powheg"+syear+schannel,1.));
  }else if(mode==3){
    samples["data"]=MakeSample("data",Sample::Type::SUM,kBlack,make_tuple(schannel+syear,1.),make_tuple("dytt"+syear,-1.),make_tuple("vv"+syear,-1.),make_tuple("wjets"+syear,-1.),make_tuple("tt"+syear,-1.));
    samples["data"].markerstyle=20;
    samples["data"].markersize=0.7;
    samples["data"].AddPrePrefix(schannel+syear+"/m60to350");
    samples["amc"]=MakeSample("aMC@NLO",Sample::Type::SUM,kRed,make_tuple("dy"+syear,1.));
    samples["amc"].AddPrePrefix(schannel+syear+"/m60to350");
    samples["powheg"]=MakeSample("powheg",Sample::Type::SUM,kBlue,make_tuple("powheg"+syear+schannel,1.));
    samples["powheg"].AddPrePrefix(schannel+syear+"/m60to350");
    samples["powheg(GEN)"]=MakeSample("powheg(GEN)",Sample::Type::SUM,kBlue,make_tuple("powhegGEN"+syear+schannel,1.));
    samples["powheg(GEN)"].AddPrePrefix(schannel+syear+"gen");
  }    
  

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
/*
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
*/
void AFBPlotter::SetupSystematics(){
  cout<<"[SetupSystematics]"<<endl;
  if(channel==Channel::ELECTRON) systematics["RECOSF"]=MakeSystematic("RECOSF",SystematicType::ENVELOPE,"_RECOSF_up _RECOSF_down",(1<<SampleFrag::Type::SIGNAL)+(1<<SampleFrag::Type::BG));
  systematics["IDSF"]=MakeSystematic("IDSF",SystematicType::ENVELOPE,"_IDSF_up _IDSF_down",(1<<SampleFrag::Type::SIGNAL)+(1<<SampleFrag::Type::BG));
  if(channel==Channel::MUON) systematics["ISOSF"]=MakeSystematic("ISOSF",SystematicType::ENVELOPE,"_ISOSF_up _ISOSF_down",(1<<SampleFrag::Type::SIGNAL)+(1<<SampleFrag::Type::BG));
  systematics["triggerSF"]=MakeSystematic("triggerSF",SystematicType::ENVELOPE,"_triggerSF_up _triggerSF_down",(1<<SampleFrag::Type::SIGNAL)+(1<<SampleFrag::Type::BG));
  systematics["PUreweight"]=MakeSystematic("PUreweight",SystematicType::ENVELOPE,"_PUreweight_up _PUreweight_down",(1<<SampleFrag::Type::SIGNAL)+(1<<SampleFrag::Type::BG));
  systematics["prefireweight"]=MakeSystematic("prefireweight",SystematicType::ENVELOPE,"_prefireweight_up _prefireweight_down",(1<<SampleFrag::Type::SIGNAL)+(1<<SampleFrag::Type::BG));
  systematics["scale"]=MakeSystematic("scale",SystematicType::ENVELOPE,"_scale_up _scale_down",(1<<SampleFrag::Type::SIGNAL)+(1<<SampleFrag::Type::BG));
  if(channel==Channel::ELECTRON) systematics["smear"]=MakeSystematic("smear",SystematicType::ENVELOPE,"_smear_up _smear_down",(1<<SampleFrag::Type::SIGNAL)+(1<<SampleFrag::Type::BG));
  systematics["alphaS"]=MakeSystematic("alphaS",SystematicType::ENVELOPE,"_alphaS_up _alphaS_down",(1<<SampleFrag::Type::SIGNAL));
  systematics["scalevariation"]=MakeSystematic("scalevariation",SystematicType::ENVELOPE,"_scalevariation0 _scalevariation1 _scalevariation2 _scalevariation3 _scalevariation4 _scalevariation6 _scalevariation8",(1<<SampleFrag::Type::SIGNAL));

  vector<TString> prefixes;
  for(int i=0;i<100;i++) prefixes.push_back(Form("_pdf%d",i));
  if(year==2017) systematics["pdf"]=MakeSystematic("pdf",SystematicType::HESSIAN,prefixes,(1<<SampleFrag::Type::SIGNAL));
  else if(year==2016) systematics["pdf"]=MakeSystematic("pdf",SystematicType::GAUSSIAN,prefixes,(1<<SampleFrag::Type::SIGNAL));
  else cout<<"###WARNING### [SetupSystematics] wrong year"<<endl;

  if(channel==Channel::ELECTRON) systematics["noRECOSF"]=MakeSystematic("noRECOSF",SystematicType::ENVELOPE,"_noRECOSF",(1<<SampleFrag::Type::SIGNAL)+(1<<SampleFrag::Type::BG));
  systematics["noIDSF"]=MakeSystematic("noIDSF",SystematicType::ENVELOPE,"_noIDSF",(1<<SampleFrag::Type::SIGNAL)+(1<<SampleFrag::Type::BG));
  if(channel==Channel::MUON) systematics["noISOSF"]=MakeSystematic("noISOSF",SystematicType::ENVELOPE,"_noISOSF",(1<<SampleFrag::Type::SIGNAL)+(1<<SampleFrag::Type::BG));
  systematics["notriggerSF"]=MakeSystematic("notriggerSF",SystematicType::ENVELOPE,"_notriggerSF",(1<<SampleFrag::Type::SIGNAL)+(1<<SampleFrag::Type::BG));
  systematics["noPUreweight"]=MakeSystematic("noPUreweight",SystematicType::ENVELOPE,"_noPUreweight",(1<<SampleFrag::Type::SIGNAL)+(1<<SampleFrag::Type::BG));
  systematics["noprefireweight"]=MakeSystematic("noprefireweight",SystematicType::ENVELOPE,"_noprefireweight",(1<<SampleFrag::Type::SIGNAL)+(1<<SampleFrag::Type::BG));
  systematics["nozptcor"]=MakeSystematic("nozptcor",SystematicType::ENVELOPE,"_nozptcor",(1<<SampleFrag::Type::SIGNAL));
  systematics["noefficiencySF"]=MakeSystematic("noefficiencySF",SystematicType::ENVELOPE,"_noefficiencySF",(1<<SampleFrag::Type::SIGNAL)+(1<<SampleFrag::Type::BG));
  if(channel==Channel::ELECTRON) systematics["IDSF_POG"]=MakeSystematic("IDSF_POG",SystematicType::ENVELOPE,"_IDSF_POG",(1<<SampleFrag::Type::SIGNAL)+(1<<SampleFrag::Type::BG));
  if(channel==Channel::ELECTRON) systematics["selective"]=MakeSystematic("selective",SystematicType::ENVELOPE,"_selective",(1<<SampleFrag::Type::DATA)+(1<<SampleFrag::Type::SIGNAL)+(1<<SampleFrag::Type::BG));
  systematics["efficiencySF"]=MakeSystematic("efficiencySF",SystematicType::MULTI,"RECOSF IDSF ISOSF triggerSF",0);
  systematics["totalsys"]=MakeSystematic("totalsys",SystematicType::MULTI,"RECOSF IDSF ISOSF triggerSF PUreweight prefireweight scale smear alphaS scalevariation pdf nozptcor",0);
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
