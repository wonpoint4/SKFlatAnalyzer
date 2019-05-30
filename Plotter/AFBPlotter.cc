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
  TString analyzer;
  AFBPlotter();

  double GetChi2(TH1* h1,TH1* h2=NULL);
};
AFBPlotter::AFBPlotter(){
  analyzer="AFBAnalyzer";
  TObjArray* arr=gSystem->GetFromPipe("find $SKFlatOutputDir$SKFlatV/"+analyzer+"/ -type f").Tokenize("\n");
  for(int i=0;i<arr->GetEntries();i++){
    TString file=((TObjString*)arr->At(i))->String();
    TString year=file("/"+analyzer+"/[0-9]+/"); year=year(analyzer.Length()+2,4);
    TString key=file(analyzer+"_.*\\.root");
    key=key(analyzer.Length()+1,key.Index(".root")-analyzer.Length()-1)+"_"+year;
    SampleFrag frag;
    frag.files.push_back(make_tuple(file,"",1.));
    samplefrags[key]=frag;
  } 
  
  vector<TString> syears={"2016","2017","2018"};
  for(const auto& syear:syears){
    samplefrags["muon"+syear]=MakeSampleFrag("muon"+syear,SampleFrag::Type::DATA,kBlack,TRegexp("SkimTree_SMP_DoubleMuon_[A-Z].*_"+syear));
    samplefrags["electron"+syear]=MakeSampleFrag("electron"+syear,SampleFrag::Type::DATA,kBlack,TRegexp("SkimTree_SMP_.*EG.*_[A-Z].*_"+syear));

    samplefrags["amc"+syear]=MakeSampleFrag("amc"+syear,SampleFrag::Type::SIGNAL,kRed,TRegexp("SkimTree_SMP_DYJets_"+syear)); //FIXME no 2018
    samplefrags["amcPt"+syear]=MakeSampleFrag("amc"+syear,SampleFrag::Type::SIGNAL,kRed,TRegexp("SkimTree_SMP_DYJets_Pt-[0-9]+To[0-9]+_"+syear));
    samplefrags["mg"+syear]=MakeSampleFrag("mg"+syear,SampleFrag::Type::SIGNAL,kRed,TRegexp("SkimTree_SMP_DYJets_MG_"+syear));
    samplefrags["powhegmuon"+syear]=MakeSampleFrag("powhegmuon"+syear,SampleFrag::Type::SIGNAL,kBlue,TRegexp("ZToMuMu_M_[0-9]+_[0-9]+_"+syear));
    samplefrags["powhegelectron"+syear]=MakeSampleFrag("powhegelectron"+syear,SampleFrag::Type::SIGNAL,kBlue,TRegexp("ZToEE_M_[0-9]+_[0-9]+_"+syear));

    samplefrags["genamc"+syear]=MakeSampleFrag("amc"+syear,SampleFrag::Type::GEN,kRed,TRegexp("SkimTree_GEN_DYJets_"+syear),"gen_"); //FIXME no 2018
    samplefrags["genpowhegmuon"+syear]=MakeSampleFrag("powhegmuon"+syear,SampleFrag::Type::SIGNAL,kBlue,TRegexp("ZToMuMu_M_[0-9]+_[0-9]+_"+syear),"gen_");
    samplefrags["genpowhegelectron"+syear]=MakeSampleFrag("genpowhegelectron"+syear,SampleFrag::Type::SIGNAL,kBlue,TRegexp("ZToEE_M_[0-9]+_[0-9]+_"+syear),"gen_");

    //FIXME no 2018
    //if(syear!="2018") samplefrags["amctt"+syear]=MakeSampleFrag("#gamma*/Z#rightarrow#tau#tau",SampleFrag::Type::BG,EColor::kGreen,"SkimTree_SMP_DYJets_"+syear,"tau_");
    //FIXME no 2016
    //samplefrags["mgtt"+syear]=MakeSampleFrag("#gamma*/Z#rightarrow#tau#tau",SampleFrag::Type::BG,EColor::kGreen,TRegexp("SkimTree_SMP_DYJets_MG_"+syear),"tau_");

    samplefrags["vv"+syear]=MakeSampleFrag("Diboson",SampleFrag::Type::BG,kBlue,TRegexp("SkimTree_SMP_[W-Z][W-Z]_pythia_"+syear));
    samplefrags["wjets"+syear]=MakeSampleFrag("W",SampleFrag::Type::BG,kYellow,"SkimTree_SMP_WJets_MG_"+syear);
    samplefrags["tt"+syear]=MakeSampleFrag("t#bar{t}",SampleFrag::Type::BG,kMagenta,"SkimTree_SMP_TTLL_powheg_"+syear);
  }
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

  if(DEBUG) cout<<"[Setup] nsample: "<<samples.size()<<endl;
  if(DEBUG) cout<<"[Setup] nsys: "<<systematics.size()<<endl;
  if(DEBUG) cout<<"[Setup] nplot: "<<plots.size()<<endl;

  return 1;
}

void AFBPlotter::SetupSamples(){
  if(DEBUG)  cout<<"[AFBPlotter::SetupSamples]"<<endl;
  if(mode==0){
    samples["data"]=MakeSample("data",Sample::Type::SUM,kBlack,make_tuple(schannel+syear,1.));
    samples["data"].markerstyle=20;
    samples["data"].markersize=0.7;
    if(syear=="2018") samples["sim"]=MakeSample("sim",Sample::Type::STACK,kBlue,make_tuple("mg"+syear,1.),make_tuple("mgtt"+syear,1.),make_tuple("vv"+syear,1.),make_tuple("wjets"+syear,1.),make_tuple("tt"+syear,1.));
    else samples["sim"]=MakeSample("sim",Sample::Type::STACK,kBlue,make_tuple("amc"+syear,1.),make_tuple("amctt"+syear,1.),make_tuple("vv"+syear,1.),make_tuple("wjets"+syear,1.),make_tuple("tt"+syear,1.));
    get<0>(samples["sim"].frags[0]).title=(channel==Channel::MUON)?"#gamma*/Z#rightarrow#mu#mu":"#gamma*/Z#rightarrowee";
  }else if(mode==1){
    samples["data"]=MakeSample("data",Sample::Type::SUM,kBlack,make_tuple(schannel+syear,1.),make_tuple("amctt"+syear,-1.),make_tuple("vv"+syear,-1.),make_tuple("wjets"+syear,-1.),make_tuple("tt"+syear,-1.));
    samples["data"].markerstyle=20;
    samples["data"].markersize=0.7;
    samples["amc"]=MakeSample((channel==Channel::MUON)?"#gamma*/Z#rightarrow#mu#mu":"#gamma*/Z#rightarrowee",Sample::Type::SUM,kRed,make_tuple("amc"+syear,1.));
  }else if(mode==2){
    samples["data"]=MakeSample("data",Sample::Type::SUM,kBlack,make_tuple(schannel+syear,1.),make_tuple("amctt"+syear,-1.),make_tuple("vv"+syear,-1.),make_tuple("wjets"+syear,-1.),make_tuple("tt"+syear,-1.));
    samples["data"].markerstyle=20;
    samples["data"].markersize=0.7;
    samples["amc"]=MakeSample("aMC@NLO",Sample::Type::SUM,kRed,make_tuple("amc"+syear,1.));
    samples["genamc"]=MakeSample("aMC@NLO(Gen)",Sample::Type::SUM,kMagenta,make_tuple("genamc"+syear,1.));
  }else if(mode==3){
    samples["data"]=MakeSample("data",Sample::Type::SUM,kBlack,make_tuple(schannel+syear,1.),make_tuple("amctt"+syear,-1.),make_tuple("vv"+syear,-1.),make_tuple("wjets"+syear,-1.),make_tuple("tt"+syear,-1.));
    samples["data"].markerstyle=20;
    samples["data"].markersize=0.7;
    samples["amc"]=MakeSample("aMC@NLO",Sample::Type::SUM,kRed,make_tuple("amc"+syear,1.));
    samples["powheg"]=MakeSample("powheg",Sample::Type::SUM,kBlue,make_tuple("powheg"+schannel+syear,1.));
  }else if(mode==4){
    samples["data"]=MakeSample("data",Sample::Type::SUM,kBlack,make_tuple(schannel+syear,1.),make_tuple("amctt"+syear,-1.),make_tuple("vv"+syear,-1.),make_tuple("wjets"+syear,-1.),make_tuple("tt"+syear,-1.));
    samples["data"].markerstyle=20;
    samples["data"].markersize=0.7;
    samples["amc"]=MakeSample("aMC@NLO",Sample::Type::SUM,kRed,make_tuple("amc"+syear,1.));
    samples["genamc"]=MakeSample("aMC@NLO(Gen)",Sample::Type::SUM,kMagenta,make_tuple("genamc"+syear,1.));
    samples["powheg"]=MakeSample("powheg",Sample::Type::SUM,kBlue,make_tuple("powheg"+schannel+syear,1.));
    samples["powheg(GEN)"]=MakeSample("powheg(GEN)",Sample::Type::SUM,kGreen,make_tuple("genpowheg"+schannel+syear,1.));
  }else if(mode==10){ // with and without EfficiencySF
    samples["data"]=MakeSample("data",Sample::Type::SUM,kBlack,make_tuple(schannel+syear,1.));
    samples["data"].markerstyle=20;
    samples["data"].markersize=0.7;

    if(syear=="2018") samples["sim"]=MakeSample("sim",Sample::Type::SUM,kRed,make_tuple("mg"+syear,1.),make_tuple("mgtt"+syear,1.),make_tuple("vv"+syear,1.),make_tuple("wjets"+syear,1.),make_tuple("tt"+syear,1.));
    else samples["sim"]=MakeSample("sim",Sample::Type::SUM,kRed,make_tuple("amc"+syear,1.),make_tuple("amctt"+syear,1.),make_tuple("vv"+syear,1.),make_tuple("wjets"+syear,1.),make_tuple("tt"+syear,1.));

    if(syear=="2018") samples["sim2"]=MakeSample("sim (w/o EffiSF)",Sample::Type::SYS,kBlue,make_tuple("mg"+syear,1.),make_tuple("mgtt"+syear,1.),make_tuple("vv"+syear,1.),make_tuple("wjets"+syear,1.),make_tuple("tt"+syear,1.));
    else samples["sim2"]=MakeSample("sim (w/o EffiSF)",Sample::Type::SYS,kBlue,make_tuple("amc"+syear,1.),make_tuple("amctt"+syear,1.),make_tuple("vv"+syear,1.),make_tuple("wjets"+syear,1.),make_tuple("tt"+syear,1.));
  }else if(mode==11){ // majority vs selective
    samples["data"]=MakeSample("data (majority)",Sample::Type::SUM,kBlack,make_tuple(schannel+syear,1.));
    samples["data"].markerstyle=20;
    samples["data"].markersize=0.7;

    samples["data2"]=MakeSample("data (selective)",Sample::Type::SYS,kRed,make_tuple(schannel+syear,1.));
    samples["data2"].markerstyle=23;
    samples["data2"].markersize=0.7;
  }else if(mode==12){ // PFISO vs trkISO
    samples["data"]=MakeSample("data (PFISO)",Sample::Type::SUM,kBlack,make_tuple(schannel+syear,1.));
    samples["data"].markerstyle=20;
    samples["data"].markersize=0.7;

    samples["data2"]=MakeSample("data (TrkISO)",Sample::Type::SYS,kRed,make_tuple(schannel+syear,1.));
    samples["data2"].markerstyle=23;
    samples["data2"].markersize=0.7;
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
  if(DEBUG)  cout<<"[SetupSystematics]"<<endl;
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
double AFBPlotter::GetChi2(TH1* h1,TH1* h2){
  double chi2=0;
  for(int i=h1->GetXaxis()->GetFirst();i<h1->GetXaxis()->GetLast()+1;i++){
    double x1=h1->GetBinContent(i);
    double ex1=h1->GetBinError(i);
    double x2=h2?h2->GetBinContent(i):0.;
    double ex2=h2?h2->GetBinError(i):0.;
    chi2+=pow((x1-x2)/(ex1-ex2),2);
  }
  chi2/=h1->GetXaxis()->GetLast()-h1->GetXaxis()->GetFirst()+1;
  return chi2;
}
    
  
