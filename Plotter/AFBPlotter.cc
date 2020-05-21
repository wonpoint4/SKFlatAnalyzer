#ifndef __AFBPLOTTER_CC__
#define __AFBPLOTTER_CC__
#include"Plotter.cc"
class AFBPlotter:public Plotter{
public:
  void SetupEntries();
  void SetupSystematics();
  int Setup(TString mode_);
  TString mode;
  TString analyzer;
  AFBPlotter(TString mode_="data *amc+bg");

  double GetChi2(TH1* h1,TH1* h2=NULL);
  void SetupTH4D();
  pair<double,double> GetRange(TString histname,TString axisname);
  pair<double,double> GetAFBAndError(TH1* costhetaCS);

  using Plotter::GetHist;
  virtual TH1* GetHist(const Sample& sample,Plot plot,TString additional_option="");
  virtual void GetHistActionForAdditionalClass(TObject*& obj,Plot plot);
  virtual void AddPlot(TString plotkey,TString option="");

  TH1* GetHistWeightedAFB(TH1* hist_forward_num,TH1* hist_forward_den,TH1* hist_backward_num,TH1* hist_backward_den);
  TH1* GetHistAFB(TH1* hist_forward,TH1* hist_backward);

  TCanvas* DrawPlot(TString plotkey,TString option="");
  map<TString,TString> plot_axisnames;
};
TCanvas* AFBPlotter::DrawPlot(TString plotkey,TString option=""){
  TCanvas* c=Plotter::DrawPlot(plotkey,option);
  if(c){
    Plot plot=plots[plotkey];
    plot.SetOption(option);
    for(auto [histname,axisname]:plot_axisnames){
      if(plot.name.Contains(TRegexp("/[^/]*"+histname))){
	if(plot.IsMultiPad()){
	  c->GetPad(2)->cd();
	}else{
	  c->cd();
	}
	TH1* axisowner=GetAxisParent(gPad);
	if(axisowner){
	  if(plot.histname.Contains(TRegexp("^mm[0-9]"))) axisname.ReplaceAll("ll","#mu#mu");
	  else if(plot.histname.Contains(TRegexp("^ee[0-9]"))) axisname.ReplaceAll("ll","ee");
	  axisowner->GetXaxis()->SetTitle(axisname);
	  gPad->Update();
	  gPad->Modified();
	}
      }
    }
    if(plot.IsMultiPad()){
      c->GetPad(1)->cd();
    }else{
      c->cd();
    }
    
    TH1* axisowner=GetAxisParent(gPad);
    if(axisowner){
      if(plot.name.Contains("AFB")) axisowner->GetYaxis()->SetTitle("A_{FB}"); 
      else if(plot.option.Contains("norm")) axisowner->GetYaxis()->SetTitle("Normalized");
      else if(plot.option.Contains("widthweight")) axisowner->GetYaxis()->SetTitle("Events / 1 GeV");
      //else axisowner->GetYaxis()->SetTitle(Form("Events / %g GeV",axisowner->GetBinWidth(1)));
      else axisowner->GetYaxis()->SetTitle("Events");
      
      if(gPad->GetPrimitive("TPave")){
	TLegend* leg=(TLegend*)gPad->GetPrimitive("TPave");
	for(auto obj:*leg->GetListOfPrimitives()){
	  TLegendEntry* entry=(TLegendEntry*)obj;
	  TString label=entry->GetLabel();
	  if(label.Contains("#rightarrowll")){
	    if(plot.histname.Contains(TRegexp("^mm[0-9]"))) label.ReplaceAll("ll","#mu#mu");
	    else if(plot.histname.Contains(TRegexp("^ee[0-9]"))) label.ReplaceAll("ll","ee");
	    entry->SetLabel(label);
	    break;
	  }
	}
      }

      //axisowner->SetTitle("");
      TLatex latex;
      latex.SetTextSize(0.07);
      latex.SetNDC();
      if(!(plot.IsMultiPad())) latex.SetTextSize(0.04);
      //latex.DrawLatex(0.17,0.92,"CMS #bf{#it{Preliminary}}");
      if(plot.histname.Contains(TRegexp("^[a-z]*2016/"))){
	//latex.DrawLatex(0.6,0.92,"35.92 fb^{-1} (13 TeV)");
      }else if(plot.histname.Contains(TRegexp("^[a-z]*2017/"))){
	//latex.DrawLatex(0.6,0.92,"41.53 fb^{-1} (13 TeV)");
      }else if(plot.histname.Contains(TRegexp("^[a-z]*2018/"))){
	//latex.DrawLatex(0.6,0.92,"59.74 fb^{-1} (13 TeV)");
      }	
      
      gPad->Update();
      gPad->Modified();
    }
  }
  return c;
}

pair<double,double> AFBPlotter::GetRange(TString histname,TString axisname){
  TString rangestring=histname(TRegexp(axisname+"[[-+0-9.]*,[-+0-9.]*]"));
  TString first=rangestring(axisname.Length()+1,rangestring.Index(',')-axisname.Length()-1);
  TString second=rangestring(rangestring.Index(',')+1,rangestring.Length()-rangestring.Index(',')-2);
  return make_pair(first.Atof(),second.Atof());
}  

TH1* AFBPlotter::GetHist(const Sample& sample,Plot plot,TString additional_option){
  TH1* hist=NULL;
  plot.SetOption(additional_option);
  if(plot.histname.Contains("weightedAFB")){
    TH1* hist_forward_num=GetHist(sample,plot,"histname:"+Replace(plot.histname,"weightedAFB","costhetaCS_num")+" umin:0 umax:1");
    TH1* hist_forward_den=GetHist(sample,plot,"histname:"+Replace(plot.histname,"weightedAFB","costhetaCS_den")+" umin:0 umax:1");
    TH1* hist_backward_num=GetHist(sample,plot,"histname:"+Replace(plot.histname,"weightedAFB","costhetaCS_num")+" umin:-1 umax:0");
    TH1* hist_backward_den=GetHist(sample,plot,"histname:"+Replace(plot.histname,"weightedAFB","costhetaCS_den")+" umin:-1 umax:0");
    hist=GetHistWeightedAFB(hist_forward_num,hist_forward_den,hist_backward_num,hist_backward_den);
  }else if(plot.histname.Contains("AFB")){
    TH1* hist_forward=GetHist(sample,plot,"histname:"+Replace(plot.histname,"AFB","costhetaCS")+" umin:0 umax:1");
    TH1* hist_backward=GetHist(sample,plot,"histname:"+Replace(plot.histname,"AFB","costhetaCS")+" umin:-1 umax:0");
    hist=GetHistAFB(hist_forward,hist_backward);
  }
  if(hist){
    sample.ApplyStyle(hist);
    hist->SetTitle(plot.name+plot.suffix);
    hist->SetDirectory(pdir);
    return hist;
  }else{
    return Plotter::GetHist(sample,plot);
  }
}
void AFBPlotter::GetHistActionForAdditionalClass(TObject*& obj,Plot plot){
  if(DEBUG>3) std::cout<<"###DEBUG### [AFBPlotter::GetHistActionForAdditionalClass(TObject*& obj,Plot plot)"<<endl;
  if(strstr(obj->ClassName(),"TH4D")!=NULL){
    TH4D* hist4d=(TH4D*)obj;
    int ixmin=0,iymin=0,izmin=0,iumin=0;
    int ixmax=-1,iymax=-1,izmax=-1,iumax=-1;
    if(plot.xmin||plot.xmax){
      ixmin=hist4d->GetXaxis()->FindBin(plot.xmin);
      ixmax=hist4d->GetXaxis()->FindBin(plot.xmax-0.00001);
    }
    if(plot.ymin||plot.ymax){
      iymin=hist4d->GetYaxis()->FindBin(plot.ymin);
      iymax=hist4d->GetYaxis()->FindBin(plot.ymax-0.00001);
    }
    if(plot.zmin||plot.zmax){
      izmin=hist4d->GetZaxis()->FindBin(plot.zmin);
      izmax=hist4d->GetZaxis()->FindBin(plot.zmax-0.00001);
    }
    if(plot.umin||plot.umax){
      iumin=hist4d->GetUaxis()->FindBin(plot.umin);
      iumax=hist4d->GetUaxis()->FindBin(plot.umax-0.00001);
    }
    TString axisstring=plot.histname(TRegexp("([x-zu]*)$"));
    if(axisstring=="") obj=(TObject*)hist4d->ProjectionU("_pu",ixmin,ixmax,iymin,iymax,izmin,izmax);
    else if(axisstring=="(x)") obj=(TObject*)hist4d->ProjectionX("_px",iymin,iymax,izmin,izmax,iumin,iumax);
    else if(axisstring=="(y)") obj=(TObject*)hist4d->ProjectionY("_py",ixmin,ixmax,izmin,izmax,iumin,iumax);
    else if(axisstring=="(z)") obj=(TObject*)hist4d->ProjectionZ("_pz",ixmin,ixmax,iymin,iymax,iumin,iumax);
    else if(axisstring=="(u)") obj=(TObject*)hist4d->ProjectionU("_pu",ixmin,ixmax,iymin,iymax,izmin,izmax);
    else{
      if(DEBUG>0) std::cout<<"###ERROR### [Plotter::GetHist] wrong axisstring or classname"<<endl;
    }
    delete hist4d;
  }else{
    if(DEBUG>1) std::cout<<"###WARNING## [AFBPlotter::GetHistActionForAdditionalClass(TObject*& obj,Plot plot)] Unsupported class name "<<obj->ClassName()<<endl;
  }
}
void AFBPlotter::AddPlot(TString plotkey,TString option){
  Plotter::AddPlot(plotkey,option);
  Plot& plot=plots[plotkey];
  if(plot.histname.Contains(TRegexp("m[[-+0-9.]*,[-+0-9.]*]/"))){
    auto range=GetRange(plot.histname,"m");
    plot.SetOption(Form("xmin:%f xmax:%f",range.first,range.second));
    plot.histname=Replace(plot.histname,"m[[-+0-9.]*,[-+0-9.]*]/","");
  }
  if(plot.histname.Contains(TRegexp("y[[-+0-9.]*,[-+0-9.]*]/"))){
    auto range=GetRange(plot.histname,"y");
    plot.SetOption(Form("ymin:%f ymax:%f",range.first,range.second));
    plot.histname=Replace(plot.histname,"y[[-+0-9.]*,[-+0-9.]*]/","");
  }
  if(plot.histname.Contains(TRegexp("pt[[-+0-9.]*,[-+0-9.]*]/"))){
    auto range=GetRange(plot.histname,"pt");
    plot.SetOption(Form("zmin:%f zmax:%f",range.first,range.second));
    plot.histname=Replace(plot.histname,"pt[[-+0-9.]*,[-+0-9.]*]/","");
  }
  plot.histname=Replace(plot.histname,"(m)","(x)");
  plot.histname=Replace(plot.histname,"(pt)","(z)");
  //if(plot.histname.Contains("costhetaCS")){plot.histname+="(u)";}
  if(plot.histname.Contains("dimass")){plot.histname=Replace(plot.histname,"dimass","costhetaCS");plot.histname+="(x)";}
  if(plot.histname.Contains("dirap")){plot.histname=Replace(plot.histname,"dirap","costhetaCS");plot.histname+="(y)";}
  if(plot.histname.Contains("dipt")){plot.histname=Replace(plot.histname,"dipt","costhetaCS");plot.histname+="(z)";}
  if(plot.histname.Contains("AFB")&&entries.size()>1) plot.SetOption(Form("type:%d",Plot::Type::CompareAndSig));
}

pair<double,double> AFBPlotter::GetAFBAndError(TH1* hist){
  if(!hist) return make_pair(0,0);
  TH1* to_delete=NULL;
  if(strstr(hist->ClassName(),"THStack")!=NULL){
    hist=GetTH1(hist);    
    to_delete=hist;
  }
  double ef,eb;
  int zerobin=hist->FindBin(0);
  double valf=hist->IntegralAndError(zerobin,hist->GetNbinsX()+2,ef);
  double valb=hist->IntegralAndError(0,zerobin-1,eb);
  if(to_delete) delete to_delete;
  return make_pair((valf-valb)/(valf+valb),2*sqrt(ef*ef*valb*valb+eb*eb*valf*valf)/pow(valf+valb,2));
}

TH1* AFBPlotter::GetHistWeightedAFB(TH1* hist_forward_num,TH1* hist_forward_den,TH1* hist_backward_num,TH1* hist_backward_den){
  if(DEBUG>3) std::cout<<"###DEBUG### [AFBPlotter::GetHistWeightedAFB(TH1* hist_forward_num,TH1* hist_forward_den,TH1* hist_backward_num,TH1* hist_backward_den)]"<<endl;
  if(!hist_forward_num||!hist_forward_den||!hist_backward_num||!hist_backward_den){
    if(DEBUG>1) std::cout<<"###WARING### [AFBPlotter::GetHistweightedAFB] It is null TH1*"<<endl;
    return NULL;
  }
  hist_forward_num=GetTH1(hist_forward_num,true);
  hist_forward_den=GetTH1(hist_forward_den,true);
  hist_backward_num=GetTH1(hist_backward_num,true);
  hist_backward_den=GetTH1(hist_backward_den,true);
  TH1* hist=(TH1*)hist_forward_num->Clone();
  hist->SetDirectory(pdir);
  hist->Reset();
  for(int i=0;i<hist->GetNbinsX()+2;i++){
    double valfn=hist_forward_num->GetBinContent(i);
    double valbn=hist_backward_num->GetBinContent(i);
    double valfd=hist_forward_den->GetBinContent(i);
    double efd=hist_forward_den->GetBinError(i);
    double valbd=hist_backward_den->GetBinContent(i);
    double ebd=hist_backward_den->GetBinError(i);
    hist->SetBinContent(i,3./8.*(valfd-valbd)/(valfn+valbn));
    hist->SetBinError(i,3./8.*(valbn*valfd+valfn*valbd)/pow(valfn+valbn,2)*sqrt(pow(efd/valfd,2)+pow(ebd/valbd,2)));
  }
  TString this_title=hist_forward_num->GetTitle();
  hist->SetTitle(this_title.ReplaceAll("forward_num","weightedAFB"));
  return hist;
}
TH1* AFBPlotter::GetHistAFB(TH1* hist_forward,TH1* hist_backward){
  if(DEBUG>3) std::cout<<"###DEBUG### [AFBPlotter::GetHistAFB(TH1* hist_forward,TH1* hist_backward)]"<<endl;
  if(!hist_forward||!hist_backward){
    if(DEBUG>1) std::cout<<"###WARING### [AFBPlotter::GetHistAFB] It is null TH1*"<<endl;
    return NULL;
  }
  hist_forward=GetTH1(hist_forward,true);
  hist_backward=GetTH1(hist_backward,true);
  TH1* hist=(TH1*)hist_forward->Clone();
  hist->SetDirectory(pdir);
  hist->Reset();
  for(int i=0;i<hist->GetNbinsX()+2;i++){
    double valf=hist_forward->GetBinContent(i);
    double ef=hist_forward->GetBinError(i);
    double valb=hist_backward->GetBinContent(i);
    double eb=hist_backward->GetBinError(i);
    if(valf+valb){
      hist->SetBinContent(i,(valf-valb)/(valf+valb));
      hist->SetBinError(i,2*sqrt(ef*ef*valb*valb+eb*eb*valf*valf)/pow(valf+valb,2));
    }
  }
  TString this_title=hist_forward->GetTitle();
  hist->SetTitle(this_title.ReplaceAll("forward","AFB"));
  return hist;
}  

AFBPlotter::AFBPlotter(TString mode_="dab"){
  SetupTH4D();
  vector<TString> files=Split(gSystem->GetFromPipe("find $SKFlatOutputDir$SKFlatV/AFBAnalyzer/ -type f"),"\n");
  for(int i=0;i<files.size();i++){
    TString file=files[i];
    TString key=Replace(file,TString()+getenv("SKFlatOutputDir")+getenv("SKFlatV")+"/AFBAnalyzer/","");
    Sample frag(key,Sample::Type::UNDEF,i%8+2);
    frag.files.push_back(make_tuple(file,1.,"",""));
    samples[key]=frag;
  } 
  
  //############# Testing z0 without requiring trigger firing
  //samples["data"]=Sample("data",Sample::Type::DATA,kBlack,20)+TRegexp("/DATA/AFBAnalyzer_SingleElectron_[A-Z]")+TRegexp("/DATA/AFBAnalyzer_EGamma_[A-Z]"); // For testing SingleElectron
  //samples["data"]=Sample("data",Sample::Type::DATA,kBlack,20)+TRegexp("/DATA/AFBAnalyzer_SingleMuon_[A-Z]"); // For testing SingleMuon
  //samples["data"]=Sample("data",Sample::Type::DATA,kBlack,20)+TRegexp("/DATA/AFBAnalyzer_DoubleEG_[A-Z]")+TRegexp("/DATA/AFBAnalyzer_EGamma_[A-Z]"); // For testing DoubleElectron
  //samples["data"]=Sample("data",Sample::Type::DATA,kBlack,20)+TRegexp("/DATA/AFBAnalyzer_DoubleMuon_[A-Z]"); // For testing DoubleMuon
  //samples["data"]=Sample("data",Sample::Type::DATA,kBlack,20)+TRegexp("/DATA/AFBAnalyzer_MuonEG_[A-Z]"); // For testing MuonEG

  //samples["data"]=Sample("data",Sample::Type::DATA,kBlack,20)+TRegexp("/DATA/AFBAnalyzer_SingleMuon_[A-Z]")+TRegexp("/DATA/AFBAnalyzer_SingleElectron_[A-Z]");
  samples["data"]=Sample("data",Sample::Type::DATA,kBlack,20)+TRegexp("/DATA/AFBAnalyzer_SkimTree_Dilepton_DoubleMuon_[A-Z]")+TRegexp("/DATA/AFBAnalyzer_SkimTree_Dilepton_DoubleEG_[A-Z]")+TRegexp("/DATA/AFBAnalyzer_SkimTree_Dilepton_MuonEG_[A-Z]");
  samples["mm2016"]=Sample("data (#mu#mu2016)",Sample::Type::DATA,kBlack,20)+TRegexp("2016.*/DATA/AFBAnalyzer_SkimTree_Dilepton_DoubleMuon_[A-Z]");
  samples["mm2017"]=Sample("data (#mu#mu2017)",Sample::Type::DATA,kRed,20)+TRegexp("2017.*/DATA/AFBAnalyzer_SkimTree_Dilepton_DoubleMuon_[A-Z]");
  samples["mm2018"]=Sample("data (#mu#mu2018)",Sample::Type::DATA,kBlue,20)+TRegexp("2018.*/DATA/AFBAnalyzer_SkimTree_Dilepton_DoubleMuon_[A-Z]");
  samples["ee2016"]=Sample("data (ee2016)",Sample::Type::DATA,kBlack,22)+TRegexp("2016.*/DATA/AFBAnalyzer_SkimTree_Dilepton_DoubleEG_[A-Z]");
  samples["ee2017"]=Sample("data (ee2017)",Sample::Type::DATA,kRed,22)+TRegexp("2017.*/DATA/AFBAnalyzer_SkimTree_Dilepton_DoubleEG_[A-Z]");
  samples["ee2018"]=Sample("data (ee2018)",Sample::Type::DATA,kBlue,22)+TRegexp("2018.*/DATA/AFBAnalyzer_SkimTree_Dilepton_EGamma_[A-Z]");
  samples["em2016"]=Sample("data (e#mu2016)",Sample::Type::DATA,kBlack,22)+TRegexp("2016.*/DATA/AFBAnalyzer_SkimTree_Dilepton_MuonEG_[A-Z]");
  samples["em2017"]=Sample("data (e#mu2017)",Sample::Type::DATA,kRed,22)+TRegexp("2017.*/DATA/AFBAnalyzer_SkimTree_Dilepton_MuonEG_[A-Z]");
  samples["em2018"]=Sample("data (e#mu2018)",Sample::Type::DATA,kBlue,22)+TRegexp("2018.*/DATA/AFBAnalyzer_SkimTree_Dilepton_MuonEG_[A-Z]");
  samples["mm"]=Sample("data (#mu#mu)",Sample::Type::DATA,kBlack,20)+TRegexp("/DATA/AFBAnalyzer_SkimTree_Dilepton_DoubleMuon_[A-Z]");
  samples["ee"]=Sample("data (ee)",Sample::Type::DATA,kRed,22)+TRegexp("/DATA/AFBAnalyzer_SkimTree_Dilepton_DoubleEG_[A-Z]")+TRegexp("/DATA/AFBAnalyzer_SkimTree_Dilepton_EGamma_[A-Z]");
  samples["em"]=Sample("data (e#mu)",Sample::Type::DATA,kRed,22)+TRegexp("/DATA/AFBAnalyzer_SkimTree_Dilepton_MuonEG_[A-Z]");

  samples["MC"]=Sample("MCtot",Sample::Type::SIGNAL,kRed)+TRegexp("/AFBAnalyzer_SkimTree_Dilepton_DYJets.root")+TRegexp("/AFBAnalyzer_SkimTree_Dilepton_[W-Z][W-Z]_pythia")+TRegexp("/AFBAnalyzer_SkimTree_Dilepton_TTLL_powheg")+TRegexp("/AFBAnalyzer_SkimTree_Dilepton_WJets_MG");
  samples["ww"]=Sample("WW",Sample::Type::BG,kBlue)+TRegexp("/AFBAnalyzer_SkimTree_Dilepton_WW_pythia");
  samples["wz"]=Sample("WZ",Sample::Type::BG,kGreen)+TRegexp("/AFBAnalyzer_SkimTree_Dilepton_WZ_pythia");
  samples["zz"]=Sample("ZZ",Sample::Type::BG,kCyan)+TRegexp("/AFBAnalyzer_SkimTree_Dilepton_ZZ_pythia");
  samples["vv"]=Sample("Diboson",Sample::Type::BG,kBlue)+TRegexp("/AFBAnalyzer_SkimTree_Dilepton_[W-Z][W-Z]_pythia");
  samples["wjets"]=Sample("W",Sample::Type::BG,kYellow)+TRegexp("/AFBAnalyzer_SkimTree_Dilepton_WJets_MG");
  samples["tt"]=Sample("t#bar{t}",Sample::Type::BG,kMagenta)+TRegexp("/AFBAnalyzer_SkimTree_Dilepton_TTLL_powheg");

  samples["amc"]=Sample("#gamma*/Z#rightarrowll",Sample::Type::SIGNAL,kRed)+TRegexp("/AFBAnalyzer_SkimTree_Dilepton_DYJets.root");
  samples["amcJet"]=Sample("#gamma*/Z#rightarrowll",Sample::Type::SIGNAL,kRed)+TRegexp("/AFBAnalyzer_DY[0-9]Jets.root");
  samples["amcPt"]=Sample("#gamma*/Z#rightarrowll",Sample::Type::SIGNAL,kRed)+TRegexp("/AFBAnalyzer_DYJets_Pt-[0-9]*To[0-9Inf]*.root");
  samples["amcM"]=Sample("#gamma*/Z#rightarrowll",Sample::Type::SIGNAL,kRed)+TRegexp("/AFBAnalyzer_DYJets_M-[0-9]*to[0-9Inf]*.root");
  samples["mg"]=Sample("#gamma*/Z#rightarrowll",Sample::Type::SIGNAL,kRed)+TRegexp("/AFBAnalyzer_DYJets_MG.root");
  TString dysamples[]={"amc","amcJet","amcPt","amcM","mg"};
  for(auto dysample:dysamples){
    samples["tau_"+dysample]="tau_"%(Sample("#gamma*/Z#rightarrow#tau#tau",Sample::Type::BG,kGreen)+dysample);
    samples["lhe_"+dysample]="lhe_"%(Sample("#gamma*/Z#rightarrowll (LHE)",Sample::Type::SIGNAL,kBlue)+dysample);
    samples["gen_"+dysample]="gen_"%(Sample("#gamma*/Z#rightarrowll (GEN)",Sample::Type::SIGNAL,kGreen)+dysample);
    samples["genfid_"+dysample]="genfid_"%(Sample("#gamma*/Z#rightarrowll (GEN fiducial)",Sample::Type::SIGNAL,kMagenta)+dysample);
    samples["truth_"+dysample]="truth_"%(Sample("#gamma*/Z#rightarrowll (truth)",Sample::Type::SIGNAL,kCyan)+dysample);
  }
    
  samples["amcPt_stack"]=Sample("DY Pt-binned",Sample::Type::STACK,kBlue)+TRegexp("/AFBAnalyzer_.*DYJets_Pt-[0-9]*To[0-9Inf]*.root");
  for(auto& sub:samples["amcPt_stack"].subs) sub.title=sub.title(TRegexp("Pt-[0-9]*To[0-9Inf]*"));
  samples["amcM_stack"]=Sample("DY M-binned",Sample::Type::STACK,kBlue)+TRegexp("/AFBAnalyzer_.*DYJets_M-[0-9]*to[0-9Inf]*.root");
  for(auto& sub:samples["amcM_stack"].subs) sub.title=sub.title(TRegexp("M-[0-9]*to[0-9Inf]*"));
  
  plot_axisnames["dimass"]="m(ll) [GeV]";
  plot_axisnames["dirap"]="y(ll)";
  plot_axisnames["dipt"]="p_{T}(ll) [GeV]";
  plot_axisnames["costhetaCS"]="cos#theta_{CS}";
  plot_axisnames["AFB(m)"]="m(ll) [GeV]";
  plot_axisnames["AFB(y)"]="y(ll)";
  plot_axisnames["AFB(pt)"]="p_{T}(ll) [GeV]";

  Setup(mode_);
}

int AFBPlotter::Setup(TString mode_){
  Reset();
  if(mode_=="dab"||mode_=="data:amc+bg") mode="data:amc+bg";
  else if(mode_=="dabs"||mode_=="data:amc+bg+ss") mode="data:amc+bg+ss";
  else if(mode_=="dba"||mode_=="data-bg:amc") mode="data-bg:amc";
  else if(mode_=="dajb"||mode_=="data:amcJet+bg") mode="data:amcJet+bg";
  else if(mode_=="dajbs"||mode_=="data:amcJet+bg+ss") mode="data:amcJet+bg+ss";
  else if(mode_=="dbaj"||mode_=="data-bg:amcJet") mode="data-bg:amcJet";
  else if(mode_=="damb"||mode_=="data:amcM+bg") mode="data:amcM+bg";
  else if(mode_=="m"||mode_=="muon") mode="muon";
  else if(mode_=="e"||mode_=="electron") mode="electron";
  else if(mode_=="ae"||mode_=="amc_electron") mode="amc_electron";
  else if(mode_=="aje"||mode_=="amcJet_electron") mode="amcJet_electron"; 
  else if(mode_=="am"||mode_=="amc_muon") mode="amc_muon";
  else if(mode_=="ajm"||mode_=="amcJet_muon") mode="amcJet_muon";
  else if(mode_=="gae"||mode_=="genamc_electron") mode="genamc_electron";
  else if(mode_=="gam"||mode_=="genamc_muon") mode="genamc_muon";
  else if(mode_=="aga"||mode_=="amc:genamc") mode="amc:genamc";
  else if(mode_=="wr") mode="wr";
  else if(mode_=="appamm") mode="amc_pp:amc_mm";
  else if(mode_=="mc") mode="mc";
  else {
    mode=mode_;
    //cout<<"Invalid mode "<<mode_<<endl;
    //exit(1);
  }

  SetupEntries();
  SetupSystematics();
  SetupPlots("plot_AFBAnalyzer/"+mode+"/plot.dat");

  if(DEBUG) std::cout<<"[AFBPlotter::Setup] nentries: "<<entries.size()<<endl;
  if(DEBUG) std::cout<<"[AFBPlotter::Setup] nsys: "<<systematics.size()<<endl;
  if(DEBUG) std::cout<<"[AFBPlotter::Setup] nplot: "<<plots.size()<<endl;

  return 1;
}
 
void AFBPlotter::SetupEntries(){
  if(DEBUG) std::cout<<"[AFBPlotter::SetupEntries] mode="<<mode<<endl;
  if(mode.Contains(TRegexp("muon$"))||mode.Contains(TRegexp("electron$"))){
    TString schannel;
    if(mode.Contains(TRegexp("muon$"))) schannel="muon";
    else if(mode.Contains(TRegexp("electron$"))) schannel="electron";
    vector<TString> syears={"2016","2017","2018"};
    for(int i=0;i<syears.size();i++){
      if(mode.Contains(TRegexp("^amc_"))) entries.push_back(Sample("amc"+syears[i],Sample::Type::SIGNAL,i+2)+TRegexp("^amc"+syears[i]+"$"));
      else if(mode.Contains(TRegexp("^genamc_"))) entries.push_back(Sample("genamc"+syears[i],Sample::Type::SIGNAL,i+2)+TRegexp("^genamc"+syears[i]+"$"));
      else if(mode.Contains(TRegexp("^amcJet_"))) entries.push_back(Sample("amcJet"+syears[i],Sample::Type::SIGNAL,i+2)+TRegexp("^amcJet"+syears[i]+"$"));
      else entries.push_back(Sample(schannel+syears[i],Sample::Type::DATA,i+2)+TRegexp("^"+schannel+syears[i]+"$"));
      entries.back().style.markerstyle=20+i;
      //entries.back().prefix=schannel+syears[i]+"/";
    }
  }else if(mode=="wr"){
    entries.push_back(Sample("W_{R}#rightarrowlN#rightarrowlljj",Sample::Type::B,kCyan+1)+TRegexp("WRtoNLtoLLJJ_WR1600_N200"));
    entries.back().style.fillcolor=0;
    entries.back().style.linewidth=2;
    entries.back().style.drawoption="hist e";
    entries.push_back(Sample("standard model",Sample::Type::STACK,kRed)+TRegexp("^amcJet201[6-8]$")+TRegexp("^amcJettt201[6-8]$")+TRegexp("^vv201[6-8]$")+TRegexp("^wjets201[6-8]$")+TRegexp("^tt201[6-8]$"));
    //entries.push_back(Sample("data",Sample::Type::DATA,kBlack)+TRegexp("^muon201[6-8]$")+TRegexp("^electron201[6-8]$"));
  }else if(mode=="amc_pp:amc_mm"){
    entries.push_back("pp_"%(Sample("e^{+}e^{+}",Sample::Type::SIGNAL,kRed)+TRegexp("^amcJet201[6-8]$")));
    entries.back().style.linewidth=2;
    entries.back().style.fillcolor=0;
    entries.back().style.drawoption="hist e";
    entries.push_back("mm_"%(Sample("e^{-}e^{-}",Sample::Type::SIGNAL,kBlue)+TRegexp("^amcJet201[6-8]$")));
    entries.back().style.linewidth=2;
    entries.back().style.fillcolor=0;
    entries.back().style.drawoption="hist e";
  }else if(mode=="mc"){
    entries.push_back(samples["amc2016"]+samples["amc2017"]+samples["amc2018"]);
    //entries.push_back(samples["amctt2016"]+samples["amctt2017"]+samples["amctt2018"]);
    entries.push_back(samples["tt2016"]+samples["tt2017"]+samples["tt2018"]);
    //entries.push_back(samples["vv2016"]+samples["vv2017"]+samples["vv2018"]);
    entries.push_back(samples["ww2016"]+samples["ww2017"]+samples["ww2018"]);
    //entries.push_back(samples["wz2016"]+samples["wz2017"]+samples["wz2018"]);
    //entries.push_back(samples["zz2016"]+samples["zz2017"]+samples["zz2018"]);
    //entries.push_back(samples["wjets2016"]+samples["wjets2017"]+samples["wjets2018"]);
  }else{
    vector<TString> entry_keys=Split(mode," ");
    for(auto entry_key:entry_keys){
      Sample entry;
      if(entry_key.BeginsWith("*")) entry=Sample("simulation",Sample::Type::STACK);
      else entry=Sample("sum",Sample::Type::SUM);
      TPRegexp("^\\*").Substitute(entry_key,"");
      TPRegexp("([+-])").Substitute(entry_key," $1","g");
      vector<TString> sample_keys=Split(entry_key," ");
      for(int i=0;i<(int)sample_keys.size();i++){
	TString sample_key=Replace(sample_keys[i],"[+-]","");
	Sample sample;
	if(samples.find(sample_key)!=samples.end()){
	  sample=samples[sample_key];
	}else if(sample_key=="bg"){
	  sample=Sample("bg",Sample::Type::STACK,kBlue);
	  if(mode.Contains("amcPt")) sample=sample+"tau_amcPt";
	  else if(mode.Contains("amcM")) sample=sample+"tau_amcM";
	  else if(mode.Contains("amc")) sample=sample+"tau_amc";
	  sample=sample+"vv"+"wjets"+"tt";
	}else if(sample_key=="ss"){
	  sample=Sample("QCD dijet",Sample::Type::BG,kCyan)+"data";
	  if(mode.Contains("amcPt")) sample=sample-"tau_amcPt";
	  else if(mode.Contains("amcM")) sample=sample-"tau_amcM";
	  else if(mode.Contains("amc")) sample=sample-"tau_amc";
	  sample="ss_"%(sample-"vv"-"wjets"-"tt");
	}else{
	  if(DEBUG>0) std::cout<<"###ERROR### [AFBPlotter::SetupEntries] No "<<entry_key<<" in samples"<<endl;
	  exit(1);
	}
	if(i==0){
	  if(entry.type==Sample::Type::SUM) entry.title=sample.title;
	  entry.style=sample.style;
	  entry=entry+sample;
	}else if(sample_keys[i].BeginsWith("+")) entry=entry+sample;
	else if(sample_keys[i].BeginsWith("-")) entry=entry-sample;
	else{
	  if(DEBUG>0) std::cout<<"###ERROR### [AFBPlotter::SetupEntries] Wrong syntax "<<entry_key<<endl;
	  exit(1);
	}
      }
      entries.push_back(entry);
    }
  }
}
  
  
void AFBPlotter::SetupSystematics(){
  if(DEBUG>0)  std::cout<<"[AFBPlotter::SetupSystematics]"<<endl;
  systematics["RECOSF"]=MakeSystematic("RECOSF",Systematic::Type::ENVELOPE,(1<<Sample::Type::SIGNAL)+(1<<Sample::Type::BG),"_RECOSF_up _RECOSF_down");
  systematics["IDSF"]=MakeSystematic("IDSF",Systematic::Type::ENVELOPE,(1<<Sample::Type::SIGNAL)+(1<<Sample::Type::BG),"_IDSF_up _IDSF_down");
  systematics["ISOSF"]=MakeSystematic("ISOSF",Systematic::Type::ENVELOPE,(1<<Sample::Type::SIGNAL)+(1<<Sample::Type::BG),"_ISOSF_up _ISOSF_down");
  systematics["triggerSF"]=MakeSystematic("triggerSF",Systematic::Type::ENVELOPE,(1<<Sample::Type::SIGNAL)+(1<<Sample::Type::BG),"_triggerSF_up _triggerSF_down");
  systematics["PUreweight"]=MakeSystematic("PUreweight",Systematic::Type::ENVELOPE,(1<<Sample::Type::SIGNAL)+(1<<Sample::Type::BG),"_PUreweight_up _PUreweight_down");
  systematics["prefireweight"]=MakeSystematic("prefireweight",Systematic::Type::ENVELOPE,(1<<Sample::Type::SIGNAL)+(1<<Sample::Type::BG),"_prefireweight_up _prefireweight_down");
  systematics["scale"]=MakeSystematic("scale",Systematic::Type::ENVELOPE,(1<<Sample::Type::SIGNAL)+(1<<Sample::Type::BG),"_scale_up _scale_down");
  systematics["smear"]=MakeSystematic("smear",Systematic::Type::ENVELOPE,(1<<Sample::Type::SIGNAL)+(1<<Sample::Type::BG),"_smear_up _smear_down");
  systematics["alphaS"]=MakeSystematic("alphaS",Systematic::Type::ENVELOPE,(1<<Sample::Type::SIGNAL),"_alphaS_up _alphaS_down");
  systematics["scalevariation"]=MakeSystematic("scalevariation",Systematic::Type::ENVELOPE,(1<<Sample::Type::SIGNAL),"_scalevariation0 _scalevariation1 _scalevariation2 _scalevariation3 _scalevariation4 _scalevariation6 _scalevariation8");

  vector<TString> prefixes;
  for(int i=0;i<100;i++) prefixes.push_back(Form("_pdf%d",i));
  systematics["pdfh"]=MakeSystematic("pdfh",Systematic::Type::HESSIAN,(1<<Sample::Type::SIGNAL),prefixes);
  systematics["pdfg"]=MakeSystematic("pdfg",Systematic::Type::GAUSSIAN,(1<<Sample::Type::SIGNAL),prefixes);

  systematics["noRECOSF"]=MakeSystematic("noRECOSF",Systematic::Type::ENVELOPE,(1<<Sample::Type::SIGNAL)+(1<<Sample::Type::BG),"_noRECOSF");
  systematics["noIDSF"]=MakeSystematic("noIDSF",Systematic::Type::ENVELOPE,(1<<Sample::Type::SIGNAL)+(1<<Sample::Type::BG),"_noIDSF");
  systematics["noISOSF"]=MakeSystematic("noISOSF",Systematic::Type::ENVELOPE,(1<<Sample::Type::SIGNAL)+(1<<Sample::Type::BG),"_noISOSF");
  systematics["notriggerSF"]=MakeSystematic("notriggerSF",Systematic::Type::ENVELOPE,(1<<Sample::Type::SIGNAL)+(1<<Sample::Type::BG),"_notriggerSF");
  systematics["noPUreweight"]=MakeSystematic("noPUreweight",Systematic::Type::ENVELOPE,(1<<Sample::Type::SIGNAL)+(1<<Sample::Type::BG),"_noPUreweight");
  systematics["noprefireweight"]=MakeSystematic("noprefireweight",Systematic::Type::ENVELOPE,(1<<Sample::Type::SIGNAL)+(1<<Sample::Type::BG),"_noprefireweight");
  systematics["nozptcor"]=MakeSystematic("nozptcor",Systematic::Type::ENVELOPE,(1<<Sample::Type::SIGNAL),"_nozptcor");
  systematics["noefficiencySF"]=MakeSystematic("noefficiencySF",Systematic::Type::ENVELOPE,(1<<Sample::Type::SIGNAL)+(1<<Sample::Type::BG),"_noefficiencySF");
  systematics["IDSF_POG"]=MakeSystematic("IDSF_POG",Systematic::Type::ENVELOPE,(1<<Sample::Type::SIGNAL)+(1<<Sample::Type::BG),"_IDSF_POG");
  systematics["selective"]=MakeSystematic("selective",Systematic::Type::ENVELOPE,(1<<Sample::Type::DATA)+(1<<Sample::Type::SIGNAL)+(1<<Sample::Type::BG),"_selective");
  systematics["efficiencySF"]=MakeSystematic("efficiencySF",Systematic::Type::MULTI,0,"RECOSF IDSF ISOSF triggerSF");
  systematics["totalsys"]=MakeSystematic("totalsys",Systematic::Type::MULTI,0,"RECOSF IDSF ISOSF triggerSF PUreweight prefireweight scale smear alphaS scalevariation pdf nozptcor");
}

void AFBPlotter::SetupTH4D(){
  if(TClass::GetClass("TH4D")==NULL){
    gROOT->ProcessLine(".L TH4D/TH4D.cxx+");
  }
}


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
    
  
#endif
