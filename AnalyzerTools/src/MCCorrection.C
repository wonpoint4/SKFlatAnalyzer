#include "MCCorrection.h"

MCCorrection::MCCorrection() : 
IgnoreNoHist(false)
{

  histDir = TDirectoryHelper::GetTempDirectory("MCCorrection");

}

void MCCorrection::ReadHistograms(){

  TString datapath = getenv("DATA_DIR");

  TDirectory* origDir = gDirectory;

  //==== ID/Trigger
  TString IDpath = datapath+"/"+TString::Itoa(DataYear,10)+"/ID/";

  string elline;
  ifstream in(IDpath+"/Electron/histmap.txt");
  while(getline(in,elline)){
    std::istringstream is( elline );

    TString tstring_elline = elline;
    if(tstring_elline.Contains("#")) continue;

    TString a,b,c,d,e,f;
    is >> a; // ID,RERCO
    is >> b; // Eff,SF
    is >> c; // <WPnames>
    is >> d; // <rootfilename>
    is >> e; // <histname>
    is >> f; // Class
    TFile *file = new TFile(IDpath+"/Electron/"+d);

    if(f=="TH2F"){
      histDir->cd();
      map_hist_Electron[a+"_"+b+"_"+c] = (TH2F *)file->Get(e)->Clone();
    }
    else if(f=="TGraphAsymmErrors"){
      histDir->cd();
      map_graph_Electron[a+"_"+b+"_"+c] = (TGraphAsymmErrors *)file->Get(e)->Clone();
    }
    else{
      cout << "[MCCorrection::MCCorrection] Wrong class type : " << elline << endl;
    }
    file->Close();
    delete file;
    origDir->cd();
  }

  cout << "[MCCorrection::MCCorrection] map_hist_Electron :" << endl;
  for(std::map< TString, TH2F* >::iterator it=map_hist_Electron.begin(); it!=map_hist_Electron.end(); it++){
    cout << it->first << endl;
  }
  cout << "[MCCorrection::MCCorrection] map_graph_Electron :" << endl;
  for(std::map< TString, TGraphAsymmErrors* >::iterator it=map_graph_Electron.begin(); it!=map_graph_Electron.end(); it++){
    cout << it->first << endl;
  }


  string elline2;
  ifstream in2(IDpath+"/Muon/histmap.txt");
  while(getline(in2,elline2)){
    std::istringstream is( elline2 );

    TString tstring_elline = elline2;
    if(tstring_elline.Contains("#")) continue;

    TString a,b,c,d,e;
    is >> a; // ID,RERCO
    is >> b; // Eff,SF
    is >> c; // <WPnames>
    is >> d; // <rootfilename>
    is >> e; // <histname>
    TFile *file = new TFile(IDpath+"/Muon/"+d);
    histDir->cd();
    map_hist_Muon[a+"_"+b+"_"+c] = (TH2F *)file->Get(e)->Clone();
    file->Close();
    delete file;
    origDir->cd();
  }

  cout << "[MCCorrection::MCCorrection] map_hist_Muon :" << endl;
  for(std::map< TString, TH2F* >::iterator it=map_hist_Muon.begin(); it!=map_hist_Muon.end(); it++){
    cout << it->first << endl;
  }


  // == Get Prefiring maps
  TString PrefirePath  = datapath+"/"+TString::Itoa(DataYear,10)+"/Prefire/";

  string elline3;
  ifstream in3(PrefirePath+"/histmap.txt");
  while(getline(in3,elline3)){
    std::istringstream is( elline3 );

    TString tstring_elline = elline3;
    if(tstring_elline.Contains("#")) continue;

    TString a,b,c;
    is >> a; // Jet, Photon
    is >> b; // <rootfilename>
    is >> c; // <histname>
    
    TFile *file = new TFile(PrefirePath+b);
    histDir->cd();
    map_hist_prefire[a + "_prefire"] = (TH2F *)file->Get(c)->Clone();
    file->Close();
    delete file;
    origDir->cd();
  }


  // == Get Pileup Reweight maps
  TString PUReweightPath = datapath+"/"+TString::Itoa(DataYear,10)+"/PileUp/";

  string elline4;
  ifstream  in4(PUReweightPath+"/histmap.txt");
  while(getline(in4,elline4)){
    std::istringstream is( elline4 );

    TString tstring_elline = elline4;
    if(tstring_elline.Contains("#")) continue;

    TString a,b,c;
    is >> a; // sample name
    is >> b; // syst
    is >> c; // rootfile name

    if(DataYear == 2017 && a!=MCSample) continue;
    
    TFile *file = new TFile(PUReweightPath+c);
    if( (TH1D *)file->Get(a+"_"+b) ){
      histDir->cd();
      map_hist_pileup[a+"_"+b+"_pileup"] = (TH1D *)file->Get(a+"_"+b)->Clone();
    }
    else{
      cout << "[MCCorrection::ReadHistograms] No : " << a + "_" + b << endl;
    }
    file->Close();
    delete file;
    origDir->cd();
  }
/*
  cout << "[MCCorrection::MCCorrection] map_hist_pileup :" << endl;
  for(std::map< TString, TH1D* >::iterator it=map_hist_pileup.begin(); it!=map_hist_pileup.end(); it++){
    cout << it->first << endl;
  }
*/

}

MCCorrection::~MCCorrection(){

}

void MCCorrection::SetMCSample(TString s){
  MCSample = s;
}
void MCCorrection::SetDataYear(int i){
  DataYear = i;
}

double MCCorrection::MuonID_SF(TString ID, double eta, double pt, int sys){

  if(ID=="Default") return 1.;

  //cout << "[MCCorrection::MuonID_SF] ID = " << ID << endl;
  //cout << "[MCCorrection::MuonID_SF] eta = " << eta << ", pt = " << pt << endl;

  double value = 1.;
  double error = 0.;

  if(DataYear!=2016){
    eta = fabs(eta);
  }

  if(ID=="NUM_TightID_DEN_genTracks" || ID=="NUM_HighPtID_DEN_genTracks"){
    //==== boundaries
    if(pt<20.) pt = 20.;
    if(pt>=120.) pt = 119.;
    if(eta>=2.4) eta = 2.39;
    if(eta<-2.4) eta = -2.4;
  }

  TH2F *this_hist = map_hist_Muon["ID_SF_"+ID];
  if(!this_hist){
    if(IgnoreNoHist) return 1.;
    else{
      cout << "[MCCorrection::MuonID_SF] No "<<"ID_SF_"+ID<<endl;
      exit(EXIT_FAILURE);
    }
  }

  int this_bin(-999);

  if(DataYear==2016){
    this_bin = this_hist->FindBin(eta,pt);
  }
  else{
    this_bin = this_hist->FindBin(pt,eta);
  }

  value = this_hist->GetBinContent(this_bin);
  error = this_hist->GetBinError(this_bin);

  //cout << "[MCCorrection::MuonID_SF] value = " << value << endl;

  return value+double(sys)*error;

}

double MCCorrection::MuonISO_SF(TString ID, double eta, double pt, int sys){

  if(ID=="Default") return 1.;

  //cout << "[MCCorrection::MuonISO_SF] eta = " << eta << ", pt = " << pt << endl;

  double value = 1.;
  double error = 0.;

  if(DataYear!=2016){
    eta = fabs(eta);
  }

  if(ID=="NUM_TightRelIso_DEN_TightIDandIPCut" || ID=="NUM_LooseRelTkIso_DEN_HighPtIDandIPCut"){
    //==== boundaries
    if(pt<20.) pt = 20.;
    if(pt>=120.) pt = 119.;
    if(eta>=2.4) eta = 2.39;
    if(eta<-2.4) eta = -2.4;
  }

  TH2F *this_hist = map_hist_Muon["ISO_SF_"+ID];
  if(!this_hist){
    if(IgnoreNoHist) return 1.;
    else{
      cout << "[MCCorrection::MuonISO_SF] No "<<"ISO_SF_"+ID<<endl;
      exit(EXIT_FAILURE);
    }
  }

  int this_bin(-999);

  if(DataYear==2016){
    this_bin = this_hist->FindBin(eta,pt);
  }
  else{
    this_bin = this_hist->FindBin(pt,eta);
  }

  value = this_hist->GetBinContent(this_bin);
  error = this_hist->GetBinError(this_bin);

  //cout << "[MCCorrection::MuonISO_SF] value = " << value << endl;

  return value+double(sys)*error;

}

double MCCorrection::MuonTrigger_Eff(TString ID, TString trig, int DataOrMC, double eta, double pt, int sys){

  if(ID=="Default") return 1.;
  if(trig=="Default") return 1.;

  //cout << "[MCCorrection::MuonTrigger_Eff] ID = " << ID << "\t" << "trig = " << trig << endl;
  //cout << "[MCCorrection::MuonTrigger_Eff] DataOrMC = " << DataOrMC << endl;
  //cout << "[MCCorrection::MuonTrigger_Eff] eta = " << eta << ", pt = " << pt << endl;

  double value = 1.;
  double error = 0.;

  eta = fabs(eta);

  //==== 2016
  if(DataYear==2016){
    if(trig=="IsoMu24"){
      if(pt<26.) return 1.; //FIXME
      if(eta>=2.4) eta = 2.39;

      if(pt>500.) pt = 499.;
    }
    else if(trig=="Mu50"){
      if(pt<52.) return 1.; //FIXME
      if(eta>=2.4) eta = 2.39;

      if(pt>800.) pt = 799.;
    }
    else{

    }
  }
  else if(DataYear==2017){
    if(trig=="IsoMu27"){
      //==== FIXME MiniAODPt Pt
      //==== FIXME 28.9918  29.0363
      //==== FIXME This event pass pt>29GeV cut, but MiniAOD pt < 29 GeV
      //==== FIXME So when I return 0., SF goes nan.. let's return 1 for now..
      if(pt<29.) return 1.; //FIXME
      if(eta>=2.4) eta = 2.39;

      if(pt>1200.) pt = 1199.;
    }
    else if(trig=="Mu50"){
      if(pt<52.) return 1.; //FIXME
      if(eta>=2.4) eta = 2.39;

      if(pt>1200.) pt = 1199.;
    }
    else{

    }
  }
  else if(DataYear==2018){
    if(trig=="IsoMu24"){
      if(pt<26.) return 1.; //FIXME
      if(eta>=2.4) eta = 2.39;

      if(pt>1200.) pt = 1199.;
    }
    else if(trig=="Mu50"){
      if(pt<52.) return 1.; //FIXME
      if(eta>=2.4) eta = 2.39;

      if(pt>1200.) pt = 1199.;
    }
    else{

    }
  }
  else{
    cout << "[MCCorrection::MuonTrigger_Eff] Wrong year : " << DataYear << endl;
    exit(EXIT_FAILURE);
  }

  TString histkey = "Trigger_Eff_DATA_"+trig+"_"+ID;
  if(DataOrMC==1) histkey = "Trigger_Eff_MC_"+trig+"_"+ID;
  //cout << "[MCCorrection::MuonTrigger_Eff] histkey = " << histkey << endl;
  TH2F *this_hist = map_hist_Muon[histkey];
  if(!this_hist){
    if(IgnoreNoHist) return 1.;
    else{
      cout << "[MCCorrection::MuonTrigger_Eff] No "<<histkey<<endl;
      exit(EXIT_FAILURE);
    }
  }

  int this_bin = this_hist->FindBin(pt,eta);

  value = this_hist->GetBinContent(this_bin);
  error = this_hist->GetBinError(this_bin);

  //cout << "[MCCorrection::MuonTrigger_Eff] value = " << value << endl;

  return value+double(sys)*error;


}

double MCCorrection::MuonTrigger_SF(TString ID, TString trig, const std::vector<Muon>& muons, int sys){

  if(ID=="Default") return 1.;
  if(trig=="Default") return 1.;

  double value = 1.;

  if(trig=="IsoMu24" || trig=="IsoMu27" || trig=="Mu50"){

    double eff_DATA = 1.;
    double eff_MC = 1.;

    for(unsigned int i=0; i<muons.size(); i++){
      eff_DATA *= ( 1.-MuonTrigger_Eff(ID, trig, 0, muons.at(i).Eta(), muons.at(i).MiniAODPt(), sys) );
      eff_MC   *= ( 1.-MuonTrigger_Eff(ID, trig, 1, muons.at(i).Eta(), muons.at(i).MiniAODPt(), -sys) );
    }

    eff_DATA = 1.-eff_DATA;
    eff_MC = 1.-eff_MC;

    value = eff_DATA/eff_MC;

/*
    if(eff_DATA==0||eff_MC==0){
      cout << "==== Zero Trigger Eff ====" << endl;
      for(unsigned int i=0;i<muons.size();i++){
        cout << muons.at(i).MiniAODPt() << "\t" << muons.at(i).Pt() << endl;
      }
    }
*/

  }

  return value;

}

double MCCorrection::MuonTrigger_SF(TString ID, TString trig, const std::vector<Muon *>& muons, int sys){

  std::vector<Muon> muvec;
  for(unsigned int i=0; i<muons.size(); i++){
    Muon this_muon = *(muons.at(i));
    muvec.push_back( this_muon );
  }

  return MuonTrigger_SF(ID, trig, muvec, sys);

}

double MCCorrection::ElectronID_SF(TString ID, double sceta, double pt, int sys){

  if(ID=="Default") return 1.;

  double value = 1.;
  double error = 0.;

  if(pt<10.) pt = 10.;
  if(pt>=500.) pt = 499.;
  if(sceta>=2.5) sceta = 2.49;
  if(sceta<-2.5) sceta = -2.5;

  if( ID.Contains("HEEP") ){

    //==== https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaRunIIRecommendations#HEEPV7_0
    //==== summary: a simple robust ID designed to be safe for high electrons.
    //==== The Et evolution of this ID must be well described in the MC therefore this ID is designed so its scale factor is flat vs Et.
    //==== As a result the HEEP differs in that it provides just a single number for the barrel and a single number for the endcap.
    //==== * note there almost certainly will have to be a retune for 2018 due to HCAL data/MC disagreements
    //==== * 2018 prompt: expected Dec 2018

    bool IsBarrel = fabs(sceta) < 1.479;
    double this_SF(1.);
    double this_SF_staterr(0.); // absolute value
    double this_SF_systerr(0.); // absolute value
    double this_SF_err(0.);

    if(DataYear==2016){
      this_SF         = (IsBarrel ? 0.971 : 0.983);
      this_SF_staterr = (IsBarrel ? 0.001 : 0.001);

      if(IsBarrel) this_SF_systerr = (pt<90. ? 0.01 : min(1.+(pt-90.)*0.0022,3.)*0.01) * this_SF;
      else         this_SF_systerr = (pt<90. ? 0.01 : min(1.+(pt-90.)*0.0143,4.)*0.01) * this_SF;

      this_SF_err = sqrt(this_SF_staterr*this_SF_staterr+this_SF_systerr*this_SF_systerr);
    }
    else if(DataYear==2017){ 
      this_SF         = (IsBarrel ? 0.967 : 0.973);
      this_SF_staterr = (IsBarrel ? 0.001 : 0.002);
      
      if(IsBarrel) this_SF_systerr = (pt<90. ? 0.01 : min(1.+(pt-90.)*0.0022,3.)*0.01) * this_SF;
      else         this_SF_systerr = (pt<90. ? 0.02 : min(1.+(pt-90.)*0.0143,5.)*0.01) * this_SF;
      
      this_SF_err = sqrt(this_SF_staterr*this_SF_staterr+this_SF_systerr*this_SF_systerr);
    }
    else if(DataYear==2018){
      //==== TODO not yet supported
      //==== copying 2017
      this_SF         = (IsBarrel ? 0.967 : 0.973);
      this_SF_staterr = (IsBarrel ? 0.001 : 0.002);

      if(IsBarrel) this_SF_systerr = (pt<90. ? 0.01 : min(1.+(pt-90.)*0.0022,3.)*0.01) * this_SF;
      else         this_SF_systerr = (pt<90. ? 0.02 : min(1.+(pt-90.)*0.0143,5.)*0.01) * this_SF;

      this_SF_err = sqrt(this_SF_staterr*this_SF_staterr+this_SF_systerr*this_SF_systerr);
    }
    else{
      cout << "[MCCorrection::ElectronID_SF] (Hist) Wrong year "<< DataYear << endl;
      exit(EXIT_FAILURE);
    }

    return this_SF+double(sys)*this_SF_err;

  }
  else{

    TH2F *this_hist = map_hist_Electron["ID_SF_"+ID];
    if(!this_hist){
      if(IgnoreNoHist) return 1.;
      else{
        cout << "[MCCorrection::ElectronID_SF] (Hist) No "<<"ID_SF_"+ID<<endl;
        exit(EXIT_FAILURE);
      }
    }

    int this_bin = this_hist->FindBin(sceta,pt);
    value = this_hist->GetBinContent(this_bin);
    error = this_hist->GetBinError(this_bin);

    return value+double(sys)*error;

  }

}

double MCCorrection::ElectronReco_SF(double sceta, double pt, int sys){

  double value = 1.;
  double error = 0.;

  TString ptrange = "ptgt20";
  if(pt<20.) ptrange = "ptlt20";

  if(pt<10.) pt = 10.;
  if(pt>=500.) pt = 499.;
  if(sceta>=2.5) sceta = 2.49;
  if(sceta<-2.5) sceta = -2.5;

  TH2F *this_hist = map_hist_Electron["RECO_SF_"+ptrange];
  if(!this_hist){
    if(IgnoreNoHist) return 1.;
    else{
      cout << "[MCCorrection::ElectronReco_SF] No "<<"RECO_SF_"+ptrange<<endl;
      exit(EXIT_FAILURE);
    }
  }

  //cout << "[MCCorrection::ElectronReco_SF] " << this_hist->GetBinContent(1,1) << endl;

  int this_bin = this_hist->FindBin(sceta,pt);
  value = this_hist->GetBinContent(this_bin);
  error = this_hist->GetBinError(this_bin);

  return value+double(sys)*error;

}

double MCCorrection::ElectronTrigger_Eff(TString ID, TString trig, int DataOrMC, double sceta, double pt, int sys){

  if(ID=="Default") return 1.;
  if(trig=="Default") return 1.;

  //cout << "[MCCorrection::ElectronTrigger_Eff] ID = " << ID << "\t" << "trig = " << trig << endl;
  //cout << "[MCCorrection::ElectronTrigger_Eff] DataOrMC = " << DataOrMC << endl;
  //cout << "[MCCorrection::ElectronTrigger_Eff] sceta = " << sceta << ", pt = " << pt << endl;

  double value = 1.;
  double error = 0.;

  //==== XXX If you have min pt, apply it here
  if(pt<50.) pt = 50.;
  if(pt>=500.) pt = 499.;

  if(sceta<-2.5) sceta = -2.5;
  if(sceta>=2.5) sceta = 2.49;

  TString histkey = "Trigger_Eff_DATA_"+trig+"_"+ID;
  if(DataOrMC==1) histkey = "Trigger_Eff_MC_"+trig+"_"+ID;
  //cout << "[MCCorrection::ElectronTrigger_Eff] histkey = " << histkey << endl;
  TH2F *this_hist = map_hist_Electron[histkey];
  if(!this_hist){
    if(IgnoreNoHist) return 1.;
    else{
      cout << "[MCCorrection::ElectronTrigger_Eff] No "<<histkey<<endl;
      exit(EXIT_FAILURE);
    }
  }

  int this_bin = this_hist->FindBin(sceta, pt);

  value = this_hist->GetBinContent(this_bin);
  error = this_hist->GetBinError(this_bin);

  //cout << "[MCCorrection::ElectronTrigger_Eff] value = " << value << endl;

  return value+double(sys)*error;


}

double MCCorrection::ElectronTrigger_SF(TString ID, TString trig, const std::vector<Electron>& electrons, int sys){

  if(ID=="Default") return 1.;
  if(trig=="Default") return 1.;

  double value = 1.;

  if(trig=="Ele27_WPTight_Gsf" || trig=="Ele35_WPTight_Gsf" || trig=="Ele32_WPTight_Gsf"){

    double eff_DATA = 1.;
    double eff_MC = 1.;

    for(unsigned int i=0; i<electrons.size(); i++){
      eff_DATA *= ( 1.-ElectronTrigger_Eff(ID, trig, 0, electrons.at(i).scEta(), electrons.at(i).Pt(), sys) );
      eff_MC   *= ( 1.-ElectronTrigger_Eff(ID, trig, 1, electrons.at(i).scEta(), electrons.at(i).Pt(), -sys) );
    }

    eff_DATA = 1.-eff_DATA;
    eff_MC = 1.-eff_MC;

    value = eff_DATA/eff_MC;

/*
    if(eff_DATA==0||eff_MC==0){
      cout << "==== Zero Trigger Eff ====" << endl;
      for(unsigned int i=0;i<electrons.size();i++){
        cout << electrons.at(i).Pt() << "\t" << electrons.at(i).Pt() << endl;
      }
    }
*/

  }

  return value;

}

double MCCorrection::ElectronTrigger_SF(TString ID, TString trig, const std::vector<Electron *>& electrons, int sys){

  std::vector<Electron> muvec;
  for(unsigned int i=0; i<electrons.size(); i++){
    Electron this_electron = *(electrons.at(i));
    muvec.push_back( this_electron );
  }

  return ElectronTrigger_SF(ID, trig, muvec, sys);

}

double MCCorrection::GetPrefireWeight(const std::vector<Photon>& photons, const std::vector<Jet>& jets, int sys){

  double photon_weight = 1.;
  double jet_weight = 1.;
  
  TH2F *photon_hist = map_hist_prefire["Photon_prefire"];
  TH2F *jet_hist = map_hist_prefire["Jet_prefire"];

  
  for(unsigned int i_pho = 0; i_pho < photons.size(); i_pho++){
    Photon current_photon = photons.at(i_pho);
    double eta = current_photon.scEta();
    double pt = current_photon.Pt();
    
    int this_bin = photon_hist->FindBin(eta, pt);

    double this_eff = photon_hist->GetBinContent(this_bin);
    double this_efferr = photon_hist->GetBinError(this_bin);

    double current_weight = 1. - (this_eff + (double)sys * this_efferr );;
    photon_weight = photon_weight * current_weight;
  }
  
  for(unsigned int i_jet = 0; i_jet < jets.size(); i_jet++){
    Jet current_jet = jets.at(i_jet);
    double eta = current_jet.Eta();
    double pt = current_jet.Pt();
    
    int this_bin = jet_hist->FindBin(eta, pt);

    double this_eff = jet_hist->GetBinContent(this_bin);
    double this_efferr = jet_hist->GetBinError(this_bin);

    double current_weight = 1.- (this_eff + (double)sys * this_efferr);
    jet_weight = jet_weight * current_weight;
  }

  return photon_weight * jet_weight;;

}


double MCCorrection::GetPileUpWeightBySampleName(int N_pileup, int syst){
  
  int this_bin = N_pileup+1;
  if(N_pileup >= 100) this_bin=100;

  TString this_histname = MCSample;
  if(syst == 0){
    this_histname += "_central_pileup";
  }
  else if(syst == -1){
    this_histname += "_sig_down_pileup";
  }
  else if(syst == 1){
    this_histname += "_sig_up_pileup";
  }
  else{
    cout << "[MCCorrection::GetPileUpWeightBySampleName] syst should be 0, -1, or +1" << endl;
    exit(EXIT_FAILURE);
  }

  TH1D *this_hist = map_hist_pileup[this_histname];
  if(!this_hist){
    cout << "[MCCorrection::GetPileUpWeightBySampleName] No " << this_histname << endl;
    exit(EXIT_FAILURE);
  }

  return this_hist->GetBinContent(this_bin);

}

double MCCorrection::GetPileUpWeight(int N_pileup, int syst){

  int this_bin = N_pileup+1;
  if(N_pileup >= 100) this_bin=100;

  TString this_histname = "MC_" + TString::Itoa(DataYear,10);
  if(syst == 0){
    this_histname += "_central_pileup";
  }
  else if(syst == -1){
    this_histname += "_sig_down_pileup";
  }
  else if(syst == 1){
    this_histname += "_sig_up_pileup";
  }
  else{
    cout << "[MCCorrection::GetPileUpWeightBySampleName] syst should be 0, -1, or +1" << endl;
    exit(EXIT_FAILURE);
  }

  TH1D *this_hist = map_hist_pileup[this_histname];
  if(!this_hist){
    cout << "[MCCorrection::GetPileUpWeightBySampleName] No " << this_histname << endl;
    exit(EXIT_FAILURE);
  }

  return this_hist->GetBinContent(this_bin);

}

double MCCorrection::GetTopPtReweight(const std::vector<Gen>& gens){
  //==== ref: https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting2017
  //==== Only top quarks in SM ttbar events must be reweighted, 
  //==== not single tops or tops from BSM production mechanisms.
  if(!MCSample.Contains("TT") || !MCSample.Contains("powheg")){
    return 1.;
  }
  //==== initialize with large number
  double toppt1=10000, toppt2=10000;
  bool found_top = false, found_atop = false;

  for(vector<Gen>::const_iterator genit=gens.begin(); genit!=gens.end(); genit++){
    
    if(genit->Status() == 22){
      if(genit->PID() == 6){
        toppt1= genit->Pt();
        found_top = true;
      }
      else if(genit->PID() == -6){
        toppt2= genit->Pt();
        found_atop = true;
      }
    }
    //==== after we found top pair, break the loop
    if(found_top && found_atop) break;
  }
  double pt_reweight = 1.;
  //==== if top pair is not found, return 1.
  //==== the measurement covers only the range pt(top)<=800GeV, otherwise, return 1.
  if(toppt1<=800 && toppt2 <=800){
    pt_reweight*=exp(0.0615-0.0005*toppt1);
    pt_reweight*=exp(0.0615-0.0005*toppt2);
    pt_reweight = sqrt(pt_reweight);
  }
  return pt_reweight;
}

