#include "Event.h"

ClassImp(Event)

Event::Event(){
  j_MCweight = 1.;
  j_nPV = -999;
  j_DataYear = -999;
}

Event::~Event(){

}

void Event::SetMCweight(double genw){
  if(genw>0) j_MCweight=1.;
  else j_MCweight=-1.;
}

void Event::SetnPV(double n){
  j_nPV = n;
}

void Event::SetTrigger(vector<string> HLT_TriggerName){
  j_HLT_TriggerName = HLT_TriggerName;
}

bool Event::PassTrigger(TString trig){
  std::vector<TString> tmp_vec;
  tmp_vec.push_back(trig);
  return PassTrigger(tmp_vec);
}

bool Event::PassTrigger(std::vector<TString> trigs){

  for(unsigned int i=0; i<trigs.size(); i++){
    TString this_check_trig = trigs.at(i);

    //cout << this_check_trig << endl;
    for(unsigned int j=0; j<j_HLT_TriggerName.size(); j++){

      //cout << j_HLT_TriggerName.at(j) << endl;
      if( j_HLT_TriggerName.at(j).find(this_check_trig)!=std::string::npos ){
        return true;
      }
    }

  }
  return false;
}

double Event::GetTriggerLumi(TString trig){

  if(j_DataYear==2016){
    if(trig=="Full"){
      return 35918.219;
    }
    else if(trig=="HLT_Mu20_v"){
      return 139.671;
    }
    else if(trig=="HLT_Mu27_v"){
      return 250.508;
    }
    else if(trig=="HLT_Mu50_v"){
      return 35918.219;
    }
    else if(trig=="HLT_Mu8_v"){
      return 3.937;
    }
    else if(trig=="HLT_Mu17_v"){
      return 282.781;
    }
    else if(trig=="HLT_Mu3_PFJet40_v"){
      return 7.408;
    }
    else if(trig=="HLT_IsoMu24_v"){
      return 35918.219;
    }
    else if(trig=="HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v"){
      return 6.988;
    }
    else if(trig=="HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v"){
      return 14.851;
    }
    else if(trig=="HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v"){
      return 62.808;
    }
    else if(trig=="HLT_DoublePhoton60_v"){
      return 35918.219;
    }
    else if(trig=="HLT_Photon25_v"){
      return 0.000;
    }
    else if(trig=="HLT_Photon33_v"){
      return 0.000;
    }
    else if(trig=="HLT_Photon50_v"){
      return 31.037;
    }
    else if(trig=="HLT_Photon75_v"){
      return 134.568;
    }
    else if(trig=="HLT_Photon90_v"){
      return 264.074;
    }
    else if(trig=="HLT_Photon120_v"){
      return 537.005;
    }
    else if(trig=="HLT_Photon150_v"){
      return 0.000;
    }
    else if(trig=="HLT_Photon175_v"){
      return 35918.219;
    }
    else if(trig=="HLT_Photon200_v"){
      return 0.000;
    }
  }
  else if(j_DataYear==2017){

    if(trig=="Full"){
      return 41527.540;
    }
    else if(trig=="HLT_Mu3_PFJet40_v"){
      return 4.612;
    }
    else if(trig=="HLT_Mu8_v"){
      //==== 180829 Norm Checked : sf = 1.41036
      //==== 181206 Norm Checked : sf = 1.33461
      //return 2.605;;
      return 2.605*1.33461;
    }
    else if(trig=="HLT_Mu17_v"){
      return 70.039;
    }
    else if(trig=="HLT_Mu19_v"){
      return 46.205;
    }
    else if(trig=="HLT_Mu20_v"){
      return 574.102;
    }
    else if(trig=="HLT_Mu27_v"){
      return 184.945;
    }
    else if(trig=="HLT_Mu50_v"){
      return 41527.540;
    }
    else if(trig=="HLT_Mu55_v"){
      return 27122.162;
    }
    else if(trig=="HLT_DoublePhoton70_v"){
      return 41527.540;
    }
    else if(trig=="HLT_DoublePhoton85_v"){
      return 41527.540;
    }
    else if(trig=="HLT_Photon25_v"){
      return 1.339;
    }
    else if(trig=="HLT_Photon33_v"){
      return 14.264;
    }
    else if(trig=="HLT_Photon50_v"){
      return 22.189;
    }
    else if(trig=="HLT_Photon75_v"){
      return 101.436;
    }
    else if(trig=="HLT_Photon90_v"){
      return 122.189;
    }
    else if(trig=="HLT_Photon120_v"){
      return 490.296;
    }
    else if(trig=="HLT_Photon150_v"){
      return 689.904;
    }
    else if(trig=="HLT_Photon175_v"){
      return 1245.822;
    }
    else if(trig=="HLT_Photon200_v"){
      return 41527.540;
    }
    else if(trig=="HLT_IsoMu24_v"){
      return 38046.662;
    }
    else if(trig=="HLT_IsoMu27_v"){
      return 41527.540;
    }
    else if(trig=="HLT_IsoMu30_v"){
      return 27122.162;
    }
    else if(trig=="HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v"){
      return 3.973;
    }
    else if(trig=="HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v"){
      return 27.699;
    }
    else if(trig=="HLT_Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30_v"){
      return 0.001;
    }
    else if(trig=="HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v"){
      return 43.468;
    }
    else if(trig=="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v"){
      return 7449.301;
    }
    else if(trig=="HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_v"){
      return 2613.503;
    }
    else if(trig=="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v"){
      return 41527.540;
    }
    else if(trig=="HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8_v"){
      return 27122.162;
    }
    else if(trig=="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"){
      return 36733.558;
    }
    else if(trig=="HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8_v"){
      return 27122.162;
    }
    else{
      cout << "[Event::GetTriggerLumi] Wrong trigger (DataYear=" << j_DataYear << ", trig=" << trig << ")" << endl;
      exit(EXIT_FAILURE);
      return -1.;
    }

  }
  else if(j_DataYear==2018){

    if(trig=="Full"){
      return 59735.969;
    }
    else{
      cout << "[Event::GetTriggerLumi] Wrong trigger (DataYear=" << j_DataYear << ", trig=" << trig << ")" << endl;
      exit(EXIT_FAILURE);
      return -1.;
    }

  }
  else{
    cout << "[Event::GetTriggerLumi] Wrong trigger (DataYear=" << j_DataYear << ", trig=" << trig << ")" << endl;
    exit(EXIT_FAILURE);
    return -1.;
  }

  cout << "[Event::GetTriggerLumi] Wrong trigger (DataYear=" << j_DataYear << ", trig=" << trig << ")" << endl;
  exit(EXIT_FAILURE);
  return -1.;

}

bool Event::IsPDForTrigger(TString trig, TString PD){

  if(j_DataYear==2017){
    if(trig=="HLT_Mu3_PFJet40_v"){
      if(PD=="SingleMuon") return true;
      else return false;
    }
    else if(trig=="HLT_Mu8_v"){
      if(PD=="DoubleMuon") return true;
      else return false;
    }
    else if(trig=="HLT_Mu17_v"){
      if(PD=="DoubleMuon") return true;
      else return false;
    }
    else if(trig=="HLT_Mu20_v"){
      if(PD=="SingleMuon") return true;
      else return false;
    }
    else if(trig=="HLT_Mu27_v"){
      if(PD=="SingleMuon") return true;
      else return false;
    }
    else if(trig=="HLT_Mu50_v"){
      if(PD=="SingleMuon") return true;
      else return false;
    }
    else if(trig=="HLT_IsoMu24_v"){
      if(PD=="SingleMuon") return true;
      else return false;
    }
    else if(trig=="HLT_IsoMu27_v"){
      if(PD=="SingleMuon") return true;
      else return false;
    }
    else if(trig=="HLT_IsoMu30_v"){
      if(PD=="SingleMuon") return true;
      else return false;
    }
    else if(trig=="HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v"){
      if(PD=="SingleElectron") return true;
      else return false;
    }
    else if(trig=="HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v"){
      if(PD=="SingleElectron") return true;
      else return false;
    }
    else if(trig=="HLT_Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30_v"){
      if(PD=="SingleElectron") return true;
      else return false;
    }
    else if(trig=="HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v"){
      if(PD=="SingleElectron") return true;
      else return false;
    }
    else if(trig=="HLT_Photon25_v"){
      if(PD=="SinglePhoton") return true;
      else return false;
    }
    else if(trig=="HLT_Photon33_v"){
      if(PD=="SinglePhoton") return true;
      else return false;
    }
    else if(trig=="HLT_Photon50_v"){
      if(PD=="SinglePhoton") return true;
      else return false;
    }
    else if(trig=="HLT_Photon75_v"){
      if(PD=="SinglePhoton") return true;
      else return false;
    }
    else if(trig=="HLT_Photon90_v"){
      if(PD=="SinglePhoton") return true;
      else return false;
    }
    else if(trig=="HLT_Photon120_v"){
      if(PD=="SinglePhoton") return true;
      else return false;
    }
    else if(trig=="HLT_Photon150_v"){
      if(PD=="SinglePhoton") return true;
      else return false;
    }
    else if(trig=="HLT_Photon175_v"){
      if(PD=="SinglePhoton") return true;
      else return false;
    }
    else if(trig=="HLT_Photon200_v"){
      if(PD=="SinglePhoton") return true;
      else return false;
    }
    else{
      cout << "[Event::IsPDForTrigger] trig = " << trig << endl;
      cout << "[Event::IsPDForTrigger] PD = " << PD << endl;
      exit(EXIT_FAILURE);
      return false;
    }
  }

  cout << "[Event::IsPDForTrigger] trig = " << trig << endl;
  cout << "[Event::IsPDForTrigger] PD = " << PD << endl;
  exit(EXIT_FAILURE);
  return false;

}

void Event::SetMET(double pt, double phi){
  j_METVector.SetPtEtaPhiE(pt, 0, phi, pt);
}

void Event::SetDataYear(int y){
  j_DataYear = y;
}
