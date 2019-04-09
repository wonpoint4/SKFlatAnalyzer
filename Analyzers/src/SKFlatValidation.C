#include "SKFlatValidation.h"

void SKFlatValidation::initializeAnalyzer(){


  //==== Triggers

  Triggers_POG_Electron.clear();
  Triggers_POG_Muon.clear();
  Triggers_POGHighPt_Electron.clear();
  Triggers_POGHighPt_Muon.clear();

  if(DataYear==2016){

    Triggers_POG_Electron = {
      "HLT_Ele27_WPTight_Gsf_v",
    };
    Triggers_POG_Muon = {
      "HLT_IsoMu24_v",
      "HLT_IsoTkMu24_v",
    };
    TriggerNameForSF_POG_Electron = "Ele27_WPTight_Gsf";
    TriggerNameForSF_POG_Muon = "IsoMu24";
    TriggerSafePt_POG_Electron = 30.;
    TriggerSafePt_POG_Muon = 26.;

    Triggers_POGHighPt_Electron = {
      "HLT_Ele27_WPTight_Gsf_v",
    };
    Triggers_POGHighPt_Muon = {
      "HLT_Mu50",
      "HLT_TkMu50_v",
    };
    TriggerNameForSF_POGHighPt_Electron = "Ele27_WPTight_Gsf";
    TriggerNameForSF_POGHighPt_Muon = "Mu50";
    TriggerSafePt_POGHighPt_Electron = 30.;
    TriggerSafePt_POGHighPt_Muon = 52.;
    DoublePhotonSafePtCut = 65.;
  }
  else if(DataYear==2017){

    Triggers_POG_Electron = {
      "HLT_Ele35_WPTight_Gsf_v",
    };
    Triggers_POG_Muon = {
      "HLT_IsoMu27_v",
    };
    TriggerNameForSF_POG_Electron = "Ele35_WPTight_Gsf";
    TriggerNameForSF_POG_Muon = "IsoMu27";
    TriggerSafePt_POG_Electron = 38.;
    TriggerSafePt_POG_Muon = 29.;

    Triggers_POGHighPt_Electron = {
      "HLT_Ele35_WPTight_Gsf_v",
    };
    Triggers_POGHighPt_Muon = {
      "HLT_Mu50_v",
    };
    TriggerNameForSF_POGHighPt_Electron = "Ele35_WPTight_Gsf";
    TriggerNameForSF_POGHighPt_Muon = "Mu50";
    TriggerSafePt_POGHighPt_Electron = 38.;
    TriggerSafePt_POGHighPt_Muon = 52.;
    DoublePhotonSafePtCut = 75.;
  }
  else if(DataYear==2018){

    Triggers_POG_Electron = {
      "HLT_Ele32_WPTight_Gsf_v",
    };
    Triggers_POG_Muon = {
      "HLT_IsoMu24_v",
    };
    TriggerNameForSF_POG_Electron = "Ele32_WPTight_Gsf";
    TriggerNameForSF_POG_Muon = "IsoMu24";
    TriggerSafePt_POG_Electron = 35.;
    TriggerSafePt_POG_Muon = 26.;

    Triggers_POGHighPt_Electron = {
      "HLT_Ele32_WPTight_Gsf_v",
    };
    Triggers_POGHighPt_Muon = {
      "HLT_Mu50_v",
      "HLT_OldMu100_v",
      "HLT_TkMu100_v",
    };
    TriggerNameForSF_POGHighPt_Electron = "Ele35_WPTight_Gsf";
    TriggerNameForSF_POGHighPt_Muon = "Mu50";
    TriggerNameForSF_POG_Electron = "Default";
    TriggerNameForSF_POG_Muon = "Default";
    TriggerSafePt_POGHighPt_Electron = 35.;
    TriggerSafePt_POGHighPt_Muon = 52.;
    DoublePhotonSafePtCut = 75.;
  }

  else{

  }

  //==== B-Tagging
  //==== add taggers and WP that you want to use in analysis
  std::vector<Jet::Tagger> vtaggers;
  vtaggers.push_back(Jet::DeepCSV);

  std::vector<Jet::WP> v_wps;
  v_wps.push_back(Jet::Medium);

  //=== list of taggers, WP, setup systematics, use period SFs
  SetupBTagger(vtaggers,v_wps, true, true);

}

void SKFlatValidation::executeEvent(){

  //==== Prefire reweight

  weight_Prefire = GetPrefireWeight(0);

  //==== AnalyzerParameter

  AnalyzerParameter param;

  //==== POG IDs

  param.Name = "POG";

  param.Electron_Tight_ID = "passMediumID";
  param.Electron_ID_SF_Key = "passMediumID";

  param.Muon_Tight_ID = "POGTightWithTightIso";
  param.Muon_ID_SF_Key = "NUM_TightID_DEN_genTracks";
  param.Muon_ISO_SF_Key = "NUM_TightRelIso_DEN_TightIDandIPCut";
  param.Muon_Trigger_SF_Key = "POGTight";

  param.Jet_ID = "HN";

  executeEventFromParameter(param);

  //==== POG HighPt

  param.Clear();

  param.Name = "POGHighPt";

  param.Electron_Tight_ID = "passHEEPID";
  param.Electron_ID_SF_Key = "HEEP";

  param.Muon_Tight_ID = "POGHighPtWithLooseTrkIso";
  param.Muon_ID_SF_Key = "NUM_HighPtID_DEN_genTracks";
  param.Muon_ISO_SF_Key = "NUM_LooseRelTkIso_DEN_HighPtIDandIPCut";
  param.Muon_Trigger_SF_Key = "POGHighPtLooseTrkIso";
  param.Muon_UseTuneP = true;

  param.Jet_ID = "HN";

  executeEventFromParameter(param);

}

void SKFlatValidation::executeEventFromParameter(AnalyzerParameter param){

  if(!PassMETFilter()) return;

  Event ev = GetEvent();

  bool PassSingleElectron, PassSingleMuon;
  TString TriggerNameForSF_Electron, TriggerNameForSF_Muon;
  int TriggerSafePt_Electron, TriggerSafePt_Muon;
  double MinLeptonPt = 20.;
  if(param.Name=="POG"){
    TriggerNameForSF_Electron = TriggerNameForSF_POG_Electron;
    TriggerNameForSF_Muon = TriggerNameForSF_POG_Muon;
    TriggerSafePt_Electron = TriggerSafePt_POG_Electron;
    TriggerSafePt_Muon = TriggerSafePt_POG_Muon;
    PassSingleElectron = ev.PassTrigger(Triggers_POG_Electron);
    PassSingleMuon = ev.PassTrigger(Triggers_POG_Muon);
  }
  else if(param.Name=="POGHighPt"){
    TriggerNameForSF_Electron = TriggerNameForSF_POGHighPt_Electron;
    TriggerNameForSF_Muon = TriggerNameForSF_POGHighPt_Muon;
    TriggerSafePt_Electron = TriggerSafePt_POGHighPt_Electron;
    TriggerSafePt_Muon = TriggerSafePt_POGHighPt_Muon;
    PassSingleElectron = ev.PassTrigger(Triggers_POGHighPt_Electron);
    PassSingleMuon = ev.PassTrigger(Triggers_POGHighPt_Muon);
  }
  else{
    cout << "[SKFlatValidation::executeEventFromParameter] Wrong param.Name = " << param.Name << endl;
    exit(EXIT_FAILURE);
  }

  std::vector<Muon> muons = GetMuons(param.Muon_Tight_ID, MinLeptonPt, 2.4);
  if(param.Muon_UseTuneP) muons = UseTunePMuon(muons);
  std::vector<Electron> electrons = GetElectrons(param.Electron_Tight_ID, MinLeptonPt, 2.5);

  std::vector<Jet> myjets = JetsVetoLeptonInside( GetJets("tight", 30., 2.4), electrons, muons);
  int NBJets=0;
  double HT=0;
  for(unsigned int i=0; i<myjets.size(); i++){
    Jet this_jet = myjets.at(i);
    HT += this_jet.Pt();
    if(IsBTagged(this_jet, Jet::DeepCSV, Jet::Medium,true,0)) NBJets++;
  }

  //==== Based on which trigger is fired
  std::vector< TString > Suffixs = {
    "SingleMuon",
    "SingleElectron",
  };
  std::vector< bool > PassTriggers = {
    PassSingleMuon && (muons.size()>=1) && (electrons.size()==0),
    PassSingleElectron && (electrons.size()>=1) && (muons.size()==0),
  };

  for(unsigned int i=0; i<Suffixs.size(); i++){

    TString Suffix = Suffixs.at(i);
    if( !PassTriggers.at(i) ) continue;

    if(Suffix.Contains("SingleMuon")){
      if( muons.at(0).Pt() < TriggerSafePt_Muon ) continue;
    }
    else if(Suffix.Contains("SingleElectron")){
      if( electrons.at(0).Pt() < TriggerSafePt_Electron ) continue;
    }
    else{

    }

    if(this->IsDATA){
      if(this->DataStream == "SingleMuon"){
        if(!( Suffix.Contains("SingleMuon") )) continue;
      }
      else if(this->DataStream == "SingleElectron"){
        if(!( Suffix.Contains("SingleElectron") )) continue;
      }
      else{

      }
    }
    
    std::vector<Lepton *> leps;
    if(Suffix.Contains("SingleMuon")){
      leps = MakeLeptonPointerVector(muons);
    }
    else if(Suffix.Contains("SingleElectron")){
      leps = MakeLeptonPointerVector(electrons);
    }
    else{

    }

    Particle METv = ev.GetMETVector();

    int n_lepton = leps.size();
    //==== DiLepton variables
    bool IsOS = false;
    Particle Z;
    if(n_lepton==2){
      IsOS = (leps[0]->Charge() != leps[1]->Charge());
      Z = (*leps[0])+(*leps[1]);
    }
    //==== SingleLepton variables
    double this_MT=-999.;
    if(n_lepton==1){
      this_MT = MT( METv, (*leps[0]) );
    }

    double weight = 1.;
    if(!IsDATA){
      //cout << "weight_norm_1invpb = " << weight_norm_1invpb << endl;
      //cout << "GetTriggerLumi = " << ev.GetTriggerLumi("Full") << endl;
      //cout << "MCweight = " << ev.MCweight() << endl;
      //cout << "weight_Prefire = " << weight_Prefire << endl;
      weight *= weight_norm_1invpb*ev.GetTriggerLumi("Full")*ev.MCweight()*weight_Prefire;

      mcCorr->IgnoreNoHist = param.MCCorrrectionIgnoreNoHist;

      //==== FIXME add third lepton veto later
      if(Suffix.Contains("SingleMuon")){
        for(unsigned int i=0; i<muons.size(); i++){

          double this_pt = muons.at(i).MiniAODPt();
          double this_eta = muons.at(i).Eta();

          double this_idsf  = mcCorr->MuonID_SF (param.Muon_ID_SF_Key,  this_eta, this_pt);
          double this_isosf = mcCorr->MuonISO_SF(param.Muon_ISO_SF_Key, this_eta, this_pt);
          double this_trigsf = mcCorr->MuonTrigger_SF(param.Muon_Trigger_SF_Key, TriggerNameForSF_Muon, muons);

          weight *= this_idsf*this_isosf*this_trigsf;

        }
      }
      if( Suffix.Contains("DiElectron") || Suffix.Contains("SingleElectron") || Suffix.Contains("DiPhoton") ){
        for(unsigned int i=0; i<electrons.size(); i++){
          double this_recosf = mcCorr->ElectronReco_SF(electrons.at(i).scEta(),electrons.at(i).Pt());
          double this_idsf = mcCorr->ElectronID_SF(param.Electron_ID_SF_Key, electrons.at(i).scEta(), electrons.at(i).Pt());
          weight *= this_recosf*this_idsf;
        }
      }

    }

    std::map<TString, bool> map_bool_To_Region;
    if(n_lepton==2){

      bool BaseDiLepSelection = (Z.M() > 15.);
      bool HigherDiLeptonPtCut = (leps[0]->Pt() > DoublePhotonSafePtCut) && (leps[1]->Pt() > DoublePhotonSafePtCut);

      if(BaseDiLepSelection){

        if( IsOS ){
          //==== generic two OS lepton
          map_bool_To_Region["OS"] = true;
          //==== OnZ event
          map_bool_To_Region["OnZ_OS"] = IsOnZ(Z.M(), 15.);
          //==== High ZMass event
          map_bool_To_Region["HigherDiLeptonPtCut_OS"] = HigherDiLeptonPtCut;
          //==== With B-jet, MET > 30 for dilepton ttbar
          map_bool_To_Region["WithBJet_METgt30_OS"] = (NBJets>0) && (METv.Pt()>30.);
        }
        else{
          //==== generic two SS lepton
          map_bool_To_Region["SS"] = true;
          //==== OnZ event
          map_bool_To_Region["OnZ_SS"] = IsOnZ(Z.M(), 15.);
          //==== High ZMass event
          map_bool_To_Region["HigherDiLeptonPtCut_SS"] = HigherDiLeptonPtCut;
          //==== With B-jet, MET > 30 for dilepton ttbar
          map_bool_To_Region["WithBJet_METgt30_SS"] = (NBJets>0) && (METv.Pt()>30.);
        }

      }

    }
    if(n_lepton==1){
      map_bool_To_Region["W_CR"] = ( METv.Pt() > 30. ) && ( this_MT > 30. );
    }

    int nPUForReweight = nPV;
    if(!IsDATA) nPUForReweight = nPileUp;

    for(std::map<TString, bool>::iterator it_map = map_bool_To_Region.begin(); it_map != map_bool_To_Region.end(); it_map++){

      TString this_region = it_map->first;
      this_region = (param.Name)+"_"+Suffix+"_"+this_region;

      if(it_map->second){

        JSFillHist(this_region, "NEvent_"+this_region, 0., weight, 1, 0., 1.);

        JSFillHist(this_region, "nPileUp_"+this_region, nPileUp, weight, 200., 0., 200.);
        JSFillHist(this_region, "nPV_"+this_region, nPV, weight, 200., 0., 200.);
        JSFillHist(this_region, "nPUForReweight_"+this_region, nPUForReweight, weight, 200., 0., 200.);

        JSFillHist(this_region, "MET_"+this_region, METv.Pt(), weight, 500, 0., 500.);
        JSFillHist(this_region, "METphi_"+this_region, METv.Phi(), weight, 60, -3., 3.);

        JSFillHist(this_region, "HT_"+this_region, HT, weight, 1000, 0., 1000.);

        JSFillHist(this_region, "Jet_Size_"+this_region, myjets.size(), weight, 10, 0., 10.);
        JSFillHist(this_region, "NBJets_"+this_region, NBJets, weight, 10, 0., 10.);

        FillLeptonPlots(leps, this_region, weight);

        if(n_lepton==1){
          JSFillHist(this_region, "MT_"+this_region, this_MT, weight, 500, 0., 500.);
        }

        if(n_lepton>=2){

          JSFillHist(this_region, "ZCand_Mass_"+this_region, Z.M(), weight, 7000, 0., 7000.);
          JSFillHist(this_region, "ZCand_Pt_"+this_region, Z.Pt(), weight, 500, 0., 500.);
          JSFillHist(this_region, "ZCand_Eta_"+this_region, Z.Eta(), weight, 60, -3., 3.);

        }



      } // END if(pass Region)

    } // END Loop mapRegionBool

  } // END Loop Suffixs

}
