#include "controlPlotDY.h"

void controlPlotDY::initializeAnalyzer(){

  MuonIDISOSFs = {"NoIDISOSF", "IDISOSF"};
  MuonTrigSFs = {"NoTrigSF", "TrigSF"};
  //For Data, it doesn't need Scale Factors.
  if(IsDATA){
    MuonIDISOSFs = {"NoIDISOSF"}; MuonTrigSFs = {"NoTrigSF"};
  }

  IsDYsamples = false;
  if(!IsDATA && (MCSample == "DYJets" || MCSample == "DYJets_MG")) IsDYsamples = true;

  IsSingleMuonData = false;
  IsDoubleMuonData = false;
  Samples = {"SingleMuon", "DoubleMuon"};

  if(IsDATA){
    if(DataStream=="SingleMuon"){
      IsSingleMuonData= true;
      Samples = {"SingleMuon"};
    }
    else{
      IsDoubleMuonData = true;
      Samples = {"DoubleMuon"};
    }
  }

  SetupZPtWeight();
}

void controlPlotDY::executeEvent(){

  AllMuons = GetAllMuons();
  gens = GetGens();
  IsFromTau = "";

  //HyonSan's ZpT correction part
    zptcor=1.;
  if(IsDYsamples){
    int parton0,parton1,hardl0,hardl1,l0,l1;
    vector<int> photons;
    GetGenIndex(gens,parton0,parton1,hardl0,hardl1,l0,l1,photons);
    int hardj0=0,hardj1=0,hardj2=0;
    Gen genhardj0,genhardj1,genhardj2;
    for(unsigned int i=parton1+1;i<gens.size();i++){
      if(gens[i].isHardProcess()){
        if(gens[i].PID()==21||abs(gens[i].PID())<7){
	  if(!hardj0){
	    hardj0=i;
	    genhardj0=gens[i];
	  }else if(!hardj1){
	    hardj1=i;
	    genhardj1=gens[i];
	  }else if(!hardj2){
	    hardj2=i;
	    genhardj2=gens[i];
	    break;
	  }
        }
      }
    }
    if(abs(gens[hardl0].PID())!=15){
      Gen genparton0=gens[parton0],genparton1=gens[parton1],genhardl0=gens[hardl0],genhardl1=gens[hardl1],genl0=gens[l0],genl1=gens[l1],genphotons;
      genparton0.SetPxPyPzE(0,0,genWeight_X1*6500.,genWeight_X1*6500.);
      genparton1.SetPxPyPzE(0,0,-genWeight_X2*6500.,genWeight_X2*6500.);
      for(unsigned int i=0;i<photons.size();i++) genphotons+=gens[photons[i]];
      TLorentzVector genZ=(genl0+genl1+genphotons);
      TLorentzVector gendilepton=genl0+genl1;
      TLorentzVector genharddilepton=genhardl0+genhardl1;
      //FIXME need to get Zptcor for 2018
      //if(DataYear!=2018) 
      zptcor*=GetZPtWeight(genZ.Pt(),genZ.Rapidity(),abs(genhardl0.PID())==13?Lepton::Flavour::MUON:Lepton::Flavour::ELECTRON);
      /*
      TString slepton=abs(genhardl0.PID())==13?"muon":"electron";
      if(genparton0.PID()==21){
        if(genparton1.PID()==21) hardprefix="gg_";
        else if(genparton1.PID()>0) hardprefix="gq_";
        else if(genparton1.PID()<0) hardprefix="gqbar_";
      }else if(genparton0.PID()>0){
        if(genparton1.PID()==21) hardprefix="qg_";
        else if(genparton1.PID()>0) hardprefix="qq_";
        else if(genparton1.PID()<0) hardprefix="qqbar_";
      }else if(genparton0.PID()<0){
        if(genparton1.PID()==21) hardprefix="qbarg_";
        else if(genparton1.PID()>0) hardprefix="qbarq_";
        else if(genparton1.PID()<0) hardprefix="qbarqbar_";
      }
      */
    }
    else IsFromTau = "tau/";
  }

  weight_Prefire = GetPrefireWeight(0);

  AnalyzerParameter param; param.Clear();
  param.Muon_Tight_ID = "POGMediumWithLooseTrkIso";

  for(unsigned int it_Sample=0; it_Sample < Samples.size(); it_Sample++){

    TString Sample = Samples.at(it_Sample);
    IsSingleMuon = false;  IsDoubleMuon = false;
    if(Sample == "SingleMuon") IsSingleMuon = true; // Go to make SingleMuon Control plots (Data, MC)
    else IsDoubleMuon = true;                       // Go to make DoubleMuon Control plots (Data, MC)

    for(unsigned int it_MuonIDISOSF=0; it_MuonIDISOSF < MuonIDISOSFs.size(); it_MuonIDISOSF++){

      TString MuonIDISOSFKey = MuonIDISOSFs.at(it_MuonIDISOSF);

      if(MuonIDISOSFKey == "IDISOSF"){
        param.Muon_ID_SF_Key = "Default";
        param.Muon_ISO_SF_Key = "NUM_MediumID_RelTrkIso010_DEN_TrackerMu";
      }
      else{
        param.Muon_ID_SF_Key = "Default";
        param.Muon_ISO_SF_Key = "Default";
      }
   
      for(unsigned int it_MuonTrigSF=0; it_MuonTrigSF < MuonTrigSFs.size(); it_MuonTrigSF++){

        TString MuonTrigSFKey = MuonTrigSFs.at(it_MuonTrigSF);
        param.Name = Sample+"/"+MuonIDISOSFKey + "_" + MuonTrigSFKey;

        if(MuonTrigSFKey == "TrigSF") param.Muon_Trigger_SF_Key = "MediumID_RelTrkIso010";
        else param.Muon_Trigger_SF_Key = "Default";

        if(!IsDATA && DataYear == 2016){
          vector<TString> periods = {"2016BF", "2016GH"};
          TString a = "", b = "", c = "", d = ""; // buffer
          a = param.Name; b = param.Muon_ID_SF_Key; c = param.Muon_ISO_SF_Key; d = param.Muon_Trigger_SF_Key;
      
          for(unsigned int i=0; i < periods.size(); i++){
	    param.Name = periods.at(i) + "_" + a;
	    param.Muon_ID_SF_Key = b + "_" + periods.at(i);
	    param.Muon_ISO_SF_Key = c + "_" + periods.at(i);
	    param.Muon_Trigger_SF_Key = d + "_" + periods.at(i);

            executeEventFromParameter(param);
            param.Name = a;    param.Muon_ID_SF_Key = b;    param.Muon_ISO_SF_Key = c;    param.Muon_Trigger_SF_Key = d;
	  }
	}
        else executeEventFromParameter(param);
      }
    }
  }

}

void controlPlotDY::executeEventFromParameter(AnalyzerParameter param){

  bool PassSingleTrigger = false;
  bool PassDoubleTrigger = false;
  Event ev = GetEvent();
  double SingleTrigger_pTcut = 999;

  if(DataYear == 2016){
    if(IsSingleMuon){
      bool IsoMu24 = ev.PassTrigger("HLT_IsoMu24_v");
      bool IsoTkMu24 = ev.PassTrigger("HLT_IsoTkMu24_v");
      PassSingleTrigger = IsoMu24 || IsoTkMu24;
      SingleTrigger_pTcut = 26;
    }
    if(IsDoubleMuon){
      bool Mu17Mu8 = ev.PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");
      bool Mu17Mu8dz = ev.PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
      bool Mu17TkMu8 = ev.PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v");
      bool Mu17TkMu8dz = ev.PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");
      PassDoubleTrigger = Mu17Mu8 || Mu17Mu8dz || Mu17TkMu8 || Mu17TkMu8dz;
    }
  }
  else if(DataYear == 2017){
    if(IsSingleMuon){
      PassSingleTrigger = ev.PassTrigger("HLT_IsoMu27_v");
      SingleTrigger_pTcut = 29;
    }
    if(IsDoubleMuon){
      PassDoubleTrigger = ev.PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v");
    }
  }
  else{ //DataYear == 2018
    if(IsSingleMuon){
      PassSingleTrigger = ev.PassTrigger("HLT_IsoMu24_v");
      SingleTrigger_pTcut = 26;
    }
    if(IsDoubleMuon){
      PassDoubleTrigger = ev.PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v");
    }
  }

  //===== Event Selection ============================
  vector<Muon> muons = SelectMuons(AllMuons, param.Muon_Tight_ID, 10, 2.4);
  if(muons.size() !=2) return;
  if(muons.at(0).Charge() * muons.at(1).Charge() > 0.) return;

  if(IsSingleMuon){
    if(!PassSingleTrigger) return;
    if(muons.at(0).Pt() <= SingleTrigger_pTcut ) return;
  }
  else{
    if(!PassDoubleTrigger) return;
    if(muons.at(0).Pt() < 20 || muons.at(1).Pt() <10) return;
  }

  Particle ZCand = muons.at(0) + muons.at(1);
  if(ZCand.M() < 60. || ZCand.M() > 120.) return;

  if(IsDYsamples && muons != MuonPromptOnly(muons, gens)) return;

  //===== Weight to MC ===============================
  double weight = 1.;

  if(!IsDATA){

    weight *= weight_norm_1invpb*ev.GetTriggerLumi("Full");
    weight *= ev.MCweight();
    weight *= weight_Prefire;

    for(unsigned int i=0; i<muons.size(); i++){
      
      TString charge = "";
      if(muons.at(i).Charge() > 0) charge = "_plus";
      else charge = "_minus";
      double this_idsf = 1., this_isosf = 1.;

      //IDISO Scale Factor
      this_isosf = mcCorr->MuonISO_SF(param.Muon_ISO_SF_Key+charge, muons.at(i).Eta(), muons.at(i).MiniAODPt());

      weight *= this_idsf * this_isosf;    
    }
    //DoubleMuon or Single trigger SF (leading -> Mu17, subleading -> Mu8)
    TString trig = "";
    double this_trigsf = 1.;
    if(IsSingleMuon){
      if(DataYear == 2017) trig = "IsoMu27";
      else trig = "IsoMu24";
    }
    else if(IsDoubleMuon) trig = "Mu17";
    this_trigsf = mcCorr->MuonTrigger_SF(param.Muon_Trigger_SF_Key, trig, muons, 0);
    weight *= this_trigsf;

    if(DataYear == 2016){
      if(param.Name.Contains("2016BF")) weight *= 20.135171588 / 36.446609816;
      else if(param.Name.Contains("2016GH")) weight *= 16.311438228/ 36.446609816;
    }
    weight *= zptcor;
  }

  FillHist(IsFromTau+param.Name+"_mass", ZCand.M(), weight, 60, 60., 120.);
  FillHist(IsFromTau+param.Name+"_ZpT", ZCand.Pt(), weight, 120, 0., 120.);
  //FillHist(IsFromTau+param.Name+"_Zeta", ZCand.Eta(), weight, 100, -5.0, 5.0);
  FillHist(IsFromTau+param.Name+"_Zrapidity", ZCand.Rapidity(), weight, 60, -3.0, 3.0);
  FillHist(IsFromTau+param.Name+"_leading_pT", muons.at(0).Pt(), weight, 120, 0., 120.);
  FillHist(IsFromTau+param.Name+"_leading_eta", muons.at(0).Eta(), weight, 48, -2.4, 2.4);
  FillHist(IsFromTau+param.Name+"_subleading_pT", muons.at(1).Pt(), weight, 120, 0., 120.);
  FillHist(IsFromTau+param.Name+"_subleading_eta", muons.at(1).Eta(), weight, 48, -2.4, 2.4);

  if(muons.at(0).Charge() > 0 ){
    FillHist(IsFromTau+param.Name+"_leading_pT_P", muons.at(0).Pt(), weight, 120, 0., 120.);
    FillHist(IsFromTau+param.Name+"_leading_eta_P", muons.at(0).Eta(), weight, 48, -2.4, 2.4);
    FillHist(IsFromTau+param.Name+"_subleading_pT_M", muons.at(1).Pt(), weight, 120, 0., 120.);
    FillHist(IsFromTau+param.Name+"_subleading_eta_M", muons.at(1).Eta(), weight, 48, -2.4, 2.4);
  }
  else{
    FillHist(IsFromTau+param.Name+"_leading_pT_M", muons.at(0).Pt(), weight, 120, 0., 120.);
    FillHist(IsFromTau+param.Name+"_leading_eta_M", muons.at(0).Eta(), weight, 48, -2.4, 2.4);
    FillHist(IsFromTau+param.Name+"_subleading_pT_P", muons.at(1).Pt(), weight, 120, 0., 120.);
    FillHist(IsFromTau+param.Name+"_subleading_eta_P", muons.at(1).Eta(), weight, 48, -2.4, 2.4);
  }

  if(IsDATA || IsDYsamples){
    vector<double> massbin = {60., 70., 78., 84., 87., 89., 91., 93., 95., 98., 104., 112., 120. };
    vector<double> rapiditybin = {0., 0.4, 0.8, 1.2, 1.6, 2.0, 2.4 };
    double cosCS = GetCosineCS(muons.at(0),muons.at(1));

    for(unsigned int x=0; x<massbin.size()-1; x++){
      if(massbin.at(x) < ZCand.M() && ZCand.M() < massbin.at(x+1)){
        for(unsigned int y=0; y<rapiditybin.size()-1; y++){
          if(rapiditybin.at(y) < fabs(ZCand.Rapidity()) && fabs(ZCand.Rapidity()) < rapiditybin.at(y+1)){
	    FillHist(IsFromTau+param.Name+"_mass"+TString::Itoa(x,10)+"_rap"+TString::Itoa(y,10)+"_CosineCS", cosCS, weight, 20, -1.0, 1.0);
	  }
        }
      }
    }
    FillHist(IsFromTau+param.Name+"_CosineCS", cosCS, weight, 60, -1.0, 1.0);
  }

  //For ZpT Correction SF.root -Before iteration
  if(IsDATA || param.Name.Contains("Hybrid_TrigSF")){
    double etabins[] = {-2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0., 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};
    double ptbins[] = {0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10., 11.,12.,13.,14.,15.,16.,17.,18.,19.,20., 22.,24.,26.,28.,30., 32.,34.,36.,38.,40.,44.,48.,52.,56.,60.};

    FillHist(IsFromTau+param.Name+"_ZpTEta", ZCand.Pt(), ZCand.Eta(), weight, sizeof(ptbins)/sizeof(ptbins[0]) - 1, ptbins, sizeof(etabins)/sizeof(etabins[0]) - 1, etabins);
    FillHist(IsFromTau+param.Name+"_ZpTRapidity", ZCand.Pt(), ZCand.Rapidity(), weight, sizeof(ptbins)/sizeof(ptbins[0]) - 1, ptbins, sizeof(etabins)/sizeof(etabins[0]) - 1, etabins);
    }
}

controlPlotDY::controlPlotDY(){

}

controlPlotDY::~controlPlotDY(){

}

double controlPlotDY::GetCosineCS(Muon mu0, Muon mu1){
  TLorentzVector lep0, lep1;
  if(mu0.Charge() < 0 && mu1.Charge() > 0){
    lep0 = mu0;
    lep1 = mu1;
  }
  else if(mu0.Charge() > 0 && mu1.Charge() < 0){
    lep0 = mu1;
    lep1 = mu0;
  }

  TLorentzVector dilepton = lep0 + lep1;
  double lep0plus = (lep0.E()+lep0.Pz())/sqrt(2);
  double lep0minus = (lep0.E()-lep0.Pz())/sqrt(2);
  double lep1plus = (lep1.E()+lep1.Pz())/sqrt(2);
  double lep1minus = (lep1.E()-lep1.Pz())/sqrt(2);
  double dimass = dilepton.M();
  double dipt = dilepton.Pt();
  double direction = dilepton.Pz()/fabs(dilepton.Pz());

  return direction*2*(lep0plus*lep1minus-lep0minus*lep1plus)/sqrt(dimass*dimass*(dimass*dimass+dipt*dipt));
}

