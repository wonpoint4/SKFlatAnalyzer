void run(){

  double lumi_periodB = 5750.490644035;
  double lumi_periodC = 2572.903488748;
  double lumi_periodD = 4242.291556970;
  double lumi_periodE = 4025.228136967;
  double lumi_periodF = 3104.509131800;
  double lumi_periodG = 7575.824256098;
  double lumi_periodH = 8650.628380028;
  double total_lumi = (lumi_periodB+lumi_periodC+lumi_periodD+lumi_periodE+lumi_periodF+lumi_periodG+lumi_periodH);

  double WeightBtoF = (lumi_periodB+lumi_periodC+lumi_periodD+lumi_periodE+lumi_periodF)/total_lumi;
  double WeightGtoH = (lumi_periodG+lumi_periodH)/total_lumi;

  //==== ID
  TFile *file_ID_BtoF = new TFile("original/RunBCDEF_SF_ID.root");
  TFile *file_ID_GtoH = new TFile("original/RunGH_SF_ID.root");

  vector<TString> hists_ID = {
"NUM_LooseID_DEN_genTracks_eta_pt",
"NUM_HighPtID_DEN_genTracks_eta_pair_newTuneP_probe_pt",
"NUM_TightID_DEN_genTracks_eta_pt",
"NUM_MediumID_DEN_genTracks_eta_pt"
  };

  TFile *file_ID_new = new TFile("output/RunAveraged_SF_ID.root","RECREATE");
  for(unsigned int i=0; i<hists_ID.size(); i++){
    TString histname = hists_ID.at(i);

    TH2D *hist_BtoF = (TH2D *)file_ID_BtoF->Get(histname);
    hist_BtoF->Scale(WeightBtoF);
    TH2D *hist_GtoH = (TH2D *)file_ID_GtoH->Get(histname);
    hist_GtoH->Scale(WeightGtoH);

    TH2D *hist_out = (TH2D *)hist_BtoF->Clone();
    hist_out->Add(hist_GtoH);
    file_ID_new->cd();
    hist_out->Write();

  }
  file_ID_new->Close();

  //==== ISO
  TFile *file_ISO_BtoF = new TFile("original/RunBCDEF_SF_ISO.root");
  TFile *file_ISO_GtoH = new TFile("original/RunGH_SF_ISO.root");

  vector<TString> hists_ISO = {
"NUM_LooseRelTkIso_DEN_HighPtIDandIPCut_eta_pair_newTuneP_probe_pt",
"NUM_TightRelIso_DEN_MediumID_eta_pt",
"NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt",
"NUM_LooseRelIso_DEN_MediumID_eta_pt",
"NUM_LooseRelIso_DEN_LooseID_eta_pt",
"NUM_LooseRelIso_DEN_TightIDandIPCut_eta_pt",
  };

  TFile *file_ISO_new = new TFile("output/RunAveraged_SF_ISO.root","RECREATE");
  for(unsigned int i=0; i<hists_ISO.size(); i++){
    TString histname = hists_ISO.at(i);

    TH2D *hist_BtoF = (TH2D *)file_ISO_BtoF->Get(histname);
    hist_BtoF->Scale(WeightBtoF);
    TH2D *hist_GtoH = (TH2D *)file_ISO_GtoH->Get(histname);
    hist_GtoH->Scale(WeightGtoH);

    TH2D *hist_out = (TH2D *)hist_BtoF->Clone();
    hist_out->Add(hist_GtoH);
    file_ISO_new->cd();
    hist_out->Write();

  }
  file_ISO_new->Close();

  //==== Trigger
  TFile *file_Trigger_BtoF = new TFile("original/EfficienciesAndSF_RunBtoF.root");
  TFile *file_Trigger_GtoH = new TFile("original/EfficienciesAndSF_Period4.root");

  vector<TString> hists_Trigger = {
"IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesDATA/pt_abseta_DATA",
"IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesMC/pt_abseta_MC",
"Mu50_OR_TkMu50_PtEtaBins/efficienciesDATA/pt_abseta_DATA",
"Mu50_OR_TkMu50_PtEtaBins/efficienciesMC/pt_abseta_MC",
  };
  vector<TString> hists_Trigger_outDirName = {
"IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesDATA/",
"IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesMC/",
"Mu50_OR_TkMu50_PtEtaBins/efficienciesDATA/",
"Mu50_OR_TkMu50_PtEtaBins/efficienciesMC/",
  };
  vector<TString> hists_Trigger_outHistName = {
"pt_abseta_DATA",
"pt_abseta_MC",
"pt_abseta_DATA",
"pt_abseta_MC",
  };

  TFile *file_Trigger_new = new TFile("output/RunAveraged_Eff_Trigger.root","RECREATE");
  for(unsigned int i=0; i<hists_Trigger.size(); i++){
    TString histname = hists_Trigger.at(i);

    TH2F *hist_BtoF = (TH2F *) file_Trigger_BtoF->Get(histname);
    hist_BtoF->Scale(WeightBtoF);
    TH2F *hist_GtoH = (TH2F *) file_Trigger_GtoH->Get(histname);
    hist_GtoH->Scale(WeightGtoH);

    TH2F *hist_out = (TH2F *)hist_BtoF->Clone();
    hist_out->Add(hist_GtoH);
    hist_out->SetName(hists_Trigger_outHistName.at(i));

    file_Trigger_new->cd();
    file_Trigger_new->mkdir( hists_Trigger_outDirName.at(i) );
    file_Trigger_new->cd( hists_Trigger_outDirName.at(i) );
    hist_out->Write();
    file_Trigger_new->cd();

  }
  file_Trigger_new->Close();

}
