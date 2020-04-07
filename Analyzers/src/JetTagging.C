#include "JetTagging.h"

JetTagging::JetTagging(){

  }

void JetTagging::initializeAnalyzer(){


}

JetTagging::~JetTagging(){

}

void JetTagging::executeEvent(){


  //==== Get L1Prefire reweight
  //==== If data, 1.;
  //==== If MC && DataYear > 2017, 1.;
  //==== If MC && DataYear <= 2017, we have to reweight the event with this value
  //==== I defined "double weight_Prefire;" in Analyzers/include/JetTagging.h
  weight_Prefire = GetPrefireWeight(0);

  //==== Declare AnalyzerParameter

  AnalyzerParameter param;

  //==== clear parameter set                                                                                                                                                                                                                                                  
  param.Clear();

  //==== set which systematic sources you want to run this time                                                                                                                                                                                                               
  //==== default syst_ is AnalyzerParameter::Central                                                                                                                                                                                                                          
  param.syst_ = AnalyzerParameter::Central;


  //==== And, Jet ID                                                                                                                                                                                                                                                          
  param.Jet_ID = "tight";
  param.FatJet_ID = "tight";


  AllJets = GetAllJets();



  //==== Now, all parameters are set. Run executeEventFromParameter() with this parameter set                                                                                                                                                                                 
  executeEventFromParameter(param);



  

}

void JetTagging::executeEventFromParameter(AnalyzerParameter param){

  //========================
  //==== MET Filter
  //========================

  if(!PassMETFilter()) return;

  Event ev = GetEvent();
  Particle METv = ev.GetMETVector();

  vector<Jet> this_AllJets = AllJets;

  /// use same ptmin / eta that CSV files use for SF
  vector<Jet> jets = SelectJets(this_AllJets, param.Jet_ID, 20., 2.4);
  std::sort(jets.begin(), jets.end(), PtComparing);

  double etabins [5] = {0.,0.6, 1.2, 1.8,2.4};
  double ptbins [7] = {20., 40., 60., 80., 100., 120., 3000.}; 
  TString yeartag= "2017";
  if(DataYear == 2016) yeartag= "2016";
  if(DataYear == 2018) yeartag= "2018";

  // === code to measure btag efficiencies in TT MC
  for(unsigned int ij = 0 ; ij < jets.size(); ij++){

    //-=== Get Efficiency for btag SF following https://github.com/rappoccio/usercode/blob/Dev_53x/EDSHyFT/plugins/BTaggingEffAnalyzer.cc
    if(DataYear == 2016){

      // taggers and values from https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation2016Legacy

      TString flav= "b";
      if(fabs(jets.at(ij).hadronFlavour()) == 4)flav= "c";
      if(fabs(jets.at(ij).hadronFlavour()) == 0)flav= "udsg";
      

      //== Fill denominator, which is all jets in TT bar sample
      JSFillHist(param.Name, "Jet_"+yeartag+"_eff_"+flav+"_denom_"+param.Name, fabs(jets.at(ij).Eta()), jets.at(ij).Pt(), ev.MCweight(), 4, etabins, 6, ptbins);
      

      if(jets.at(ij).GetTaggerResult(Jet::DeepCSV) > 0.2217)
	JSFillHist(param.Name, "Jet_"+yeartag+"DeepCSV_Loose_eff_"+flav+"_num_"+param.Name, fabs(jets.at(ij).Eta()), jets.at(ij).Pt(), ev.MCweight(), 4, etabins, 6, ptbins);
      if(jets.at(ij).GetTaggerResult(Jet::DeepCSV) > 0.6321)
	JSFillHist(param.Name, "Jet_"+yeartag+"DeepCSV_Medium_eff_"+flav+"_num_"+param.Name, fabs(jets.at(ij).Eta()), jets.at(ij).Pt(), ev.MCweight(), 4, etabins, 6, ptbins);
      if(jets.at(ij).GetTaggerResult(Jet::DeepCSV) > 0.8953)
	JSFillHist(param.Name, "Jet_"+yeartag+"DeepCSV_Tight_eff_"+flav+"_num_"+param.Name, fabs(jets.at(ij).Eta()), jets.at(ij).Pt(), ev.MCweight(), 4, etabins, 6, ptbins);
      
      
      if(jets.at(ij).GetTaggerResult(Jet::DeepJet) > 0.0614)
        JSFillHist(param.Name, "Jet_"+yeartag+"DeepJet_Loose_eff_"+flav+"_num_"+param.Name, fabs(jets.at(ij).Eta()), jets.at(ij).Pt(), ev.MCweight(), 4, etabins, 6, ptbins);
      if(jets.at(ij).GetTaggerResult(Jet::DeepJet) > 0.3093)
        JSFillHist(param.Name, "Jet_"+yeartag+"DeepJet_Medium_eff_"+flav+"_num_"+param.Name, fabs(jets.at(ij).Eta()), jets.at(ij).Pt(), ev.MCweight(), 4, etabins, 6, ptbins);
      if(jets.at(ij).GetTaggerResult(Jet::DeepJet) > 0.7221)
        JSFillHist(param.Name, "Jet_"+yeartag+"DeepJet_Tight_eff_"+flav+"_num_"+param.Name, fabs(jets.at(ij).Eta()), jets.at(ij).Pt(), ev.MCweight(), 4, etabins, 6, ptbins);

    }// 2016      
    
    
    if(DataYear == 2017){
      
      // Tagger and values from https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
      
      TString flav= "b";
      if(fabs(jets.at(ij).hadronFlavour()) == 4)flav= "c";
      if(fabs(jets.at(ij).hadronFlavour()) == 0)flav= "udsg";
      
      
      JSFillHist(param.Name, "Jet_"+yeartag+"_eff_"+flav+"_denom_"+param.Name, fabs(jets.at(ij).Eta()), jets.at(ij).Pt(), ev.MCweight(), 4, etabins, 6, ptbins);
      
      if(jets.at(ij).GetTaggerResult(Jet::DeepCSV) > 0.1522)        
	JSFillHist(param.Name, "Jet_"+yeartag+"DeepCSV_Loose_eff_"+flav+"_num_"+param.Name, fabs(jets.at(ij).Eta()), jets.at(ij).Pt(), ev.MCweight(), 4, etabins, 6, ptbins);
      if(jets.at(ij).GetTaggerResult(Jet::DeepCSV) > 0.4941)     
	JSFillHist(param.Name, "Jet_"+yeartag+"DeepCSV_Medium_eff_"+flav+"_num_"+param.Name, fabs(jets.at(ij).Eta()), jets.at(ij).Pt(), ev.MCweight(), 4, etabins, 6, ptbins);
      if(jets.at(ij).GetTaggerResult(Jet::DeepCSV) > 0.8001)        
	JSFillHist(param.Name, "Jet_"+yeartag+"DeepCSV_Tight_eff_"+flav+"_num_"+param.Name, fabs(jets.at(ij).Eta()), jets.at(ij).Pt(), ev.MCweight(), 4, etabins, 6, ptbins);
      
      if(jets.at(ij).GetTaggerResult(Jet::CSVv2) > 0.5803)
	JSFillHist(param.Name, "Jet_"+yeartag+"CSVv2_Loose_eff_"+flav+"_num_"+param.Name, fabs(jets.at(ij).Eta()), jets.at(ij).Pt(), ev.MCweight(), 4, etabins, 6, ptbins);
      if(jets.at(ij).GetTaggerResult(Jet::CSVv2) > 0.8838)
	JSFillHist(param.Name, "Jet_"+yeartag+"CSVv2_Medium_eff_"+flav+"_num_"+param.Name, fabs(jets.at(ij).Eta()), jets.at(ij).Pt(), ev.MCweight(), 4, etabins, 6, ptbins);
      if(jets.at(ij).GetTaggerResult(Jet::CSVv2) > 0.9693)
	JSFillHist(param.Name, "Jet_"+yeartag+"CSVv2_Tight_eff_"+flav+"_num_"+param.Name, fabs(jets.at(ij).Eta()), jets.at(ij).Pt(), ev.MCweight(), 4, etabins, 6, ptbins);

      if(jets.at(ij).GetTaggerResult(Jet::DeepJet) > 0.0521)
	JSFillHist(param.Name, "Jet_"+yeartag+"DeepJet_Loose_eff_"+flav+"_num_"+param.Name, fabs(jets.at(ij).Eta()), jets.at(ij).Pt(), ev.MCweight(), 4, etabins, 6, ptbins);
      if(jets.at(ij).GetTaggerResult(Jet::DeepJet) > 0.3033)
        JSFillHist(param.Name, "Jet_"+yeartag+"DeepJet_Medium_eff_"+flav+"_num_"+param.Name, fabs(jets.at(ij).Eta()), jets.at(ij).Pt(), ev.MCweight(), 4, etabins, 6, ptbins);
      if(jets.at(ij).GetTaggerResult(Jet::DeepJet) > 0.7489)
        JSFillHist(param.Name, "Jet_"+yeartag+"DeepJet_Tight_eff_"+flav+"_num_"+param.Name, fabs(jets.at(ij).Eta()), jets.at(ij).Pt(), ev.MCweight(), 4, etabins, 6, ptbins);
     
      
    }// 2017

    if(DataYear == 2018){

      // Tagger and values from https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X

      TString flav= "b";
      if(fabs(jets.at(ij).hadronFlavour()) == 4)flav= "c";
      if(fabs(jets.at(ij).hadronFlavour()) == 0)flav= "udsg";


      JSFillHist(param.Name, "Jet_"+yeartag+"_eff_"+flav+"_denom_"+param.Name, fabs(jets.at(ij).Eta()), jets.at(ij).Pt(), ev.MCweight(), 4, etabins, 6, ptbins);

      if(jets.at(ij).GetTaggerResult(Jet::DeepCSV) > 0.1241)
        JSFillHist(param.Name, "Jet_"+yeartag+"DeepCSV_Loose_eff_"+flav+"_num_"+param.Name, fabs(jets.at(ij).Eta()), jets.at(ij).Pt(), ev.MCweight(), 4, etabins, 6, ptbins);
      if(jets.at(ij).GetTaggerResult(Jet::DeepCSV) > 0.4184)
        JSFillHist(param.Name, "Jet_"+yeartag+"DeepCSV_Medium_eff_"+flav+"_num_"+param.Name, fabs(jets.at(ij).Eta()), jets.at(ij).Pt(), ev.MCweight(), 4, etabins, 6, ptbins);
      if(jets.at(ij).GetTaggerResult(Jet::DeepCSV) > 0.7527)
        JSFillHist(param.Name, "Jet_"+yeartag+"DeepCSV_Tight_eff_"+flav+"_num_"+param.Name, fabs(jets.at(ij).Eta()), jets.at(ij).Pt(), ev.MCweight(), 4, etabins, 6, ptbins);


      if(jets.at(ij).GetTaggerResult(Jet::DeepJet) > 0.0494)
        JSFillHist(param.Name, "Jet_"+yeartag+"DeepJet_Loose_eff_"+flav+"_num_"+param.Name, fabs(jets.at(ij).Eta()), jets.at(ij).Pt(), ev.MCweight(), 4, etabins, 6, ptbins);
      if(jets.at(ij).GetTaggerResult(Jet::DeepJet) > 0.2770)
        JSFillHist(param.Name, "Jet_"+yeartag+"DeepJet_Medium_eff_"+flav+"_num_"+param.Name, fabs(jets.at(ij).Eta()), jets.at(ij).Pt(), ev.MCweight(), 4, etabins, 6, ptbins);
      if(jets.at(ij).GetTaggerResult(Jet::DeepJet) > 0.7264)
        JSFillHist(param.Name, "Jet_"+yeartag+"DeepJet_Tight_eff_"+flav+"_num_"+param.Name, fabs(jets.at(ij).Eta()), jets.at(ij).Pt(), ev.MCweight(), 4, etabins, 6, ptbins);

    }// 2018                                                                                                                                                                                                                                                                    


  }


}



