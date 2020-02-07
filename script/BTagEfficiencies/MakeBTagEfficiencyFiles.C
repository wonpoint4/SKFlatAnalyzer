#include <string.h>
#include "TFile.h"
#include "TH2.h"
#include <iostream>
#include "TString.h"
#include <sstream>
#include <iostream>
#include <fstream>

void MakeBTagEfficiencyFiles();
void MakeBTagEfficiencyFiles(TString year, bool set_verbose=false);

void MakeBTagEfficiencyFiles(){
  cout << "Will run 2016/2017/2018 samples. To run individual years use MakeBTagEfficiencyFiles(year)" << endl;
  
  MakeBTagEfficiencyFiles("16");
  MakeBTagEfficiencyFiles("17");
  MakeBTagEfficiencyFiles("18");
  
}

void MakeBTagEfficiencyFiles(TString year, bool set_verbose=false){

  if(year != "16" && year != "17" && year != "18" ){
    cout << "Year is set wrong, run with root -l -q -b 'MakeBTagEfficiencyFiles.C(year)' where year is a string equal to 16,17 or 18 as input" << endl;
    return;
  }

  //=== For now have TT LL sample only, can have full list of efficiencies for each background

  TString sample = "JetTagging_TTLL_powheg.root";

  //=== list of hist names for numerator
  /*
    KEY: TH2DJet_2017CSVv2_Loose_eff_b_num_;1
    KEY: TH2DJet_2017CSVv2_Loose_eff_c_num_;1
    KEY: TH2DJet_2017CSVv2_Loose_eff_udsg_num_;1
    KEY: TH2DJet_2017CSVv2_Medium_eff_b_num_;1
    KEY: TH2DJet_2017CSVv2_Medium_eff_c_num_;1
    KEY: TH2DJet_2017CSVv2_Medium_eff_udsg_num_;1
    KEY: TH2DJet_2017CSVv2_Tight_eff_b_num_;1
    KEY: TH2DJet_2017CSVv2_Tight_eff_c_num_;1
    KEY: TH2DJet_2017CSVv2_Tight_eff_udsg_num_;1
    KEY: TH2DJet_2017DeepCSV_Loose_eff_b_num_;1
    KEY: TH2DJet_2017DeepCSV_Loose_eff_c_num_;1
    KEY: TH2DJet_2017DeepCSV_Loose_eff_udsg_num_;1
    KEY: TH2DJet_2017DeepCSV_Medium_eff_b_num_;1
    KEY: TH2DJet_2017DeepCSV_Medium_eff_c_num_;1
    KEY: TH2DJet_2017DeepCSV_Medium_eff_udsg_num_;1
    KEY: TH2DJet_2017DeepCSV_Tight_eff_b_num_;1
    KEY: TH2DJet_2017DeepCSV_Tight_eff_c_num_;1
    KEY: TH2DJet_2017DeepCSV_Tight_eff_udsg_num_;1
  */
  
  TString skout_dir = getenv("SKFlatOutputDir");
  TString path= skout_dir+"/Run2Legacy_v3/JetTagging/20"+year+"/" + sample ;


  TFile * f = new TFile(path );

  ofstream lumi_file;
  string lfile =  "BTagEfficienciesTTbar"+string(year)+".C";
  
  lumi_file.open(lfile.c_str());
  lumi_file.setf(ios::fixed,ios::floatfield);
  
  lumi_file.precision(5);
  lumi_file <<  "/* These Efficiencies have been derived for 20" << year << " "<< sample <<" events and should " << endl;
  lumi_file <<  "be used only for the same MC samples or for events with similar topology */  " << endl;

  lumi_file <<  "\n" << endl;

  vector <TString> flavour;
  flavour.push_back("B");
  flavour.push_back("C");
  flavour.push_back("Light");

  vector <TString> fileflavour;
  fileflavour.push_back("b");
  fileflavour.push_back("c");
  fileflavour.push_back("udsg");
  
  for(unsigned int _a = 0; _a < fileflavour.size(); _a++){

    lumi_file <<  "float BTagSFUtil::TagEfficiency"<< flavour.at(_a)<<"_20"<< year<<"(float JetPt, float JetEta) {" <<endl;
    lumi_file <<  "\n" <<endl;

    vector <TString> taggers;
    taggers.push_back("CSVv2_Loose");
    taggers.push_back("CSVv2_Medium");
    taggers.push_back("CSVv2_Tight");
    taggers.push_back("DeepCSV_Loose");
    taggers.push_back("DeepCSV_Medium");
    taggers.push_back("DeepCSV_Tight");
    taggers.push_back("DeepJet_Loose");
    taggers.push_back("DeepJet_Medium");
    taggers.push_back("DeepJet_Tight");

    TH2D* denominator  = (TH2D*)f->Get(("Jet_20"+year+"_eff_"+fileflavour.at(_a)+"_denom_"));
    
    for(unsigned int it = 0 ; it < taggers.size() ; it++){

      TH2D* numerator  = (TH2D*)f->Get(("Jet_20"+year+taggers.at(it)+"_eff_"+fileflavour.at(_a)+"_num_"));

      if(!numerator) continue;

      lumi_file <<  "\t if (TaggerOP==\""<<taggers.at(it)<<"\") { " <<endl;
      if(set_verbose) cout << ("Jet_20"+year+taggers.at(it)+"_eff_"+fileflavour.at(_a)+"_num_") << endl;
      numerator->Divide(denominator);

      std::vector<float> jpt1;
      jpt1.push_back(20.);
      jpt1.push_back(40.);
      jpt1.push_back(60.);
      jpt1.push_back(80.);
      jpt1.push_back(100.);
      jpt1.push_back(120.);
      std::vector<float> jpt2;
      jpt2.push_back(40.);
      jpt2.push_back(60.);
      jpt2.push_back(80.);
      jpt2.push_back(100.);
      jpt2.push_back(120.);
      jpt2.push_back(3000.);
     
      std::vector<float> jeta1;
      jeta1.push_back(0.);
      jeta1.push_back(0.6);
      jeta1.push_back(1.2);
      jeta1.push_back(1.8);
      std::vector<float> jeta2;
      jeta2.push_back(0.6);
      jeta2.push_back(1.2);
      jeta2.push_back(1.8);
      jeta2.push_back(2.4);
      
      for(unsigned int x = 0; x < jpt1.size(); x++){
	if(x == 0)
	  {  lumi_file <<  "\t\t if (JetPt > "<< jpt1.at(x)<<" && JetPt <= "<<  jpt2.at(x) <<"){  " << endl;
	    if(set_verbose)cout <<  "   if (JetPt > "<< jpt1.at(x)<<" && JetPt <= "<<  jpt2.at(x) <<"){  " << endl;
	  }
	else {
	  lumi_file <<  "\t\t else if (JetPt > "<< jpt1.at(x)<<" && JetPt <= "<<  jpt2.at(x) <<"){  " << endl;
	  if(set_verbose)cout <<  "   else if (JetPt > "<< jpt1.at(x)<<" && JetPt <= "<<  jpt2.at(x) <<"){  " << endl;
	  
	}

	for(unsigned int y = 0; y < jeta1.size(); y++){
	  if(y == 0){
	    lumi_file <<  "\t\t\t if      (fabs(JetEta) > " <<jeta1.at(y) <<"   && fabs(JetEta) <= "<<  jeta2.at(y) << ") return " << numerator->GetBinContent(y+1,x+1)<< " ;  " << endl;
	    if(set_verbose)cout <<  "     if      (fabs(JetEta) > " <<jeta1.at(y) <<"   && fabs(JetEta) <= "<<  jeta2.at(y) << ") return " << numerator->GetBinContent(y+1,x+1)<< " ;  " << endl;

	  }
	  else {
	    lumi_file <<  "\t\t\t else if (fabs(JetEta) > " <<jeta1.at(y) <<" && fabs(JetEta) <= "<<  jeta2.at(y) << ") return " << numerator->GetBinContent(y+1,x+1)<< ";  " << endl;
	    if(set_verbose) cout <<  "     else if (fabs(JetEta) > " <<jeta1.at(y) <<" && fabs(JetEta) <= "<<  jeta2.at(y) << ") return " << numerator->GetBinContent(y+1,x+1)<< ";  " << endl;
	  }
	}
	lumi_file <<  "\t\t\t } " << endl;
      }
      lumi_file <<  "\t\t } " << endl;
      
    }
    lumi_file <<  "\t return 1.; " << endl;
    lumi_file <<  "} " << endl;

  }

  cout << "**************************" << endl;
  cout << "Output file " << lfile << " created" << endl;
  cout << "**************************" << endl;
}
