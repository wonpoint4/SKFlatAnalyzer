#include "SkimTree_GEN.h"

void SkimTree_GEN::initializeAnalyzer(){

  outfile->cd();
  cout << "[SkimTree_GEN::initializeAnalyzer()] gDirectory = " << gDirectory->GetName() << endl;
  newtree = fChain->CloneTree(0);

  newtree->SetBranchStatus("*",0);
  newtree->SetBranchStatus("gen*",1);
}

void SkimTree_GEN::executeEvent(){

  newtree->Fill();

}

void SkimTree_GEN::executeEventFromParameter(AnalyzerParameter param){

}

SkimTree_GEN::SkimTree_GEN(){

}

SkimTree_GEN::~SkimTree_GEN(){

}

void SkimTree_GEN::WriteHist(){

  outfile->mkdir("recoTree");
  outfile->cd("recoTree");
  newtree->Write();
  outfile->cd();

}


