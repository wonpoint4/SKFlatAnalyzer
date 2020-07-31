#include "SkimTree_SS2lOR3l.h"

void SkimTree_SS2lOR3l::initializeAnalyzer(){

  outfile->cd();
  cout << "[SkimTree_SS2lOR3l::initializeAnalyzer()] gDirectory = " << gDirectory->GetName() << endl;
  newtree = fChain->CloneTree(0);

}

void SkimTree_SS2lOR3l::executeEvent(){

  Event ev;

  std::vector<Muon>     muonPreColl     = GetMuons("NOCUT", 8., 2.4);
  std::vector<Electron> electronPreColl = GetElectrons("NOCUT", 8., 2.5);

  std::sort(muonPreColl.begin(), muonPreColl.end(), PtComparing);
  std::sort(electronPreColl.begin(), electronPreColl.end(), PtComparing);

  int NEl  = electronPreColl.size();
  int NMu  = muonPreColl.size();
  int NLep = NEl+NMu;
  bool HasSS2lOR3l = false;
  bool LeadLepPt17 = false;
  if      ( NLep >= 3 ){ HasSS2lOR3l = true; }
  else if ( NLep == 2 ){
    int QTot = 0;
    for(unsigned int it_m=0; it_m<muonPreColl.size(); it_m++){ QTot+=muonPreColl.at(it_m).Charge(); }
    for(unsigned int it_e=0; it_e<electronPreColl.size(); it_e++){ QTot+=electronPreColl.at(it_e).Charge(); }
    if(abs(QTot)==2) HasSS2lOR3l = true;
  }
  if(NMu>0 && muonPreColl.at(0).Pt()>17.    ) LeadLepPt17 = true;
  if(NEl>0 && electronPreColl.at(0).Pt()>17.) LeadLepPt17 = true;

  if( !(HasSS2lOR3l && LeadLepPt17) ) return;

  newtree->Fill();

}

void SkimTree_SS2lOR3l::executeEventFromParameter(AnalyzerParameter param){

}

SkimTree_SS2lOR3l::SkimTree_SS2lOR3l(){

  newtree = NULL;

}

SkimTree_SS2lOR3l::~SkimTree_SS2lOR3l(){

}

void SkimTree_SS2lOR3l::WriteHist(){

  outfile->mkdir("recoTree");
  outfile->cd("recoTree");
  newtree->Write();
  outfile->cd();

}

