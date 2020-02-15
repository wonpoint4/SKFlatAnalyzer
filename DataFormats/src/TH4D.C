#include "TH4D.h"

ClassImp(TH4D);

TH4D::TH4D(){
  fXaxis.SetName("xaxis");
}
TH4D::~TH4D(){
  for(auto hist:hists)
    if(hist) delete hist;
  hists.clear();
}
TH4D::TH4D(const char *name,const char *title,Int_t nbinsx,Double_t xlow,Double_t xup
	   ,Int_t nbinsy,Double_t ylow,Double_t yup
	   ,Int_t nbinsz,Double_t zlow,Double_t zup
	   ,Int_t nbinsu,Double_t ulow,Double_t uup){
  SetNameTitle(name,title);
  fXaxis.SetName("xaxis");
  if (nbinsx <= 0) {
    Warning("TH4D","nbinsx is <=0 - set to nbinsx = 1");
    nbinsx = 1;
  }
  fXaxis.Set(nbinsx,xlow,xup);
  for(int i=0;i<nbinsx+2;i++){
    hists.push_back(new TH3D(Form("%d",i),Form("%d",i),nbinsy,ylow,yup,nbinsz,zlow,zup,nbinsu,ulow,uup));
    hists.back()->Sumw2();
    hists.back()->SetDirectory(0);
  }
}
TH4D::TH4D(const char *name,const char *title,Int_t nbinsx,const Double_t *xbins
	   ,Int_t nbinsy,const Double_t *ybins
	   ,Int_t nbinsz,const Double_t *zbins
	   ,Int_t nbinsu,const Double_t *ubins){
  SetNameTitle(name,title);
  fXaxis.SetName("xaxis");
  if (nbinsx <= 0) {
    Warning("TH4D","nbinsx is <=0 - set to nbinsx = 1");
    nbinsx = 1;
  }
  fXaxis.Set(nbinsx,xbins);
  for(int i=0;i<nbinsx+2;i++){
    hists.push_back(new TH3D(Form("%d",i),Form("%d",i),nbinsy,ybins,nbinsz,zbins,nbinsu,ubins));
    hists.back()->Sumw2();
    hists.back()->SetDirectory(0);
  }
}

Int_t TH4D::Fill(Double_t x, Double_t y, Double_t z, Double_t u, Double_t w){
  return hists.at(fXaxis.FindBin(x))->Fill(y,z,u,w);
}
Int_t TH4D::Fill(Double_t x, Double_t y, Double_t z, Double_t u){
  return hists.at(fXaxis.FindBin(x))->Fill(y,z,u);
}

Double_t TH4D::IntegralAndError(Int_t binx1, Int_t binx2, Int_t biny1, Int_t biny2, Int_t binz1, Int_t binz2, Int_t binu1, Int_t binu2, Double_t & error, Option_t *option) const {
  Double_t integral=0;
  Double_t error2sum=0;
  for(int i=binx1;i<=binx2;i++){
    Double_t this_error=0;
    integral+=hists[i]->IntegralAndError(biny1,biny2,binz1,binz2,binu1,binu2,this_error,option);
    error2sum=this_error*this_error;
  }
  error=sqrt(error2sum);
  return integral;
}
Double_t TH4D::Integral(Int_t binx1, Int_t binx2, Int_t biny1, Int_t biny2, Int_t binz1, Int_t binz2, Int_t binu1, Int_t binu2, Option_t *option) const {
  Double_t dummy=0;
  return IntegralAndError(binx1,binx2,biny1,biny2,binz1,binz2,binu1,binu2,dummy,option);
}
Double_t TH4D::Integral(Option_t *option) const {
  return Integral(fXaxis.GetFirst(),fXaxis.GetLast(),
		  hists.at(0)->GetXaxis()->GetFirst(),hists.at(0)->GetXaxis()->GetLast(),
		  hists.at(0)->GetYaxis()->GetFirst(),hists.at(0)->GetYaxis()->GetLast(),
		  hists.at(0)->GetZaxis()->GetFirst(),hists.at(0)->GetZaxis()->GetLast(),option);
}
  
bool TH4D::CheckConsistency() const {
  int nbinsx=fXaxis.GetNbins();
  bool consistent=true;
  for(int i=1;i<nbinsx+2;i++){
    consistent&=TH1::CheckConsistency((TH1*)hists.at(0),(TH1*)hists.at(1));
  }
  if(!consistent){
    Error("TH4D","not consistent");
    return false;
  }
  return true;
}
bool TH4D::CheckConsistency(const TH4D* h1,const TH4D* h2){
  if(!h1->CheckConsistency()) return false;
  if(!h2->CheckConsistency()) return false;
  if(!CheckEqualAxes(h1->GetXaxis(),h2->GetXaxis())) return false; 
  return true;
}
