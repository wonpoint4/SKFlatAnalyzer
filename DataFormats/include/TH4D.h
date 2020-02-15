#ifndef HS_TH4D
#define HS_TH4D
#include <vector>
#include "TH3.h"
#include "TFile.h"

class TH4D : public TH1{
 public:
  std::vector<TH3D*> hists;
  
  TH4D();
  ~TH4D();
  TH4D(const char *name,const char *title,Int_t nbinsx,Double_t xlow,Double_t xup
      ,Int_t nbinsy,Double_t ylow,Double_t yup
      ,Int_t nbinsz,Double_t zlow,Double_t zup
      ,Int_t nbinsu,Double_t ulow,Double_t uup);
  TH4D(const char *name,const char *title,Int_t nbinsx,const Double_t *xbins
      ,Int_t nbinsy,const Double_t *ybins
      ,Int_t nbinsz,const Double_t *zbins
      ,Int_t nbinsu,const Double_t *ubins);
  //TH4D(TString path);
  
  virtual Int_t Fill(Double_t){return -1;} // MayNotUse
  virtual Int_t Fill(Double_t, Double_t){return -1;} // MayNotUse
  virtual Int_t Fill(const char*, Double_t){return -1;} // MayNotUse
  virtual Double_t IntegralAndError(Int_t binx1, Int_t binx2, Double_t & err, Option_t *option="") const{return -1;} // MayNotUse
  virtual Double_t Integral(Int_t, Int_t, Option_t*) const{return -1;} // MayNotUse

  virtual Int_t Fill(Double_t x, Double_t y, Double_t z, Double_t u, Double_t w);
  virtual Int_t Fill(Double_t x, Double_t y, Double_t z, Double_t u);
  virtual Double_t IntegralAndError(Int_t binx1, Int_t binx2, Int_t biny1, Int_t biny2, Int_t binz1, Int_t binz2, Int_t binu1, Int_t binu2, Double_t & error, Option_t *option="") const;
  virtual Double_t Integral(Int_t binx1, Int_t binx2, Int_t biny1, Int_t biny2, Int_t binz1, Int_t binz2, Int_t binu1, Int_t binu2, Option_t *option="") const;
  virtual Double_t Integral(Option_t *option="") const;

  virtual const TAxis* GetXaxis() const {return &fXaxis;}
  virtual const TAxis* GetYaxis() const {if(hists.size()) return hists[0]->GetXaxis();else return NULL;}
  virtual const TAxis* GetZaxis() const {if(hists.size()) return hists[0]->GetYaxis();else return NULL;}
  virtual const TAxis* GetUaxis() const {if(hists.size()) return hists[0]->GetZaxis();else return NULL;}

  virtual bool CheckConsistency() const;
  static bool CheckConsistency(const TH4D* h1,const TH4D* h2);

  //virtual Int_t WriteFile(TString path=".");
  //virtual Int_t ReadFile(TString path);

  ClassDef(TH4D,1);
};
  
#endif
