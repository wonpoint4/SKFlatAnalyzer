#ifndef ElectroWeakAnalysis_RocelecoR_H
#define ElectroWeakAnalysis_RocelecoR_H

#include "RoccoR_common.h"

class RocelecoR{

 private:
  enum TYPE{MC, DT};
  enum TVAR{Default, Replica, Symhes};

  static const double MPHI; 

  int NETA;
  int NPHI; 
  double DPHI;
  std::vector<double> etabin;

  struct CorParams{double M; double A;};

  struct RocOne{
    RocRes RR;
    std::vector<std::vector<CorParams>> CP[2];
  };

  int nset;
  std::vector<int> nmem;
  std::vector<int> tvar;
  std::vector<std::vector<RocOne>> RC;
  int etaBin(double eta) const;
  int phiBin(double phi) const;
  template <typename T> double error(T f) const;

 public:
  RocelecoR(); 
  RocelecoR(std::string filename); 
  void init(std::string filename);
  void reset();

  const RocRes& getRes(int s=0, int m=0) const {return RC[s][m].RR;}
  double getM(int T, int H, int F, int s=0, int m=0) const{return RC[s][m].CP[T][H][F].M;}
  double getA(int T, int H, int F, int s=0, int m=0) const{return RC[s][m].CP[T][H][F].A;}
  double getK(int T, int H, int s=0, int m=0)        const{return RC[s][m].RR.resol[H].kRes[T];}
  double kGenSmear(double pt, double eta, double v, double u, RocRes::TYPE TT=RocRes::Data, int s=0, int m=0) const;
  double kScaleMC(int Q, double pt, double eta, double phi, int s=0, int m=0) const;

  double kScaleDT(int Q, double pt, double eta, double phi, int s=0, int m=0) const;
  double kSpreadMC(int Q, double pt, double eta, double phi, double gt, int s=0, int m=0) const;
  double kSmearMC(int Q, double pt, double eta, double phi, int n, double u, int s=0, int m=0) const;

  double kScaleDTerror(int Q, double pt, double eta, double phi) const;
  double kSpreadMCerror(int Q, double pt, double eta, double phi, double gt) const;
  double kSmearMCerror(int Q, double pt, double eta, double phi, int n, double u) const;

  //old, should only be used with 2017v0
  double kScaleFromGenMC(int Q, double pt, double eta, double phi, int n, double gt, double w, int s=0, int m=0) const; 
  double kScaleAndSmearMC(int Q, double pt, double eta, double phi, int n, double u, double w, int s=0, int m=0) const;  
  double kScaleFromGenMCerror(int Q, double pt, double eta, double phi, int n, double gt, double w) const; 
  double kScaleAndSmearMCerror(int Q, double pt, double eta, double phi, int n, double u, double w) const;  
};

#endif

