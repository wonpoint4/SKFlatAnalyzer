#include "RoccoR_common.h"

const double CrystalBall::pi = 3.14159;
const double CrystalBall::sqrtPiOver2 = sqrt(CrystalBall::pi/2.0);
const double CrystalBall::sqrt2 = sqrt(2.0);

RocRes::RocRes(){
  reset();
}

void RocRes::reset(){
  NETA=0;
  NTRK=0;
  NMIN=0;
  std::vector<ResParams>().swap(resol);
}

int RocRes::etaBin(double eta) const{
  double abseta=fabs(eta);
  for(int i=0; i<NETA-1; ++i) if(abseta<resol[i+1].eta) return i;
  return NETA-1;
}

int RocRes::trkBin(double x, int h, TYPE T) const{
  for(int i=0; i<NTRK-1; ++i) if(x<resol[h].nTrk[T][i+1]) return i;
  return NTRK-1;
}

double RocRes::Sigma(double pt, int H, int F) const{
  double dpt=pt-45;
  const ResParams &rp = resol[H];
  return rp.rsPar[0][F] + rp.rsPar[1][F]*dpt + rp.rsPar[2][F]*dpt*dpt;
}

double RocRes::rndm(int H, int F, double w) const{
  const ResParams &rp = resol[H];
  return rp.nTrk[MC][F]+(rp.nTrk[MC][F+1]-rp.nTrk[MC][F])*w;
}

double RocRes::kSpread(double gpt, double rpt, double eta, int n, double w) const{
  int H = etaBin(fabs(eta));
  int F = n>NMIN ? n-NMIN : 0;
  double v = rndm(H, F, w);
  int D = trkBin(v, H, Data);
  double kold = gpt / rpt;
  const ResParams &rp = resol[H];
  double u = rp.cb[F].cdf( (kold-1.0)/rp.kRes[MC]/Sigma(gpt,H,F) );
  double knew = 1.0 + rp.kRes[Data]*Sigma(gpt,H,D)*rp.cb[D].invcdf(u);

  if(knew<0) return 1.0;
  return kold/knew;
}


double RocRes::kSpread(double gpt, double rpt, double eta) const{
  int H = etaBin(fabs(eta));
  const auto &k = resol[H].kRes;
  double x = gpt/rpt;
  return x / (1.0 + (x-1.0)*k[Data]/k[MC]);
}

double RocRes::kSmear(double pt, double eta, TYPE type, double v, double u) const{
  int H = etaBin(fabs(eta));
  int F = trkBin(v, H);
  const ResParams &rp = resol[H];
  double x = rp.kRes[type] * Sigma(pt, H, F) * rp.cb[F].invcdf(u);
  return 1.0/(1.0+x);
}

double RocRes::kSmear(double pt, double eta, TYPE type, double w, double u, int n) const{
  int H = etaBin(fabs(eta));
  int F = n-NMIN;
  if(type==Data) F = trkBin(rndm(H, F, w), H, Data);
  const ResParams &rp = resol[H];
  double x = rp.kRes[type] * Sigma(pt, H, F) * rp.cb[F].invcdf(u);
  return 1.0/(1.0+x);
}

double RocRes::kExtra(double pt, double eta, int n, double u, double w) const{
  int H = etaBin(fabs(eta));
  int F = n>NMIN ? n-NMIN : 0;
  const ResParams &rp = resol[H];
  double v = rp.nTrk[MC][F]+(rp.nTrk[MC][F+1]-rp.nTrk[MC][F])*w;
  int D = trkBin(v, H, Data);
  double RD = rp.kRes[Data]*Sigma(pt, H, D);
  double RM = rp.kRes[MC]*Sigma(pt, H, F);
  double x = RD>RM ? sqrt(RD*RD-RM*RM)*rp.cb[F].invcdf(u) : 0;
  if(x<=-1) return 1.0;
  return 1.0/(1.0 + x);
}

double RocRes::kExtra(double pt, double eta, int n, double u) const{
  int H = etaBin(fabs(eta));
  int F = n>NMIN ? n-NMIN : 0;
  const ResParams &rp = resol[H];
  double d = rp.kRes[Data];
  double m = rp.kRes[MC];
  double x = d>m ? sqrt(d*d-m*m) * Sigma(pt, H, F) * rp.cb[F].invcdf(u) : 0;
  if(x<=-1) return 1.0;
  return 1.0/(1.0 + x);
}
