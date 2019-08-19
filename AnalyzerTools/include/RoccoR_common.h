#ifndef ElectroWeakAnalysis_RoccoR_common_H
#define ElectroWeakAnalysis_RoccoR_common_H

#include <boost/math/special_functions/erf.hpp>
struct CrystalBall{
  static const double pi;
  static const double sqrtPiOver2;
  static const double sqrt2;

  double m;
  double s;
  double a;
  double n;

  double B;
  double C;
  double D;
  double N;

  double NA;
  double Ns;
  double NC;
  double F;
  double G;
  double k;

  double cdfMa;
  double cdfPa;

CrystalBall():m(0),s(1),a(10),n(10){
  init();
}

  void init(){
    double fa = fabs(a);
    double ex = exp(-fa*fa/2);
    double A  = pow(n/fa, n) * ex;
    double C1 = n/fa/(n-1) * ex;
    double D1 = 2 * sqrtPiOver2 * erf(fa/sqrt2);

    B = n/fa-fa;
    C = (D1+2*C1)/C1;
    D = (D1+2*C1)/2;

    N = 1.0/s/(D1+2*C1);
    k = 1.0/(n-1);

    NA = N*A;
    Ns = N*s;
    NC = Ns*C1;
    F = 1-fa*fa/n;
    G = s*n/fa;

    cdfMa = cdf(m-a*s);
    cdfPa = cdf(m+a*s);
  }

  double pdf(double x) const{
    double d=(x-m)/s;
    if(d<-a) return NA*pow(B-d, -n);
    if(d>a) return NA*pow(B+d, -n);
    return N*exp(-d*d/2);
  }

  double pdf(double x, double ks, double dm) const{
    double d=(x-m-dm)/(s*ks);
    if(d<-a) return NA/ks*pow(B-d, -n);
    if(d>a) return NA/ks*pow(B+d, -n);
    return N/ks*exp(-d*d/2);
  }

  double cdf(double x) const{
    double d = (x-m)/s;
    if(d<-a) return NC / pow(F-s*d/G, n-1);
    if(d>a) return NC * (C - pow(F+s*d/G, 1-n) );
    return Ns * (D - sqrtPiOver2 * erf(-d/sqrt2));
  }

  double invcdf(double u) const{
    if(u<cdfMa) return m + G*(F - pow(NC/u, k));
    if(u>cdfPa) return m - G*(F - pow(C-u/NC, -k) );
    return m - sqrt2 * s * boost::math::erf_inv((D - u/Ns )/sqrtPiOver2);
  }
};


struct RocRes{
  enum TYPE {MC, Data, Extra};

  struct ResParams{
    double eta;
    double kRes[2];
    std::vector<double> nTrk[2];
    std::vector<double> rsPar[3];
    std::vector<CrystalBall> cb;
  ResParams():eta(0){for(auto& k: kRes) k=1;}
  };

  int NETA;
  int NTRK;
  int NMIN;

  std::vector<ResParams> resol;

  RocRes();

  int etaBin(double x) const;
  int trkBin(double x, int h, TYPE T=MC) const;
  void reset();

  double rndm(int H, int F, double v) const;
  double Sigma(double pt, int H, int F) const;
  double kSpread(double gpt, double rpt, double eta, int nlayers, double w) const;
  double kSpread(double gpt, double rpt, double eta) const;
  double kSmear(double pt, double eta, TYPE type, double v, double u) const;
  double kSmear(double pt, double eta, TYPE type, double v, double u, int n) const;
  double kExtra(double pt, double eta, int nlayers, double u, double w) const;
  double kExtra(double pt, double eta, int nlayers, double u) const;
};
#endif
