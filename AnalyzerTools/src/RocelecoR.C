#ifndef ElectroWeakAnalysis_RocelecoR
#define ElectroWeakAnalysis_RocelecoR

#include <fstream>
#include <sstream>
#include <stdexcept>
#include "RocelecoR.h"

RocelecoR::RocelecoR(){}

RocelecoR::RocelecoR(std::string filename){
    init(filename);
}

void RocelecoR::reset(){
    NETA=0;
    NPHI=0;
    std::vector<double>().swap(etabin);
    nset=0;
    std::vector<int>().swap(nmem);
    std::vector<std::vector<RocOne>>().swap(RC);

}


void RocelecoR::init(std::string filename){
    std::ifstream in(filename.c_str());
    if(in.fail()) throw std::invalid_argument("RocelecoR::init could not open file " + filename);

    int RMIN(0), RTRK(0), RETA(0);
    std::vector<double> BETA;

    std::string tag;
    int type, sys, mem, var, bin;	
    std::string s;
    while(std::getline(in, s)){
	std::stringstream ss(s); 
	std::string first4=s.substr(0,4);
	if(first4=="NSET"){
	    ss >> tag >> nset;
	    nmem.resize(nset);
	    tvar.resize(nset);
	    RC.resize(nset);
	}
	else if(first4=="NMEM") {
	    ss >> tag;
	    for(int i=0; i<nset; ++i) {
		ss >> nmem[i];
		RC[i].resize(nmem[i]);
	    }
	}
	else if(first4=="TVAR") {
	    ss >> tag;
	    for(int i=0; i<nset; ++i) ss >> tvar[i];
	}
	else if(first4=="RMIN") ss >> tag >> RMIN;
	else if(first4=="RTRK") ss >> tag >> RTRK;
	else if(first4=="RETA") {
	    ss >> tag >> RETA;
	    BETA.resize(RETA+1);
	    for(auto &h: BETA) ss >> h;

	}
	else if(first4=="CPHI") {
	    ss >> tag >> NPHI; 
	    DPHI=2*CrystalBall::pi/NPHI;
	}
	else if(first4=="CETA")  {
	    ss >> tag >> NETA;
	    etabin.resize(NETA+1);
	    for(auto& h: etabin) ss >> h;
	}
	else{ 
	    ss >> sys >> mem >> tag;
	    auto &rc = RC[sys][mem]; 
	    rc.RR.NETA=RETA;
	    rc.RR.NTRK=RTRK;
	    rc.RR.NMIN=RMIN;
	    auto &resol = rc.RR.resol;
	    if(resol.empty()){
		resol.resize(RETA);
		for(size_t ir=0; ir<resol.size(); ++ir){
		    auto &r = resol[ir];
		    r.eta = BETA[ir];
		    r.cb.resize(RTRK);
		    for(auto i:{0,1})r.nTrk[i].resize(RTRK+1);
		    for(auto i:{0,1,2})r.rsPar[i].resize(RTRK);
		}
	    }
	    auto &cp = rc.CP;
	    for(TYPE T:{MC,DT}){
		if(cp[T].empty()){
		    cp[T].resize(NETA);
		    for(auto &i: cp[T]) i.resize(NPHI);
		}
	    }

	    if(tag=="R"){
		ss >> var >> bin; 
		for(int i=0; i<RTRK; ++i) {
		    switch(var){
			case 0: ss >> resol[bin].rsPar[var][i]; break;
			case 1: ss >> resol[bin].rsPar[var][i]; break;
			case 2: ss >> resol[bin].rsPar[var][i]; resol[bin].rsPar[var][i]/=100; break; 
			case 3: ss >> resol[bin].cb[i].s; break; 
			case 4: ss >> resol[bin].cb[i].a; break; 
			case 5: ss >> resol[bin].cb[i].n; break; 
			default: break;
		    }
		}
	    }
	    else if(tag=="T") {
		ss >> type >> bin; 
		for(int i=0; i<RTRK+1; ++i) ss >> resol[bin].nTrk[type][i];
	    }
	    else if(tag=="F") {
		ss >> type; 
		for(int i=0; i<RETA; ++i) ss >> resol[i].kRes[type];

	    }
	    else if(tag=="C") {
		ss >> type >> var >> bin; 
		for(int i=0; i<NPHI; ++i){
		    auto &x = cp[type][bin][i];
		    if(var==0) { ss >> x.M; /*x.M = 1.0+x.M/100;*/}
		    else if(var==1){ ss >> x.A; /*x.A/=100; */}
		}
	    }
	}
    }

    for(auto &rcs: RC)
	for(auto &rcm: rcs)
	    for(auto &r: rcm.RR.resol)
		for(auto &i: r.cb) i.init();

    in.close();
}

const double RocelecoR::MPHI=-CrystalBall::pi;

int RocelecoR::etaBin(double x) const{
    for(int i=0; i<NETA-1; ++i) if(x<etabin[i+1]) return i;
    return NETA-1;
}

int RocelecoR::phiBin(double x) const{
    int ibin=(x-MPHI)/DPHI;
    if(ibin<0) return 0; 
    if(ibin>=NPHI) return NPHI-1;
    return ibin;
}

double RocelecoR::kScaleDT(int Q, double pt, double eta, double phi, int s, int m) const{
    int H = etaBin(eta);
    int F = phiBin(phi);
    return 2.0/RC[s][m].CP[DT][H][F].M; // + Q*RC[s][m].CP[DT][H][F].A*pt);
}

double RocelecoR::kScaleMC(int Q, double pt, double eta, double phi, int s, int m) const{
    int H = etaBin(eta);
    int F = phiBin(phi);
    return 2.0/RC[s][m].CP[MC][H][F].M;// + Q*RC[s][m].CP[MC][H][F].A*pt);
}

double RocelecoR::kSpreadMC(int Q, double pt, double eta, double phi, double gt, int s, int m) const{
    const auto& rc=RC[s][m];
    int H = etaBin(eta);
    int F = phiBin(phi);
    double k=1.0/(rc.CP[MC][H][F].M + Q*rc.CP[MC][H][F].A*pt);
    return k*rc.RR.kSpread(gt, k*pt, eta);
}

double RocelecoR::kSmearMC(int Q, double pt, double eta, double phi, int n, double u, int s, int m) const{
    const auto& rc=RC[s][m];
    int H = etaBin(eta);
    int F = phiBin(phi);
    double k=1.0/(rc.CP[MC][H][F].M + Q*rc.CP[MC][H][F].A*pt);
    return k*rc.RR.kExtra(k*pt, eta, n, u);
}


double RocelecoR::kScaleFromGenMC(int Q, double pt, double eta, double phi, int n, double gt, double w, int s, int m) const{
    const auto& rc=RC[s][m];
    int H = etaBin(eta);
    int F = phiBin(phi);
    double k=1.0/(rc.CP[MC][H][F].M + Q*rc.CP[MC][H][F].A*pt);
    return k*rc.RR.kSpread(gt, k*pt, eta, n, w);
}

double RocelecoR::kScaleAndSmearMC(int Q, double pt, double eta, double phi, int n, double u, double w, int s, int m) const{
    const auto& rc=RC[s][m];
    int H = etaBin(eta);
    int F = phiBin(phi);
    double k=1.0/(rc.CP[MC][H][F].M + Q*rc.CP[MC][H][F].A*pt);
    return k*rc.RR.kExtra(k*pt, eta, n, u, w);
}

double RocelecoR::kGenSmear(double pt, double eta, double v, double u, RocRes::TYPE TT, int s, int m) const{
    return RC[s][m].RR.kSmear(pt, eta, TT, v, u);
}

template <typename T>
double RocelecoR::error(T f) const{
    double sum=0;
    for(int s=0; s<nset; ++s){
	for(int i=0; i<nmem[s]; ++i) {
	    double d = f(s,i) - f(0,0); 
	    sum += d*d/nmem[s];
	}
    }
    return sqrt(sum);
}

double RocelecoR::kScaleDTerror(int Q, double pt, double eta, double phi) const{
    return error([this, Q, pt, eta, phi](int s, int m) {return kScaleDT(Q, pt, eta, phi, s, m);});
}

double RocelecoR::kSpreadMCerror(int Q, double pt, double eta, double phi, double gt) const{
    return error([this, Q, pt, eta, phi, gt](int s, int m){return kSpreadMC(Q, pt, eta, phi, gt, s, m);});
}

double RocelecoR::kSmearMCerror(int Q, double pt, double eta, double phi, int n, double u) const{
    return error([this, Q, pt, eta, phi, n, u](int s, int m){return kSmearMC(Q, pt, eta, phi, n, u, s, m);});
}

double RocelecoR::kScaleFromGenMCerror(int Q, double pt, double eta, double phi, int n, double gt, double w) const{
    return error([this, Q, pt, eta, phi, n, gt, w](int s, int m) {return kScaleFromGenMC(Q, pt, eta, phi, n, gt, w, s, m);});
}

double RocelecoR::kScaleAndSmearMCerror(int Q, double pt, double eta, double phi, int n, double u, double w) const{
    return error([this, Q, pt, eta, phi, n, u, w](int s, int m) {return kScaleAndSmearMC(Q, pt, eta, phi, n, u, w, s, m);});
}

#endif


