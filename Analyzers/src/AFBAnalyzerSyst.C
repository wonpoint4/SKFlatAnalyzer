#include "AFBAnalyzerSyst.h"

void AFBAnalyzerSyst::initializeAnalyzer(){
  SMPAnalyzerCore::initializeAnalyzer(); //setup zpt roc z0 
  IsSkimmed=GetSkimName()!="" ? true : false;
}
void AFBAnalyzerSyst::executeEvent(){
  //cout<<"Event:"<<event<<endl;
  //// FIXME some events of DYJets has nan PDF weights. I don't know why...
  if(MCSample=="DYJets"&&!isnormal(weight_Scale->at(0))) return;

  ///////////////// GEN level /////////////////////
  genfid_b0=NULL;
  for(auto& lhe:lhes){
    if(lhe.Status()!=1) continue;
    if(abs(lhe.ID())!=5) continue;
    if(lhe.Pt()<40) continue;
    if(fabs(lhe.Eta())>2.4) continue;
    if(genfid_b0&&genfid_b0->Pt()>lhe.Pt()) continue;
    genfid_b0=&lhe;
  }  
  executeEventGen();

  ///////////////// RECO level /////////////////////
  if(!IsDATA||DataStream.Contains("DoubleMuon")){
    executeEventWithParameter(MakeParameter("mm"));
    executeEventWithParameter(MakeParameter("MM","fake"));
  }
  if(!IsDATA||DataStream.Contains("DoubleEG")||DataStream.Contains("EGamma")){
    executeEventWithParameter(MakeParameter("ee"));
    executeEventWithParameter(MakeParameter("EE","fake"));
  }
}
SMPAnalyzerCore::Parameter AFBAnalyzerSyst::MakeParameter(TString key,TString option){
  Parameter p=SMPAnalyzerCore::MakeParameter(key,option);
  p.variationbits=NominalWeight|SystematicWeight|EfficiencyWeight|LeptonCorrection;
  if((IsDYSample&&!p.hprefix.Contains("tau_"))||IsTTLLSample){
    p.variationbits|=PDFWeight;
  }
  return p;
}
void AFBAnalyzerSyst::executeEventGen(){
  if(IsSkimmed) return;
  if(!IsDYSample&&!IsTTLLSample) return;
  Parameter p;
  double letacut=2.4;
  if(abs(lhe_l0.ID())==11&&abs(lhe_l1.ID())==11){
    p=MakeParameter("ee");
    letacut=2.5;
  }else if(abs(lhe_l0.ID())==13&&abs(lhe_l1.ID())==13){
    p=MakeParameter("mm");
  }else return;
  
  bool genfid=false;
  if(TMath::Max(gen_l0_dressed.Pt(),gen_l1_dressed.Pt())>p.c.lepton0pt){
    if(TMath::Min(gen_l0_dressed.Pt(),gen_l1_dressed.Pt())>p.c.lepton1pt){
      if(fabs(gen_l0_dressed.Eta())<letacut&&fabs(gen_l1_dressed.Eta())<letacut){
	genfid=true;
      }
    }
  }

  TString region="0bjet/";
  if(genfid_b0) region="nbjet/";
  
  TLorentzVector dilepton=gen_l0_dressed+gen_l1_dressed;
  TLorentzVector lhe_dilepton=lhe_l0+lhe_l1;
  double dimass=dilepton.M();
  
  double costCS=GetCosThetaCS(&gen_l0_dressed,&gen_l1_dressed);
  double lhe_costCS=GetCosThetaCS(&lhe_l0,&lhe_l1);
  double correct_costCS=-999;
  if(gen_p0.PID()==21||gen_p0.PID()==22){
    if(gen_p1.PID()==21||gen_p1.PID()==22) correct_costCS=GetCosThetaCS(&gen_l0_dressed,&gen_l1_dressed,0);
    else if(gen_p1.PID()>0) correct_costCS=GetCosThetaCS(&gen_l0_dressed,&gen_l1_dressed,-1);
    else if(gen_p1.PID()<0) correct_costCS=GetCosThetaCS(&gen_l0_dressed,&gen_l1_dressed,1);
  }else if(gen_p0.PID()>0){
    if(gen_p1.PID()==21||gen_p1.PID()==22) correct_costCS=GetCosThetaCS(&gen_l0_dressed,&gen_l1_dressed,1);
    else if(gen_p1.PID()>0) correct_costCS=GetCosThetaCS(&gen_l0_dressed,&gen_l1_dressed,0);
    else if(gen_p1.PID()<0) correct_costCS=GetCosThetaCS(&gen_l0_dressed,&gen_l1_dressed,1);
  }else if(gen_p0.PID()<0){
    if(gen_p1.PID()==21||gen_p1.PID()==22) correct_costCS=GetCosThetaCS(&gen_l0_dressed,&gen_l1_dressed,-1);
    else if(gen_p1.PID()>0) correct_costCS=GetCosThetaCS(&gen_l0_dressed,&gen_l1_dressed,-1);
    else if(gen_p1.PID()<0) correct_costCS=GetCosThetaCS(&gen_l0_dressed,&gen_l1_dressed,0);
  }
  if(correct_costCS==-999){
    cout<<"wrong pid for parton: "<<gen_p0.PID()<<" "<<gen_p1.PID()<<endl;
    exit(EXIT_FAILURE);
  }


  map<TString,double> weightmap;
  weightmap[""]=p.w.lumiweight*p.w.zptweight*p.w.weakweight*p.w.topptweight;
  weightmap["_noweight"]=p.w.lumiweight;
  weightmap["_nozptweight"]=p.w.lumiweight*p.w.weakweight*p.w.topptweight;
  weightmap["_noweakweight"]=p.w.lumiweight*p.w.zptweight*p.w.topptweight;
  weightmap["_notopptweight"]=p.w.lumiweight*p.w.weakweight;

  for(auto [vsuffix,weight]:weightmap){
    TString pre=p.prefix+region+p.hprefix;
    TString suf=p.suffix+vsuffix;
    FillHist(pre+"correct_dimassCS"+suf,dimass,correct_costCS,weight,afb_mbinnum,afb_mbin,grid_costbinnum,grid_costbin);
    FillHist(pre+"lhe_dimassCS"+suf,lhe_dilepton.M(),lhe_costCS,weight,afb_mbinnum,afb_mbin,grid_costbinnum,grid_costbin);
    FillHist(pre+"gen_dimassCS"+suf,dimass,costCS,weight,afb_mbinnum,afb_mbin,grid_costbinnum,grid_costbin);
    if(genfid) FillHist(pre+"genfid_dimassCS"+suf,dimass,costCS,weight,afb_mbinnum,afb_mbin,grid_costbinnum,grid_costbin);
    
    if(genfid_b0){
      double costRecoil=GetCosThetaRecoil(&gen_l0_dressed,&gen_l1_dressed,genfid_b0);
      FillHist(pre+"gen_dimassRecoil"+suf,dimass,costRecoil,weight,afb_mbinnum,afb_mbin,grid_costbinnum,grid_costbin);
      if(genfid) FillHist(pre+"genfid_dimassRecoil"+suf,dimass,costRecoil,weight,afb_mbinnum,afb_mbin,grid_costbinnum,grid_costbin);
    }
  }

  // for AN PUweight section
  weightmap[""]=p.w.lumiweight*p.w.zptweight*p.w.weakweight*p.w.topptweight*p.w.PUweight;
  weightmap["_noPUweight"]=p.w.lumiweight*p.w.zptweight*p.w.weakweight*p.w.topptweight;
  weightmap["_PUweight_up"]=p.w.lumiweight*p.w.zptweight*p.w.weakweight*p.w.topptweight*p.w.PUweight_up;
  weightmap["_PUweight_down"]=p.w.lumiweight*p.w.zptweight*p.w.weakweight*p.w.topptweight*p.w.PUweight_down;

  for(auto [vsuffix,weight]:weightmap){
    TString pre=p.prefix+region+p.hprefix;
    TString suf=p.suffix+vsuffix;
    FillHist(pre+"gen_nPU"+suf,nPileUp,weight,100,0,100);    
  }

}

void AFBAnalyzerSyst::FillHistsUnfold(Parameter& preco,Parameter& pgen){
  TString gen_channel="";
  if(abs(lhe_l0.ID())==11&&abs(lhe_l1.ID())==11) gen_channel="ee";
  else if(abs(lhe_l0.ID())==13&&abs(lhe_l1.ID())==13) gen_channel="mm";
  TLorentzVector gen_ll_dressed=gen_l0_dressed+gen_l1_dressed;
  TString gen_region="0bjet/";
  if(genfid_b0) gen_region="nbjet/";

  for(TString region:{"0bjet/","nbjet/"}){
    double genweight=0;
    double genm=-1;
    double geny=-100;
    double genpt=-1;
    double gencost=0;
    if(gen_channel==preco.channel){
      if(gen_region==region||region=="0bjet/"){
	if(TMath::Max(gen_l0_dressed.Pt(),gen_l1_dressed.Pt())>25){
	  if(TMath::Min(gen_l0_dressed.Pt(),gen_l1_dressed.Pt())>15){
	    if(fabs(gen_l0_dressed.Eta())<2.4&&fabs(gen_l1_dressed.Eta())<2.4){
	      if(gen_ll_dressed.M()>52){
		genweight=pgen.weight;
		genm=gen_ll_dressed.M();
		geny=gen_ll_dressed.Rapidity();
		genpt=gen_ll_dressed.Pt();
		if(region=="0bjet/")
		  gencost=GetCosThetaCS(&gen_l0_dressed,&gen_l1_dressed);
		else
		  gencost=GetCosThetaRecoil(&gen_l0_dressed,&gen_l1_dressed,genfid_b0);
	      }
	    }
	  }
	}
      }
    }
    
    double recoweight=0;
    double recom=-1;
    double recoy=-100;
    double recopt=-1;
    double recocost=0;
    if(PassSelection(preco)){
      if(region=="0bjet/"&&(preco.bjets.size()==0||preco.bjets.at(0).Pt()<preco.c.jetpt)){
	recoweight=preco.weight;
	recom=(*preco.lepton0+*preco.lepton1).M();
	recoy=(*preco.lepton0+*preco.lepton1).Rapidity();
	recopt=(*preco.lepton0+*preco.lepton1).Pt();
	recocost=GetCosThetaCS(preco.lepton0,preco.lepton1);
      }else if(region=="nbjet/"&&preco.bjets.size()&&preco.bjets.at(0).Pt()>preco.c.jetpt){
	recoweight=preco.weight;
	recom=(*preco.lepton0+*preco.lepton1).M();
	recoy=(*preco.lepton0+*preco.lepton1).Rapidity();
	recopt=(*preco.lepton0+*preco.lepton1).Pt();
	recocost=GetCosThetaRecoil(preco.lepton0,preco.lepton1,&preco.bjets.at(0));
      }
    }
    TString pre=preco.prefix+region+preco.hprefix;
    TString suf=preco.suffix+preco.vsuffix;
   
    int imbin=GetUnfoldBin(afb_mbinnum,afb_mbin,genm,gencost);
    int jmbin=GetUnfoldBin(afb_mbinnum,afb_mbin,recom,recocost);
    int iybin=GetUnfoldBin(afb_ybinnum,afb_ybin,geny,gencost);
    int jybin=GetUnfoldBin(afb_ybinnum,afb_ybin,recoy,recocost);
    int iptbin=GetUnfoldBin(afb_ptbinnum,afb_ptbin,genpt,gencost);
    int jptbin=GetUnfoldBin(afb_ptbinnum,afb_ptbin,recopt,recocost);
    
    if(imbin||jmbin){
      FillHist(pre+"response_afbm"+suf,imbin,jmbin,recoweight,2*afb_mbinnum,1,2*afb_mbinnum+1,2*afb_mbinnum,1,2*afb_mbinnum+1);
      FillHist(pre+"response_afbm"+suf,imbin,0,genweight-recoweight,2*afb_mbinnum,1,2*afb_mbinnum+1,2*afb_mbinnum,1,2*afb_mbinnum+1);
    }
    
    for(int k=0;k<grid_mbinnum;k++){
      {
	int ibin=(genm>=grid_mbin[k]&&genm<grid_mbin[k+1]) ? iybin : 0;
	int jbin=(recom>=grid_mbin[k]&&recom<grid_mbin[k+1]) ? jybin : 0;
	if(ibin||jbin){
	  FillHist(pre+Form("response_afby_m%d",k)+suf,ibin,jbin,recoweight,2*afb_ybinnum,1,2*afb_ybinnum+1,2*afb_ybinnum,1,2*afb_ybinnum+1);
	  FillHist(pre+Form("response_afby_m%d",k)+suf,ibin,0,genweight-recoweight,2*afb_ybinnum,1,2*afb_ybinnum+1,2*afb_ybinnum,1,2*afb_ybinnum+1);
	}
      }
      {
	int ibin=(genm>=grid_mbin[k]&&genm<grid_mbin[k+1]) ? iptbin : 0;
	int jbin=(recom>=grid_mbin[k]&&recom<grid_mbin[k+1]) ? jptbin : 0;
	if(ibin||jbin){
	  FillHist(pre+Form("response_afbpt_m%d",k)+suf,ibin,jbin,recoweight,2*afb_ptbinnum,1,2*afb_ptbinnum+1,2*afb_ptbinnum,1,2*afb_ptbinnum+1);
	  FillHist(pre+Form("response_afbpt_m%d",k)+suf,ibin,0,genweight-recoweight,2*afb_ptbinnum,1,2*afb_ptbinnum+1,2*afb_ptbinnum,1,2*afb_ptbinnum+1);
	}
      }
    }
  }  
}
AFBAnalyzerSyst::Variations AFBAnalyzerSyst::MakeVariations(const Parameter& p){
  Variations v=SMPAnalyzerCore::MakeVariations(p);
  //EvalVariationsBcharge(p,v);
  return v;
}
void AFBAnalyzerSyst::FillHistsSyst(Parameter p,Variations& vs){
  p.SetLeptons();
  Parameter pgen=p.Clone();
  ResetRecoWeights(pgen);
  Variations genvariations=MakeVariations(pgen);
  for(auto& [vsuf,variation]:vs){
    Apply(p,vsuf,variation);
    Apply(pgen,vsuf,genvariations[vsuf]);
    if(!IsSkimmed){
      FillHistsUnfold(p,pgen);
      continue;
    }

    if(!PassSelection(p,p.vsuffix=="")) continue;
    TLorentzVector dilepton=*p.lepton0+*p.lepton1;
    double dimass=dilepton.M();
    double dirap=dilepton.Rapidity();
    double dipt=dilepton.Pt();
    TString region="0bjet/";
    if(p.bjets.size()&&p.bjets.at(0).Pt()>p.c.jetpt) region="nbjet/";
    double cost=0;
    double costCS=0;
    double costRecoil=0;
    if(region=="0bjet/"){
      costCS=GetCosThetaCS(p.lepton0,p.lepton1);
      cost=costCS;
    }else{
      costCS=GetCosThetaCS(p.lepton0,p.lepton1);
      costRecoil=GetCosThetaRecoil(p.lepton0,p.lepton1,&p.bjets.at(0));
      cost=costRecoil;
    }
    
    TString pre=p.prefix+region+p.hprefix;
    TString suf=p.suffix+p.vsuffix;
    
    if(region=="0bjet/"){
      FillHist(pre+"dimass"+suf,dimass,cost,p.weight,afb_mbinnum,afb_mbin,grid_costbinnum,grid_costbin);
      FillHist(pre+"dirap"+suf,dirap,cost,p.weight,afb_ybinnum,afb_ybin,grid_costbinnum,grid_costbin);
      FillHist(pre+"dipt"+suf,dipt,cost,p.weight,afb_ptbinnum,afb_ptbin,grid_costbinnum,grid_costbin);
      FillHist(pre+"cost"+suf,cost,p.weight,20,-1,1);
    }else{
      FillHist(pre+"dimass"+suf,dimass,cost,p.weight,afb_mbinnum,afb_mbin,grid_costbinnum,grid_costbin);
      FillHist(pre+"dirap"+suf,dirap,cost,p.weight,afb_ybinnum,afb_ybin,grid_costbinnum,grid_costbin);
      FillHist(pre+"dipt"+suf,dipt,cost,p.weight,afb_ptbinnum,afb_ptbin,grid_costbinnum,grid_costbin);
      FillHist(pre+"cost"+suf,cost,p.weight,20,-1,1);      
      if(vsuf==""){
	FillHist(pre+"dimassCS"+suf,dimass,costCS,p.weight,afb_mbinnum,afb_mbin,grid_costbinnum,grid_costbin);
	FillHist(pre+"dirapCS"+suf,dirap,costCS,p.weight,afb_ybinnum,afb_ybin,grid_costbinnum,grid_costbin);
	FillHist(pre+"diptCS"+suf,dipt,costCS,p.weight,afb_ptbinnum,afb_ptbin,grid_costbinnum,grid_costbin);
	FillHist(pre+"costCS"+suf,costCS,p.weight,20,-1,1);      
      }	
    }
    for(int im=0;im<grid_mbinnum;im++){
      if(dimass>=grid_mbin[im]&&dimass<grid_mbin[im+1]){
	if(region=="0bjet/"){
	  FillHist(pre+Form("dirap_m%d",im)+suf,dirap,cost,p.weight,afb_ybinnum,afb_ybin,grid_costbinnum,grid_costbin);
	  FillHist(pre+Form("dipt_m%d",im)+suf,dipt,cost,p.weight,afb_ptbinnum,afb_ptbin,grid_costbinnum,grid_costbin);
	}else{
	  FillHist(pre+Form("dirap_m%d",im)+suf,dirap,cost,p.weight,afb_ybinnum,afb_ybin,grid_costbinnum,grid_costbin);
	  FillHist(pre+Form("dipt_m%d",im)+suf,dipt,cost,p.weight,afb_ptbinnum,afb_ptbin,grid_costbinnum,grid_costbin);
	  if(vsuf==""){
	    FillHist(pre+Form("dirapCS_m%d",im)+suf,dirap,costCS,p.weight,afb_ybinnum,afb_ybin,grid_costbinnum,grid_costbin);
	    FillHist(pre+Form("diptCS_m%d",im)+suf,dipt,costCS,p.weight,afb_ptbinnum,afb_ptbin,grid_costbinnum,grid_costbin);
	  }
	}
      }
    }
    FillHist(pre+"lpt"+suf,p.lepton0->Pt(),p.weight,lptbinnum,lptbin);
    FillHist(pre+"lpt"+suf,p.lepton1->Pt(),p.weight,lptbinnum,lptbin);
    FillHist(pre+"leta"+suf,p.lepton0->Eta(),p.weight,50,-2.5,2.5);
    FillHist(pre+"leta"+suf,p.lepton1->Eta(),p.weight,50,-2.5,2.5);

    int njet=count_if(p.jets.begin(),p.jets.end(),[&p](auto& jet){return jet.Pt()>p.c.jetpt;});
    FillHist(pre+"jets"+suf,njet,p.weight,10,0,10);
    if(njet){
      FillHist(pre+"j0pt"+suf,p.jets.at(0).Pt(),p.weight,100,0,200);
    }
    int nbjet=count_if(p.bjets.begin(),p.bjets.end(),[&p](auto& jet){return jet.Pt()>p.c.jetpt;});
    FillHist(pre+"bjets"+suf,nbjet,p.weight,10,0,10);
    if(nbjet){
      FillHist(pre+"b0pt"+suf,p.bjets.at(0).Pt(),p.weight,100,0,200);
      FillHist(pre+"b0charge"+suf,p.bjets.at(0).userFloat["AFBCharge"],p.weight,100,-5,5);
    }
    if(vsuf==""||vsuf.Contains("z0weight")){
      FillHist(pre+"z0"+suf,vertex_Z,p.weight,120,-15,15);
    }
    if(vsuf==""||vsuf.Contains("PUweight")){
      FillHist(pre+"nPV"+suf,nPV,p.weight,100,0,100);
      FillHist(pre+"rho"+suf,Rho,p.weight,50,0,50);
      FillHist(pre+"met"+suf,pfMET_Type1_pt,p.weight,100,0,200);
    }
  }
}

AFBAnalyzerSyst::AFBAnalyzerSyst(){}
AFBAnalyzerSyst::~AFBAnalyzerSyst(){}

bool AFBAnalyzerSyst::PassSelection(Parameter& p,bool cutflow){
  if(!SMPAnalyzerCore::PassSelection(p,cutflow)) return false;
  if((*p.lepton0+*p.lepton1).M()<52) return false;
  return true;
}
