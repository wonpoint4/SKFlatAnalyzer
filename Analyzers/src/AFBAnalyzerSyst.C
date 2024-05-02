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
      if(gen_region==region){
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
    if(region=="0bjet/"){
      cost=GetCosThetaCS(p.lepton0,p.lepton1);
    }else{
      cost=GetCosThetaRecoil(p.lepton0,p.lepton1,&p.bjets.at(0));
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
    }
    for(int im=0;im<grid_mbinnum;im++){
      if(dimass>=grid_mbin[im]&&dimass<grid_mbin[im+1]){
	if(region=="0bjet/"){
	  FillHist(pre+Form("dirap_m%d",im)+suf,dirap,cost,p.weight,afb_ybinnum,afb_ybin,grid_costbinnum,grid_costbin);
	  FillHist(pre+Form("dipt_m%d",im)+suf,dipt,cost,p.weight,afb_ptbinnum,afb_ptbin,grid_costbinnum,grid_costbin);
	}else{
	  FillHist(pre+Form("dirap_m%d",im)+suf,dirap,cost,p.weight,afb_ybinnum,afb_ybin,grid_costbinnum,grid_costbin);
	  FillHist(pre+Form("dipt_m%d",im)+suf,dipt,cost,p.weight,afb_ptbinnum,afb_ptbin,grid_costbinnum,grid_costbin);
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
  }
}

AFBAnalyzerSyst::AFBAnalyzerSyst(){}
AFBAnalyzerSyst::~AFBAnalyzerSyst(){}

bool AFBAnalyzerSyst::PassSelection(Parameter& p,bool cutflow){
  if(!SMPAnalyzerCore::PassSelection(p,cutflow)) return false;
  if((*p.lepton0+*p.lepton1).M()<52) return false;
  return true;
}
