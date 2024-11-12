#include "ZptWeight.h"
ZptWeight::ZptWeight(){
}
ZptWeight::~ZptWeight(){
}
void ZptWeight::executeEvent(){
  ////////////////////////check genlevel//////////////////
  if(IsDYSample){
    if(abs(lhe_l0.ID())!=15){
      Parameter p;
      if(abs(lhe_l0.ID())==11) p=MakeParameter("ee");
      else if(abs(lhe_l1.ID())==13) p=MakeParameter("mm");
      TLorentzVector genZ=(gen_l0+gen_l1);
      FillHist(p.prefix+"gen_mypt_nozptweight",genZ.M(),fabs(genZ.Rapidity()),genZ.Pt(),p.w.lumiweight*p.w.PUweight*p.w.weakweight,massbinnum,massbin,ybinnum,ybin,ptbinnum,ptbin);
      FillHist(p.prefix+"gen_mypt",genZ.M(),fabs(genZ.Rapidity()),genZ.Pt(),p.w.lumiweight*p.w.PUweight*p.w.weakweight*p.w.zptweight,massbinnum,massbin,ybinnum,ybin,ptbinnum,ptbin);
    }
  }

  if(!IsDATA||DataStream.Contains("DoubleMuon")) 
    executeEventWithParameter(MakeParameter("mm"));
  if(!IsDATA||DataStream.Contains("DoubleEG")||DataStream.Contains("EGamma")) 
    executeEventWithParameter(MakeParameter("ee"));
  // if(!IsDATA||DataStream.Contains("SingleMuon")) 
  //   executeEventWithParameter(MakeParameter("mu"));
  // if(!IsDATA||DataStream.Contains("SingleElectron")||DataStream.Contains("EGamma")) 
  //   executeEventWithParameter(MakeParameter("el"));
}
SMPAnalyzerCore::Parameter ZptWeight::MakeParameter(TString key,TString option){
  Parameter p=SMPAnalyzerCore::MakeParameter(key,option);
  p.variationbits|=EfficiencyWeight;
  return p;
}

SMPAnalyzerCore::Variations ZptWeight::MakeVariations(const Parameter& p){
  Variations v;
  TLorentzVector genZ=(gen_l0+gen_l1);
  AddVariationWeight(v,"",p.default_weight);
  if(!IsDATA){
    if(IsDYSample){
      AddVariationWeight(v,"_zptg",p.default_weight/p.w.zptweight*fZptCorrection->GetZptWeight(genZ.Pt()));
      AddVariationWeight(v,"_zptgy",p.default_weight/p.w.zptweight*fZptCorrection->GetZptWeight(genZ.Pt(),genZ.Rapidity()));
      AddVariationWeight(v,"_zptgym",p.default_weight/p.w.zptweight*fZptCorrection->GetZptWeight(genZ.Pt(),genZ.Rapidity(),genZ.M()));
    }
    AddVariationWeight(v,"_nozptweight",p.default_weight/p.w.zptweight);
    if(p.channel[0]=='e'){
      int set=p.w.electronIDSF_sys.size()-1;
      AddVariationWeight(v,"_effresidual",p.default_weight/p.w.electronIDSF*p.w.electronIDSF_sys[set][0]);
      AddVariationWeight(v,"_nozptweight_effresidual",p.default_weight/p.w.zptweight/p.w.electronIDSF*p.w.electronIDSF_sys[set][0]);
    }else if(p.channel[0]=='m'){
      int set=p.w.muonIDSF_sys.size()-1;
      AddVariationWeight(v,"_effresidual",p.default_weight/p.w.muonIDSF*p.w.muonIDSF_sys[set][0]);
      AddVariationWeight(v,"_nozptweight_effresidual",p.default_weight/p.w.zptweight/p.w.muonIDSF*p.w.muonIDSF_sys[set][0]);
    }
    if(IsDYSample){
      // for(unsigned int i=0;i<weight_Scale->size();i++){
      // 	AddVariationWeight(v,Form("_scalevariation%d",i),p.default_weight*weight_Scale->at(i));
      // }
      for(unsigned int i=0;i<weight_Scale->size();i++){
	double scale=1.;
	if(isnormal(weight_Scale->at(i))) scale=weight_Scale->at(i);
	AddVariationWeight(v,Form("_nozptweight_scalevariation%d",i),p.default_weight/p.w.zptweight*scale);
      }
    }
  }
  return v;
}
void ZptWeight::FillHists(Parameter& p){
  TLorentzVector dilepton=(*p.lepton0)+(*p.lepton1);
  double dimass=dilepton.M();
  double dipt=dilepton.Pt();
  double dirap=fabs(dilepton.Rapidity());
  TLorentzVector genZ=(gen_l0+gen_l1);
 
  TString pre=p.prefix+p.hprefix;
  TString suf=p.suffix+p.vsuffix;

  FillHist(pre+"myptgpt"+suf,dimass,dirap,dipt,genZ.Pt(),p.weight,massbinnum,massbin,ybinnum,ybin,ptbinnum,ptbin,ptbinnum,ptbin);

  if(dimass>77&&dimass<106){
    FillHist(pre+"pt"+suf,dipt,p.weight,ptbinnum,ptbin);
    FillHist(pre+"lpt"+suf,p.lepton0->Pt(),p.weight,ptbinnum,ptbin);
    FillHist(pre+"lpt"+suf,p.lepton1->Pt(),p.weight,ptbinnum,ptbin);
    FillHist(pre+"leta"+suf,p.lepton0->Eta(),p.weight,50,-2.5,2.5);
    FillHist(pre+"leta"+suf,p.lepton1->Eta(),p.weight,50,-2.5,2.5);    
  }
}
