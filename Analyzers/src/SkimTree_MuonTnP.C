#include "SkimTree_MuonTnP.h"

SkimTree_MuonTnP::SkimTree_MuonTnP(){
}
SkimTree_MuonTnP::~SkimTree_MuonTnP(){
}
void SkimTree_MuonTnP::initializeAnalyzer(){
  SMPAnalyzerCore::initializeAnalyzer(); //setup zpt roc z0
  outfile->cd();

  newtree=new TTree("Events","Events");

  newtree->Branch("run",&run);
  newtree->Branch("event",&event);
  newtree->Branch("lumi",&lumi);
  if(!IsDATA){
    newtree->Branch("weight",&weight);
    newtree->Branch("PUweight",&PUweight);
    newtree->Branch("PUweight_up",&PUweight_up);
    newtree->Branch("PUweight_down",&PUweight_down);
    newtree->Branch("prefireweight",&prefireweight);
    newtree->Branch("prefireweight_up",&prefireweight_up);
    newtree->Branch("prefireweight_down",&prefireweight_down);
    newtree->Branch("zptweight",&zptweight);
    newtree->Branch("z0weight",&z0weight);
  }

  newtree->Branch("probe_isTracker",&probe_isTracker);
  newtree->Branch("probe_isGlobal",&probe_isGlobal);
  newtree->Branch("probe_isSA",&probe_isSA);
  newtree->Branch("probe_isSA_unique",&probe_isSA_unique);
  newtree->Branch("probe_isTight",&probe_isTight);
  newtree->Branch("probe_isMedium",&probe_isMedium);
  newtree->Branch("probe_isMedium2016a",&probe_isMedium2016a);
  newtree->Branch("probe_TkIsoLoose",&probe_TkIsoLoose);
  newtree->Branch("probe_PFIsoTight",&probe_PFIsoTight);
  newtree->Branch("probe_IsoMu24",&probe_IsoMu24);
  newtree->Branch("probe_IsoMu27",&probe_IsoMu27);
  newtree->Branch("probe_Mu17Leg1",&probe_Mu17Leg1);
  newtree->Branch("probe_Mu8Leg2",&probe_Mu8Leg2);  

  newtree->Branch("probe_pt",&probe_pt);
  newtree->Branch("probe_pt_cor",&probe_pt_cor);
  newtree->Branch("probe_eta",&probe_eta);
  newtree->Branch("probe_phi",&probe_phi);
  newtree->Branch("probe_q",&probe_q);

  newtree->Branch("tag_IsoMu24",&tag_IsoMu24);
  newtree->Branch("tag_IsoMu27",&tag_IsoMu27);
  newtree->Branch("tag_isTight",&tag_isTight);
  newtree->Branch("tag_isMedium",&tag_isMedium);
  newtree->Branch("tag_isMedium2016a",&tag_isMedium2016a);
  newtree->Branch("tag_TkIsoLoose",&tag_TkIsoLoose);
  newtree->Branch("tag_PFIsoTight",&tag_PFIsoTight);

  newtree->Branch("tag_pt",&tag_pt);
  newtree->Branch("tag_pt_cor",&tag_pt_cor);
  newtree->Branch("tag_eta",&tag_eta);
  newtree->Branch("tag_phi",&tag_phi);
  newtree->Branch("tag_q",&tag_q);
  
  newtree->Branch("pair_mass",&pair_mass);
  newtree->Branch("pair_mass_cor",&pair_mass_cor);
  newtree->Branch("pair_pt",&pair_pt);
  newtree->Branch("pair_pt_cor",&pair_pt_cor);
  newtree->Branch("pair_EMTF",&pair_EMTF);
  
  if(!IsDATA){
    newtree->Branch("pair_gen_matched",&pair_gen_matched);
    newtree->Branch("pair_gen_mass",&pair_gen_mass);  
    newtree->Branch("probe_gen_pt",&probe_gen_pt);
    newtree->Branch("probe_gen_eta",&probe_gen_eta);
    newtree->Branch("probe_gen_phi",&probe_gen_phi);
    newtree->Branch("probe_gen_dR",&probe_gen_dR);
    newtree->Branch("probe_gen_reldpt",&probe_gen_reldpt);
  }
  
}
void SkimTree_MuonTnP::executeEvent(){
  if(!IsDATA||DataStream.Contains("SingleMuon")){
    Parameter p=MakeParameter("mu");
    p.triggers={"HLT_IsoMu24_v","HLT_IsoTkMu24_v","HLT_IsoMu27_v"};
    p.SetMuonKeys("Default","Default",{"Default"});
    p.SetMuons(MuonMomentumCorrection(GetAllMuons(),0,0,0));
    executeEventWithParameter(p);
  }
}
bool SkimTree_MuonTnP::PassSelection(Parameter& p){
  if(!PassMETFilter()) return false;
  return true;
}
void SkimTree_MuonTnP::FillHists(Parameter& p){
  map<Muon*,Gen*> genmatching;
  if(!IsDATA){
    for(Gen* gen:{&gen_l0_bare,&gen_l1_bare}){
      vector<Muon*> cands={};
      for(Muon& muon:p.muons)
	if(gen->DeltaR(muon)<0.2) cands.push_back(&muon);
      double mindpt=1000.;
      Muon *matched=NULL;
      for(Muon* cand:cands){
	double dpt=fabs((cand->Pt()-gen->Pt())/gen->Pt());
	if(dpt<mindpt&&genmatching.find(cand)==genmatching.end()){
	  mindpt=dpt;
	  matched=cand;
	}
      }
      if(matched){
	genmatching[matched]=gen;
      }
    }

    weight=p.w.lumiweight;
    PUweight=p.w.PUweight;
    PUweight_up=p.w.PUweight_up;
    PUweight_down=p.w.PUweight_down;
    prefireweight=p.w.prefireweight;
    prefireweight_up=p.w.prefireweight_up;
    prefireweight_down=p.w.prefireweight_down;
    zptweight=p.w.zptweight;
    z0weight=p.w.z0weight;
  }

  for(Muon& tag:p.muons){
    for(Muon& probe:p.muons){
      if(&tag==&probe) continue;
      if((tag+probe).M()<40) continue;

      probe_isTracker=probe.IsType(Muon::Type::TrackerMuon);
      probe_isGlobal=probe.IsType(Muon::Type::GlobalMuon);
      probe_isSA=probe.IsType(Muon::Type::StandAloneMuon);
      probe_isSA_unique=probe_isSA;
      if(probe_isSA){
	for(Muon& dup:p.muons){
	  if(&dup==&probe) continue;
	  if(!dup.IsType(Muon::Type::GlobalMuon)) continue;
	  if(probe.DeltaR(dup)>0.2) continue;
	  if(fabs(probe.Pt()/dup.Pt()-1)>0.3) continue;
	  probe_isSA_unique=false;
	}
      }
	
      probe_isTight=probe.isPOGTight();
      probe_isMedium=probe.isPOGMedium_nohip();
      probe_isMedium2016a=probe.isPOGMedium_hip();
      probe_TkIsoLoose=probe.PassSelector(Muon::Selector::TkIsoLoose);
      probe_PFIsoTight=probe.PassSelector(Muon::Selector::PFIsoTight);
      probe_IsoMu24=PassSLT1(&probe);
      probe_IsoMu27=PassSLT2(&probe);
      probe_Mu17Leg1=PassDLT1(&probe);
      probe_Mu8Leg2=PassDLT2(&probe);
      probe_pt=probe.MiniAODPt();
      probe_pt_cor=probe.Pt();
      probe_eta=probe.Eta();
      probe_phi=probe.Phi();
      probe_q=probe.Charge();
      
      tag_IsoMu24=PassSLT1(&tag);
      tag_IsoMu27=PassSLT2(&tag);
      if(!(tag_IsoMu24||tag_IsoMu27)) continue;
      tag_isTight=tag.isPOGTight();
      tag_isMedium=tag.isPOGMedium_nohip();
      tag_isMedium2016a=tag.isPOGMedium_hip();
      tag_TkIsoLoose=tag.PassSelector(Muon::Selector::TkIsoLoose);
      tag_PFIsoTight=tag.PassSelector(Muon::Selector::PFIsoTight);
      tag_pt=tag.MiniAODPt();
      tag_pt_cor=tag.Pt();
      tag_eta=tag.Eta();
      tag_phi=tag.Phi();
      tag_q=tag.Charge();
	
      TLorentzVector pair_cor=tag+probe;
      TLorentzVector pair=tag.MiniAODPt()/tag.Pt()*tag+probe.MiniAODPt()/probe.Pt()*probe;
      pair_mass=pair.M();
      pair_mass_cor=pair_cor.M();
      pair_pt=pair.Pt();
      pair_pt_cor=pair_cor.Pt();
      pair_EMTF=tag_eta*probe_eta > 0 && fabs(tag_eta) > 0.9 && fabs(probe_eta) > 0.9 && fabs(tag_phi-probe_phi) < 70/180.*3.141592;
      
      if(!IsDATA){
	Gen *probe_gen=NULL,*tag_gen=NULL;
	if(genmatching.find(&probe)!=genmatching.end())
	  probe_gen=genmatching[&probe];
	if(genmatching.find(&tag)!=genmatching.end())
	  tag_gen=genmatching[&tag];
	
	if(probe_gen){
	  probe_gen_pt=probe_gen->Pt();
	  probe_gen_eta=probe_gen->Eta();
	  probe_gen_phi=probe_gen->Phi();
	  probe_gen_dR=probe_gen->DeltaR(probe);
	  probe_gen_reldpt=(probe.Pt()-probe_gen->Pt())/probe_gen->Pt();
	}else{
	  probe_gen_pt=0.;
	  probe_gen_eta=0.;
	  probe_gen_phi=0.;
	  probe_gen_dR=0.;
	  probe_gen_reldpt=0.;
	}
	if(probe_gen&&tag_gen){
	  pair_gen_matched=true;
	  pair_gen_mass=(*probe_gen+*tag_gen).M();
	}else{
	  pair_gen_matched=false;
	  pair_gen_mass=0.;
	}
      }
      newtree->Fill();
    }
  }
}

void SkimTree_MuonTnP::WriteHist(){
  outfile->mkdir("muon");
  outfile->cd("muon");
  newtree->Write();
  outfile->cd();
}
