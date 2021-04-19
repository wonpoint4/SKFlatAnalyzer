#include "AFBAnalyzer.h"

void AFBAnalyzer::initializeAnalyzer(){
  SMPAnalyzerCore::initializeAnalyzer(); //setup zpt roc z0 PUJet 
  //SetupToy(100);
  SetupCosThetaWeight();
  IsNominalRun=!HasFlag("SYS")&&!HasFlag("PDFSYS");
  
  vector<JetTagging::Parameters> jtps={JetTagging::Parameters(JetTagging::DeepCSV,JetTagging::Tight,JetTagging::mujets,JetTagging::mujets),
				       JetTagging::Parameters(JetTagging::DeepCSV,JetTagging::Loose,JetTagging::mujets,JetTagging::mujets),
				       JetTagging::Parameters(JetTagging::DeepJet,JetTagging::Tight,JetTagging::mujets,JetTagging::mujets),
				       JetTagging::Parameters(JetTagging::DeepJet,JetTagging::Loose,JetTagging::mujets,JetTagging::mujets)};
  mcCorr->SetJetTaggingParameters(jtps);

  if(fChain->GetListOfFiles()->GetEntries()){
    TString filename=fChain->GetListOfFiles()->At(0)->GetTitle();
    if(filename.Contains("SkimTree_")) IsSkimmed=true;
    else IsSkimmed=false;
  }else{
    cout<<"[AFBAnalyzer::initializeAnalyzer] no input file"<<endl;
    exit(EXIT_FAILURE);
  }
}
void AFBAnalyzer::executeEvent(){
  //GetToyWeight();
  costhetaweight=1.;
  costhetaweight_up=1.;
  costhetaweight_down=1.;

  Muons.clear();  Electrons.clear();  jets.clear();
  Muons=SMPGetMuons("POGMediumWithLooseTrkIso",0.0,2.4);
  Electrons=SMPGetElectrons("passMediumID",0.0,2.4);
  jets=GetJets("tightLepVeto",20,2.4);
  std::sort(jets.begin(),jets.end(),PtComparing);
  pujetweight = GetPUJetWeight(jets, 0);

  realjets.clear();
  for(unsigned int l=0; l<jets.size();l++){
    //jets passing PileUp ID Medium
    if(jets.at(l).Pt()<30){
      if(jets.at(l).PileupJetId() <0.18) continue;
    }else if(jets.at(l).Pt()<50){
      if(jets.at(l).PileupJetId() <0.61) continue;
    }
    //jets spliting from (sub)leading electrons, and muons - jet cleaning
    if(Muons.size()>0 && Muons.at(0).Pt()>20 && jets.at(l).DeltaR(Muons.at(0)) <0.4) continue;
    if(Muons.size()>1 && Muons.at(1).Pt()>10 && jets.at(l).DeltaR(Muons.at(1)) <0.4) continue;
    if(Electrons.size()>0 && Electrons.at(0).Pt()>25 && jets.at(l).DeltaR(Electrons.at(0)) <0.4) continue;
    if(Electrons.size()>1 && Electrons.at(1).Pt()>15 && jets.at(l).DeltaR(Electrons.at(1)) <0.4) continue;

    realjets.push_back(jets.at(l));
  }

  vector <TString> btag = {"", "DeepJET"};
  for(unsigned int b=0; b<btag.size(); b++){
    //bjet setting
    n_bjet=0;
    n_powerbjet=0;
    bjet_charge=-3.5;
    bjets.clear();
    JetTagging::Parameters jtpT = JetTagging::Parameters(JetTagging::DeepCSV,JetTagging::Tight,JetTagging::mujets,JetTagging::mujets);
    JetTagging::Parameters jtpL = JetTagging::Parameters(JetTagging::DeepCSV,JetTagging::Loose,JetTagging::mujets,JetTagging::mujets);
    if(btag.at(b)=="DeepJET"){
      jtpT = JetTagging::Parameters(JetTagging::DeepJet,JetTagging::Tight,JetTagging::mujets,JetTagging::mujets);
      jtpL = JetTagging::Parameters(JetTagging::DeepJet,JetTagging::Loose,JetTagging::mujets,JetTagging::mujets);
    }

    /*
      for(const auto& jet:realjets){
      // For leading b, ID:Tight, pT>30
      if(mcCorr->IsBTagged_2a(jtp,jet) && jet.Pt()>30){
      //n_bjet++;
      bjets.push_back(jet);
      }
      // For 2nd b-veto, ID:Loose, pT>20
      else if(mcCorr->IsBTagged_2a(jtp2,jet)){
      n_bjet++;
      if(jet.Pt()>30) n_powerbjet++;
      }
      }
    */

    for(const auto& jet:realjets){
      // For leading b, ID:Tight, pT>30
      if(jet.GetTaggerResult(jtpT.j_Tagger) > mcCorr->GetJetTaggingCutValue(jtpT.j_Tagger, jtpT.j_WP) && jet.Pt()>30){
	bjets.push_back(jet);
      }
      // For 2nd b-veto, ID:Loose, pT>20
      else if(jet.GetTaggerResult(jtpL.j_Tagger) > mcCorr->GetJetTaggingCutValue(jtpL.j_Tagger, jtpL.j_WP)){
	n_bjet++;
	if(jet.Pt()>30) n_powerbjet++;
      }
    }
    btagweight = GetBTaggingReweight_1a_2WP(realjets, jtpT, jtpL, "central");

    softmus.clear();
    softmus=GetMuons("POGLoose",0,2.4);
    std::sort(softmus.begin(),softmus.end(),PtComparing);
    bmuon.clear();
    for(unsigned int l=0; l<softmus.size(); l++){
      if(bjets.size()==0) break;
      if(softmus.at(l).TrkIso()/softmus.at(l).Pt() <0.1) continue;
      if(softmus.at(l).P()*sin(softmus.at(l).Angle(bjets.at(0).Vect())) <1.0) continue;
      if(abs(softmus.at(l).IP3D())/softmus.at(l).IP3Derr() <2.5) continue;
      if(bjets.at(0).DeltaR(softmus.at(l))<0.4) bmuon.push_back(softmus.at(l));
    }

    if(IsDYSample && abs(lhe_l0.ID())!=15){
      //////////////////////// Check LHE /////////////////////////
      Parameter p;
      double letacut=2.4;
      if(abs(lhe_l0.ID())==11){
	p=MakeParameter("ee"+btag.at(b));
	p.c.lepton0pt=25;
	p.c.lepton1pt=15;
      }else if(abs(lhe_l0.ID())==13){
	p=MakeParameter("mm"+btag.at(b));
	p.c.lepton0pt=20;
	p.c.lepton1pt=10;
      }else{
	cout<<"[AFBAnalyzer::executeEvent()] something is wrong l0.ID="<<abs(lhe_l0.ID())<<endl;
	vector<LHE> lhes=GetLHEs();
	for(auto& lhe:lhes) lhe.Print();
	exit(EXIT_FAILURE);
      }

      //////////////////////// GEN /////////////////////////
      TLorentzVector gen_Z=gen_l0+gen_l1;
      double gen_Zmass=gen_Z.M();
      double gen_Zrap=gen_Z.Rapidity();
      double gen_Zpt=gen_Z.Pt();
      double gen_cost_correct=-999;
      if(gen_p0.PID()==21){
	if(gen_p1.PID()==21) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,0);
	else if(gen_p1.PID()>0) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,-1);
	else if(gen_p1.PID()<0) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,1);
      }else if(gen_p0.PID()>0){
	if(gen_p1.PID()==21) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,1);
	else if(gen_p1.PID()>0) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,0);
	else if(gen_p1.PID()<0) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,1);
      }else if(gen_p0.PID()<0){
	if(gen_p1.PID()==21) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,-1);
	else if(gen_p1.PID()>0) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,-1);
	else if(gen_p1.PID()<0) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,0);
      }
      if(gen_cost_correct==-999){
	cout<<"wrong pid for parton"<<endl;
	exit(EXIT_FAILURE);
      }
      costhetaweight=GetCosThetaWeight(gen_Zmass,gen_Zpt,gen_cost_correct,"_pdg");
      costhetaweight_up=GetCosThetaWeight(gen_Zmass,gen_Zpt,gen_cost_correct,"_up");
      costhetaweight_down=GetCosThetaWeight(gen_Zmass,gen_Zpt,gen_cost_correct,"_down");
      
      map<TString,double> map_weight;
      map_weight[""]=p.w.lumiweight;
      //map_weight[""]=p.w.lumiweight*p.w.zptweight*costhetaweight;
      //map_weight["_noweight"]=p.w.lumiweight;
      //map_weight["_nozptweight"]=p.w.lumiweight*costhetaweight;
      //map_weight["_nocosthetaweight"]=p.w.lumiweight*p.w.zptweight;

      /// Matching Study
      if(bjets.size()>0 && bjets.at(0).DeltaR(gen_j0) <0.3){
	//if(bmuon.size() >0){
	  vector<Gen> genmu;
	  genmu.clear();
	  for(unsigned int k=0; k<gens.size(); k++){
	    if(gens.at(k).Status() != 1) continue;
	    if(gens.at(k).MotherIndex() < 0) continue;
	    int motherIndex = gens.at(k).MotherIndex();
	    int motherPID = gens.at(motherIndex).PID();
	    while(abs(motherPID) == 13){
	      if(gens.at(motherIndex).MotherIndex() <0) break;
	      motherIndex = gens.at(motherIndex).MotherIndex();
	      motherPID = gens.at(motherIndex).PID();
	    }
	    if(motherPID == 23) continue;
	    if((gens.at(k).DeltaR(gen_j0) <0.4 || gens.at(k).DeltaR(bjets.at(0)) <0.4) && abs(gens.at(k).PID())==13) genmu.push_back(gens.at(k));
	  }
	  std::sort(genmu.begin(),genmu.end(),PtComparing);
	  FillHist(p.prefix+p.hprefix+"genmuNum",genmu.size(),map_weight,7,0,7);

	  /*for(unsigned int l=0; l<bmuon.size(); l++){
	    TString charge = "";
	    if(bmuon.at(l).Charge() > 0) charge = "P";
	    else charge = "M";
	    FillHist(p.prefix+p.hprefix+"bmuon"+Form("%d",l)+charge+"_PFiso",bmuon.at(l).RelIso(),map_weight,50,0,1);
	    FillHist(p.prefix+p.hprefix+"bmuon"+Form("%d",l)+charge+"_Trkiso",bmuon.at(l).TrkIso()/bmuon.at(l).Pt(),map_weight,50,0,1);
	    FillHist(p.prefix+p.hprefix+"bmuon"+Form("%d",l)+"_pT",(-1)*bmuon.at(l).Charge()*bmuon.at(l).Pt(),map_weight,80,-100,100);
	    FillHist(p.prefix+p.hprefix+"bmuon"+Form("%d",l)+"_p",(-1)*bmuon.at(l).Charge()*bmuon.at(l).P(),map_weight,80,-100,100);
            FillHist(p.prefix+p.hprefix+"bmuon"+Form("%d",l)+"_pT2",(-1)*bmuon.at(l).Charge()*bmuon.at(l).P()*sin(bmuon.at(l).DeltaR(bjets.at(0))),map_weight,80,-10,10);
	    FillHist(p.prefix+p.hprefix+"bmuon"+Form("%d",l)+"_pT3",(-1)*bmuon.at(l).Charge()*bmuon.at(l).P()*sin(bmuon.at(l).Angle(bjets.at(0).Vect())),map_weight,80,-10,10);
            FillHist(p.prefix+p.hprefix+"bmuon"+Form("%d",l)+"_pT23",bmuon.at(l).P()*(sin(bmuon.at(l).DeltaR(bjets.at(0)))-sin(bmuon.at(l).Angle(bjets.at(0).Vect()))),map_weight,200,-10,10);
	    FillHist(p.prefix+p.hprefix+"bmuon"+Form("%d",l)+"_pTRatio",(-1)*bmuon.at(l).Charge()*(bmuon.at(l).Pt())/(bjets.at(0).Pt()),map_weight,40,-1,1);
	    //FillHist(p.prefix+p.hprefix+"bmuon"+Form("%d",l)+"_charge",bmuon.at(l).Charge(),map_weight,4,-2,2);
	    FillHist(p.prefix+p.hprefix+"bmuon"+Form("%d",l)+"_deltaR",(-1)*bmuon.at(l).Charge()*bmuon.at(l).DeltaR(bjets.at(0)),map_weight,40,-1,1);
	    FillHist(p.prefix+p.hprefix+"bmuon"+Form("%d",l)+"_IP3D",(-1)*bmuon.at(l).Charge()*abs(bmuon.at(l).IP3D()),map_weight,200,-5,5);
	    FillHist(p.prefix+p.hprefix+"bmuon"+Form("%d",l)+"_SIP3D",(-1)*bmuon.at(l).Charge()*abs(bmuon.at(l).IP3D())/bmuon.at(l).IP3Derr(),map_weight,40,-20,20);
	    FillHist(p.prefix+p.hprefix+"bmuon"+Form("%d",l)+charge+"_PassMedium",bmuon.at(l).PassID("POGMedium"),map_weight,2,0,2);
	    FillHist(p.prefix+p.hprefix+"bmuon"+Form("%d",l)+charge+"_PassTight",bmuon.at(l).PassID("POGTight"),map_weight,2,0,2);
	    if(l==0) FillHist(p.prefix+p.hprefix+"bmuon"+Form("%d",l)+charge+"_jetcharge",bjets.at(0).Charge(),map_weight,200,-2,2);
	    */
	  if(genmu.size() !=0){
	      for(unsigned int k=0; k<genmu.size(); k++){
		int motherflavor = 0;
		int motherIndex = genmu.at(k).MotherIndex();
		int motherPID = gens.at(motherIndex).PID();
		while(abs(motherPID) == 13){
		  if(gens.at(motherIndex).MotherIndex() <0) break;
		  motherIndex = gens.at(motherIndex).MotherIndex();
		  motherPID = gens.at(motherIndex).PID();
		}
		if(abs(motherPID) < 100 ) motherflavor = abs(motherPID);
		else if(abs(motherPID) > 1000 && abs(motherPID) < 10000) motherflavor = abs(motherPID)/1000;
		else motherflavor = (abs(motherPID) % 1000)/100;
		unsigned int l=0;
		if(l==0) FillHist(p.prefix+p.hprefix+"genmu"+Form("%d",k)+"_deltaR_genj0",genmu.at(k).DeltaR(gen_j0)*genmu.at(k).PID()/13,map_weight,40,-1,1);
		if(l==0) FillHist(p.prefix+p.hprefix+"genmu"+Form("%d",k)+"_deltaR_bjet",genmu.at(k).DeltaR(bjets.at(0))*genmu.at(k).PID()/13,map_weight,40,-1,1);
		if(l==0) FillHist(p.prefix+p.hprefix+"genmu"+Form("%d",k)+"_mother_PID",motherPID,map_weight,1200,-600,600);
		if(l==0) FillHist(p.prefix+p.hprefix+"genmu"+Form("%d",k)+"_PID",genmu.at(k).PID(),map_weight,50,-25,25);
		if(l==0) FillHist(p.prefix+p.hprefix+"genmu"+Form("%d",k)+"_mother_flavor",motherflavor*genmu.at(k).PID()/13,map_weight,100,-50,50);
		if(l==0) FillHist(p.prefix+p.hprefix+"genmu"+Form("%d",k)+"_pT",genmu.at(k).Pt()*genmu.at(k).PID()/13,map_weight,80,-100,100);
		if(l==0) FillHist(p.prefix+p.hprefix+"genmu"+Form("%d",k)+"_p",genmu.at(k).P()*genmu.at(k).PID()/13,map_weight,80,-100,100);
		if(l==0) FillHist(p.prefix+p.hprefix+"genmu"+Form("%d",k)+"_pT2",genmu.at(k).P()*sin(genmu.at(k).DeltaR(bjets.at(0)))*genmu.at(k).PID()/13,map_weight,80,-10,10);
                if(l==0) FillHist(p.prefix+p.hprefix+"genmu"+Form("%d",k)+"_pT3",genmu.at(k).P()*sin(genmu.at(k).Angle(bjets.at(0).Vect()))*genmu.at(k).PID()/13,map_weight,80,-10,10);
		if(l==0) FillHist(p.prefix+p.hprefix+"genmu"+Form("%d",k)+"_pT23",genmu.at(k).P()*(sin(genmu.at(k).DeltaR(bjets.at(0)))-sin(genmu.at(k).Angle(bjets.at(0).Vect()))),map_weight,200,-10,10);
		if(l==0) FillHist(p.prefix+p.hprefix+"genmu"+Form("%d",k)+"_pT4",genmu.at(k).P()*sin(genmu.at(k).DeltaR(gen_j0))*genmu.at(k).PID()/13,map_weight,80,-10,10);
                if(l==0) FillHist(p.prefix+p.hprefix+"genmu"+Form("%d",k)+"_pT5",genmu.at(k).P()*sin(genmu.at(k).Angle(gen_j0.Vect()))*genmu.at(k).PID()/13,map_weight,80,-10,10);
                if(l==0) FillHist(p.prefix+p.hprefix+"genmu"+Form("%d",k)+"_pT45",genmu.at(k).P()*(sin(genmu.at(k).DeltaR(gen_j0))-sin(genmu.at(k).Angle(gen_j0.Vect()))),map_weight,200,-10,10);
		if(l==0) FillHist(p.prefix+p.hprefix+"genmu"+Form("%d",k)+"_pTRatio",genmu.at(k).Pt()/bjets.at(0).Pt()*genmu.at(k).PID()/13,map_weight,20,-1,1);
	      }
	    }
      }
      //FillHist(p.prefix+p.hprefix+"bmuon"+Form("%d",l)+charge+"_genmu"+Form("%d",k)+"_dCharge",bmuon.at(l).Charge()+(genmu.at(k).PID()/13),map_weight,6,-3,3);
      //FillHist(p.prefix+p.hprefix+"bmuon"+Form("%d",l)+charge+"_genmu"+Form("%d",k)+"_deltaR",genmu.at(k).DeltaR(bmuon.at(l)),map_weight,20,0,1);
      //FillHist(p.prefix+p.hprefix+"bmuon"+Form("%d",l)+charge+"_genmu"+Form("%d",k)+"_dpT",genmu.at(k).Pt()-bmuon.at(l).Pt(),map_weight,100,-50,50);
      //FillHist(p.prefix+p.hprefix+"bmuon"+Form("%d",l)+charge+"_genmu"+Form("%d",k)+"_dp",genmu.at(k).P()-bmuon.at(l).P(),map_weight,100,-50,50);

      //////////////// Fill LHE,Gen hists //////////////////////
      if(IsNominalRun){
	vector<TString> accp = {""};//, "accep_"};
	for(unsigned int k=0;k<accp.size();k++){
	  if(accp.at(k)=="accep_"){
	    if(abs(lhe_l0.Rapidity())>2.4 || abs(lhe_l1.Rapidity())>2.4 || abs(lhe_j0.Rapidity())>2.4) continue;
	  }
	  FillHist(p.prefix+p.hprefix+accp.at(k)+"lhe_jetID",lhe_j0.ID(),map_weight,50,-25,25);
	  FillHist(p.prefix+p.hprefix+accp.at(k)+"lhe_jetpT",lhe_j0.Pt(),map_weight,200,0,1000);
	  FillHist(p.prefix+p.hprefix+accp.at(k)+"lhe_jeteta",lhe_j0.Eta(),map_weight,200,-10,10);
	  FillHist(p.prefix+p.hprefix+accp.at(k)+"lhe_jet_Zdy",(lhe_l0+lhe_l1).Rapidity()-lhe_j0.Rapidity(),map_weight,200,-10,10);
	  FillHist(p.prefix+p.hprefix+accp.at(k)+"lhe_jet_Zdphi",(lhe_l0+lhe_l1).DeltaPhi(lhe_j0),map_weight,100,-5,5);
	  FillHist(p.prefix+p.hprefix+accp.at(k)+"lhe_jet_Zy",(lhe_l0+lhe_l1).Rapidity(),map_weight,200,-10,10);
          FillHist(p.prefix+p.hprefix+accp.at(k)+"lhe_jet_Zby",(lhe_l0+lhe_l1+lhe_j0).Rapidity(),map_weight,200,-10,10);
	  FillHist(p.prefix+p.hprefix+accp.at(k)+"lhe_jet_Zy_2D",(lhe_l0+lhe_l1).Rapidity(),lhe_j0.Rapidity(),map_weight,200,-10,10,200,-10,10);
	  /*
	  FillHist(p.prefix+p.hprefix+accp.at(k)+"lhe_jet_Zdytest0",((lhe_l0+lhe_l1).Rapidity()-lhe_j0.Rapidity())>0?1:0,map_weight,10,-5,5);
	  FillHist(p.prefix+p.hprefix+accp.at(k)+"lhe_jet_Zdytest1",(((lhe_l0+lhe_l1).Rapidity()-lhe_j0.Rapidity()) * ((lhe_l0+lhe_l1).Rapidity()+lhe_j0.Rapidity()))>0?1:0,map_weight,10,-5,5);
	  FillHist(p.prefix+p.hprefix+accp.at(k)+"lhe_jet_Zdytest2",((lhe_l0+lhe_l1).Rapidity()*lhe_j0.Rapidity())>0?1:0,map_weight,10,-5,5);
	  FillHist(p.prefix+p.hprefix+accp.at(k)+"lhe_jet_Zdytest3",(lhe_j0.Rapidity())>0?1:0,map_weight,10,-5,5);
	  FillHist(p.prefix+p.hprefix+accp.at(k)+"lhe_jet_Zdytest4",(((lhe_l0+lhe_l1).Rapidity()-lhe_j0.Rapidity()) * (lhe_j0.Rapidity()))>0?1:0,map_weight,10,-5,5);
	  FillHist(p.prefix+p.hprefix+accp.at(k)+"lhe_jet_Zdytest5",(((lhe_l0+lhe_l1).Rapidity()-lhe_j0.Rapidity()) * ((lhe_l0+lhe_l1).Rapidity()+lhe_j0.Rapidity()) * (lhe_l0+lhe_l1).Rapidity()*lhe_j0.Rapidity())>0?1:0,map_weight,10,-5,5);
	  FillHist(p.prefix+p.hprefix+accp.at(k)+"lhe_jet_Zdytest6",(0.9*(lhe_l0+lhe_l1).Rapidity()-lhe_j0.Rapidity())>0?1:0,map_weight,10,-5,5);
	  FillHist(p.prefix+p.hprefix+accp.at(k)+"lhe_jet_Zdytest7",(0.8*(lhe_l0+lhe_l1).Rapidity()-lhe_j0.Rapidity())>0?1:0,map_weight,10,-5,5);
	  FillHist(p.prefix+p.hprefix+accp.at(k)+"lhe_jet_Zdytest8",(0.7*(lhe_l0+lhe_l1).Rapidity()-lhe_j0.Rapidity())>0?1:0,map_weight,10,-5,5);
	  FillHist(p.prefix+p.hprefix+accp.at(k)+"lhe_jet_Zdytest9",(0.6*(lhe_l0+lhe_l1).Rapidity()-lhe_j0.Rapidity())>0?1:0,map_weight,10,-5,5);
	  FillHist(p.prefix+p.hprefix+accp.at(k)+"lhe_jet_Zdytest10",(0.5*(lhe_l0+lhe_l1).Rapidity()-lhe_j0.Rapidity())>0?1:0,map_weight,10,-5,5);
	  FillHist(p.prefix+p.hprefix+accp.at(k)+"lhe_jet_Zdytest11",(0.4*(lhe_l0+lhe_l1).Rapidity()-lhe_j0.Rapidity())>0?1:0,map_weight,10,-5,5);
	  FillHist(p.prefix+p.hprefix+accp.at(k)+"gen_jet_Zdytest0",((gen_l0+gen_l1).Rapidity()-gen_j0.Rapidity())>0?1:0,map_weight,10,-5,5);
	  FillHist(p.prefix+p.hprefix+accp.at(k)+"gen_jet_Zdytest1",(((gen_l0+gen_l1).Rapidity()-gen_j0.Rapidity()) * ((gen_l0+gen_l1).Rapidity()+gen_j0.Rapidity()))>0?1:0,map_weight,10,-5,5);
	  FillHist(p.prefix+p.hprefix+accp.at(k)+"gen_jet_Zdytest2",((gen_l0+gen_l1).Rapidity()*gen_j0.Rapidity())>0?1:0,map_weight,10,-5,5);
	  FillHist(p.prefix+p.hprefix+accp.at(k)+"gen_jet_Zdytest3",(gen_j0.Rapidity())>0?1:0,map_weight,10,-5,5);
	  FillHist(p.prefix+p.hprefix+accp.at(k)+"gen_jet_Zdytest4",(((gen_l0+gen_l1).Rapidity()-gen_j0.Rapidity()) * (gen_j0.Rapidity()))>0?1:0,map_weight,10,-5,5);
	  FillHist(p.prefix+p.hprefix+accp.at(k)+"gen_jet_Zdytest5",(((gen_l0+gen_l1).Rapidity()-gen_j0.Rapidity()) * ((gen_l0+gen_l1).Rapidity()+gen_j0.Rapidity()) * (gen_l0+gen_l1).Rapidity()*gen_j0.Rapidity())>0?1:0,map_weight,10,-5,5);
	  FillHist(p.prefix+p.hprefix+accp.at(k)+"gen_jet_Zdytest6",(0.9*(gen_l0+gen_l1).Rapidity()-gen_j0.Rapidity())>0?1:0,map_weight,10,-5,5);
	  FillHist(p.prefix+p.hprefix+accp.at(k)+"gen_jet_Zdytest7",(0.8*(gen_l0+gen_l1).Rapidity()-gen_j0.Rapidity())>0?1:0,map_weight,10,-5,5);
	  FillHist(p.prefix+p.hprefix+accp.at(k)+"gen_jet_Zdytest8",(0.7*(gen_l0+gen_l1).Rapidity()-gen_j0.Rapidity())>0?1:0,map_weight,10,-5,5);
	  FillHist(p.prefix+p.hprefix+accp.at(k)+"gen_jet_Zdytest9",(0.6*(gen_l0+gen_l1).Rapidity()-gen_j0.Rapidity())>0?1:0,map_weight,10,-5,5);
	  FillHist(p.prefix+p.hprefix+accp.at(k)+"gen_jet_Zdytest10",(0.5*(gen_l0+gen_l1).Rapidity()-gen_j0.Rapidity())>0?1:0,map_weight,10,-5,5);
	  FillHist(p.prefix+p.hprefix+accp.at(k)+"gen_jet_Zdytest11",(0.4*(gen_l0+gen_l1).Rapidity()-gen_j0.Rapidity())>0?1:0,map_weight,10,-5,5);
	  */

	  FillHist(p.prefix+p.hprefix+accp.at(k)+"gen_jetID",gen_j0.PID(),map_weight,50,-25,25);
	  FillHist(p.prefix+p.hprefix+accp.at(k)+"gen_jetpT",gen_j0.Pt(),map_weight,200,0,1000);
	  FillHist(p.prefix+p.hprefix+accp.at(k)+"gen_jeteta",gen_j0.Eta(),map_weight,200,-10,10);
	  FillHist(p.prefix+p.hprefix+accp.at(k)+"gen_jet_Zdy",(gen_l0+gen_l1).Rapidity()-gen_j0.Rapidity(),map_weight,200,-10,10);
	  FillHist(p.prefix+p.hprefix+accp.at(k)+"gen_jet_Zdphi",(gen_l0+gen_l1).DeltaPhi(gen_j0),map_weight,100,-5,5);
	  FillHist(p.prefix+p.hprefix+accp.at(k)+"gen_jet_Zy",(gen_l0+gen_l1).Rapidity(),map_weight,200,-10,10);
          FillHist(p.prefix+p.hprefix+accp.at(k)+"gen_jet_Zby",(gen_l0+gen_l1+gen_j0).Rapidity(),map_weight,200,-10,10);
	  //FillHist(p.prefix+p.hprefix+"gen_jet_Zdy",((gen_l0+gen_l1).Rapidity()-gen_j0.Rapidity())>0?1:0,map_weight,10,-5,5);
	  FillHist(p.prefix+p.hprefix+accp.at(k)+"gen_jet_Zy_2D",(gen_l0+gen_l1).Rapidity(),gen_j0.Rapidity(),map_weight,200,-10,10,200,-10,10);
	  FillHist(p.prefix+p.hprefix+accp.at(k)+"lhe_gen_dpT",lhe_j0.Pt()-gen_j0.Pt(),map_weight,200,-50,50);
	  FillHist(p.prefix+p.hprefix+accp.at(k)+"lhe_gen_dR",lhe_j0.DeltaR(gen_j0),map_weight,100,0,5);
	  FillHist(p.prefix+p.hprefix+accp.at(k)+"lhe_gen_dID",lhe_j0.ID()-gen_j0.PID(),map_weight,100,-50,50);

	  if(bjets.size()>0){
	    FillHist(p.prefix+p.hprefix+accp.at(k)+"lhe_bjetdpT",lhe_j0.Pt()-bjets.at(0).Pt(),map_weight,200,-100,100);
	    FillHist(p.prefix+p.hprefix+accp.at(k)+"lhe_bjetdR",lhe_j0.DeltaR(bjets.at(0)),map_weight,100,0,5);
	    FillHist(p.prefix+p.hprefix+accp.at(k)+"gen_bjetdpT",gen_j0.Pt()-bjets.at(0).Pt(),map_weight,200,-100,100);
	    FillHist(p.prefix+p.hprefix+accp.at(k)+"gen_bjetdR",gen_j0.DeltaR(bjets.at(0)),map_weight,100,0,5);
	  }
	}

	if(gen_j0.PID()==0 || gen_j0.PID()==21) bjet_charge *= gRandom->Rndm()>0.5?1:-1;
	else bjet_charge *= gen_j0.PID()>0?1:-1; // bjet_charge is -3.5 when lhe_j0 is quark(b), and is +3.5 when lhe_j0 is anti-quark(bbar) ->This makes difference in direction of GetCosThetaCS(Recoil)

	bjet = lhe_j0;
	FillHistsAFB(p.prefix,"lhe_","",(Particle*)&lhe_l0,(Particle*)&lhe_l1,map_weight);
	bjet = gen_j0;
	FillHistsAFB(p.prefix,"gen_","",(Particle*)&gen_l0,(Particle*)&gen_l1,map_weight);
	//FillHistsAFB(p.prefix,"gen_","_dressed",(Particle*)&gen_l0_dressed,(Particle*)&gen_l1_dressed,map_weight);
	//FillHistsAFB(p.prefix,"gen_","_bare",(Particle*)&gen_l0_bare,(Particle*)&gen_l1_bare,map_weight);
	if(gen_l0.Pt()>p.c.lepton0pt&&gen_l1.Pt()>p.c.lepton1pt&&fabs(gen_l0.Eta())<letacut&&fabs(gen_l1.Eta())<letacut){
	  FillHistsAFB(p.prefix,"genfid_","",(Particle*)&gen_l0,(Particle*)&gen_l1,map_weight);
	}
	//if(gen_l0_dressed.Pt()>p.c.lepton0pt&&gen_l1_dressed.Pt()>p.c.lepton1pt&&fabs(gen_l0_dressed.Eta())<letacut&&fabs(gen_l1_dressed.Eta())<letacut){
	  //FillHistsAFB(p.prefix,"genfid_","_dressed",(Particle*)&gen_l0_dressed,(Particle*)&gen_l1_dressed,map_weight);
	//}
	//if(gen_l0_bare.Pt()>p.c.lepton0pt&&gen_l1_bare.Pt()>p.c.lepton1pt&&fabs(gen_l0_bare.Eta())<letacut&&fabs(gen_l1_bare.Eta())<letacut){
	  //FillHistsAFB(p.prefix,"genfid_","_bare",(Particle*)&gen_l0_bare,(Particle*)&gen_l1_bare,map_weight);
	//}
	FillHist(p.prefix+"gen_costhetaCS_correct",gen_Zmass,gen_Zrap,gen_Zpt,gen_cost_correct,map_weight,afb_mbinnum,(double*)afb_mbin,afb_ybinnum,(double*)afb_ybin,afb_ptbinnum,(double*)afb_ptbin,20,-1,1);
      }
    }

    ///////////////// RECO level /////////////////////
    if(!IsDATA||DataStream.Contains("SingleMuon")){
      //executeEventWithParameter(MakeParameter("me"));
      //executeEventWithParameter(MakeParameter("mM"));
      //executeEventWithParameter(MakeParameter("mu"));
    }
    if(!IsDATA||DataStream.Contains("DoubleMuon")){
      executeEventWithParameter(MakeParameter("mm"+btag.at(b)));
      //executeEventWithParameter(MakeParameter("MM"));
    }
    if(!IsDATA||DataStream.Contains("SingleElectron")||DataStream.Contains("EGamma")){
      //executeEventWithParameter(MakeParameter("em"));
      //executeEventWithParameter(MakeParameter("eE"));
      //executeEventWithParameter(MakeParameter("el"));
    }
    if(!IsDATA||DataStream.Contains("DoubleEG")||DataStream.Contains("EGamma")){
      executeEventWithParameter(MakeParameter("ee"+btag.at(b)));
      //executeEventWithParameter(MakeParameter("EE"));
    }
    //executeEventWithChannelName(prefix+"NoTrig"+Form("%d",DataYear)); //Not require firing trigger (For Data, just triggers in Stream)
    //if(event.PassTrigger(mmtrigger)){
    //  if(!IsDATA||DataStream.Contains("DoubleMuon")) executeEventWithChannelName(prefix+"DilepTrig"+Form("%d",DataYear));
    //}
    //if(event.PassTrigger(eetrigger)&&!event.PassTrigger(mmtrigger)){
    //  if(!IsDATA||DataStream.Contains("DoubleEG")||DataStream.Contains("EGamma")) executeEventWithChannelName(prefix+"DilepTrig"+Form("%d",DataYear));
    //}
    //if(event.PassTrigger(emtrigger)&&!(event.PassTrigger(mmtrigger) || event.PassTrigger(eetrigger))){
    //  if(!IsDATA||DataStream.Contains("MuonEG")) executeEventWithChannelName(prefix+"DilepTrig"+Form("%d",DataYear));
    //}
  }
}
SMPAnalyzerCore::Parameter AFBAnalyzer::MakeParameter(TString key){
  Parameter p=SMPAnalyzerCore::MakeParameter(key);

  if(IsNominalRun) p.weightbit|=NominalWeight;
  if(HasFlag("SYS")&&!IsDATA) p.weightbit|=SystematicWeight;
  if(HasFlag("PDFSYS")&&!IsDATA) p.weightbit|=PDFWeight;

  if(HasFlag("DYbStudy")) p.prefix+="DYbStudy/";
  if(HasFlag("highmet")) p.prefix+="highmet/";

  return p;
}
bool AFBAnalyzer::PassSelection(Parameter& p){
  if(p.prefix.Contains("highmet")){
    if(pfMET_Type1_pt<60) return false;
    if(IsNominalRun) FillCutflow(p.prefix+p.hprefix+"cutflow","METCut",p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight);
  }
  /*
  int n_bjet=0;
  if(p.prefix.Contains("bjet")||p.prefix.Contains("nobjet")){
    std::vector<Jet> jets=GetJets("tightLepVeto",30,2.7);
    std::sort(jets.begin(),jets.end(),PtComparing);

    JetTagging::Parameters jtp = JetTagging::Parameters(JetTagging::DeepCSV,JetTagging::Medium,JetTagging::mujets,JetTagging::mujets);
    for(const auto& jet:jets)
      if(mcCorr->IsBTagged_2a(jtp,jet))
	n_bjet++;
    
    if(p.prefix.Contains("bjet")&&!n_bjet) return false;
    if(p.prefix.Contains("nobjet")&&n_bjet) return false;
    if(IsNominalRun) FillCutflow(p.prefix+p.hprefix+"cutflow","BJetCut",p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight);
    }*/

  if(!SMPAnalyzerCore::PassSelection(p)) return false;  
  return true;
}

void AFBAnalyzer::FillHists(Parameter& p){
  ///////////////////////map_weight//////////////////
  map<TString,double> map_weight;
  if(p.weightbit&NominalWeight){
    map_weight[""]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF*p.w.CFSF*pujetweight*btagweight;
    map_weight["_nozpt"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF*p.w.CFSF*pujetweight*btagweight;
    map_weight["_nopujet"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF*p.w.CFSF*btagweight;
    map_weight["_nobtag"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF*p.w.CFSF*pujetweight;
    //map_weight["_noCFSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF;
    //map_weight["_CFSF_up"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF*p.w.CFSF_up;
    //map_weight["_CFSF_down"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF*p.w.CFSF_down;
  }
  if(p.weightbit&SystematicWeight){
    map_weight["_noPUweight"]=p.w.lumiweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF;
    map_weight["_PUweight_up"]=p.w.lumiweight*p.w.PUweight_up*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF;
    map_weight["_PUweight_down"]=p.w.lumiweight*p.w.PUweight_down*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF;
    
    map_weight["_noprefireweight"]=p.w.lumiweight*p.w.PUweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF;
    map_weight["_prefireweight_up"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight_up*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF;
    map_weight["_prefireweight_down"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight_down*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF;
    
    map_weight["_nozptweight"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF;

    map_weight["_noz0weight"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF;

    map_weight["_nocosthetaweight"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF;
    map_weight["_costhetaweight_up"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight_up*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF;
    map_weight["_costhetaweight_down"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight_down*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF;
    
    map_weight["_noefficiencySF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight;
    
    map_weight["_noRECOSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.IDSF*p.w.ISOSF*p.w.triggerSF;
    map_weight["_RECOSF_up"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF_up*p.w.IDSF*p.w.ISOSF*p.w.triggerSF;
    map_weight["_RECOSF_down"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF_down*p.w.IDSF*p.w.ISOSF*p.w.triggerSF;
    
    map_weight["_noIDSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.ISOSF*p.w.triggerSF;
    map_weight["_IDSF_up"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF_up*p.w.ISOSF*p.w.triggerSF;
    map_weight["_IDSF_down"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF_down*p.w.ISOSF*p.w.triggerSF;
    
    map_weight["_noISOSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.triggerSF;
    map_weight["_ISOSF_up"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF_up*p.w.triggerSF;
    map_weight["_ISOSF_down"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF_down*p.w.triggerSF;
    
    map_weight["_notriggerSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF;
    map_weight["_triggerSF_up"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF_up;
    map_weight["_triggerSF_down"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF_down;
    
  }
  if(p.weightbit&PDFWeight){
    for(unsigned int i=0;i<PDFWeights_Scale->size();i++){
      map_weight[Form("_scalevariation%d",i)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF*PDFWeights_Scale->at(i);
    }
    for(unsigned int i=0;i<PDFWeights_Error->size();i++){
      map_weight[Form("_pdf%d",i)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF*PDFWeights_Error->at(i);
    }
    if(PDFWeights_AlphaS->size()==2){
      map_weight["_alphaS_up"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF*PDFWeights_AlphaS->at(0);
      map_weight["_alphaS_down"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*costhetaweight*p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF*PDFWeights_AlphaS->at(1);
    }
  }
 
  ///////////////////////DY+b Specific Event Selection///////////////////////
  double eventweight = p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.zptweight*costhetaweight;

  if(p.weightbit&NominalWeight) FillHist(p.prefix+p.hprefix+"nlepton"+p.suffix,(p.muons.size()+p.electrons.size()),eventweight,10,0,10);
  if(p.weightbit&NominalWeight) FillHist(p.prefix+p.hprefix+"nleptonMu"+p.suffix,p.muons.size(),eventweight,10,0,10);
  if(p.weightbit&NominalWeight) FillHist(p.prefix+p.hprefix+"nleptonEl"+p.suffix,p.electrons.size(),eventweight,10,0,10);

  eventweight *= p.w.RECOSF*p.w.IDSF*p.w.ISOSF*p.w.triggerSF*pujetweight*btagweight;
  FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"PUJetSF_bef",eventweight/pujetweight/btagweight);
  FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"PUJetSF",eventweight/btagweight);
  FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"btagSF",eventweight);

  TLorentzVector *l0,*l1;
  l0=p.lepton0;
  l1=p.lepton1;

  if((*l0+*l1).M() <52 ) return;
  //if(p.weightbit&NominalWeight) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"mass52Cut",eventweight);

  FillHist(p.prefix+p.hprefix+"OneTightb_bef"+p.suffix,bjets.size(),eventweight,5,0,5);
  FillHist(p.prefix+p.hprefix+"OneTightb_bef_jet"+p.suffix,realjets.size(),eventweight,10,0,10);
  FillHist(p.prefix+p.hprefix+"OneTightb_bef_normjet"+p.suffix,jets.size(),eventweight,10,0,10);
  FillHist(p.prefix+p.hprefix+"OneTightb_bef_jet_nobSF"+p.suffix,realjets.size(),eventweight/btagweight,10,0,10);
  FillHist(p.prefix+p.hprefix+"OneTightb_bef_normjet_nobSF"+p.suffix,jets.size(),eventweight/btagweight,10,0,10);
  if(bjets.size() !=1) return;
  if(p.weightbit&NominalWeight) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"oneTightb",eventweight);
  bjet_charge = bjets.at(0).Charge();
  bjet = bjets.at(0);

  FillHist(p.prefix+p.hprefix+"2bveto_bef"+p.suffix,n_bjet,eventweight,5,0,5);
  FillHist(p.prefix+p.hprefix+"2bveto_bef_powerbjet"+p.suffix,n_powerbjet,eventweight,5,0,5);
  if(n_bjet !=0) return;
  if(p.weightbit&NominalWeight) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"2bveto",eventweight);

  FillHist(p.prefix+p.hprefix+"2jveto_bef"+p.suffix,realjets.size(),eventweight,7,0,7);
  FillHist(p.prefix+p.hprefix+"2jveto_bef_normjet"+p.suffix,jets.size(),eventweight,10,0,10);
  int n_powerjet = 0;
  for(unsigned int k=0; k<realjets.size(); k++){
    if(realjets.at(k).Pt() > 30) n_powerjet += 1;
  }
  FillHist(p.prefix+p.hprefix+"2jveto_bef_powerjet"+p.suffix,n_powerjet,eventweight,7,0,7);
  //if(realjets.size() >1) return;
  if(n_powerjet >1) return;
  if(p.weightbit&NominalWeight) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"2jveto",eventweight);

  FillHist(p.prefix+p.hprefix+"MET90Cut_bef"+p.suffix,pfMET_Type1_pt,eventweight,150,0,150);
  if(pfMET_Type1_pt >90) return;
  if(p.weightbit&NominalWeight) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"MET90Cut",eventweight);

  FillHist(p.prefix+p.hprefix+"ZbdPhiCut_bef"+p.suffix,(*l0+*l1).DeltaPhi(bjets.at(0)),eventweight,140,-3.5,3.5);
  if(abs((*l0+*l1).DeltaPhi(bjets.at(0))) <2) return;
  if(p.weightbit&NominalWeight) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"ZbdPhiCut",eventweight);

  FillHist(p.prefix+p.hprefix+"ZbpT50Cut_bef"+p.suffix,(*l0+*l1+bjets.at(0)).Pt(),eventweight,200,0,200);
  if((*l0+*l1+bjets.at(0)).Pt() >50) return;
  if(p.weightbit&NominalWeight) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"ZbpT50Cut",eventweight);

  FillHist(p.prefix+p.hprefix+"ZpT20Cut_bef"+p.suffix,(*l0+*l1).Pt(),eventweight,200,0,200);
  if((*l0+*l1).Pt() <20) return;
  if(p.weightbit&NominalWeight) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"ZpT20Cut",eventweight);

  if(bmuon.size() >0){
    for(unsigned int l=0; l<bmuon.size(); l++){
      TString charge = "";
      if(bmuon.at(l).Charge() > 0) charge = "P";
      else charge = "M";
      FillHist(p.prefix+p.hprefix+"bmuon"+Form("%d",l)+charge+"_PFiso",bmuon.at(l).RelIso(),eventweight,50,0,1);
      FillHist(p.prefix+p.hprefix+"bmuon"+Form("%d",l)+charge+"_Trkiso",bmuon.at(l).TrkIso()/bmuon.at(l).Pt(),eventweight,50,0,1);
      FillHist(p.prefix+p.hprefix+"bmuon"+Form("%d",l)+charge+"_pT",bmuon.at(l).Pt(),eventweight,50,0,100);
      FillHist(p.prefix+p.hprefix+"bmuon"+Form("%d",l)+charge+"_p",bmuon.at(l).P(),eventweight,50,0,100);
      FillHist(p.prefix+p.hprefix+"bmuon"+Form("%d",l)+charge+"_pT2",bmuon.at(l).P()*sin(bmuon.at(l).DeltaR(bjets.at(0))),eventweight,50,0,10);
      FillHist(p.prefix+p.hprefix+"bmuon"+Form("%d",l)+charge+"_pT3",bmuon.at(l).P()*sin(bmuon.at(l).Angle(bjets.at(0).Vect())),eventweight,100,0,10);
      FillHist(p.prefix+p.hprefix+"bmuon"+Form("%d",l)+charge+"_pT4",bmuon.at(l).P()*bmuon.at(l).DeltaR(bjets.at(0)),eventweight,50,0,10);
      FillHist(p.prefix+p.hprefix+"bmuon"+Form("%d",l)+charge+"_pTRatio",bmuon.at(l).Pt()/bjets.at(0).Pt(),eventweight,20,0,1);
      //FillHist(p.prefix+p.hprefix+"bmuon"+Form("%d",l)+"_charge",bmuon.at(l).Charge(),eventweight,4,-2,2);
      FillHist(p.prefix+p.hprefix+"bmuon"+Form("%d",l)+charge+"_deltaRlep",min(p.lepton0->DeltaR(bmuon.at(l)),p.lepton1->DeltaR(bmuon.at(l))),eventweight,40,0,2);
      FillHist(p.prefix+p.hprefix+"bmuon"+Form("%d",l)+charge+"_deltaR",bmuon.at(l).DeltaR(bjets.at(0)),eventweight,40,0,2);
      FillHist(p.prefix+p.hprefix+"bmuon"+Form("%d",l)+charge+"_PassMedium",bmuon.at(l).PassID("POGMedium"),eventweight,2,0,2);
      FillHist(p.prefix+p.hprefix+"bmuon"+Form("%d",l)+charge+"_PassTight",bmuon.at(l).PassID("POGTight"),eventweight,2,0,2);
      FillHist(p.prefix+p.hprefix+"bmuon"+Form("%d",l)+charge+"_IP3D",abs(bmuon.at(l).IP3D()),eventweight,100,0,5);
      FillHist(p.prefix+p.hprefix+"bmuon"+Form("%d",l)+charge+"_SIP3D",abs(bmuon.at(l).IP3D())/bmuon.at(l).IP3Derr(),eventweight,40,0,20);
      TVector3 bboost = bjets.at(0).BoostVector();
      bmuon.at(l).Boost(-bboost);
      FillHist(p.prefix+p.hprefix+"bmuon"+Form("%d",l)+charge+"_p2",bmuon.at(l).P(),eventweight,100,0,10);
      if(bmuon.at(l).P() >1.3){
	bmuon.at(l).Boost(bboost);
	FillHist(p.prefix+p.hprefix+"bmuon"+Form("%d",l)+charge+"_pT3_p2_13cut",bmuon.at(l).P()*sin(bmuon.at(l).Angle(bjets.at(0).Vect())),eventweight,100,0,10);
      }
      else bmuon.at(l).Boost(bboost);
      if(bmuon.at(l).P()*sin(bmuon.at(l).Angle(bjets.at(0).Vect())) >0.9){
	FillHist(p.prefix+p.hprefix+"bmuon"+Form("%d",l)+charge+"_pT3_pT3_09cut",bmuon.at(l).P()*sin(bmuon.at(l).Angle(bjets.at(0).Vect())),eventweight,50,0,5);
	bmuon.at(l).Boost(-bboost);
	FillHist(p.prefix+p.hprefix+"bmuon"+Form("%d",l)+charge+"_p2_pT3_09cut",bmuon.at(l).P(),eventweight,50,0,5);
	if(bmuon.at(l).P() >1.0){
	  FillHist(p.prefix+p.hprefix+"bmuon"+Form("%d",l)+charge+"_p2_pT3_09_p2_10cut",bmuon.at(l).P(),eventweight,50,0,5);
	  bmuon.at(l).Boost(bboost);
	  FillHist(p.prefix+p.hprefix+"bmuon"+Form("%d",l)+charge+"_pT3_pT3_09_p2_10cut",bmuon.at(l).P()*sin(bmuon.at(l).Angle(bjets.at(0).Vect())),eventweight,50,0,5);
	}
      }
      if(l==0) FillHist(p.prefix+p.hprefix+"bmuon"+Form("%d",l)+charge+"_jetcharge",bjets.at(0).Charge(),eventweight,200,-2,2);
    }
  }

  FillHist(p.prefix+p.hprefix+"bjetCharge030Cut_bef",bjet_charge,eventweight,100,-1,1);
  if(bmuon.size()!=0){
    //p.suffix += "_SemiLeptonic";
    if(bmuon.at(0).Charge()>0) bjet_charge += 2.;
    else if(bmuon.at(0).Charge()<0) bjet_charge -= 2.;
    else bjet_charge = 3.5;
  }
  FillHist(p.prefix+p.hprefix+"bjetCharge030Cut_bef2",bjet_charge,eventweight,300,-3,3);
  if((bjet_charge < -1.0 || bjet_charge > 1.0)) FillHist(p.prefix+p.hprefix+"bjetCharge030Cut_bef3",bjet_charge,eventweight,300,-3,3);
  else FillHist(p.prefix+p.hprefix+"bjetCharge030Cut_bef4",bjet_charge,eventweight,300,-3,3);
  if(bjet_charge > -0.3 && bjet_charge < 0.3) return;
  if(p.weightbit&NominalWeight) FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"bjetcharge030Cut",eventweight);

  /*
    if(p.weightbit&NominalWeight){
    FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"RECO",eventweight*RECOSF);
    FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"ID",eventweight*RECOSF*IDSF);
    FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"ISO",eventweight*RECOSF*IDSF*ISOSF);
    FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"trigger",eventweight*RECOSF*IDSF*ISOSF*triggerSF);
    FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"PUjet",eventweight*RECOSF*IDSF*ISOSF*triggerSF*pujetweight);
    FillCutflow(p.prefix+p.hprefix+"cutflow"+p.suffix,"BTagging",eventweight*RECOSF*IDSF*ISOSF*triggerSF*pujetweight*btagweight);
    }
  */

  ///////////////////////fill hists///////////////////////
  if(HasFlag("TOY")) FillHistsToy(p.prefix,p.hprefix,p.suffix,(Particle*)p.lepton0,(Particle*)p.lepton1,map_weight);
  else{
    FillHistsAFB(p.prefix,p.hprefix,p.suffix,(Particle*)p.lepton0,(Particle*)p.lepton1,map_weight);

    FillHist(p.prefix+p.hprefix+"bjetpT",bjets.at(0).Pt(),map_weight,200,0,1000);
    FillHist(p.prefix+p.hprefix+"bjeteta",bjets.at(0).Eta(),map_weight,100,-5,5);
    FillHist(p.prefix+p.hprefix+"bjetM",bjets.at(0).M(),map_weight,200,0,20);
    FillHist(p.prefix+p.hprefix+"bjetCharge",bjet_charge,map_weight,400,-4,4);
    FillHist(p.prefix+p.hprefix+"bjetPUID",bjets.at(0).PileupJetId(),map_weight,200,-2,2);

    FillHist(p.prefix+p.hprefix+"yZ",(*l0+*l1).Rapidity(),map_weight,60,-3,3);
    FillHist(p.prefix+p.hprefix+"yb",bjets.at(0).Rapidity(),map_weight,60,-3,3);
    FillHist(p.prefix+p.hprefix+"Zb_y",(*l0+*l1+bjets.at(0)).Rapidity(),map_weight,100,-5,5);
    FillHist(p.prefix+p.hprefix+"Zb_pT",(*l0+*l1+bjets.at(0)).Pt(),map_weight,100,0,100);
    FillHist(p.prefix+p.hprefix+"Zb_dy",(*l0+*l1).Rapidity()-bjets.at(0).Rapidity(),map_weight,100,-5,5);
    FillHist(p.prefix+p.hprefix+"Zb_dphi",(*l0+*l1).DeltaPhi(bjets.at(0)),map_weight,100,-5,5);
    FillHist(p.prefix+p.hprefix+"Zb_y2D",(*l0+*l1).Rapidity(),bjets.at(0).Rapidity(),map_weight,30,-3,3,30,-3,3);

    FillHist(p.prefix+p.hprefix+"mll",(*l0+*l1).M(),map_weight,60,60,120);
    FillHist(p.prefix+p.hprefix+"pTll",(*l0+*l1).Pt(),map_weight,100,0,100);
    FillHist(p.prefix+p.hprefix+"yll",(*l0+*l1).Rapidity(),map_weight,50,-2.5,2.5);
    /*
      for(unsigned int j=0;j<jets.size();j++){
      FillHist(p.prefix+p.hprefix+"jet_"+Form("%i",j)+"_pT",jets.at(j).Pt(),map_weight,200,0,1000);
      FillHist(p.prefix+p.hprefix+"jet_"+Form("%i",j)+"_dpT",jets.at(j).Pt()-bjets.at(0).Pt(),map_weight,200,-100,100);
      FillHist(p.prefix+p.hprefix+"jet_"+Form("%i",j)+"_PUID",jets.at(j).GetPileupJetId(),map_weight,200,-2,2);
      }
      for(unsigned int j=0;j<realjets.size();j++){
      FillHist(p.prefix+p.hprefix+"realjet_"+Form("%i",j)+"_pT",realjets.at(j).Pt(),map_weight,200,0,1000);
      FillHist(p.prefix+p.hprefix+"realjet_"+Form("%i",j)+"_dpT",realjets.at(j).Pt()-bjets.at(0).Pt(),map_weight,200,-100,100);
      FillHist(p.prefix+p.hprefix+"realjet_"+Form("%i",j)+"_PUID",realjets.at(j).GetPileupJetId(),map_weight,200,-2,2);
      }*/

    //if(IsDYSample&&p.hprefix==""&&IsNominalRun){
    if(IsDYSample&&IsNominalRun){
      vector<Gen> gens=GetGens();
      /*
      if(p.lepton0->LeptonFlavour()!=Lepton::ELECTRON && p.lepton0->LeptonFlavour()!=Lepton::MUON){
	PrintGens(gens);
	if(p.lepton0) cout<<"p.lepton0 true"<<endl;
	else cout<<"p.lepton0 false"<<endl;
	if(p.lepton1) cout<<"p.lepton1 true"<<endl;
	else cout<<"p.lepton1 false"<<endl;
	cout<<p.lepton0->Pt()<<", "<<p.lepton1->Pt()<<endl;
	cout<<"p.leptons size : "<<p.leptons.size()<<", p.muons size : "<<p.muons.size()<<", p.electrons size : "<<p.electrons.size()<<endl;
	cout<<"p.c.lepton0 pt : "<<p.c.lepton0pt<<", p.c.lepton1 pt : "<<p.c.lepton1pt<<", p.c.muon0 pt : "<<p.c.muon0pt<<", p.c.muon1 pt : "<<p.c.muon1pt<<endl;
	cout<<"p.leptons.at(0) pt : "<<p.leptons.at(0)->Pt()<<", p.leptons.at(1) pt : "<<p.leptons.at(1)->Pt()<<", p.muons.at(0) pt : "<<p.muons.at(0).Pt()<<", p.muons.at(1) pt : "<<p.muons.at(1).Pt()<<endl;
	cout<<"&p.muons.at(0) pt : "<<(&p.muons.at(0))->Pt()<<", &p.muons.at(1) pt : "<<(&p.muons.at(1))->Pt()<<endl;
	if(p.leptons.at(0)->Pt() > p.c.lepton0pt) cout<<"p.leptons.at(0) pt > p.c.lepton0 pt"<<endl;
	else if(p.leptons.at(0)->Pt() > p.c.lepton0pt) cout<<"p.leptons.at(0) pt < p.c.lepton0 pt"<<endl;
	else if(p.leptons.at(0)->Pt() == p.c.lepton0pt) cout<<"p.leptons.at(0) pt == p.c.lepton0 pt"<<endl;
	else cout<<"p.leptons.at(0) pt ?? p.c.lepton0 pt"<<endl;
	cout<<p.channel<<", "<<p.prefix<<", "<<p.hprefix<<", "<<p.suffix<<endl;
	p.leptons = {};
	p.leptons.push_back(&p.muons.at(0));
	p.leptons.push_back(&p.muons.at(1));
	cout<<"After Reset p.leptons. p.leptons size : "<<p.leptons.size()<<", p.muons size : "<<p.muons.size()<<", p.electrons size : "<<p.electrons.size()<<endl;
	cout<<"p.leptons.at(0) pt : "<<p.leptons.at(0)->Pt()<<", p.leptons.at(1) pt : "<<p.leptons.at(1)->Pt()<<endl;
      }
      if(p.lepton1->LeptonFlavour()!=Lepton::ELECTRON && p.lepton1->LeptonFlavour()!=Lepton::MUON){
	PrintGens(gens);
	cout<<p.lepton0<<", "<<p.lepton1;
      }
      */
      Gen truth_l0=GetGenMatchedLepton(*p.lepton0,gens);
      Gen truth_l1=GetGenMatchedLepton(*p.lepton1,gens);
      if(!truth_l0.IsEmpty()&&!truth_l1.IsEmpty()) 
	FillHistsAFB(p.prefix,"truth_",p.suffix,(Particle*)&truth_l0,(Particle*)&truth_l1,map_weight);
      //else cout<<"no matching"<<endl;
    }
  }
}

AFBAnalyzer::AFBAnalyzer(){}
AFBAnalyzer::~AFBAnalyzer(){
  //DeleteToy();
  DeleteCosThetaWeight();
}
double AFBAnalyzer::GetCosThetaCS(const Particle *p0,const Particle *p1,int direction){
  const TLorentzVector *l0,*l1;
  if(p0->Charge()<0&&p1->Charge()>0){
    l0=p0;
    l1=p1;
  }else if(p0->Charge()>0&&p1->Charge()<0){
    l0=p1;
    l1=p0;
  }else if(strcmp(p0->ClassName(),"LHE")==0){ 
    if(((LHE*)p0)->ID()>0&&((LHE*)p1)->ID()<0){
      l0=p0;
      l1=p1;
    }else if(((LHE*)p0)->ID()<0&&((LHE*)p1)->ID()>0){
      l0=p1;
      l1=p0;
    }else{
      if(gRandom->Rndm()<0.5){
	l0=p0;
	l1=p1;
      }else{
	l0=p1;
	l1=p0;
      }      
    } 
  }else{
    if(gRandom->Rndm()<0.5){
      l0=p0;
      l1=p1;
    }else{
      l0=p1;
      l1=p0;
    }      
  }

  TLorentzVector dilepton=*l0+*l1;
  double l0pp=(l0->E()+l0->Pz())/sqrt(2);
  double l0pm=(l0->E()-l0->Pz())/sqrt(2);
  double l1pp=(l1->E()+l1->Pz())/sqrt(2);
  double l1pm=(l1->E()-l1->Pz())/sqrt(2);
  double dimass=dilepton.M();
  double dipt=dilepton.Pt();
  if(bjet*bjet==0){
    if(direction==0) direction=dilepton.Pz()>0?1:-1;
  }else{
    if(direction==0) direction=(0.75*dilepton.Rapidity()-bjet.Rapidity())>0?1:-1;
    if(bjet_charge>0) direction *= -1.;
  }
  return direction*2*(l0pp*l1pm-l0pm*l1pp)/sqrt(dimass*dimass*(dimass*dimass+dipt*dipt));
} 
double AFBAnalyzer::GetCosThetaRecoil(const Particle *p0,const Particle *p1,int direction){
  const TLorentzVector *lm,*lp;
  if(p0->Charge()<0&&p1->Charge()>0){
    lm=p0;
    lp=p1;
  }else if(p0->Charge()>0&&p1->Charge()<0){
    lm=p1;
    lp=p0;
  }else if(strcmp(p0->ClassName(),"LHE")==0){ 
    if(((LHE*)p0)->ID()>0&&((LHE*)p1)->ID()<0){
      lm=p0;
      lp=p1;
    }else if(((LHE*)p0)->ID()<0&&((LHE*)p1)->ID()>0){
      lm=p1;
      lp=p0;
    }else{
      if(gRandom->Rndm()<0.5){
	lm=p0;
	lp=p1;
      }else{
	lm=p1;
	lp=p0;
      }      
    } 
  }else{
    if(gRandom->Rndm()<0.5){
      lm=p0;
      lp=p1;
    }else{
      lm=p1;
      lp=p0;
    }      
  }

  if(bjet_charge>0) direction *= -1.;
  return direction*((*lm-*lp)*bjet)/((*lm+*lp)*bjet);
}
double AFBAnalyzer::GetCosTheta(const vector<Lepton*>& leps,const vector<Jet>& jets,TString option,double fcut){
  TLorentzVector *l0,*l1;
  if(leps.at(0)->Charge()<0&&leps.at(1)->Charge()>0){
    l0=leps.at(0);
    l1=leps.at(1);
  }else if(leps.at(0)->Charge()>0&&leps.at(1)->Charge()<0){
    l0=leps.at(1);
    l1=leps.at(0);
  }else{
    if(gRandom->Rndm()<0.5){
      l0=leps.at(0);
      l1=leps.at(1);
    }else{
      l0=leps.at(1);
      l1=leps.at(0);
    }      
  }
  TLorentzVector dilepton=*l0+*l1;
  TLorentzVector j0;
  if(jets.size()>0){
    if((jets[0]+dilepton).Pt()/dilepton.Pt()<fcut)
      j0=jets[0];
  }
  TLorentzVector system=dilepton+j0;
  double direction=dilepton.Pz()/fabs(dilepton.Pz());
  double systemmass=system.M();
  //double systempt=system.Pt();
  double systempz=system.Pz();
  TLorentzVector p0(0,0,(sqrt(systemmass*systemmass+systempz*systempz)+systempz)/2.,(sqrt(systemmass*systemmass+systempz*systempz)+systempz)/2.);
  TLorentzVector p1(0,0,-(sqrt(systemmass*systemmass+systempz*systempz)-systempz)/2.,(sqrt(systemmass*systemmass+systempz*systempz)-systempz)/2.);
  TLorentzVector particle=*l0;
  TVector3 b1=system.BoostVector();
  dilepton.Boost(-b1);particle.Boost(-b1);
  p0.Boost(-b1);p1.Boost(-b1);j0.Boost(-b1);
  if(j0.E()>0){
    if(p0.Angle(j0.Vect())<p1.Angle(j0.Vect())){
      if(option.Contains("C")) p0-=j0;
      if(option.Contains("B")) direction=-1;
    }else{
      if(option.Contains("C")) p1-=j0;
      if(option.Contains("B")) direction=1;
    }
  }
  TVector3 b2=dilepton.BoostVector();
  p0.Boost(-b2);p1.Boost(-b2);particle.Boost(-b2);
  return direction*cos(particle.Angle(p0.Vect().Unit()-p1.Vect().Unit()));
}
void AFBAnalyzer::FillHistsToy(TString pre,TString hpre,TString suf,Particle* l0,Particle* l1,map<TString,double> map_weight){
  int n_toy=toy_random.size();
  for(int i=0;i<n_toy;i++) FillHistsAFB(pre,hpre,suf+Form("_toy%d",i),l0,l1,Multiply(map_weight,toy_weight[i]));
}
void AFBAnalyzer::FillHistsAFB(TString pre,TString hpre,TString suf,Particle* l0,Particle* l1,map<TString,double> map_weight){
  TLorentzVector dilepton=(*l0)+(*l1);
  double dimass=dilepton.M();
  double dirap=dilepton.Rapidity();
  double dipt=dilepton.Pt();

  double cost=GetCosThetaCS(l0,l1);
  double h=0.5*pow(dipt/dimass,2)/(1+pow(dipt/dimass,2))*(1-3*cost*cost);
  double den_weight=0.5*fabs(cost)/pow(1+cost*cost+h,2);
  double num_weight=0.5*cost*cost/pow(1+cost*cost+h,3);

  FillHist(pre+hpre+"costhetaCS"+suf,dimass,dirap,dipt,cost,map_weight,afb_mbinnum,(double*)afb_mbin,afb_ybinnum,(double*)afb_ybin,afb_ptbinnum,(double*)afb_ptbin,20,-1,1);
  //FillHist(pre+hpre+"costhetaCS_den"+suf,dimass,dirap,dipt,cost,Multiply(map_weight,den_weight),afb_mbinnum,(double*)afb_mbin,afb_ybinnum,(double*)afb_ybin,afb_ptbinnum,(double*)afb_ptbin,20,-1,1);
  //FillHist(pre+hpre+"costhetaCS_num"+suf,dimass,dirap,dipt,cost,Multiply(map_weight,num_weight),afb_mbinnum,(double*)afb_mbin,afb_ybinnum,(double*)afb_ybin,afb_ptbinnum,(double*)afb_ptbin,20,-1,1);
  if(bjet*bjet!=0){
    double cosr=GetCosThetaRecoil(l0,l1);
    FillHist(pre+hpre+"costhetaRecoil"+suf,dimass,dirap,dipt,cosr,map_weight,afb_mbinnum,(double*)afb_mbin,afb_ybinnum,(double*)afb_ybin,afb_ptbinnum,(double*)afb_ptbin,20,-1,1);
  }

  if(!HasFlag("PDFSYS")){
    FillHist(pre+hpre+"l0pt"+suf,dimass,dirap,dipt,l0->Pt(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,lptbinnum,(double*)lptbin);
    FillHist(pre+hpre+"l1pt"+suf,dimass,dirap,dipt,l1->Pt(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,lptbinnum,(double*)lptbin);
    FillHist(pre+hpre+"lpt"+suf,dimass,dirap,dipt,l0->Pt(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,lptbinnum,(double*)lptbin);
    FillHist(pre+hpre+"lpt"+suf,dimass,dirap,dipt,l1->Pt(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,lptbinnum,(double*)lptbin);
    
    FillHist(pre+hpre+"l0eta"+suf,dimass,dirap,dipt,l0->Eta(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,60,-3,3);
    FillHist(pre+hpre+"l1eta"+suf,dimass,dirap,dipt,l1->Eta(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,60,-3,3);
    FillHist(pre+hpre+"leta"+suf,dimass,dirap,dipt,l0->Eta(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,60,-3,3);
    FillHist(pre+hpre+"leta"+suf,dimass,dirap,dipt,l1->Eta(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,60,-3,3);

    if(!hpre.Contains("gen")&&!hpre.Contains("lhe")&&!hpre.Contains("truth")){
      FillHist(pre+hpre+"z0"+suf,dimass,dirap,dipt,vertex_Z,map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,120,-15,15);
      FillHist(pre+hpre+"met"+suf,dimass,dirap,dipt,pfMET_Type1_pt,map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,100,0,200);  
      FillHist(pre+hpre+"nPV"+suf,dimass,dirap,dipt,nPV,map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,60,0,60);  
      FillHist(pre+hpre+"nPileUp"+suf,dimass,dirap,dipt,nPileUp,map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,60,0,60);  
      FillHist(pre+hpre+"bjetCh"+suf,dimass,dirap,dipt,bjet_charge,map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,20,-5,5);
    }
  }
}
void AFBAnalyzer::FillHardHists(TString pre,TString suf,const Gen& genparton0,const Gen& genparton1,const Gen& genhardl0,const Gen& genhardl1,const Gen& genhardj0,double w){
  Gen genhardl=genhardl0.PID()>0?genhardl0:genhardl1;
  TLorentzVector genZ=genhardl0+genhardl1;
  TLorentzVector genpp=genparton0+genparton1;

  FillHist(pre+"cos_l_p0"+suf,cos(genhardl.Angle(genparton0.Vect())),w,100,-1,1);
  FillHist(pre+"cos_l_p0_asym"+suf,cos(genhardl.Angle(genparton0.Vect())),w/2,100,-1,1);
  FillHist(pre+"cos_l_p0_asym"+suf,-1.*cos(genhardl.Angle(genparton0.Vect())),-w/2,100,-1,1);
  if(cos(genhardl.Angle(genparton0.Vect()))>0) FillHist(pre+"cos_l_p0_forward"+suf,genZ.M(),w,fine_mbinnum,(double*)fine_mbin);
  else FillHist(pre+"cos_l_p0_backward"+suf,genZ.M(),w,fine_mbinnum,(double*)fine_mbin);

  FillHist(pre+"cos_l_p1"+suf,cos(genhardl.Angle(genparton1.Vect())),w,100,-1,1);
  FillHist(pre+"cos_l_p1_asym"+suf,cos(genhardl.Angle(genparton1.Vect())),w/2,100,-1,1);
  FillHist(pre+"cos_l_p1_asym"+suf,-1.*cos(genhardl.Angle(genparton1.Vect())),-w/2,100,-1,1);
  if(cos(genhardl.Angle(genparton1.Vect()))>0) FillHist(pre+"cos_l_p1_forward"+suf,genZ.M(),w,fine_mbinnum,(double*)fine_mbin);
  else FillHist(pre+"cos_l_p1_backward"+suf,genZ.M(),w,fine_mbinnum,(double*)fine_mbin);
  
  FillHist(pre+"cos_Z_p0"+suf,cos(genZ.Angle(genparton0.Vect())),w,100,-1,1);
  FillHist(pre+"cos_Z_p0_asym"+suf,cos(genZ.Angle(genparton0.Vect())),w/2,100,-1,1);
  FillHist(pre+"cos_Z_p0_asym"+suf,-1.*cos(genZ.Angle(genparton0.Vect())),-w/2,100,-1,1);
  
  FillHist(pre+"cos_Z_p1"+suf,cos(genZ.Angle(genparton1.Vect())),w,100,-1,1);
  FillHist(pre+"cos_Z_p1_asym"+suf,cos(genZ.Angle(genparton1.Vect())),w/2,100,-1,1);
  FillHist(pre+"cos_Z_p1_asym"+suf,-1.*cos(genZ.Angle(genparton1.Vect())),-w/2,100,-1,1);
  
  
  FillHist(pre+"Zrap"+suf,genZ.Rapidity(),w,100,-5,5);
  FillHist(pre+"Zrap_asym"+suf,genZ.Rapidity(),w/2,100,-5,5);
  FillHist(pre+"Zrap_asym"+suf,-1.*genZ.Rapidity(),-w/2,100,-5,5);
  
  FillHist(pre+"pprap"+suf,genpp.Rapidity(),w,100,-5,5);
  FillHist(pre+"pprap_asym"+suf,genpp.Rapidity(),w/2,100,-5,5);
  FillHist(pre+"pprap_asym"+suf,-1.*genpp.Rapidity(),-w/2,100,-5,5);
  
  if(!genhardj0.IsEmpty()){
    FillHist(pre+"cos_l_j0"+suf,cos(genhardl.Angle(genhardj0.Vect())),w,100,-1,1);
    FillHist(pre+"cos_l_j0_asym"+suf,cos(genhardl.Angle(genhardj0.Vect())),w/2,100,-1,1);
    FillHist(pre+"cos_l_j0_asym"+suf,-1.*cos(genhardl.Angle(genhardj0.Vect())),-w/2,100,-1,1);
    if(cos(genhardl.Angle(genhardj0.Vect()))>0) FillHist(pre+"cos_l_j0_forward"+suf,genZ.M(),w,fine_mbinnum,(double*)fine_mbin);
    else FillHist(pre+"cos_l_j0_backward"+suf,genZ.M(),w,fine_mbinnum,(double*)fine_mbin);

    FillHist(pre+"cos_Z_j0"+suf,cos(genZ.Angle(genhardj0.Vect())),w,100,-1,1);
    FillHist(pre+"cos_Z_j0_asym"+suf,cos(genZ.Angle(genhardj0.Vect())),w/2,100,-1,1);
    FillHist(pre+"cos_Z_j0_asym"+suf,-1.*cos(genZ.Angle(genhardj0.Vect())),-w/2,100,-1,1);
    
    FillHist(pre+"cos_j0_p0"+suf,cos(genhardj0.Angle(genparton0.Vect())),w,100,-1,1);
    FillHist(pre+"cos_j0_p0_asym"+suf,cos(genhardj0.Angle(genparton0.Vect())),w/2,100,-1,1);
    FillHist(pre+"cos_j0_p0_asym"+suf,-1.*cos(genhardj0.Angle(genparton0.Vect())),-w/2,100,-1,1);
    
    FillHist(pre+"cos_j0_p1"+suf,cos(genhardj0.Angle(genparton1.Vect())),w,100,-1,1);
    FillHist(pre+"cos_j0_p1_asym"+suf,cos(genhardj0.Angle(genparton1.Vect())),w/2,100,-1,1);
    FillHist(pre+"cos_j0_p1_asym"+suf,-1.*cos(genhardj0.Angle(genparton1.Vect())),-w/2,100,-1,1);
    
    FillHist(pre+"jeta"+suf,genhardj0.Eta(),w,100,-5,5);
    FillHist(pre+"jeta_asym"+suf,genhardj0.Eta(),w/2,100,-5,5);
    FillHist(pre+"jeta_asym"+suf,-1.*genhardj0.Eta(),-w/2,100,-5,5);
  }
}
void AFBAnalyzer::SetupToy(int n_toy){
  DeleteToy();
  for(int i=0;i<n_toy;i++){
    toy_random.push_back(new TRandom3);
    toy_weight.push_back(-9999.);
  } 
}
void AFBAnalyzer::DeleteToy(){
  int n_toy=toy_random.size();
  for(int i=0;i<n_toy;i++) delete toy_random.at(i);
  toy_random.clear();
  toy_weight.clear();
}
void AFBAnalyzer::GetToyWeight(){
  int n_toy=toy_random.size();
  if(fChain->GetTree()->GetReadEntry()==0){
    TString filename=fChain->GetFile()->GetName();
    for(int i=0;i<n_toy;i++) toy_random[i]->SetSeed((filename+Form("%d",i)).MD5().Hash());
  }
  for(int i=0;i<n_toy;i++)
    toy_weight[i]=toy_random[i]->PoissonD(1.);
}

void AFBAnalyzer::FillHistToy(TString histname, double value, double weight, int n_bin, double x_min, double x_max){
  int n_toy=toy_random.size();
  for(int i=0;i<n_toy;i++) FillHist(histname+Form("_toy%d",i),value,weight*toy_weight[i],n_bin,x_min,x_max);
}
void AFBAnalyzer::FillHistToy(TString histname, double value, map<TString,double> weights, int n_bin, double x_min, double x_max){
  for(const auto& [suffix,weight]:weights) FillHistToy(histname+suffix,value,weight,n_bin,x_min,x_max);
}
void AFBAnalyzer::FillHistToy(TString histname, double value, double weight, int n_bin, double *xbins){ 
  int n_toy=toy_random.size();
  for(int i=0;i<n_toy;i++) FillHist(histname+Form("_toy%d",i),value,weight*toy_weight[i],n_bin,xbins);
}
void AFBAnalyzer::FillHistToy(TString histname, double value, map<TString,double> weights, int n_bin, double *xbins){ 
  for(const auto& [suffix,weight]:weights) FillHistToy(histname+suffix,value,weight,n_bin,xbins);
}
void AFBAnalyzer::FillHistToy(TString histname, double value_x, double value_y, double value_z, double weight, int n_binx, double *xbins, int n_biny, double *ybins, int n_binz, double *zbins){
  int n_toy=toy_random.size();
  for(int i=0;i<n_toy;i++) FillHist(histname+Form("_toy%d",i),value_x,value_y,value_z,weight,n_binx,xbins,n_biny,ybins,n_binz,zbins);
}
void AFBAnalyzer::FillHistToy(TString histname, double value_x, double value_y, double value_z, map<TString,double> weights, int n_binx, double *xbins, int n_biny, double *ybins, int n_binz, double *zbins){
  for(const auto& [suffix,weight]:weights) FillHistToy(histname+suffix,value_x,value_y,value_z,weight,n_binx,xbins,n_biny,ybins,n_binz,zbins);
}
void AFBAnalyzer::SetupCosThetaWeight(){
  cout<<"[AFBAnalyzer::SetupCosThetaWeight] Setup"<<endl;
  TString datapath=getenv("DATA_DIR");
  ifstream file_check(datapath+"/"+GetEra()+"/SMP/CosThetaWeight.root");
  bool isexist=file_check.is_open();
  file_check.close();
  if(!isexist){
    cout<<"[AFBAnalyzer::SetupCosThetaWeight] no CosThetaWeight.root"<<endl;
    return;
  }
  TFile fcost(datapath+"/"+GetEra()+"/SMP/CosThetaWeight.root");
  for(const auto&& key:*(fcost.GetListOfKeys())){
    TObject* obj=((TKey*)key)->ReadObj();
    if(!obj->InheritsFrom("TH3D")) continue;
    TH3D* hist=(TH3D*)obj;
    cout<<"[AFBAnalyzer::SetupCosThetaWeight] get "<<hist->GetName()<<endl;
    map_hist_cost[hist->GetName()]=hist;
    hist->SetDirectory(0);
  }
}
void AFBAnalyzer::DeleteCosThetaWeight(){
  for(auto& iter:map_hist_cost)
    if(iter.second) delete iter.second;
}
double AFBAnalyzer::GetCosThetaWeight(double mass,double pt,double cost,TString suffix){
  double val=1.;
  if(!IsDYSample) return val;
  TString MCName=MCSample;
  if(MCName.Contains(TRegexp("^DY[0-9]Jets$"))) MCName="DYJets";
  if(MCName.Contains(TRegexp("^DYJets_Pt-[0-9]*To[0-9Inf]*$"))) MCName="DYJets";
  if(MCName.Contains(TRegexp("^DYJets_M-[0-9]*to[0-9Inf]*$"))) MCName="DYJets";
  TString hname=MCName+suffix;
  auto it=map_hist_cost.find(hname);
  if(it!=map_hist_cost.end())
    val*=GetBinContentUser(it->second,mass,pt,cost,0);
  if(val==0) val=1.;
  return val;
}
