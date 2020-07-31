#include "AFBAnalyzer.h"

void AFBAnalyzer::initializeAnalyzer(){
  SMPAnalyzerCore::initializeAnalyzer(); //setup zpt roc z0 
  SetupCosThetaWeight();
  
  vector<JetTagging::Parameters> jtps={JetTagging::Parameters(JetTagging::DeepCSV,JetTagging::Tight,JetTagging::mujets,JetTagging::mujets)};
  mcCorr->SetJetTaggingParameters(jtps);

  SetupToy(100);
  IsNominalRun=!HasFlag("SYS")&&!HasFlag("PDFSYS");
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
  GetToyWeight();

  if(!HLT_TriggerName) HLT_TriggerName=new vector<string>;
  event=GetEvent();
  GetEventWeights();

  hardprefix="";
  costhetaweight=1.;
  costhetaweight_up=1.;
  costhetaweight_down=1.;

  muons.clear();  electrons.clear();  jets.clear();
  muons=SMPGetMuons("POGMediumWithLooseTrkIso",0.0,2.4);
  electrons=SMPGetElectrons("passMediumID",0.0,2.4);
  jets=GetJets("tightLepVeto",30,2.4);
  std::sort(jets.begin(),jets.end(),PtComparing);
  realjets.clear();
  for(unsigned int l=0; l<jets.size();l++){
    //jets passing PileUp ID
    if(jets.at(l).Pt()<30){
      if(jets.at(l).GetPileupJetId() <0.18) continue;
    }else if(jets.at(l).Pt()<50){
      if(jets.at(l).GetPileupJetId() <0.61) continue;
    }
    //jets spliting from (sub)leading electrons, and muons - jet cleaning
    if(muons.size()>0 && muons.at(0).Pt()>20 && jets.at(l).DeltaR(muons.at(0)) <0.4) continue;
    if(muons.size()>1 && muons.at(1).Pt()>10 && jets.at(l).DeltaR(muons.at(1)) <0.4) continue;
    if(electrons.size()>0 && electrons.at(0).Pt()>25 && jets.at(l).DeltaR(electrons.at(0)) <0.4) continue;
    if(electrons.size()>1 && electrons.at(1).Pt()>15 && jets.at(l).DeltaR(electrons.at(1)) <0.4) continue;

    realjets.push_back(jets.at(l));
  }
  //bjet setting
  n_bjet=0;
  bjet_charge=-1.5;
  bjet_rap=0;
  bjets.clear();
  JetTagging::Parameters jtp = JetTagging::Parameters(JetTagging::DeepCSV,JetTagging::Tight,JetTagging::mujets,JetTagging::mujets);
  for(const auto& jet:realjets){
    if(mcCorr->IsBTagged_2a(jtp,jet)){
      n_bjet++;
      bjets.push_back(jet);
    }
  }

  softmus.clear();
  softmus=GetMuons("POGTight",0,2.4);
  std::sort(softmus.begin(),softmus.end(),PtComparing);

  if(IsDYSample){
    //////////////////////// Check LHE /////////////////////////
    vector<LHE> lhes=GetLHEs();
    LHE lhe_l0,lhe_l1,lhe_j0;
    GetDYLHEParticles(lhes,lhe_l0,lhe_l1,lhe_j0);
    TString channelname="";
    if(abs(lhe_l0.ID())!=15){
      double l0ptcut,l1ptcut,letacut;
      if(abs(lhe_l0.ID())==11){
	channelname=Form("ee%d",DataYear);
	l0ptcut=25.;
	l1ptcut=15.;
	letacut=2.4;
      }else if(abs(lhe_l0.ID())==13){
	channelname=Form("mm%d",DataYear);
	l0ptcut=20.;
	l1ptcut=10.;
	letacut=2.4;
      }else{
	cout<<"[AFBAnalyzer::executeEvent()] something is wrong l0.ID="<<abs(lhe_l0.ID())<<endl;
	for(auto& lhe:lhes) lhe.Print();
	exit(EXIT_FAILURE);
      }

      /*
      vector<LHE> lhes_final;
      lhes_final.clear();
      for(unsigned int i=0;i<lhes.size();i++){
	if(lhes[i].Status()==1) lhes_final.push_back(lhes[i]);
      }
      if(lhes_final.size() != 3) return;
      */

      //////////////////////// GEN /////////////////////////
      vector<Gen> gens=GetGens();
      Gen gen_parton0,gen_parton1,gen_l0,gen_l1,gen_l0_dressed,gen_l1_dressed,gen_l0_bare,gen_l1_bare,gen_j0;
      GetDYGenParticles(gens,gen_parton0,gen_parton1,gen_l0,gen_l1,gen_j0,3);
      //GetDYGenParticles(gens,gen_parton0,gen_parton1,gen_l0_dressed,gen_l1_dressed,gen_j0,1);
      //GetDYGenParticles(gens,gen_parton0,gen_parton1,gen_l0_bare,gen_l1_bare,gen_j0,0);
      /*
    int hardj0=0,hardj1=0,hardj2=0;
    Gen genhardj0,genhardj1,genhardj2;
    for(unsigned int i=parton1+1;i<gens.size();i++){
      if(gens[i].isHardProcess()){
	if(gens[i].PID()==21||abs(gens[i].PID())<7){
	  if(!hardj0){
	    hardj0=i;
	    genhardj0=gens[i];
	  }else if(!hardj1){
	    hardj1=i;
	    genhardj1=gens[i];
	  }else if(!hardj2){
	    hardj2=i;
	    genhardj2=gens[i];
	    break;
	  }
	}
      }
    }
    */
      //genparton0.SetPxPyPzE(0,0,genWeight_X1*6500.,genWeight_X1*6500.);
      //genparton1.SetPxPyPzE(0,0,-genWeight_X2*6500.,genWeight_X2*6500.);
      //for(unsigned int i=0;i<photons.size();i++) genphotons+=gens[photons[i]];
      /*
      if(genparton0.PID()==21){
	if(genparton1.PID()==21) hardprefix="gg_";
	else if(genparton1.PID()>0) hardprefix="gq_";
	else if(genparton1.PID()<0) hardprefix="gqbar_";
      }else if(genparton0.PID()>0){
	if(genparton1.PID()==21) hardprefix="qg_";
	else if(genparton1.PID()>0) hardprefix="qq_";
	else if(genparton1.PID()<0) hardprefix="qqbar_";
      }else if(genparton0.PID()<0){
	if(genparton1.PID()==21) hardprefix="qbarg_";
	else if(genparton1.PID()>0) hardprefix="qbarq_";
	else if(genparton1.PID()<0) hardprefix="qbarqbar_";
      }
      int nhardjet=(hardj0?1:0)+(hardj1?1:0)+(hardj2?1:0);
      */
     
      TLorentzVector gen_dilepton=gen_l0+gen_l1;
      double gen_dimass=gen_dilepton.M();
      double gen_dirap=gen_dilepton.Rapidity();
      double gen_dipt=gen_dilepton.Pt();
      double gen_cost_correct=-999;
      if(gen_parton0.PID()==21){
	if(gen_parton1.PID()==21) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,0);
	else if(gen_parton1.PID()>0) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,-1);
	else if(gen_parton1.PID()<0) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,1);
      }else if(gen_parton0.PID()>0){
	if(gen_parton1.PID()==21) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,1);
	else if(gen_parton1.PID()>0) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,0);
	else if(gen_parton1.PID()<0) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,1);
      }else if(gen_parton0.PID()<0){
	if(gen_parton1.PID()==21) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,-1);
	else if(gen_parton1.PID()>0) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,-1);
	else if(gen_parton1.PID()<0) gen_cost_correct=GetCosThetaCS(&gen_l0,&gen_l1,0);
      }
      if(gen_cost_correct==-999){
	cout<<"wrong pid for parton"<<endl;
	exit(EXIT_FAILURE);
      }
      costhetaweight=GetCosThetaWeight(gen_dimass,gen_dipt,gen_cost_correct,"_pdg");
      costhetaweight_up=GetCosThetaWeight(gen_dimass,gen_dipt,gen_cost_correct,"_up");
      costhetaweight_down=GetCosThetaWeight(gen_dimass,gen_dipt,gen_cost_correct,"_down");
      
      map<TString,double> map_weight;
      //map_weight[""]=lumiweight*zptweight*costhetaweight;
      //map_weight["_noweight"]=lumiweight;
      map_weight[""]=lumiweight; 

      // This test yields nothing as expected!
      //if(lhes[0].ID()!=gen_parton0.PID()||lhes[1].ID()!=gen_parton1.PID()) tauprefix+="Weird1_";
      //if(lhes[0].Pz()<0||lhes[1].Pz()>0) tauprefix+="Weird2_";

      // Only qG, Gq collisions
      if((abs(gen_parton0.PID())<=5&&gen_parton1.PID()==21) || (gen_parton0.PID()==21&&abs(gen_parton1.PID())<=5)){
	//if((abs(gen_parton0.PID())<=5&&gen_parton1.PID()==21)) tauprefix+="qG_";
	//else tauprefix+="Gq_";

	if(gen_parton0.PID()==5||gen_parton1.PID()==5)  tauprefix+="Dyb_";
	else if(gen_parton0.PID()==-5||gen_parton1.PID()==-5)  tauprefix+="Dybbar_";
        else if(gen_parton0.PID()==4||gen_parton1.PID()==4)  tauprefix+="Dyc_";
        else if(gen_parton0.PID()==-4||gen_parton1.PID()==-4){
	  tauprefix+="Dycbar_";
	  /*
	  if(gen_j0.PID()!=-4){
	    cout<<"Dyc but gen_j0 isn't cbar, ID is "<<gen_j0.PID()<<endl;
	    PrintGens(gens);
	  }
	  if(lhe_j0.ID()!=-4){
	    cout<<"Dyc but lhe_j0 isn't cbar, ID is "<<lhe_j0.ID()<<endl;
	    for(int i=0;i<(int)lhes.size();i++){
	      cout<<lhes[i].Index()<<"\t"<<lhes[i].ID()<<"\t"<<lhes[i].Status()<<"\t"<<lhes[i].E()<<"\t"<<lhes[i].Px()<<"\t"<<lhes[i].Py()<<"\t"<<lhes[i].Pz()<<"\t"<<lhes[i].Eta()<<"\t"<<lhes[i].M()<<"\t"<<endl;
	    }
	  }
	  */
	}
	else if(abs(gen_parton0.PID())<4||abs(gen_parton1.PID())<4) tauprefix+="Dyudsg_";
	else tauprefix+="Weird_"; //this should be empty
      }
      //qq qqbar GG collsions
      else{
	/*
	if(gen_j0.PID()==5 && lhe_j0.ID()==5) tauprefix+="bBkg2_";
        else if(gen_j0.PID()==-5 && lhe_j0.ID()==-5) tauprefix+="bBkg3_";
	else if(gen_j0.PID()==4 && lhe_j0.ID()==4) tauprefix+="cBkg2_";
        else if(gen_j0.PID()==-4 && lhe_j0.ID()==-4) tauprefix+="cBkg3_";
	else if(abs(lhe_j0.ID())==5) tauprefix+="bBkg4_";
        else if(abs(gen_j0.PID())==5) tauprefix+="bBkg5_";
        else if(abs(lhe_j0.ID())==4) tauprefix+="cBkg4_";
        else if(abs(gen_j0.PID())==4) tauprefix+="cBkg5_";
	*/
      }
    
      //////////////// Fill LHE,Gen hists //////////////////////
      if(IsNominalRun){
	vector<TString> accp = {""};//, "accep_"};
	for(unsigned int k=0;k<accp.size();k++){
	  if(accp.at(k)=="accep_"){
	    if(abs(lhe_l0.Rapidity())>2.4 || abs(lhe_l1.Rapidity())>2.4 || abs(lhe_j0.Rapidity())>2.4) continue;
	  }
	  FillHist(channelname+"/"+tauprefix+accp.at(k)+"lhe_jetID",lhe_j0.ID(),map_weight,50,-25,25);
	  FillHist(channelname+"/"+tauprefix+accp.at(k)+"lhe_jetpT",lhe_j0.Pt(),map_weight,200,0,1000);
	  FillHist(channelname+"/"+tauprefix+accp.at(k)+"lhe_jeteta",lhe_j0.Eta(),map_weight,200,-10,10);
	  FillHist(channelname+"/"+tauprefix+accp.at(k)+"lhe_jet_Zdy",(lhe_l0+lhe_l1).Rapidity()-lhe_j0.Rapidity(),map_weight,200,-10,10);
	  FillHist(channelname+"/"+tauprefix+accp.at(k)+"lhe_jet_Zdphi",(lhe_l0+lhe_l1).DeltaPhi(lhe_j0),map_weight,100,-5,5);
	  FillHist(channelname+"/"+tauprefix+accp.at(k)+"lhe_jet_Zy",(lhe_l0+lhe_l1).Rapidity(),map_weight,200,-10,10);
	  FillHist(channelname+"/"+tauprefix+accp.at(k)+"lhe_jet_Zy_2D",(lhe_l0+lhe_l1).Rapidity(),lhe_j0.Rapidity(),map_weight,200,-10,10,200,-10,10);
	  /*FillHist(channelname+"/"+tauprefix+accp.at(k)+"lhe_jet_Zdytest0",((lhe_l0+lhe_l1).Rapidity()-lhe_j0.Rapidity())>0?1:0,map_weight,10,-5,5);
	  FillHist(channelname+"/"+tauprefix+accp.at(k)+"lhe_jet_Zdytest1",(((lhe_l0+lhe_l1).Rapidity()-lhe_j0.Rapidity()) * ((lhe_l0+lhe_l1).Rapidity()+lhe_j0.Rapidity()))>0?1:0,map_weight,10,-5,5);
	  FillHist(channelname+"/"+tauprefix+accp.at(k)+"lhe_jet_Zdytest2",((lhe_l0+lhe_l1).Rapidity()*lhe_j0.Rapidity())>0?1:0,map_weight,10,-5,5);
	  FillHist(channelname+"/"+tauprefix+accp.at(k)+"lhe_jet_Zdytest3",(lhe_j0.Rapidity())>0?1:0,map_weight,10,-5,5);
	  FillHist(channelname+"/"+tauprefix+accp.at(k)+"lhe_jet_Zdytest4",(((lhe_l0+lhe_l1).Rapidity()-lhe_j0.Rapidity()) * (lhe_j0.Rapidity()))>0?1:0,map_weight,10,-5,5);
	  FillHist(channelname+"/"+tauprefix+accp.at(k)+"lhe_jet_Zdytest5",(((lhe_l0+lhe_l1).Rapidity()-lhe_j0.Rapidity()) * ((lhe_l0+lhe_l1).Rapidity()+lhe_j0.Rapidity()) * (lhe_l0+lhe_l1).Rapidity()*lhe_j0.Rapidity())>0?1:0,map_weight,10,-5,5);
	  FillHist(channelname+"/"+tauprefix+accp.at(k)+"lhe_jet_Zdytest6",(0.9*(lhe_l0+lhe_l1).Rapidity()-lhe_j0.Rapidity())>0?1:0,map_weight,10,-5,5);
	  FillHist(channelname+"/"+tauprefix+accp.at(k)+"lhe_jet_Zdytest7",(0.8*(lhe_l0+lhe_l1).Rapidity()-lhe_j0.Rapidity())>0?1:0,map_weight,10,-5,5);
	  FillHist(channelname+"/"+tauprefix+accp.at(k)+"lhe_jet_Zdytest8",(0.7*(lhe_l0+lhe_l1).Rapidity()-lhe_j0.Rapidity())>0?1:0,map_weight,10,-5,5);
	  FillHist(channelname+"/"+tauprefix+accp.at(k)+"lhe_jet_Zdytest9",(0.6*(lhe_l0+lhe_l1).Rapidity()-lhe_j0.Rapidity())>0?1:0,map_weight,10,-5,5);
	  FillHist(channelname+"/"+tauprefix+accp.at(k)+"lhe_jet_Zdytest10",(0.5*(lhe_l0+lhe_l1).Rapidity()-lhe_j0.Rapidity())>0?1:0,map_weight,10,-5,5);
	  FillHist(channelname+"/"+tauprefix+accp.at(k)+"lhe_jet_Zdytest11",(0.4*(lhe_l0+lhe_l1).Rapidity()-lhe_j0.Rapidity())>0?1:0,map_weight,10,-5,5);
	  FillHist(channelname+"/"+tauprefix+accp.at(k)+"gen_jet_Zdytest0",((gen_l0+gen_l1).Rapidity()-gen_j0.Rapidity())>0?1:0,map_weight,10,-5,5);
	  FillHist(channelname+"/"+tauprefix+accp.at(k)+"gen_jet_Zdytest1",(((gen_l0+gen_l1).Rapidity()-gen_j0.Rapidity()) * ((gen_l0+gen_l1).Rapidity()+gen_j0.Rapidity()))>0?1:0,map_weight,10,-5,5);
	  FillHist(channelname+"/"+tauprefix+accp.at(k)+"gen_jet_Zdytest2",((gen_l0+gen_l1).Rapidity()*gen_j0.Rapidity())>0?1:0,map_weight,10,-5,5);
	  FillHist(channelname+"/"+tauprefix+accp.at(k)+"gen_jet_Zdytest3",(gen_j0.Rapidity())>0?1:0,map_weight,10,-5,5);
	  FillHist(channelname+"/"+tauprefix+accp.at(k)+"gen_jet_Zdytest4",(((gen_l0+gen_l1).Rapidity()-gen_j0.Rapidity()) * (gen_j0.Rapidity()))>0?1:0,map_weight,10,-5,5);
	  FillHist(channelname+"/"+tauprefix+accp.at(k)+"gen_jet_Zdytest5",(((gen_l0+gen_l1).Rapidity()-gen_j0.Rapidity()) * ((gen_l0+gen_l1).Rapidity()+gen_j0.Rapidity()) * (gen_l0+gen_l1).Rapidity()*gen_j0.Rapidity())>0?1:0,map_weight,10,-5,5);
	  FillHist(channelname+"/"+tauprefix+accp.at(k)+"gen_jet_Zdytest6",(0.9*(gen_l0+gen_l1).Rapidity()-gen_j0.Rapidity())>0?1:0,map_weight,10,-5,5);
	  FillHist(channelname+"/"+tauprefix+accp.at(k)+"gen_jet_Zdytest7",(0.8*(gen_l0+gen_l1).Rapidity()-gen_j0.Rapidity())>0?1:0,map_weight,10,-5,5);
	  FillHist(channelname+"/"+tauprefix+accp.at(k)+"gen_jet_Zdytest8",(0.7*(gen_l0+gen_l1).Rapidity()-gen_j0.Rapidity())>0?1:0,map_weight,10,-5,5);
	  FillHist(channelname+"/"+tauprefix+accp.at(k)+"gen_jet_Zdytest9",(0.6*(gen_l0+gen_l1).Rapidity()-gen_j0.Rapidity())>0?1:0,map_weight,10,-5,5);
	  FillHist(channelname+"/"+tauprefix+accp.at(k)+"gen_jet_Zdytest10",(0.5*(gen_l0+gen_l1).Rapidity()-gen_j0.Rapidity())>0?1:0,map_weight,10,-5,5);
	  FillHist(channelname+"/"+tauprefix+accp.at(k)+"gen_jet_Zdytest11",(0.4*(gen_l0+gen_l1).Rapidity()-gen_j0.Rapidity())>0?1:0,map_weight,10,-5,5);
	  */
	  FillHist(channelname+"/"+tauprefix+accp.at(k)+"gen_jetID",gen_j0.PID(),map_weight,50,-25,25);
	  FillHist(channelname+"/"+tauprefix+accp.at(k)+"gen_jetpT",gen_j0.Pt(),map_weight,200,0,1000);
	  FillHist(channelname+"/"+tauprefix+accp.at(k)+"gen_jeteta",gen_j0.Eta(),map_weight,200,-10,10);
	  FillHist(channelname+"/"+tauprefix+accp.at(k)+"gen_jet_Zdy",(gen_l0+gen_l1).Rapidity()-gen_j0.Rapidity(),map_weight,200,-10,10);
	  FillHist(channelname+"/"+tauprefix+accp.at(k)+"gen_jet_Zdphi",(gen_l0+gen_l1).DeltaPhi(gen_j0),map_weight,100,-5,5);
	  FillHist(channelname+"/"+tauprefix+accp.at(k)+"gen_jet_Zy",(gen_l0+gen_l1).Rapidity(),map_weight,200,-10,10);
	  //FillHist(channelname+"/"+tauprefix+"gen_jet_Zdy",((gen_l0+gen_l1).Rapidity()-gen_j0.Rapidity())>0?1:0,map_weight,10,-5,5);
	  FillHist(channelname+"/"+tauprefix+accp.at(k)+"gen_jet_Zy_2D",(gen_l0+gen_l1).Rapidity(),gen_j0.Rapidity(),map_weight,200,-10,10,200,-10,10);

	  FillHist(channelname+"/"+tauprefix+accp.at(k)+"lhe_gen_dpT",lhe_j0.Pt()-gen_j0.Pt(),map_weight,200,-50,50);
	  FillHist(channelname+"/"+tauprefix+accp.at(k)+"lhe_gen_dR",lhe_j0.DeltaR(gen_j0),map_weight,100,0,5);
	  FillHist(channelname+"/"+tauprefix+accp.at(k)+"lhe_gen_dID",lhe_j0.ID()-gen_j0.PID(),map_weight,100,-50,50);

	  if(bjets.size()>0){
	    FillHist(channelname+"/"+tauprefix+accp.at(k)+"lhe_bjetdpT",lhe_j0.Pt()-bjets.at(0).Pt(),map_weight,200,-100,100);
	    FillHist(channelname+"/"+tauprefix+accp.at(k)+"lhe_bjetdR",lhe_j0.DeltaR(bjets.at(0)),map_weight,100,0,5);
	    FillHist(channelname+"/"+tauprefix+accp.at(k)+"gen_bjetdpT",gen_j0.Pt()-bjets.at(0).Pt(),map_weight,200,-100,100);
	    FillHist(channelname+"/"+tauprefix+accp.at(k)+"gen_bjetdR",gen_j0.DeltaR(bjets.at(0)),map_weight,100,0,5);
	  }
	}
        FillHists(channelname,"lhe_","",(Particle*)&lhe_l0,(Particle*)&lhe_l1,map_weight);
	FillHists(channelname,"gen_","",(Particle*)&gen_l0,(Particle*)&gen_l1,map_weight);
	FillHists(channelname,"gen_","_dressed",(Particle*)&gen_l0_dressed,(Particle*)&gen_l1_dressed,map_weight);
	FillHists(channelname,"gen_","_bare",(Particle*)&gen_l0_bare,(Particle*)&gen_l1_bare,map_weight);
	if(gen_l0.Pt()>l0ptcut&&gen_l1.Pt()>l1ptcut&&fabs(gen_l0.Eta())<letacut&&fabs(gen_l1.Eta())<letacut){
	  FillHists(channelname,"genfid_","",(Particle*)&gen_l0,(Particle*)&gen_l1,map_weight);
	}
	FillHist(channelname+"/"+tauprefix+"gen_costhetaCS_correct",gen_dimass,gen_dirap,gen_dipt,gen_cost_correct,map_weight,afb_mbinnum,(double*)afb_mbin,afb_ybinnum,(double*)afb_ybin,afb_ptbinnum,(double*)afb_ptbin,20,-1,1);
      }
    }
  }
  
  if(!PassMETFilter()) return;
  TString prefix="";
			
  ///////////////// cutflow ///////////////////
  if(IsNominalRun){
    FillCutflow(prefix+tauprefix+"cutflow","lumi",lumiweight);
    FillCutflow(prefix+tauprefix+"cutflow","PU",lumiweight*PUweight);
    FillCutflow(prefix+tauprefix+"cutflow","prefire",lumiweight*PUweight*prefireweight);
    FillCutflow(prefix+tauprefix+"cutflow","zpt",lumiweight*PUweight*prefireweight*zptweight);
    FillCutflow(prefix+tauprefix+"cutflow","z0",lumiweight*PUweight*prefireweight*zptweight*z0weight);
    FillCutflow(prefix+tauprefix+"cutflow","costheta",lumiweight*PUweight*prefireweight*zptweight*z0weight*costhetaweight);
  }

  vector<TString> mmtrigger, eetrigger, emtrigger;
  if(DataYear==2016){
    mmtrigger = {
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v",
      "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
      "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"
    };
    eetrigger = {
      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
    };
    emtrigger = {
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"
    };
  }
  else if(DataYear==2017){
    mmtrigger = {
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v",
    };
    eetrigger = {
      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v",
    };
    emtrigger = {
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"
    };
  }
  else if(DataYear==2018){
    mmtrigger = {
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v",
    };
    eetrigger = {
      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v",
    };
    emtrigger = {
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"
    };
  }

  if(event.PassTrigger(mmtrigger)){
    if(!IsDATA||DataStream.Contains("DoubleMuon")) executeEventWithChannelName(prefix+"mm"+Form("%d",DataYear));
  }
  if(event.PassTrigger(eetrigger)&&!event.PassTrigger(mmtrigger)){
    if(!IsDATA||DataStream.Contains("DoubleEG")||DataStream.Contains("EGamma")) executeEventWithChannelName(prefix+"ee"+Form("%d",DataYear));
  }
  /*
  //if(event.PassTrigger(emtrigger)){
  //  if(!IsDATA||DataStream.Contains("MuonEG")) executeEventWithChannelName(prefix+"em"+Form("%d",DataYear));
  //}

  //executeEventWithChannelName(prefix+"NoTrig"+Form("%d",DataYear)); //Not require firing trigger (For Data, just triggers in Stream)

  if(event.PassTrigger(mmtrigger)){
    if(!IsDATA||DataStream.Contains("DoubleMuon")) executeEventWithChannelName(prefix+"DilepTrig"+Form("%d",DataYear));
  }
  if(event.PassTrigger(eetrigger)&&!event.PassTrigger(mmtrigger)){
    if(!IsDATA||DataStream.Contains("DoubleEG")||DataStream.Contains("EGamma")) executeEventWithChannelName(prefix+"DilepTrig"+Form("%d",DataYear));
  }
  if(event.PassTrigger(emtrigger)&&!(event.PassTrigger(mmtrigger) || event.PassTrigger(eetrigger))){
    if(!IsDATA||DataStream.Contains("MuonEG")) executeEventWithChannelName(prefix+"DilepTrig"+Form("%d",DataYear));
  }
  */
}

void AFBAnalyzer::executeEventWithChannelName(TString channelname){
  map<TString,vector<Muon>> map_muons;
  map<TString,vector<Electron>> map_electrons;
  map<TString,Parameter> map_parameter;
  
  if(channelname.Contains(TRegexp("mm20[0-9][0-9]"))){
    Parameter p("IDISO_SF_MediumID_trkIsoLoose_Q","",{"Mu17Leg1_MediumID_trkIsoLoose_Q","Mu8Leg2_MediumID_trkIsoLoose_Q"},20.,10.);
    
    map_muons[""]=MuonMomentumCorrection(muons,0);
    map_parameter[""]=p.Clone(MakeLeptonPointerVector(map_muons[""]),
			      (IsNominalRun?NominalWeight:0)
			      +(HasFlag("SYS")&&!IsDATA?SystematicWeight:0)
			      +(HasFlag("PDFSYS")&&!IsDATA?PDFWeight:0)
			      );
    
    if(HasFlag("SYS")){
      map_muons["_scale_up"]=MuonMomentumCorrection(map_muons[""],+1);
      map_parameter["_scale_up"]=p.Clone(MakeLeptonPointerVector(map_muons["_scale_up"]));
					 
      map_muons["_scale_down"]=MuonMomentumCorrection(map_muons[""],-1);
      map_parameter["_scale_down"]=p.Clone(MakeLeptonPointerVector(map_muons["_scale_down"]));
      
      map_muons["_noroccor"]=MuonMomentumCorrection(map_muons[""],0,-1);
      map_parameter["_noroccor"]=p.Clone(MakeLeptonPointerVector(map_muons["_noroccor"]));
    }
  }else if(channelname.Contains(TRegexp("ee20[0-9][0-9]"))){
    Parameter p("ID_SF_MediumID_Q",{"Ele23Leg1_MediumID_Q","Ele12Leg2_MediumID_Q"},25.,15.);
    
    map_electrons["_noroccor"]=electrons;
    map_electrons[""]=ElectronEnergyCorrection(map_electrons["_noroccor"],0,0);
    map_parameter[""]=p.Clone(MakeLeptonPointerVector(map_electrons[""]),
			      (IsNominalRun?NominalWeight:0)
			      +(HasFlag("SYS")&&!IsDATA?SystematicWeight:0)
			      +(HasFlag("PDFSYS")&&!IsDATA?PDFWeight:0)
			      );

    if(HasFlag("SYS")){
      map_electrons["_scale_up"]=ScaleElectrons(map_electrons[""],1);
      std::sort(map_electrons["_scale_up"].begin(),map_electrons["_scale_up"].end(),PtComparing);
      map_parameter["_scale_up"]=p.Clone(MakeLeptonPointerVector(map_electrons["_scale_up"]));
      
      map_electrons["_scale_down"]=ScaleElectrons(map_electrons[""],-1);
      std::sort(map_electrons["_scale_down"].begin(),map_electrons["_scale_down"].end(),PtComparing);
      map_parameter["_scale_down"]=p.Clone(MakeLeptonPointerVector(map_electrons["_scale_down"]));
      
      map_electrons["_smear_up"]=SmearElectrons(map_electrons[""],1);
      std::sort(map_electrons["_smear_up"].begin(),map_electrons["_smear_up"].end(),PtComparing);
      map_parameter["_smear_up"]=p.Clone(MakeLeptonPointerVector(map_electrons["_smear_up"]));
      
      map_electrons["_smear_down"]=SmearElectrons(map_electrons[""],-1);
      std::sort(map_electrons["_smear_down"].begin(),map_electrons["_smear_down"].end(),PtComparing);
      map_parameter["_smear_down"]=p.Clone(MakeLeptonPointerVector(map_electrons["_smear_down"]));

      map_electrons["_eta2p5"]=ElectronEnergyCorrection(SMPGetElectrons("passMediumID",0.0,2.5),0,0);
      map_parameter["_eta2p5"]=p.Clone(MakeLeptonPointerVector(map_electrons["_eta2p5"]));

      map_parameter["_noroccor"]=p.Clone(MakeLeptonPointerVector(map_electrons["_noroccor"]));
      
      map_electrons["_noEcor"]=ElectronEnergyCorrection(map_electrons[""],-1,0);
      map_parameter["_noEcor"]=p.Clone(MakeLeptonPointerVector(map_electrons["_noEcor"]));

    }
  }else if(channelname.Contains(TRegexp("em20[0-9][0-9]")) || channelname.Contains("Trig") || channelname.Contains("Dilep")){
    Parameter p;
    p.electronIDSF="ID_SF_MediumID_Q";
    p.muonIDSF="IDISO_SF_MediumID_trkIsoLoose_Q";
    p.triggerSF={"",""};
    p.lep0ptcut=25.;
    p.lep1ptcut=15.;

    map_electrons["_noroccor"]=electrons;
    map_electrons[""]=ElectronEnergyCorrection(map_electrons["_noroccor"],0,0);
    map_muons[""]=MuonMomentumCorrection(muons,0);

    std::vector<Lepton *> emu = MakeLeptonPointerVector(map_electrons[""]);
    std::vector<Lepton *> mu = MakeLeptonPointerVector(map_muons[""]);
    emu.insert(emu.end(),mu.begin(),mu.end());
    std::sort(emu.begin(),emu.end(),PtComparingPtr);

    map_parameter[""]=p.Clone(emu,
			      (IsNominalRun?NominalWeight:0)
			      +(HasFlag("SYS")&&!IsDATA?SystematicWeight:0)
			      +(HasFlag("PDFSYS")&&!IsDATA?PDFWeight:0)
			      );

  }else{
    cout<<"[AFBAnalyzer::executeEventWithPrefix] wrong channelname"<<endl;
    return;
  }
  
  ///////////////////////lepton selection///////////////////////
  for(const auto& [suffix,p]:map_parameter){
    TString prefix=tauprefix;
    double eventweight=lumiweight*PUweight*prefireweight*z0weight*zptweight*costhetaweight;

    if(p.weightbit&NominalWeight) FillHist(channelname+"/"+prefix+"nlepton"+suffix,p.leps.size(),eventweight,10,0,10);
    if(p.leps.size()>=2){
      if(HasFlag("REGION_cf")){
	if(p.leps.at(0)->Charge()>0&&p.leps.at(1)->Charge()>0) prefix="pp_"+prefix;
	else if(p.leps.at(0)->Charge()<0&&p.leps.at(1)->Charge()<0) prefix="mm_"+prefix;
	else continue;
      }else{
	if(p.leps.at(0)->Charge()*p.leps.at(1)->Charge()>0) prefix="ss_"+prefix;
      }
      if(channelname.Contains(TRegexp("em20[0-9][0-9]"))){
	if(p.leps.at(0)->LeptonFlavour() == p.leps.at(1)->LeptonFlavour()) continue;
      }
      if(p.weightbit&NominalWeight) FillCutflow(channelname+"/"+prefix+"cutflow"+suffix,"dilepton",eventweight);
      if(p.leps.at(0)->Pt()>p.lep0ptcut&&p.leps.at(1)->Pt()>p.lep1ptcut){
	if(p.weightbit&NominalWeight) FillCutflow(channelname+"/"+prefix+"cutflow"+suffix,"ptcut",eventweight);

	FillHist(channelname+"/"+prefix+"nbjet"+suffix,bjets.size(),eventweight,7,0,7);
	if(n_bjet!=1) return;
        if(p.weightbit&NominalWeight) FillCutflow(channelname+"/"+prefix+"cutflow"+suffix,"onebjet",eventweight);
	bjet_charge=bjets.at(0).Charge();
	bjet_rap=bjets.at(0).Rapidity();

	if(pfMET_Type1_pt>60) return;
	if(p.weightbit&NominalWeight) FillCutflow(channelname+"/"+prefix+"cutflow"+suffix,"MET60Cut",eventweight);
	/*
	if(bjets.at(0).Pt()<30){
	  if(bjets.at(0).GetPileupJetId() <0.18) return;
	}else if(bjets.at(0).Pt()<50){
	  if(bjets.at(0).GetPileupJetId() <0.61) return;
	}  
	if(p.weightbit&NominalWeight) FillCutflow(channelname+"/"+prefix+"cutflow"+suffix,"bjetPUIDCut",eventweight);

	if(min(p.leps.at(0)->DeltaR(bjets.at(0)),p.leps.at(1)->DeltaR(bjets.at(0)))<0.4) return;
	if(p.weightbit&NominalWeight) FillCutflow(channelname+"/"+prefix+"cutflow"+suffix,"bjetcleaning",eventweight);
	*/
	TLorentzVector *l0,*l1;
	l0=p.leps.at(0);
	l1=p.leps.at(1);
	if(abs((*l0+*l1).DeltaPhi(bjets.at(0)))<1.57) return;
        if(p.weightbit&NominalWeight) FillCutflow(channelname+"/"+prefix+"cutflow"+suffix,"ZbdPhiCut",eventweight);

        FillHist(channelname+"/"+prefix+"nrealjet",realjets.size(),eventweight,10,0,10);  
	FillHist(channelname+"/"+prefix+"njet",jets.size(),eventweight,10,0,10);

	if(realjets.size()>2) return;
	if(p.weightbit&NominalWeight) FillCutflow(channelname+"/"+prefix+"cutflow"+suffix,"2realjet",eventweight);

	bmuon.clear();
	for(unsigned int l=0; l<softmus.size(); l++){
	  if(bjets.at(0).DeltaR(softmus.at(l))<0.4 && (softmus.at(l).Pt())/(bjets.at(0).Pt())<0.6 && min(p.leps.at(0)->DeltaR(softmus.at(l)),p.leps.at(1)->DeltaR(softmus.at(l)))>0.4 ) bmuon.push_back(softmus.at(l));
	  //if(softmus.at(l).RelIso()>0.3) bmuon.push_back(softmus.at(l));
	}
	if(bmuon.size() >0){
	  for(unsigned int l=0; l<bmuon.size(); l++){
	    TString charge = "";
	    if(bmuon.at(l).Charge() > 0) charge = "P";
	    else charge = "M";
	    FillHist(channelname+"/"+prefix+"bmuon"+Form("%d",l)+charge+"_PFiso",bmuon.at(l).RelIso(),eventweight,20,0,2);
	    FillHist(channelname+"/"+prefix+"bmuon"+Form("%d",l)+charge+"_Trkiso",(bmuon.at(l).TrkIso())/(bmuon.at(l).Pt()),eventweight,20,0,2);
	    FillHist(channelname+"/"+prefix+"bmuon"+Form("%d",l)+charge+"_pt",bmuon.at(l).Pt(),eventweight,100,0,100);
	    FillHist(channelname+"/"+prefix+"bmuon"+Form("%d",l)+charge+"_ptRatio",(bmuon.at(l).Pt())/(bjets.at(0).Pt()),eventweight,20,0,1);
	    //FillHist(channelname+"/"+prefix+"bmuon"+Form("%d",l)+"_charge",bmuon.at(l).Charge(),eventweight,4,-2,2);
	    FillHist(channelname+"/"+prefix+"bmuon"+Form("%d",l)+charge+"_deltaRlep",min(p.leps.at(0)->DeltaR(bmuon.at(l)),p.leps.at(1)->DeltaR(bmuon.at(l))),eventweight,20,0,2);
	    FillHist(channelname+"/"+prefix+"bmuon"+Form("%d",l)+charge+"_deltaR",bmuon.at(l).DeltaR(bjets.at(0)),eventweight,60,0,3);
	    FillHist(channelname+"/"+prefix+"bmuon"+Form("%d",l)+charge+"_PassMedium",bmuon.at(l).PassID("POGMedium"),eventweight,2,0,2);
	    FillHist(channelname+"/"+prefix+"bmuon"+Form("%d",l)+charge+"_PassTight",bmuon.at(l).PassID("POGTight"),eventweight,2,0,2);
	    FillHist(channelname+"/"+prefix+"bmuon"+Form("%d",l)+charge+"_IP3D",abs(bmuon.at(l).IP3D()),eventweight,20,0,10);
	    FillHist(channelname+"/"+prefix+"bmuon"+Form("%d",l)+charge+"_SIP3D",abs(bmuon.at(l).IP3D())/bmuon.at(l).IP3Derr(),eventweight,20,0,10);
	    if(l==0) FillHist(channelname+"/"+prefix+"bmuon"+Form("%d",l)+charge+"_jetcharge",bjets.at(0).Charge(),eventweight,200,-2,2);
	  }
	}

        //if(bmuon.size()==0 && (bjet_charge > -0.3 && bjet_charge < 0.3)) return;
        if(bjet_charge > -0.3 && bjet_charge < 0.3) return;
	if(p.weightbit&NominalWeight) FillCutflow(channelname+"/"+prefix+"cutflow"+suffix,"bjetcharge030Cut",eventweight);

	/////////////////efficiency scale factors///////////////////
	double IDSF=1.,IDSF_up=1.,IDSF_down=1.;
	double ISOSF=1.,ISOSF_up=1.,ISOSF_down=1.;
	double RECOSF=1.,RECOSF_up=1.,RECOSF_down=1.;
	if(!IsDATA){
	  for(const auto& lep:p.leps){
	    TString LeptonIDSF_key="";
            if(lep->LeptonFlavour()==Lepton::ELECTRON){
	      LeptonIDSF_key=p.electronIDSF;

              double this_pt,this_eta;
              this_pt=((Electron*)lep)->UncorrPt();
              this_eta=((Electron*)lep)->scEta();

              double this_RECOSF=mcCorr->ElectronReco_SF(this_eta,this_pt,0);
              double this_RECOSF_up=mcCorr->ElectronReco_SF(this_eta,this_pt,1);
              double this_RECOSF_down=mcCorr->ElectronReco_SF(this_eta,this_pt,-1);
              RECOSF*=this_RECOSF; RECOSF_up*=this_RECOSF_up; RECOSF_down*=this_RECOSF_down;
	    }else if(lep->LeptonFlavour()==Lepton::MUON){
	      LeptonIDSF_key=p.muonIDSF;

	      double this_ISOSF=Lepton_SF(p.muonISOSF,lep,0);
	      double this_ISOSF_up=Lepton_SF(p.muonISOSF,lep,1);
	      double this_ISOSF_down=Lepton_SF(p.muonISOSF,lep,-1);
	      ISOSF*=this_ISOSF; ISOSF_up*=this_ISOSF_up; ISOSF_down*=this_ISOSF_down;
	    }

            double this_IDSF=Lepton_SF(LeptonIDSF_key,lep,0);
            double this_IDSF_up=Lepton_SF(LeptonIDSF_key,lep,1);
            double this_IDSF_down=Lepton_SF(LeptonIDSF_key,lep,-1);
            IDSF*=this_IDSF; IDSF_up*=this_IDSF_up; IDSF_down*=this_IDSF_down;
          }
	}
      
	double triggerSF=1.,triggerSF_up=1.,triggerSF_down=1.;
	if(!IsDATA){
          if(p.triggerSF.size()==1){
            triggerSF*=LeptonTrigger_SF(p.triggerSF[0],p.leps,0);
            triggerSF_up*=LeptonTrigger_SF(p.triggerSF[0],p.leps,1);
            triggerSF_down*=LeptonTrigger_SF(p.triggerSF[0],p.leps,-1);
          }else if(p.triggerSF.size()==2){
            triggerSF*=DileptonTrigger_SF(p.triggerSF[0],p.triggerSF[1],p.leps,0);
            triggerSF_up*=DileptonTrigger_SF(p.triggerSF[0],p.triggerSF[1],p.leps,1);
            triggerSF_down*=DileptonTrigger_SF(p.triggerSF[0],p.triggerSF[1],p.leps,-1);
          }
	}

	if(p.weightbit&NominalWeight){
	  FillCutflow(channelname+"/"+prefix+"cutflow"+suffix,"RECO",eventweight*RECOSF);
	  FillCutflow(channelname+"/"+prefix+"cutflow"+suffix,"ID",eventweight*RECOSF*IDSF);
	  FillCutflow(channelname+"/"+prefix+"cutflow"+suffix,"ISO",eventweight*RECOSF*IDSF*ISOSF);
	  FillCutflow(channelname+"/"+prefix+"cutflow"+suffix,"trigger",eventweight*RECOSF*IDSF*ISOSF*triggerSF);
	}

	///////////////////////map_weight//////////////////
	map<TString,double> map_weight;
	if(p.weightbit&NominalWeight){

	  map_weight[""]=lumiweight*PUweight*prefireweight*zptweight*z0weight*costhetaweight*RECOSF*IDSF*ISOSF*triggerSF;
	}
	if(p.weightbit&SystematicWeight){
	  map_weight["_noPUweight"]=lumiweight*prefireweight*zptweight*z0weight*costhetaweight*RECOSF*IDSF*ISOSF*triggerSF;
	  map_weight["_PUweight_up"]=lumiweight*PUweight_up*prefireweight*zptweight*z0weight*costhetaweight*RECOSF*IDSF*ISOSF*triggerSF;
	  map_weight["_PUweight_down"]=lumiweight*PUweight_down*prefireweight*zptweight*z0weight*costhetaweight*RECOSF*IDSF*ISOSF*triggerSF;
	  
	  map_weight["_noprefireweight"]=lumiweight*PUweight*zptweight*z0weight*costhetaweight*RECOSF*IDSF*ISOSF*triggerSF;
	  map_weight["_prefireweight_up"]=lumiweight*PUweight*prefireweight_up*zptweight*z0weight*costhetaweight*RECOSF*IDSF*ISOSF*triggerSF;
	  map_weight["_prefireweight_down"]=lumiweight*PUweight*prefireweight_down*zptweight*z0weight*costhetaweight*RECOSF*IDSF*ISOSF*triggerSF;

	  map_weight["_nozptweight"]=lumiweight*PUweight*prefireweight*z0weight*costhetaweight*RECOSF*IDSF*ISOSF*triggerSF;
	  
	  map_weight["_noz0weight"]=lumiweight*PUweight*prefireweight*zptweight*costhetaweight*RECOSF*IDSF*ISOSF*triggerSF;

	  map_weight["_nocosthetaweight"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF;
	  map_weight["_costhetaweight_up"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*costhetaweight_up*RECOSF*IDSF*ISOSF*triggerSF;
	  map_weight["_costhetaweight_down"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*costhetaweight_down*RECOSF*IDSF*ISOSF*triggerSF;
	  
	  map_weight["_noefficiencySF"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*costhetaweight;
	  
	  map_weight["_noRECOSF"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*costhetaweight*IDSF*ISOSF*triggerSF;
	  map_weight["_RECOSF_up"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*costhetaweight*RECOSF_up*IDSF*ISOSF*triggerSF;
	  map_weight["_RECOSF_down"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*costhetaweight*RECOSF_down*IDSF*ISOSF*triggerSF;
	  
	  map_weight["_noIDSF"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*costhetaweight*RECOSF*ISOSF*triggerSF;
	  map_weight["_IDSF_up"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*costhetaweight*RECOSF*IDSF_up*ISOSF*triggerSF;
	  map_weight["_IDSF_down"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*costhetaweight*RECOSF*IDSF_down*ISOSF*triggerSF;
	  
	  map_weight["_noISOSF"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*costhetaweight*RECOSF*IDSF*triggerSF;
	  map_weight["_ISOSF_up"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*costhetaweight*RECOSF*IDSF*ISOSF_up*triggerSF;
	  map_weight["_ISOSF_down"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*costhetaweight*RECOSF*IDSF*ISOSF_down*triggerSF;
	  
	  map_weight["_notriggerSF"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*costhetaweight*RECOSF*IDSF*ISOSF;
	  map_weight["_triggerSF_up"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*costhetaweight*RECOSF*IDSF*ISOSF*triggerSF_up;
	  map_weight["_triggerSF_down"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*costhetaweight*RECOSF*IDSF*ISOSF*triggerSF_down;
	  
	  map_weight["_nocosthetaweight"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF;
	  map_weight["_costhetaweight_up"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*costhetaweight_up*RECOSF*IDSF*ISOSF*triggerSF;
	  map_weight["_costhetaweight_down"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*costhetaweight_down*RECOSF*IDSF*ISOSF*triggerSF;

	}
	if(p.weightbit&PDFWeight){
	  for(unsigned int i=0;i<PDFWeights_Scale->size();i++){
	    map_weight[Form("_scalevariation%d",i)]=lumiweight*PUweight*prefireweight*zptweight*z0weight*costhetaweight*RECOSF*IDSF*ISOSF*triggerSF*PDFWeights_Scale->at(i);
	  }
	  for(unsigned int i=0;i<PDFWeights_Error->size();i++){
	    map_weight[Form("_pdf%d",i)]=lumiweight*PUweight*prefireweight*zptweight*z0weight*costhetaweight*RECOSF*IDSF*ISOSF*triggerSF*PDFWeights_Error->at(i);
	  }
	  if(PDFWeights_AlphaS->size()==2){
	    map_weight["_alphaS_up"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*costhetaweight*RECOSF*IDSF*ISOSF*triggerSF*PDFWeights_AlphaS->at(0);
	    map_weight["_alphaS_down"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*costhetaweight*RECOSF*IDSF*ISOSF*triggerSF*PDFWeights_AlphaS->at(1);
	  }
	}

	///////////////////////fill hists///////////////////////
        //if(bjet_charge > -0.33 && bjet_charge < 0.33) prefix += "bbmixed_";
	if(HasFlag("TOY")) FillHistsToy(channelname,prefix,suffix,(Particle*)p.leps[0],(Particle*)p.leps[1],map_weight);
	else{
	  FillHists(channelname,prefix,suffix,(Particle*)p.leps[0],(Particle*)p.leps[1],map_weight);
	  FillHist(channelname+"/"+prefix+"bjetpT",bjets.at(0).Pt(),map_weight,200,0,1000);
	  FillHist(channelname+"/"+prefix+"bjeteta",bjets.at(0).Eta(),map_weight,100,-5,5);
	  FillHist(channelname+"/"+prefix+"bjetM",bjets.at(0).M(),map_weight,200,0,20);
	  FillHist(channelname+"/"+prefix+"bjetCharge",bjet_charge,map_weight,200,-2,2);
	  FillHist(channelname+"/"+prefix+"bjetPUID",bjets.at(0).GetPileupJetId(),map_weight,200,-2,2);

          FillHist(channelname+"/"+prefix+"yZ",(*l0+*l1).Rapidity(),map_weight,60,-3,3);
          FillHist(channelname+"/"+prefix+"yb",bjets.at(0).Rapidity(),map_weight,60,-3,3);
	  FillHist(channelname+"/"+prefix+"Zb_dy",(*l0+*l1).Rapidity()-bjets.at(0).Rapidity(),map_weight,100,-5,5);
	  FillHist(channelname+"/"+prefix+"Zb_dphi",(*l0+*l1).DeltaPhi(bjets.at(0)),map_weight,100,-5,5);
          FillHist(channelname+"/"+prefix+"Zb_y2D",(*l0+*l1).Rapidity(),bjets.at(0).Rapidity(),map_weight,30,-3,3,30,-3,3);

	  /*
	  for(unsigned int j=0;j<jets.size();j++){
	    FillHist(channelname+"/"+prefix+"jet_"+Form("%i",j)+"_pT",jets.at(j).Pt(),map_weight,200,0,1000);
            FillHist(channelname+"/"+prefix+"jet_"+Form("%i",j)+"_dpT",jets.at(j).Pt()-bjets.at(0).Pt(),map_weight,200,-100,100);
            FillHist(channelname+"/"+prefix+"jet_"+Form("%i",j)+"_PUID",jets.at(j).GetPileupJetId(),map_weight,200,-2,2);
	  }
	  for(unsigned int j=0;j<realjets.size();j++){
	    FillHist(channelname+"/"+prefix+"realjet_"+Form("%i",j)+"_pT",realjets.at(j).Pt(),map_weight,200,0,1000);
            FillHist(channelname+"/"+prefix+"realjet_"+Form("%i",j)+"_dpT",realjets.at(j).Pt()-bjets.at(0).Pt(),map_weight,200,-100,100);
            FillHist(channelname+"/"+prefix+"realjet_"+Form("%i",j)+"_PUID",realjets.at(j).GetPileupJetId(),map_weight,200,-2,2);
	    }*/
	  if(IsDYSample&&prefix==""&&IsNominalRun){
	    vector<Gen> gens=GetGens();
	    Gen truth_l0=GetGenMatchedLepton(*p.leps[0],gens);
	    Gen truth_l1=GetGenMatchedLepton(*p.leps[1],gens);
	    if(!truth_l0.IsEmpty()&&!truth_l1.IsEmpty()) 
		FillHists(channelname,"truth_",suffix,(Particle*)&truth_l0,(Particle*)&truth_l1,map_weight);
	    //else cout<<"no matching"<<endl;
	  }
	}
      }
    }
  }
}
AFBAnalyzer::AFBAnalyzer(){}
AFBAnalyzer::~AFBAnalyzer(){
  DeleteToy();
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
  if(direction==0) direction=(0.75*dilepton.Rapidity()-bjet_rap)>0?1:-1;
  if(bjet_charge>0.3) direction *= -1.;
  return direction*2*(l0pp*l1pm-l0pm*l1pp)/sqrt(dimass*dimass*(dimass*dimass+dipt*dipt));
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
void AFBAnalyzer::FillHistsToy(TString channelname,TString pre,TString suf,Particle* l0,Particle* l1,map<TString,double> map_weight){
  int n_toy=toy_random.size();
  for(int i=0;i<n_toy;i++) FillHists(channelname,pre,suf+Form("_toy%d",i),l0,l1,Multiply(map_weight,toy_weight[i]));
}
void AFBAnalyzer::FillHists(TString channelname,TString pre,TString suf,Particle* l0,Particle* l1,map<TString,double> map_weight){
  TLorentzVector dilepton=(*l0)+(*l1);
  double dimass=dilepton.M();
  double dirap=dilepton.Rapidity();
  double dipt=dilepton.Pt();

  double cost=GetCosThetaCS(l0,l1);
  double h=0.5*pow(dipt/dimass,2)/(1+pow(dipt/dimass,2))*(1-3*cost*cost);
  double den_weight=0.5*fabs(cost)/pow(1+cost*cost+h,2);
  double num_weight=0.5*cost*cost/pow(1+cost*cost+h,3);
  FillHist(channelname+"/"+pre+"costhetaCS"+suf,dimass,dirap,dipt,cost,map_weight,afb_mbinnum,(double*)afb_mbin,afb_ybinnum,(double*)afb_ybin,afb_ptbinnum,(double*)afb_ptbin,20,-1,1);
  FillHist(channelname+"/"+pre+"costhetaCS_den"+suf,dimass,dirap,dipt,cost,Multiply(map_weight,den_weight),afb_mbinnum,(double*)afb_mbin,afb_ybinnum,(double*)afb_ybin,afb_ptbinnum,(double*)afb_ptbin,20,-1,1);
  FillHist(channelname+"/"+pre+"costhetaCS_num"+suf,dimass,dirap,dipt,cost,Multiply(map_weight,num_weight),afb_mbinnum,(double*)afb_mbin,afb_ybinnum,(double*)afb_ybin,afb_ptbinnum,(double*)afb_ptbin,20,-1,1);

  if(!HasFlag("PDFSYS")){
    FillHist(channelname+"/"+pre+"l0pt"+suf,dimass,dirap,dipt,l0->Pt(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,lptbinnum,(double*)lptbin);
    FillHist(channelname+"/"+pre+"l1pt"+suf,dimass,dirap,dipt,l1->Pt(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,lptbinnum,(double*)lptbin);
    FillHist(channelname+"/"+pre+"lpt"+suf,dimass,dirap,dipt,l0->Pt(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,lptbinnum,(double*)lptbin);
    FillHist(channelname+"/"+pre+"lpt"+suf,dimass,dirap,dipt,l1->Pt(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,lptbinnum,(double*)lptbin);
    
    FillHist(channelname+"/"+pre+"l0eta"+suf,dimass,dirap,dipt,l0->Eta(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,60,-3,3);
    FillHist(channelname+"/"+pre+"l1eta"+suf,dimass,dirap,dipt,l1->Eta(),map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,60,-3,3);

    if(!pre.Contains("gen")&&!pre.Contains("lhe")&&!pre.Contains("truth")){
      FillHist(channelname+"/"+pre+"z0"+suf,dimass,dirap,dipt,vertex_Z,map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,120,-15,15);
      FillHist(channelname+"/"+pre+"met"+suf,dimass,dirap,dipt,pfMET_Type1_pt,map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,100,0,200);
      FillHist(channelname+"/"+pre+"nPV"+suf,dimass,dirap,dipt,nPV,map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,60,0,60);  
      FillHist(channelname+"/"+pre+"nPileUp"+suf,dimass,dirap,dipt,nPileUp,map_weight,grid_mbinnum,(double*)grid_mbin,grid_ybinnum,(double*)grid_ybin,grid_ptbinnum,(double*)grid_ptbin,60,0,60);  
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
  ifstream file_check(datapath+"/"+TString::Itoa(DataYear,10)+"/CosTheta/CosThetaWeight.root");
  bool isexist=file_check.is_open();
  file_check.close();
  if(!isexist){
    cout<<"[AFBAnalyzer::SetupCosThetaWeight] no CosThetaWeight.root"<<endl;
    return;
  }
  TFile fcost(datapath+"/"+TString::Itoa(DataYear,10)+"/CosTheta/CosThetaWeight.root");
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
  MCName=Replace(MCName,"DY[0-9]Jets","DYJets");
  TString hname=MCName+suffix;
  auto it=map_hist_cost.find(hname);
  if(it!=map_hist_cost.end())
    val*=GetBinContentUser(it->second,mass,pt,cost,0);
  if(val==0) val=1.;
  return val;
}
