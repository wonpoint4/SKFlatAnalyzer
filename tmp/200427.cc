void get_table(){
  AFBPlotter aa;
  DEBUG=0;
  TH1D* hists[2][3];
  hists[0][0]=new TH1D("mean_data","mean #mu",9,0,9);
  hists[0][1]=new TH1D("sigma_data","sigma #sigma",9,0,9);
  hists[0][2]=new TH1D("sigma_mean_data","#sigma/#mu",9,0,9);
  hists[1][0]=new TH1D("mean_mc","mean #mu",9,0,9);
  hists[1][1]=new TH1D("sigma_mc","sigma #sigma",9,0,9);
  hists[1][2]=new TH1D("sigma_mean_mc","#sigma/#mu",9,0,9);
  
  vector<TString> years={"2016","2017","2018"};
  vector<TString> suffixes={"_noEcor","_noroccor",""};
  for(int iyear=0;iyear<3;iyear++){
    TString year=years[iyear];
    for(int i=0;i<2;i++){
      for(int j=0;j<3;j++){
	cout<<year<<" ";
	if(i==0) cout<<"data ";
	else cout<<"MC ";
	if(j==0) cout<<"mean ";
	else if(j==1) cout<<"sigma ";
	else if(j==2) cout<<"s/m ";
	for(int isuf=0;isuf<3;isuf++){
	  TString suffix=suffixes[isuf];
	  TH1* hist=aa.GetTH1(aa.GetHist(i,"ee"+year+"/m[80,100]/dimass"+suffix,"widthweight"),true);
	  hist->Fit("gaus","q","",86,96);
	  double var,err;
	  if(j<2){
	    var=hist->GetFunction("gaus")->GetParameter(j+1);
	    err=hist->GetFunction("gaus")->GetParError(j+1);
	  }else if(j==2){
	    var=hist->GetFunction("gaus")->GetParameter(2)/hist->GetFunction("gaus")->GetParameter(1);
	    err=hist->GetFunction("gaus")->GetParError(2)/hist->GetFunction("gaus")->GetParameter(1);
	  }
	  cout<<var<<" ";
	  hists[i][j]->SetBinContent(iyear*3+isuf+1,var);
	  hists[i][j]->SetBinError(iyear*3+isuf+1,err);
	  //hists[i][j]->GetXaxis()->SetBinLabel(iyear*3+isuf+1,year+" "+suffix);
	}
	cout<<endl;
      }
    }
  }
  TCanvas* c1=new TCanvas;
  hists[0][0]->SetMarkerStyle(20);
  hists[0][0]->SetMarkerSize(1);
  hists[0][0]->SetStats(0);
  hists[0][0]->GetYaxis()->SetTitle("mean #mu [GeV]");
  hists[0][0]->GetYaxis()->SetRangeUser(90,91.5);
  hists[0][0]->GetXaxis()->SetNdivisions(303,0);
  hists[0][0]->GetXaxis()->SetLabelSize(0);
  hists[0][0]->Draw();

  hists[1][0]->SetMarkerStyle(22);
  hists[1][0]->SetMarkerSize(1);
  hists[1][0]->SetMarkerColor(2);
  hists[1][0]->SetLineColor(2);
  hists[1][0]->Draw("same");
  gPad->SetGridx();


  TCanvas* c2=new TCanvas;
  hists[0][1]->SetMarkerStyle(20);
  hists[0][1]->SetMarkerSize(1);
  hists[0][1]->SetStats(0);
  hists[0][1]->GetYaxis()->SetTitle("sigma #sigma [GeV]");
  hists[0][1]->GetYaxis()->SetRangeUser(2.5,3.0);
  hists[0][1]->GetXaxis()->SetNdivisions(303,0);
  hists[0][1]->GetXaxis()->SetLabelSize(0);
  hists[0][1]->Draw();

  hists[1][1]->SetMarkerStyle(22);
  hists[1][1]->SetMarkerSize(1);
  hists[1][1]->SetMarkerColor(2);
  hists[1][1]->SetLineColor(2);
  hists[1][1]->Draw("same");
  gPad->SetGridx();


  TCanvas* c3=new TCanvas;
  hists[0][2]->SetMarkerStyle(20);
  hists[0][2]->SetMarkerSize(1);
  hists[0][2]->SetStats(0);
  hists[0][2]->GetYaxis()->SetTitle("#sigma/#mu");
  hists[0][2]->GetYaxis()->SetRangeUser(0.027,0.034);
  hists[0][2]->GetXaxis()->SetNdivisions(303,0);
  hists[0][2]->GetXaxis()->SetLabelSize(0);
  hists[0][2]->Draw();

  hists[1][2]->SetMarkerStyle(22);
  hists[1][2]->SetMarkerSize(1);
  hists[1][2]->SetMarkerColor(2);
  hists[1][2]->SetLineColor(2);
  hists[1][2]->Draw("same");
  gPad->SetGridx();

}
