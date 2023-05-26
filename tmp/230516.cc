#include"Plotter.cc"

void plot(TString sim="mi"){
  Plotter::SetupStyle();
  TH1::AddDirectory(0);
  TFile f("LTAnalyzer/"+sim+".root");
  {
    TH1* unfolded_j0_a0=(TH1*)f.Get("unfolded_j0_a0");
    TH1* unfolded_j0_a2=(TH1*)f.Get("unfolded_j0_a2");
    TH1* gen_j0_a0=(TH1*)f.Get("gen_j0_a0");
    TH1* gen_j0_a2=(TH1*)f.Get("gen_j0_a2");
    
    TCanvas *c=new TCanvas;
    unfolded_j0_a0->SetStats(0);
    unfolded_j0_a0->Draw("hist e");
    unfolded_j0_a0->GetXaxis()->SetTitle("p_{T}(ll) [GeV]");
    unfolded_j0_a0->GetYaxis()->SetRangeUser(0,1);
    unfolded_j0_a0->GetYaxis()->SetTitle("A_{i}");
    unfolded_j0_a0->SetTitle("");    
    unfolded_j0_a0->SetLineStyle(2);
    unfolded_j0_a0->SetMarkerStyle(26);
    unfolded_j0_a2->Draw("same hist e");
    unfolded_j0_a2->SetLineColor(2);
    unfolded_j0_a2->SetLineStyle(2);
    unfolded_j0_a2->SetMarkerStyle(26);
    unfolded_j0_a2->SetMarkerColor(2);
    gen_j0_a0->Draw("same hist e");
    gen_j0_a0->SetMarkerStyle(24);
    gen_j0_a2->Draw("same hist e");
    gen_j0_a2->SetLineColor(2);
    gen_j0_a2->SetMarkerStyle(24);
    gen_j0_a2->SetMarkerColor(2);
    TLegend *leg=new TLegend(0.6,0.5,0.89,0.15);
    leg->AddEntry(gen_j0_a0,"true A_{0}");
    leg->AddEntry(gen_j0_a2,"true A_{2}");
    leg->AddEntry(unfolded_j0_a0,"unfolded A_{0}");
    leg->AddEntry(unfolded_j0_a2,"unfolded A_{2}");
    leg->Draw();
    TLatex tex;
    tex.SetNDC(true);
    tex.DrawLatex(0.7,0.6,"N_{j} #leq 1");
  }
  {
    TH1* unfolded_j1_a0=(TH1*)f.Get("unfolded_j1_a0");
    TH1* unfolded_j1_a2=(TH1*)f.Get("unfolded_j1_a2");
    TH1* gen_j1_a0=(TH1*)f.Get("gen_j1_a0");
    TH1* gen_j1_a2=(TH1*)f.Get("gen_j1_a2");
    
    TCanvas *c=new TCanvas;
    unfolded_j1_a0->SetStats(0);
    unfolded_j1_a0->Draw("hist e");
    unfolded_j1_a0->GetXaxis()->SetTitle("p_{T}(ll) [GeV]");
    unfolded_j1_a0->GetYaxis()->SetRangeUser(0,1);
    unfolded_j1_a0->GetYaxis()->SetTitle("A_{i}");
    unfolded_j1_a0->SetTitle("");    
    unfolded_j1_a0->SetLineStyle(2);
    unfolded_j1_a0->SetMarkerStyle(26);
    unfolded_j1_a2->Draw("same hist e");
    unfolded_j1_a2->SetLineColor(2);
    unfolded_j1_a2->SetLineStyle(2);
    unfolded_j1_a2->SetMarkerStyle(26);
    unfolded_j1_a2->SetMarkerColor(2);
    gen_j1_a0->Draw("same hist e");
    gen_j1_a0->SetMarkerStyle(24);
    gen_j1_a2->Draw("same hist e");
    gen_j1_a2->SetLineColor(2);
    gen_j1_a2->SetMarkerStyle(24);
    gen_j1_a2->SetMarkerColor(2);
    TLegend *leg=new TLegend(0.6,0.5,0.89,0.15);
    leg->AddEntry(gen_j1_a0,"true A_{0}");
    leg->AddEntry(gen_j1_a2,"true A_{2}");
    leg->AddEntry(unfolded_j1_a0,"unfolded A_{0}");
    leg->AddEntry(unfolded_j1_a2,"unfolded A_{2}");
    leg->Draw();
    TLatex tex;
    tex.SetNDC(true);
    tex.DrawLatex(0.7,0.6,"N_{j} > 1");
  }
}
void plot2(TString sim="mi"){
  Plotter::SetupStyle();
  TH1::AddDirectory(0);
  TFile f("LTAnalyzer/"+sim+".root");
  TH1* unfolded_j0_a0=(TH1*)f.Get("unfolded_j0_a0");
  TH1* unfolded_j0_a2=(TH1*)f.Get("unfolded_j0_a2");
  TH1* gen_j0_a0=(TH1*)f.Get("gen_j0_a0");
  TH1* gen_j0_a2=(TH1*)f.Get("gen_j0_a2");
  TH1* unfolded_j1_a0=(TH1*)f.Get("unfolded_j1_a0");
  TH1* unfolded_j1_a2=(TH1*)f.Get("unfolded_j1_a2");
  TH1* gen_j1_a0=(TH1*)f.Get("gen_j1_a0");
  TH1* gen_j1_a2=(TH1*)f.Get("gen_j1_a2");
  unfolded_j0_a0->Add(unfolded_j0_a2,-1);
  unfolded_j1_a0->Add(unfolded_j1_a2,-1);
  gen_j0_a0->Add(gen_j0_a2,-1);
  gen_j1_a0->Add(gen_j1_a2,-1);
    
  TCanvas *c=new TCanvas;
  unfolded_j0_a0->SetStats(0);
  unfolded_j0_a0->Draw("hist e");
  unfolded_j0_a0->GetXaxis()->SetTitle("p_{T}(ll) [GeV]");
  unfolded_j0_a0->GetYaxis()->SetRangeUser(-0.14,0.5);
  unfolded_j0_a0->GetYaxis()->SetTitle("A_{0}-A_{2}");
  unfolded_j0_a0->SetTitle("");    
  unfolded_j0_a0->SetLineStyle(2);
  unfolded_j0_a0->SetMarkerStyle(26);
  unfolded_j1_a0->Draw("same hist e");
  unfolded_j1_a0->SetLineColor(2);
  unfolded_j1_a0->SetLineStyle(2);
  unfolded_j1_a0->SetMarkerStyle(26);
  unfolded_j1_a0->SetMarkerColor(2);
  gen_j0_a0->Draw("same hist e");
  gen_j0_a0->SetMarkerStyle(24);
  gen_j1_a0->Draw("same hist e");
  gen_j1_a0->SetLineColor(2);
  gen_j1_a0->SetMarkerStyle(24);
  gen_j1_a0->SetMarkerColor(2);
  TLegend *leg=new TLegend(0.6,0.5,0.89,0.89);
  leg->AddEntry(gen_j0_a0,"true A_{0}-A_{2}, N_{j} #leq 1");
  leg->AddEntry(gen_j1_a0,"true A_{0}-A_{2}, N_{j} > 1");
  leg->AddEntry(unfolded_j0_a0,"unfolded A_{0}-A_{2}, N_{j} #leq 1");
  leg->AddEntry(unfolded_j1_a0,"unfolded A_{0}-A_{2}, N_{j} > 1");
  leg->Draw();
}
void plot_sim_compare(){
  Plotter::SetupStyle();
  TH1::AddDirectory(0);
  TFile fmi("LTAnalyzer/mi.root");
  TFile fmg("LTAnalyzer/mg.root");
  {
    TH1* gen_a0_mi=(TH1*)fmi.Get("gen_j0_a0");
    TH1* gen_a2_mi=(TH1*)fmi.Get("gen_j0_a2");
    TH1* gen_a0_a2_mi=(TH1*)fmi.Get("gen_j0_a0");
    gen_a0_a2_mi->Add(gen_a2_mi,-1);
    TH1* gen_a0_mg=(TH1*)fmg.Get("gen_j0_a0");
    TH1* gen_a2_mg=(TH1*)fmg.Get("gen_j0_a2");
    TH1* gen_a0_a2_mg=(TH1*)fmg.Get("gen_j0_a0");
    gen_a0_a2_mg->Add(gen_a2_mg,-1);
    
    TCanvas *c=new TCanvas;
    gen_a0_mi->SetStats(0);
    gen_a0_mi->GetXaxis()->SetTitle("p_{T}(ll) [GeV]");
    gen_a0_mi->GetYaxis()->SetRangeUser(0,1);
    gen_a0_mi->GetYaxis()->SetTitle("A_{i} or A_{0}-A{2}");
    gen_a0_mi->SetTitle("");

    gen_a0_mi->Draw("hist e");
    gen_a0_mi->SetLineColor(1);
    gen_a0_mi->SetLineStyle(1);
    gen_a0_mi->SetMarkerColor(1);
    gen_a0_mi->SetMarkerStyle(24);
    
    gen_a2_mi->Draw("same hist e");
    gen_a2_mi->SetLineColor(2);
    gen_a2_mi->SetLineStyle(1);
    gen_a2_mi->SetMarkerColor(2);
    gen_a2_mi->SetMarkerStyle(24);

    gen_a0_a2_mi->Draw("same hist e");
    gen_a0_a2_mi->SetLineColor(4);
    gen_a0_a2_mi->SetLineStyle(1);
    gen_a0_a2_mi->SetMarkerColor(4);
    gen_a0_a2_mi->SetMarkerStyle(24);

    gen_a0_mg->Draw("same hist e");
    gen_a0_mg->SetLineColor(1);
    gen_a0_mg->SetLineStyle(2);
    gen_a0_mg->SetMarkerColor(1);
    gen_a0_mg->SetMarkerStyle(26);
    
    gen_a2_mg->Draw("same hist e");
    gen_a2_mg->SetLineColor(2);
    gen_a2_mg->SetLineStyle(2);
    gen_a2_mg->SetMarkerColor(2);
    gen_a2_mg->SetMarkerStyle(26);

    gen_a0_a2_mg->Draw("same hist e");
    gen_a0_a2_mg->SetLineColor(4);
    gen_a0_a2_mg->SetLineStyle(2);
    gen_a0_a2_mg->SetMarkerColor(4);
    gen_a0_a2_mg->SetMarkerStyle(26);

    TLegend *leg=new TLegend(0.6,0.5,0.89,0.15);
    leg->AddEntry(gen_a0_mi,"MiNNLO gen A_{0}");
    leg->AddEntry(gen_a2_mi,"MiNNLO gen A_{2}");
    leg->AddEntry(gen_a0_a2_mi,"MiNNLO gen A_{0}-A_{2}");
    leg->AddEntry(gen_a0_mg,"MG gen A_{0}");
    leg->AddEntry(gen_a2_mg,"MG gen A_{2}");
    leg->AddEntry(gen_a0_a2_mg,"MG gen A_{0}-A_{2}");
    leg->Draw();
    TLatex tex;
    tex.SetNDC(true);
    tex.DrawLatex(0.7,0.6,"N_{j} #leq 1");
  }
  {
    TH1* gen_a0_mi=(TH1*)fmi.Get("gen_j1_a0");
    TH1* gen_a2_mi=(TH1*)fmi.Get("gen_j1_a2");
    TH1* gen_a0_a2_mi=(TH1*)fmi.Get("gen_j1_a0");
    gen_a0_a2_mi->Add(gen_a2_mi,-1);
    TH1* gen_a0_mg=(TH1*)fmg.Get("gen_j1_a0");
    TH1* gen_a2_mg=(TH1*)fmg.Get("gen_j1_a2");
    TH1* gen_a0_a2_mg=(TH1*)fmg.Get("gen_j1_a0");
    gen_a0_a2_mg->Add(gen_a2_mg,-1);
    
    TCanvas *c=new TCanvas;
    gen_a0_mi->SetStats(0);
    gen_a0_mi->GetXaxis()->SetTitle("p_{T}(ll) [GeV]");
    gen_a0_mi->GetYaxis()->SetRangeUser(0,1);
    gen_a0_mi->GetYaxis()->SetTitle("A_{i} or A_{0}-A{2}");
    gen_a0_mi->SetTitle("");

    gen_a0_mi->Draw("hist e");
    gen_a0_mi->SetLineColor(1);
    gen_a0_mi->SetLineStyle(1);
    gen_a0_mi->SetMarkerColor(1);
    gen_a0_mi->SetMarkerStyle(24);
    
    gen_a2_mi->Draw("same hist e");
    gen_a2_mi->SetLineColor(2);
    gen_a2_mi->SetLineStyle(1);
    gen_a2_mi->SetMarkerColor(2);
    gen_a2_mi->SetMarkerStyle(24);

    gen_a0_a2_mi->Draw("same hist e");
    gen_a0_a2_mi->SetLineColor(4);
    gen_a0_a2_mi->SetLineStyle(1);
    gen_a0_a2_mi->SetMarkerColor(4);
    gen_a0_a2_mi->SetMarkerStyle(24);

    gen_a0_mg->Draw("same hist e");
    gen_a0_mg->SetLineColor(1);
    gen_a0_mg->SetLineStyle(2);
    gen_a0_mg->SetMarkerColor(1);
    gen_a0_mg->SetMarkerStyle(26);
    
    gen_a2_mg->Draw("same hist e");
    gen_a2_mg->SetLineColor(2);
    gen_a2_mg->SetLineStyle(2);
    gen_a2_mg->SetMarkerColor(2);
    gen_a2_mg->SetMarkerStyle(26);

    gen_a0_a2_mg->Draw("same hist e");
    gen_a0_a2_mg->SetLineColor(4);
    gen_a0_a2_mg->SetLineStyle(2);
    gen_a0_a2_mg->SetMarkerColor(4);
    gen_a0_a2_mg->SetMarkerStyle(26);

    TLegend *leg=new TLegend(0.6,0.5,0.89,0.15);
    leg->AddEntry(gen_a0_mi,"MiNNLO gen A_{0}");
    leg->AddEntry(gen_a2_mi,"MiNNLO gen A_{2}");
    leg->AddEntry(gen_a0_a2_mi,"MiNNLO gen A_{0}-A_{2}");
    leg->AddEntry(gen_a0_mg,"MG gen A_{0}");
    leg->AddEntry(gen_a2_mg,"MG gen A_{2}");
    leg->AddEntry(gen_a0_a2_mg,"MG gen A_{0}-A_{2}");
    leg->Draw();
    TLatex tex;
    tex.SetNDC(true);
    tex.DrawLatex(0.7,0.6,"N_{j} > 1");
  }
}
