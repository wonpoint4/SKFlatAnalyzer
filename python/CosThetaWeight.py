import sys,os,math
import ROOT
def Check(hden,hnum):
    if not hden.InheritsFrom("TH4D"):
        print("Source is not TH4D... return NULL")
        return None
    if not hnum.InheritsFrom("TH4D"):
        print("Target is not TH4D... return NULL")
        return None
    if not hden.CheckConsistency(hnum):
        print("Source and Target are not consistent... return NULL")
        return None

def GetCosThetaWeight(hden,hnum,correlation=0):
    Check(hden,hnum)
    
    nx=hden.GetXaxis().GetNbins()
    ny=hden.GetYaxis().GetNbins()
    nz=hden.GetZaxis().GetNbins()-10
    nu=hden.GetUaxis().GetNbins()

    hout_all=ROOT.TH4D("","",
                       nx,hden.GetXaxis().GetXbins().GetArray(),
                       ny,hden.GetYaxis().GetXbins().GetArray(),
                       nz,hden.GetZaxis().GetXbins().GetArray(),
                       nu,hden.GetUaxis().GetXbins().GetArray())
    for ix in range(nx+2):
        stop=False
        for iy in range(ny+2):
            for iz in range(nz+2):
                this_hden=hden.ProjectionU("d_{}_{}_{}".format(ix,iy,iz),ix,ix,iy,iy,iz,iz)
                this_hnum=hnum.ProjectionU("n_{}_{}_{}".format(ix,iy,iz),ix,ix,iy,iy,iz,iz)
                if this_hden.Integral()==0 or this_hnum.Integral()==0: continue
                this_hden.Scale(1/this_hden.Integral())
                this_hnum.Scale(1/this_hnum.Integral())
                for iu in range(1,nu+1):
                    val_num=this_hnum.GetBinContent(iu)
                    err_num=this_hnum.GetBinError(iu)
                    val_den=this_hden.GetBinContent(iu)
                    err_den=this_hden.GetBinError(iu)
                    if val_num==0 or val_den==0:
                        val=0
                        err=0
                        #if ix!=0 and iz!=0: stop=True
                    else:
                        val=val_num/val_den
                        if correlation==0:
                            err=abs(val_num/val_den*math.sqrt(pow(err_num/val_num,2)+pow(err_den/val_den,2)))
                        else:
                            err=abs(val_num/val_den*(err_num/val_num-correlation*err_den/val_den))
                    hout_all.SetBinContent(ix,iy,iz,iu,val)
                    hout_all.SetBinError(ix,iy,iz,iu,err)
                    print(ix,iy,iz,iu,val,err)
                this_hnum.Delete()
                this_hden.Delete()
                if stop: break;
            if stop: break;
        if stop: break;
            
    nx=ix-1
    hout=ROOT.TH4D("","",
                   nx,hden.GetXaxis().GetXbins().GetArray(),
                   ny,hden.GetYaxis().GetXbins().GetArray(),
                   nz,hden.GetZaxis().GetXbins().GetArray(),
                   nu,hden.GetUaxis().GetXbins().GetArray())
    for ix in range(nx+1): # range(nx+1) -> do not use overflow bin
        for iy in range(ny+2):
            for iz in range(nz+2):
                for iu in range(nu+2):
                    hout.SetBinContent(ix,iy,iz,iu,hout_all.GetBinContent(ix,iy,iz,iu))
                    hout.SetBinError(ix,iy,iz,iu,hout_all.GetBinError(ix,iy,iz,iu))
    hout_all.Delete()
    hout.SetStats(0)
    return hout

def GetCosThetaWeight3D(hden,hnum,correlation=0):
    Check(hden,hnum)
    
    nx=hden.GetXaxis().GetNbins()-12
    ny=hden.GetYaxis().GetNbins()
    nz=hden.GetZaxis().GetNbins()-3
    nu=hden.GetUaxis().GetNbins()

    hout_all=ROOT.TH3D("","",
                       nx,hden.GetXaxis().GetXbins().GetArray(),
                       nz,hden.GetZaxis().GetXbins().GetArray(),
                       nu,hden.GetUaxis().GetXbins().GetArray())
    for ix in range(nx+2):
        stop=False
        for iz in range(nz+2):
            this_hden=hden.ProjectionU("d_{}_{}".format(ix,iz),ix,ix,0,-1,iz,iz)
            this_hnum=hnum.ProjectionU("n_{}_{}".format(ix,iz),ix,ix,0,-1,iz,iz)
            if this_hden.Integral()==0 or this_hnum.Integral()==0: 
                this_hnum.Delete()
                this_hden.Delete()
                continue
            this_hden.Scale(1/this_hden.Integral())
            this_hnum.Scale(1/this_hnum.Integral())
            for iu in range(1,nu+1):
                val_num=this_hnum.GetBinContent(iu)
                err_num=this_hnum.GetBinError(iu)
                val_den=this_hden.GetBinContent(iu)
                err_den=this_hden.GetBinError(iu)
                if val_num==0 or val_den==0:
                    val=0
                    err=0
                    print((ix,hden.GetXaxis().GetBinLowEdge(ix),hden.GetXaxis().GetBinLowEdge(ix+1)),(iz,hden.GetZaxis().GetBinLowEdge(iz),hden.GetZaxis().GetBinLowEdge(iz+1)),iu,val,err)
                    #stop=True
                else:
                    val=val_num/val_den
                    if correlation==0:
                        err=abs(val_num/val_den*math.sqrt(pow(err_num/val_num,2)+pow(err_den/val_den,2)))
                    else:
                        err=abs(val_num/val_den*(err_num/val_num-correlation*err_den/val_den))
                hout_all.SetBinContent(ix,iz,iu,val)
                hout_all.SetBinError(ix,iz,iu,err)
                #print((ix,hden.GetXaxis().GetBinLowEdge(ix),hden.GetXaxis().GetBinLowEdge(ix+1)),(iz,hden.GetZaxis().GetBinLowEdge(iz),hden.GetZaxis().GetBinLowEdge(iz+1)),iu,val,err)
            this_hnum.Delete()
            this_hden.Delete()
            if stop: break;
        if stop: break;
            
    nx=ix-1
    hout=ROOT.TH3D("","",
                   nx,hden.GetXaxis().GetXbins().GetArray(),
                   nz,hden.GetZaxis().GetXbins().GetArray(),
                   nu,hden.GetUaxis().GetXbins().GetArray())
    for ix in range(nx+1): # range(nx+1) -> do not use overflow bin
        for iz in range(nz+2):
            for iu in range(nu+2):
                hout.SetBinContent(ix,iz,iu,hout_all.GetBinContent(ix,iz,iu))
                hout.SetBinError(ix,iz,iu,hout_all.GetBinError(ix,iz,iu))
    hout_all.Delete()
    hout.SetStats(0)
    return hout

def GetCosThetaWeight2D(hden,hnum,axis="X",correlation=0,maxerror=0.1):
    Check(hden,hnum)
    
    if axis=="X":
        axis1=hden.GetXaxis()
    elif axis=="Y":
        axis1=hden.GetYaxis()
    elif axis=="Z":
        axis1=hden.GetZaxis()
    else:
        print("wrong axisname")
        return None
    axis2=hden.GetUaxis()

    n1=axis1.GetNbins()
    n2=axis2.GetNbins()

    hout_all=ROOT.TH2D()
    hout_all.SetBins(n1,axis1.GetXbins().GetArray(),n2,axis2.GetXbins().GetArray())
    for i in range(n1+2):
        if axis=="X":
            this_hden=hden.ProjectionU("d_{}".format(i),i,i,0,-1,0,-1)
            this_hnum=hnum.ProjectionU("n_{}".format(i),i,i,0,-1,0,-1)
        elif axis=="Y":
            this_hden=hden.ProjectionU("d_{}".format(i),0,-1,i,i,0,-1)
            this_hnum=hnum.ProjectionU("n_{}".format(i),0,-1,i,i,0,-1)
        elif axis=="Z":
            this_hden=hden.ProjectionU("d_{}".format(i),0,-1,0,-1,i,i)
            this_hnum=hnum.ProjectionU("n_{}".format(i),0,-1,0,-1,i,i)
            
        stop=False
        for iu in range(1,n2+1):
            if this_hden.GetBinContent(iu)==0: stop=True; break;
            if this_hnum.GetBinContent(iu)==0: stop=True; break;
        if stop: 
            if i==0: continue;
            break; 
        this_hden.Scale(1/this_hden.Integral())
        this_hnum.Scale(1/this_hnum.Integral())
        for iu in range(1,n2+1):
            val_num=this_hnum.GetBinContent(iu)
            err_num=this_hnum.GetBinError(iu)
            val_den=this_hden.GetBinContent(iu)
            err_den=this_hden.GetBinError(iu)
            if val_num==0 or val_den==0:
                val=0
                err=0
            else:
                val=val_num/val_den
                if correlation==0:
                    err=abs(val_num/val_den*math.sqrt(pow(err_num/val_num,2)+pow(err_den/val_den,2)))
                else:
                    err=abs(val_num/val_den*(err_num/val_num-correlation*err_den/val_den))
            if maxerror:
                if err>maxerror: stop=True; break;
            hout_all.SetBinContent(i,iu,val)
            hout_all.SetBinError(i,iu,err)
        this_hnum.Delete()
        this_hden.Delete()
        if stop: break;
            
    n1=i-1
    hout=ROOT.TH2D()
    hout.SetBins(n1,axis1.GetXbins().GetArray(),n2,axis2.GetXbins().GetArray())
    for i in range(n1+1):
        for iu in range(n2+2):
            hout.SetBinContent(i,iu,hout_all.GetBinContent(i,iu))
            hout.SetBinError(i,iu,hout_all.GetBinError(i,iu))
    hout_all.Delete()
    hout.SetStats(0)
    return hout
                
if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(description="get CosThetaWeight histogram")
    parser.add_argument("--input",dest="Input",type=str,help="Input root file path")
    parser.add_argument("--output",dest="Output",type=str,help="Output root file name")
    parser.add_argument("--test",dest="TestMode",default=0,type=int,help="test mode (0 or 1)")
    parser.add_argument("--note",action="store_true",help="Save plots for note as PDF format")
    args=parser.parse_args()

    if args.note:
        ROOT.gInterpreter.ProcessLine(".L Plotter/AFBPlotter.cc")
        p=ROOT.AFBPlotter("")
        p.SetupPlots("fig/costhetaweight/central_private.dat")
        p.AddFile("private","/data8/Users/hsseo/GeneratorTools/Hist/DY4D_MG_DY_NLO0_aew.root")
        p.SetupEntries("amc private")
        p.entries[1].replace["20[0-9][0-9]"]=""
        p.entries[1].replace["/gen_"]="/"
        p.entries[1].replace["_noweight"]=""
        p.SavePlot("compare_central_private_m","histname:[em][em]2017/gen_costhetaCS_noweight(x) xmin:60 xmax:120 norm widthweight pdf")
        p.SavePlot("compare_central_private_y","histname:[em][em]2017/gen_costhetaCS_noweight(y) xmin:60 xmax:120 norm 1:TMleg pdf")
        p.SavePlot("compare_central_private_pt","histname:[em][em]2017/gen_costhetaCS_noweight(z) xmin:60 xmax:120 norm widthweight zmin:2 zmax:400 logx 2:widey pdf")
        p.SavePlot("compare_central_private_AFB_m","histname:[em][em]2017/gen_AFB_noweight(x) xmin:60 xmax:120 1:TLleg pdf")
        p.SavePlot("compare_central_private_AFB_y","histname:[em][em]2017/gen_AFB_noweight(y) xmin:60 xmax:120 1:TMleg pdf")
        p.SavePlot("compare_central_private_AFB_pt","histname:[em][em]2017/gen_AFB_noweight(z) xmin:60 xmax:120 zmin:2 zmax:400 logx pdf")

        p.SetupEntries("data *amc+bg")
        p.SavePlot("before_ee","histname:ee2017/AFB_nocosthetaweight(x) xmin:60 xmax:120 1:TLleg type:6 pdf")
        p.SavePlot("before_mm","histname:mm2017/AFB_nocosthetaweight(x) xmin:60 xmax:120 1:TLleg type:6 pdf")
        p.SavePlot("after_ee","histname:ee2017/AFB(x) xmin:60 xmax:120 1:TLleg type:6 pdf")
        p.SavePlot("after_mm","histname:mm2017/AFB(x) xmin:60 xmax:120 1:TLleg type:6 pdf")

        p.SetupEntries("private private")
        p.PrintEntries(1)
        p.entries[0].title="sin^{2}#theta_{W} = 0.22225"
        p.entries[1].title="sin^{2}#theta_{W} = 0.23155"
        p.entries[1].replace["(\([x-zu,]*\))$"]="_0$1";
        p.SavePlot("compare_default_pdg_AFB_m","histname:[em][em]/AFB_correct(x) xmin:60 xmax:120 1:TLleg pdf")
        p.SavePlot("compare_default_pdg_AFB_y","histname:[em][em]/AFB_correct(y) xmin:60 xmax:120 1:TMleg pdf")
        p.SavePlot("compare_default_pdg_AFB_pt","histname:[em][em]/AFB_correct(z) xmin:60 xmax:120 zmin:2 zmax:400 logx pdf")
        
        ROOT.gStyle.SetPaintTextFormat(".2f")
        fofficial=ROOT.TFile("/data6/Users/hsseo/SKFlatOutput/Run2Legacy_v4/AFBAnalyzer/2017/AFBAnalyzer_DYJets.root")
        hee_official=fofficial.Get("ee2017/gen_costhetaCS_noweight")
        hee_official_correct=fofficial.Get("ee2017/gen_costhetaCS_correct_noweight")

        fprivate=ROOT.TFile("/data8/Users/hsseo/GeneratorTools/Hist/DY4D_MG_DY_NLO0_aew.root")
        hee_private=fprivate.Get("ee/costhetaCS")
        hee_private0=fprivate.Get("ee/costhetaCS_0")
        hee_private_correct=fprivate.Get("ee/costhetaCS_correct")
        hee_private_correct0=fprivate.Get("ee/costhetaCS_correct_0")
                
        hmm_private=fprivate.Get("mm/costhetaCS")
        hmm_private0=fprivate.Get("mm/costhetaCS_0")
        hmm_private_correct=fprivate.Get("mm/costhetaCS_correct")
        hmm_private_correct0=fprivate.Get("mm/costhetaCS_correct_0")

        h_private_correct=fprivate.Get("ee/costhetaCS_correct")
        h_private_correct0=fprivate.Get("ee/costhetaCS_correct_0")
        h_private_correct.Add(hmm_private_correct)
        h_private_correct0.Add(hmm_private_correct0)
        
        c0=ROOT.TCanvas()
        c0.cd()
        hout0=GetCosThetaWeight2D(h_private_correct,h_private_correct0,correlation=1)
        hout0.SetMaximum(1.15)
        hout0.SetMinimum(0.85)
        hout0.GetYaxis().SetTitle("cos#theta")
        hout0.GetXaxis().SetTitle("m(ll)")
        hout0.GetXaxis().SetRangeUser(60,-1)
        hout0.Draw("text colz")
        ROOT.gPad.SetLogx()
        c0.SaveAs("fig/costhetaweight/costhetaweight2D_m.pdf")
        c0.Delete()

        c1=ROOT.TCanvas()
        c1.cd()
        hout1=GetCosThetaWeight2D(h_private_correct,h_private_correct0,axis="Y",correlation=1)
        hout1.GetYaxis().SetTitle("cos#theta")
        hout1.GetXaxis().SetTitle("y(ll)")
        hout1.Draw("text colz")
        c1.SaveAs("fig/costhetaweight/costhetaweight2D_y.pdf")
        c1.Delete()

        c2=ROOT.TCanvas()
        c2.cd()
        hout2=GetCosThetaWeight2D(h_private_correct,h_private_correct0,axis="Z",correlation=1)
        hout2.GetYaxis().SetTitle("cos#theta")
        hout2.GetXaxis().SetTitle("p_{T}(ll)")
        hout2.Draw("text colz")
        ROOT.gPad.SetLogx()
        c2.SaveAs("fig/costhetaweight/costhetaweight2D_pt.pdf")
        c2.Delete()

        hee_m=GetCosThetaWeight2D(hee_private_correct,hee_private_correct0,correlation=1)
        hmm_m=GetCosThetaWeight2D(hmm_private_correct,hmm_private_correct0,correlation=1)
        hratio=hee_m.Clone("ratio")
        hratio.Divide(hmm_m);

        c3=ROOT.TCanvas()
        c3.cd()
        hratio.SetMaximum(1.15)
        hratio.SetMinimum(0.85)
        hratio.GetYaxis().SetTitle("cos#theta")
        hratio.GetXaxis().SetTitle("m(ll)")
        hratio.GetXaxis().SetRangeUser(80,120)
        hratio.Draw("text colz")
        ROOT.gPad.SetLogx()
        c3.SaveAs("fig/costhetaweight/compare_ee_mm.pdf")
        c3.Delete()


        exit()
        
    if args.TestMode==1:
        ROOT.gStyle.SetPaintTextFormat(".2f")

        fofficial=ROOT.TFile("/data6/Users/hsseo/SKFlatOutput/Run2Legacy_v4/AFBAnalyzer/2017/AFBAnalyzer_DYJets.root")
        hofficial=fofficial.Get("ee2017/gen_costhetaCS_noweight")
        hofficial_correct=fofficial.Get("ee2017/gen_costhetaCS_correct_noweight")

        fprivate=ROOT.TFile("/data8/Users/hsseo/GeneratorTools/Hist/DY4D_MG_DY_NLO0_aew.root")
        hprivate=fprivate.Get("ee/costhetaCS")
        hprivate0=fprivate.Get("ee/costhetaCS_0")
        hprivate10=fprivate.Get("ee/costhetaCS_10")
        hprivate_correct=fprivate.Get("ee/costhetaCS_correct")
        hprivate_correct0=fprivate.Get("ee/costhetaCS_correct_0")
        hprivate_correct10=fprivate.Get("ee/costhetaCS_correct_10")
        
        htest=GetCosThetaWeight3D(hprivate_correct0,hprivate_correct10,correlation=1)
        
        c0=ROOT.TCanvas()
        c0.cd()
        hout0=GetCosThetaWeight2D(hofficial_correct,hprivate_correct)
        hout0.SetMaximum(1.2)
        hout0.SetMinimum(0.8)
        hout0.GetXaxis().SetRangeUser(120,-1)
        hout0.Draw("text e colz")
        
        c1=ROOT.TCanvas()
        c1.cd()
        hout1=GetCosThetaWeight2D(hprivate_correct0,hprivate_correct10,axis="X",correlation=1)
        hout1.SetMaximum(1.2)
        hout1.SetMinimum(0.8)
        hout1.GetXaxis().SetRangeUser(120,-1)
        hout1.Draw("text e colz")
        
        c2=ROOT.TCanvas()
        c2.cd()
        hout2=GetCosThetaWeight2D(hprivate_correct0,hprivate_correct10,axis="Y",correlation=1)
        hout2.SetMaximum(1.2)
        hout2.SetMinimum(0.8)
        hout2.Draw("text e colz")
        
        c3=ROOT.TCanvas()
        c3.cd()
        hout3=GetCosThetaWeight2D(hprivate_correct0,hprivate_correct10,axis="Z",correlation=1)
        hout3.SetMaximum(1.2)
        hout3.SetMinimum(0.8)
        hout3.Draw("text e colz")
        ROOT.gPad.SetLogx()
        
        #hout4=GetCosThetaWeight(hden3,hnum3,correlation=1)
        a=input()
        exit()
    elif args.TestMode==2:
        ROOT.gStyle.SetPaintTextFormat(".2f")

        fofficial=ROOT.TFile("/data6/Users/hsseo/SKFlatOutput/Run2Legacy_v4/AFBAnalyzer/2017/AFBAnalyzer_DYJets.root")
        hee_official=fofficial.Get("ee2017/gen_costhetaCS_noweight")
        hee_official_correct=fofficial.Get("ee2017/gen_costhetaCS_correct_noweight")

        fprivate=ROOT.TFile("/data8/Users/hsseo/GeneratorTools/Hist/DY4D_MG_DY_NLO0_aew.root")
        hee_private=fprivate.Get("ee/costhetaCS")
        hee_private0=fprivate.Get("ee/costhetaCS_0")
        hee_private_correct=fprivate.Get("ee/costhetaCS_correct")
        hee_private_correct0=fprivate.Get("ee/costhetaCS_correct_0")
        hmm_private=fprivate.Get("mm/costhetaCS")
        hmm_private0=fprivate.Get("mm/costhetaCS_0")
        hmm_private_correct=fprivate.Get("mm/costhetaCS_correct")
        hmm_private_correct0=fprivate.Get("mm/costhetaCS_correct_0")
        
        c0=ROOT.TCanvas()
        c0.cd()
        hout0=GetCosThetaWeight2D(hee_private_correct,hee_private_correct0,correlation=1)
        hout0.SetMaximum(1.15)
        hout0.SetMinimum(0.85)
        hout0.GetYaxis().SetTitle("cos#theta")
        hout0.GetXaxis().SetTitle("m(ll)")
        hout0.GetXaxis().SetRangeUser(60,-1)
        hout0.Draw("text colz")
        ROOT.gPad.SetLogx()

        c1=ROOT.TCanvas()
        c1.cd()
        hout1=GetCosThetaWeight2D(hee_private_correct,hee_private_correct0,axis="Y",correlation=1)
        hout1.GetYaxis().SetTitle("cos#theta")
        hout1.GetXaxis().SetTitle("y(ll)")
        hout1.Draw("text colz")

        c2=ROOT.TCanvas()
        c2.cd()
        hout2=GetCosThetaWeight2D(hee_private_correct,hee_private_correct0,axis="Z",correlation=1)
        hout2.GetYaxis().SetTitle("cos#theta")
        hout2.GetXaxis().SetTitle("p_{T}(ll)")
        hout2.Draw("text colz")
        ROOT.gPad.SetLogx()

        c3=ROOT.TCanvas()
        c3.cd()
        hout3=GetCosThetaWeight2D(hee_official_correct,hee_private_correct)
        hout3.SetMaximum(1.15)
        hout3.SetMinimum(0.85)
        hout3.GetYaxis().SetTitle("cos#theta")
        hout3.GetXaxis().SetTitle("m(ll)")
        hout3.GetXaxis().SetRangeUser(60,-1)
        hout3.Draw("text colz")
        ROOT.gPad.SetLogx()

        c5=ROOT.TCanvas()
        c5.cd()
        hout5=GetCosThetaWeight2D(hee_official_correct,hee_private_correct,axis="Z")
        hout5.GetYaxis().SetTitle("cos#theta")
        hout5.GetXaxis().SetTitle("p_{T}(ll)")
        hout5.Draw("text colz")
        ROOT.gPad.SetLogx()

        c6=ROOT.TCanvas()
        c6.cd()
        hout6=GetCosThetaWeight2D(hmm_private_correct,hee_private_correct)
        hout6.SetMaximum(1.15)
        hout6.SetMinimum(0.85)
        hout6.GetYaxis().SetTitle("cos#theta")
        hout6.GetXaxis().SetTitle("m(ll)")
        hout6.GetXaxis().SetRangeUser(60,-1)
        hout6.Draw("text colz")
        ROOT.gPad.SetLogx()

        a=input()
        exit()
        
    
    
    fin=ROOT.TFile(args.Input)
    hee_den=fin.Get("ee/costhetaCS_correct")
    hmm_den=fin.Get("mm/costhetaCS_correct")
    h_den=fin.Get("ee/costhetaCS_correct")
    h_den.Add(hmm_den)
    hists=[]
    suffixes=["_pdg","_up","_down","_aew125","_aew128","_aew129","_aew130","_aew133"]
    for i in range(len(suffixes)):
        hee_num=fin.Get("ee/costhetaCS_correct_"+str(i))
        hmm_num=fin.Get("mm/costhetaCS_correct_"+str(i))
        h_num=fin.Get("ee/costhetaCS_correct_"+str(i))
        h_num.Add(hmm_num)
        h_out=GetCosThetaWeight3D(h_den,h_num,correlation=1)
        h_out.SetName("DYJets"+suffixes[i])
        h_out.SetTitle("DYJets"+suffixes[i])
        hists+=[h_out]

    fout=ROOT.TFile(args.Output,"recreate")
    for hist in hists:
        hist.Write()

    exit()
    
