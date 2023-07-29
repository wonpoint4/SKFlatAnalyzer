## eff residual systematic set
import os,sys,re
import ROOT
ROOT.gROOT.LoadMacro("./Plotter/EfficiencyPlotter.cc")
ROOT.TH1.AddDirectory(0)
def addResidual(infilename):
    print infilename
    outfilename=infilename.replace(".root","_residual.root")

    if "2016preVFP" in infilename:
        era="2016a"
    elif "2016postVFP" in infilename:
        era="2016b"
    elif "2017" in infilename:
        era="2017"
    elif "2018" in infilename:
        era="2018"
    
    if "Electron" in infilename:
        channel="el"
    elif "Muon" in infilename:
        channel="mu"

    os.system("cp {} {}".format(infilename,outfilename));
    f=ROOT.TFile(outfilename,"update")
    hsim=f.Get("sim")
    hsf_origin=f.Get("sf")

    plotter=ROOT.EfficiencyPlotter("data mi+tau_mi+vv+wjets+tt+st+qcdss+aa")
    hsf=plotter.GetHist(0,channel+era+"/m80to100/letapt","noproject")
    hsf_denom=plotter.GetHist(1,channel+era+"/m80to100/letapt","noproject")
    scale=hsf.Integral()/hsf_denom.Integral()
    hsf_denom.Scale(scale)
    hsf.Divide(hsf_denom)
    hsf.Multiply(hsf_origin)
    hdata=hsim.Clone("data")
    hdata.Multiply(hsf)
    f.cd()
    for h in [hdata, hsim, hsf]:
        for i in range(h.GetNcells()):
            h.SetBinError(i,0)
    iset=max([int(re.match("sf_s([0-9]+)",key.GetName()).group(1)) for key in f.GetListOfKeys() if re.match("sf_s([0-9]+)",key.GetName())])+1
    hdata.SetNameTitle("data_s{}m0".format(iset),"residual")
    hsim.SetNameTitle("sim_s{}m0".format(iset),"residual")
    hsf.SetNameTitle("sf_s{}m0".format(iset),"residual")
    hdata.Write()
    hsim.Write()
    hsf.Write()
    f.Close()

if __name__=="__main__":
    if sys.argv[1]=="all":
        for era in ["2016preVFP","2016postVFP","2017","2018"]:
            erashort=era.replace("preVFP","a").replace("postVFP","b")
            addResidual("data/Run2UltraLegacy_v3/{}/ID/Electron/{}_MediumID_v15.root".format(era,erashort))
            addResidual("data/Run2UltraLegacy_v3/{}/ID/Muon/MediumID_LooseTrkIso_Inclusive_v14.root".format(era))
    else:
        addResidual(sys.argv[1])
    
