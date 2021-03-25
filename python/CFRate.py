import os,sys
import math
from array import array
import ROOT as rt

def evaluate(args):
    rt.gROOT.LoadMacro("./Plotter/ZpeakPlotter.cc")
    from ROOT import ZpeakPlotter
    plotter=ZpeakPlotter("data mg")
    rt.Verbosity=0

    ptbins=[10,30,40,50,70,100,200]
    etabins=[0.0,1.0,1.4,1.7,2.0,2.5]
    cfdata=rt.TH2D("cfdata","cfdata",len(etabins)-1,array('d',etabins),len(ptbins)-1,array('d',ptbins));
    cfmc=rt.TH2D("cfmc","cfmc",len(etabins)-1,array('d',etabins),len(ptbins)-1,array('d',ptbins));

    rt.RooMsgService.instance().setGlobalKillBelow(rt.RooFit.WARNING)
    Import=getattr(rt.RooWorkspace,"import")

    if not os.path.exists("fig/CFRate"):
        os.makedirs("fig/CFRate")

    for ip in range(len(ptbins)-1):
        for ie in range(len(etabins)-1):
            binstring="Xmin:{} Xmax:{} Ymin:{} Ymax:{} absy".format(ptbins[ip],ptbins[ip+1],etabins[ie],etabins[ie+1])
            hdataos=plotter.GetHist(0,"ee"+args.eras+"/dimass_even","xmin:52 xmax:150 prject:z "+binstring)
            hdatass=plotter.GetHist(0,"ee"+args.eras+"/ss_dimass_even","xmin:52 xmax:150 prject:z "+binstring)
            hmcos=plotter.GetHist(1,"ee"+args.eras+"/dimass_even","xmin:52 xmax:150 prject:z "+binstring)
            hmcss=plotter.GetHist(1,"ee"+args.eras+"/ss_dimass_even","xmin:52 xmax:150 prject:z "+binstring)

            hmcsscf=plotter.GetHist(1,"ee"+args.eras+"/ss_dimass_cf_even","xmin:52 xmax:150 prject:z "+binstring)

            w=rt.RooWorkspace("w")
            x=w.factory("x[52,150]")
            dataos=rt.RooDataHist("dataos","dataos",rt.RooArgList(x),rt.RooFit.Import(hdataos))
            datass=rt.RooDataHist("datass","datass",rt.RooArgList(x),rt.RooFit.Import(hdatass))
            mcos=rt.RooDataHist("mcos","mcos",rt.RooArgList(x),hmcos)
            mcospdf=rt.RooHistPdf("mcospdf","mcospdf",rt.RooArgSet(x),mcos,0)
            mcss=rt.RooDataHist("mcss","mcss",rt.RooArgList(x),hmcss)
            mcsspdf=rt.RooHistPdf("mcsspdf","mcsspdf",rt.RooArgSet(x),mcss,0)
            Import(w,mcospdf)
            Import(w,mcsspdf)
    
            bgos=w.factory("CMSShape::bgos(x,alphaos[70,30,100],betaos[0.05,0.0,0.1],gammaos[0,0.2],peak[90])")
            modelos=w.factory("SUM::modelos(fsigos[0.8,0.1,1]*mcospdf,bgos)")

            bgss=w.factory("CMSShape::bgss(x,alphass[70,30,100],betass[0.05,0.0,0.1],gammass[0,0.2],peak[90])")
            modelss=w.factory("SUM::modelss(fsigss[0.8,0.1,1]*mcsspdf,bgss)")

            cos=rt.TCanvas("cos")
            modelos.fitTo(dataos)
            plotos=x.frame(rt.RooFit.Title("os "+binstring))
            dataos.plotOn(plotos)
            modelos.plotOn(plotos)
            modelos.plotOn(plotos,rt.RooFit.Components("bgos"),rt.RooFit.LineStyle(rt.kDashed))
            plotos.Draw()
            cos.SaveAs("fig/CFRate/pt{}to{}_eta{}to{}_OS".format(ptbins[ip],ptbins[ip+1],etabins[ie],etabins[ie+1]).replace(".","p")+".png")

            css=rt.TCanvas("css")
            modelss.fitTo(datass)
            plotss=x.frame(rt.RooFit.Title("ss "+binstring)).Clone()
            datass.plotOn(plotss)
            modelss.plotOn(plotss)
            modelss.plotOn(plotss,rt.RooFit.Components("bgss"),rt.RooFit.LineStyle(rt.kDashed))
            plotss.Draw()
            css.SaveAs("fig/CFRate/pt{}to{}_eta{}to{}_SS".format(ptbins[ip],ptbins[ip+1],etabins[ie],etabins[ie+1]).replace(".","p")+".png")
            
            hdataoserr=rt.Double(0.)
            hdataosval=hdataos.IntegralAndError(hdataos.GetXaxis().GetFirst(),hdataos.GetXaxis().GetLast(),hdataoserr)
            hdatasserr=rt.Double(0.)
            hdatassval=hdatass.IntegralAndError(hdatass.GetXaxis().GetFirst(),hdatass.GetXaxis().GetLast(),hdatasserr)
            fsigosval=w.var("fsigos").getVal()
            fsigoserr=w.var("fsigos").getError()
            fsigssval=w.var("fsigss").getVal()
            fsigsserr=w.var("fsigss").getError()

            hmcoserr=rt.Double(0.)
            hmcosval=hmcos.IntegralAndError(hmcos.GetXaxis().GetFirst(),hmcos.GetXaxis().GetLast(),hmcoserr)
            hmcsserr=rt.Double(0.)
            hmcssval=hmcss.IntegralAndError(hmcss.GetXaxis().GetFirst(),hmcss.GetXaxis().GetLast(),hmcsserr)
            hmcsscferr=rt.Double(0.)
            hmcsscfval=hmcsscf.IntegralAndError(hmcsscf.GetXaxis().GetFirst(),hmcsscf.GetXaxis().GetLast(),hmcsscferr)

            this_cfdata=hdatassval*fsigssval/hdataosval/fsigosval*hmcsscfval/hmcssval
            this_cfdata_err=this_cfdata*math.sqrt((hdatasserr/hdatassval)**2+(fsigsserr/fsigssval)**2+(hdataoserr/hdataosval)**2+(fsigoserr/fsigosval)**2+(hmcsscferr/hmcsscfval)**2+(hmcsserr/hmcssval)**2)
            this_cfmc=hmcsscfval/hmcosval
            this_cfmc_err=this_cfmc*math.sqrt((hmcsscferr/hmcsscfval)**2+(hmcoserr/hmcosval)**2)
            
            if this_cfdata==this_cfdata:
                cfdata.SetBinContent(ie+1,ip+1,this_cfdata)
                cfdata.SetBinError(ie+1,ip+1,this_cfdata_err)
            if this_cfmc==this_cfmc: 
                cfmc.SetBinContent(ie+1,ip+1,this_cfmc)
                cfmc.SetBinError(ie+1,ip+1,this_cfmc_err)

            #r=raw_input()


    cfsf=cfdata.Clone("cfsf")
    cfsf.SetTitle("cfsf")
    cfsf.Divide(cfmc)
    cfsf_this=cfsf.Clone("cfsf_this")
    cfsf_this.SetTitle("cfsf_this")
    cffilepath=os.getenv("SKFlat_WD")+"/data/"+os.getenv("SKFlatV")+"/"+args.era+"/SMP/CFRate.root"
    cfsf_old=None
    if os.path.exists(cffilepath):
        fold=rt.TFile(cffilepath)
        cfsf_old=fold.Get("cfsf")
        if cfsf_old:
            cfsf_old.SetDirectory(0)
            h=cfsf_old.Clone()
            for i in range(h.GetNcells()):
                h.SetBinError(i,0)
            cfsf.Multiply(h)
        fold.Close()

    f=rt.TFile("CFRate.root","recreate")

    cfdata.SetOption("colz text e")
    cfdata.Write("cfdata")
    cfmc.SetOption("colz text e")
    cfmc.Write("cfmc")
    cfsf.SetOption("colz text e")
    cfsf.Write("cfsf")
    cfsf_this.SetOption("colz text e")
    cfsf_this.Write("cfsf_this")
    if cfsf_old:
        cfsf_old.SetOption("colz text e")
        cfsf_old.Write("cfsf_old")


    cgcf_eta=rt.TCanvas("cgcf_eta")
    for i in range(len(ptbins)-1):
        data_px=cfdata.ProjectionX("cfdata_px_{}_{}".format(ptbins[i],ptbins[i+1]),i+1,i+1)
        data_px.SetLineColor(i+1)
        data_px.SetStats(0)
        if i==0: 
            data_px.Draw()
        else:
            data_px.Draw("same")
        mc_px=cfmc.ProjectionX("cfmc_px_{}_{}".format(ptbins[i],ptbins[i+1]),i+1,i+1)
        mc_px.SetLineColor(i+1)
        mc_px.SetLineStyle(2)
        mc_px.SetStats(0)
        mc_px.Draw("same")
    cgcf_eta.Write("cf_eta")

    cgcf_pt=rt.TCanvas("cgcf_pt")
    for i in range(len(etabins)-1):
        data_py=cfdata.ProjectionY("cfdata_py_{}_{}".format(etabins[i],etabins[i+1]),i+1,i+1)
        data_py.SetLineColor(i+1)
        data_py.SetStats(0)
        if i==0: 
            data_py.Draw()
        else:
            data_py.Draw("same")
        mc_py=cfmc.ProjectionY("cfmc_py_{}_{}".format(etabins[i],etabins[i+1]),i+1,i+1)
        mc_py.SetLineColor(i+1)
        mc_py.SetLineStyle(2)
        mc_py.SetStats(0)
        mc_py.Draw("same")
    cgcf_pt.Write("cf_pt")

    cgcfsf_eta=rt.TCanvas("cgcfsf_eta")
    for i in range(len(ptbins)-1):
        sf_px=cfsf.ProjectionX("cfsf_px_{}_{}".format(ptbins[i],ptbins[i+1]),i+1,i+1)
        sf_px.SetLineColor(i+1)
        sf_px.SetStats(0)
        if i==0: 
            sf_px.Draw()
        else:
            sf_px.Draw("same")
    cgcfsf_eta.Write("cfsf_eta")

    cgcfsf_pt=rt.TCanvas("cgcfsf_pt")
    for i in range(len(etabins)-1):
        sf_py=cfsf.ProjectionY("cfsf_py_{}_{}".format(etabins[i],etabins[i+1]),i+1,i+1)
        sf_py.SetLineColor(i+1)
        sf_py.SetStats(0)
        if i==0: 
            sf_py.Draw()
        else:
            sf_py.Draw("same")
    cgcfsf_pt.Write("cfsf_pt")

    f.Close()

    #r=raw_input()

def validate(args):
    rt.gROOT.LoadMacro("./Plotter/ZpeakPlotter.cc")
    from ROOT import ZpeakPlotter
    plotter=ZpeakPlotter("data ^mg+tau_mg+vv+wjets+tttw")
    rt.Verbosity=0

    ptbins=[10,30,40,50,70,100,200]
    etabins=[ 0.4*x for x in range(7)]

    binstring=""
    hdataos=plotter.GetHist(0,"ee"+args.eras+"/dimass_odd","xmin:52 xmax:150 prject:z "+binstring)
    hmcos=plotter.GetTH1(plotter.GetHist(1,"ee"+args.eras+"/dimass_odd","xmin:52 xmax:150 prject:z "+binstring))
    normsf=hdataos.Integral()/hmcos.Integral()
    plotter.entries[1].weight=normsf
    plotter.DrawPlot("ee"+args.eras+"/dimass_odd","xmin:54 xmax:150 prject:z rebin:2 2:widewidey "+binstring)
    plotter.DrawPlot("ee"+args.eras+"/ss_dimass_odd","xmin:54 xmax:150 prject:z rebin:4 2:widewidey "+binstring)
    plotter.entries[1].weight=1.
    r=raw_input()

    for ip in range(len(ptbins)-1):
        binstring="Xmin:{} Xmax:{}".format(ptbins[ip],ptbins[ip+1])
        print binstring
        hdataos=plotter.GetHist(0,"ee"+args.eras+"/dimass_odd","xmin:52 xmax:150 prject:z "+binstring)
        hmcos=plotter.GetTH1(plotter.GetHist(1,"ee"+args.eras+"/dimass_odd","xmin:52 xmax:150 prject:z "+binstring))
        normsf=hdataos.Integral()/hmcos.Integral()
        plotter.entries[1].weight=normsf
        plotter.DrawPlot("ee"+args.eras+"/dimass_odd","xmin:54 xmax:150 prject:z rebin:2 2:widewidey "+binstring)
        plotter.DrawPlot("ee"+args.eras+"/ss_dimass_odd","xmin:54 xmax:150 prject:z rebin:4 2:widewidey "+binstring)
        plotter.entries[1].weight=1.
        r=raw_input()

    for ie in range(len(etabins)-1):
        binstring="Ymin:{} Ymax:{} absy".format(etabins[ie],etabins[ie+1])
        print binstring
        hdataos=plotter.GetHist(0,"ee"+args.eras+"/dimass_odd","xmin:52 xmax:150 prject:z "+binstring)
        hmcos=plotter.GetTH1(plotter.GetHist(1,"ee"+args.eras+"/dimass_odd","xmin:52 xmax:150 prject:z "+binstring))
        normsf=hdataos.Integral()/hmcos.Integral()
        plotter.entries[1].weight=normsf
        plotter.DrawPlot("ee"+args.eras+"/dimass_odd","xmin:54 xmax:150 prject:z rebin:2 2:widewidey "+binstring)
        plotter.DrawPlot("ee"+args.eras+"/ss_dimass_odd","xmin:54 xmax:150 prject:z rebin:4 2:widewidey "+binstring)
        plotter.entries[1].weight=1.
        r=raw_input()

    for ip in range(len(ptbins)-1):
        for ie in range(len(etabins)-1):
            binstring="Xmin:{} Xmax:{} Ymin:{} Ymax:{} absy".format(ptbins[ip],ptbins[ip+1],etabins[ie],etabins[ie+1])
            print binstring
            hdataos=plotter.GetHist(0,"ee"+args.eras+"/dimass_odd","xmin:52 xmax:150 prject:z "+binstring)
            hmcos=plotter.GetTH1(plotter.GetHist(1,"ee"+args.eras+"/dimass_odd","xmin:52 xmax:150 prject:z "+binstring))
            normsf=hdataos.Integral()/hmcos.Integral()
            plotter.entries[1].weight=normsf
            plotter.DrawPlot("ee"+args.eras+"/dimass_odd","xmin:54 xmax:150 prject:z rebin:2 2:widewidey "+binstring)
            plotter.DrawPlot("ee"+args.eras+"/ss_dimass_odd","xmin:54 xmax:150 prject:z rebin:4 2:widewidey "+binstring)
            plotter.entries[1].weight=1.

            r=raw_input()

if __name__=="__main__":
    import argparse

    parser=argparse.ArgumentParser()
    parser.add_argument("action",help="evaluate(eval),validate(val)")
    parser.add_argument("era",help="2016preVFP(2016a),2016postVFP(2016b),2017,2018")
    args=parser.parse_args()

    if args.era in ["2016preVFP","2016a"]:
        args.era="2016preVFP"
        args.eras="2016a"
    elif args.era in ["2016postVFP","2016b"]:
        args.era="2016preVFP"
        args.eras="2016a"
    elif args.era=="2017":
        args.era="2017"
        args.eras="2017"
    elif args.era=="2018":
        args.era="2018"
        args.eras="2018"
    else:
        print "unavailable era",args.era
        exit(1)

    if args.action in ["evaludate","eval"]:
        evaluate(args)
    elif args.action in ["validate","val"]:
        validate(args)
    else:
        print "unavailable action",args.action
        exit(1)
