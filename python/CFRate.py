import os,sys
import math,ctypes
from array import array
import ROOT as rt
rt.gROOT.LoadMacro("./Plotter/ZpeakPlotter.cc")
from ROOT import ZpeakPlotter

def calc_eff(valp,valf,errp=None,errf=None):
    if valp+valf==0: return 0.,0.
    eff=valp/(valp+valf)
    if errp is None and errf is None:
        return eff
    err = 1/(valp+valf)**2*math.sqrt(errp*errp*valf*valf+errf*errf*valp*valp)
    return eff,err

def evaluate(args):
    plotter=ZpeakPlotter("data-"+args.bgkey.replace("+wjets","").replace("+","-")+" "+args.dykey)
    rt.Verbosity=0

    ptbins=[10,30,40,50,70,90,200]
    #ptbins=[10,200]
    etabins=[0.0,1.0,1.5,1.7,2.0,2.5]
    #etabins=[0,1.5,2.5]
    cfdata=rt.TH2D("cfdata","cfdata",len(etabins)-1,array('d',etabins),len(ptbins)-1,array('d',ptbins));
    cfmc=rt.TH2D("cfmc","cfmc",len(etabins)-1,array('d',etabins),len(ptbins)-1,array('d',ptbins));
    cfscale=rt.TH2D("cfscale","cfscale",len(etabins)-1,array('d',etabins),len(ptbins)-1,array('d',ptbins));
    for h in [cfdata,cfmc,cfscale]:
        for i in range(h.GetNcells()):
            if h.GetBinContent(i)<0:
                h.SetBinContent(i,0)

    rt.RooMsgService.instance().setGlobalKillBelow(rt.RooFit.WARNING)
    Import=getattr(rt.RooWorkspace,"import")

    if not os.path.exists("fig/CFRate"):
        os.makedirs("fig/CFRate")
    rt.gROOT.SetBatch(True)
    for ip in range(len(ptbins)-1):
        for ie in range(len(etabins)-1):
            binstring="Xmin:{} Xmax:{} Ymin:{} Ymax:{} absX project:z ".format(etabins[ie],etabins[ie+1],ptbins[ip],ptbins[ip+1])
            hdataos=plotter.GetHist(0,"ee"+args.eras+"/dimass","xmin:76 xmax:106 prject:z "+binstring)
            hdatass=plotter.GetHist(0,"ee"+args.eras+"/ss_dimass","xmin:76 xmax:106 prject:z "+binstring)
            hmcos=plotter.GetHist(1,"ee"+args.eras+"/dimass","xmin:76 xmax:106 prject:z "+binstring)
            hmcss=plotter.GetHist(1,"ee"+args.eras+"/ss_dimass","xmin:76 xmax:106 prject:z "+binstring)

            hmcsscf=plotter.GetHist(1,"ee"+args.eras+"/ss_dimass_cf","xmin:76 xmax:106 prject:z "+binstring)

            w=rt.RooWorkspace("w")
            x=w.factory("x[70,110]")
            x.setRange("fit_range",76,106)
            dataos=rt.RooDataHist("dataos","dataos",rt.RooArgList(x),rt.RooFit.Import(hdataos))
            datass=rt.RooDataHist("datass","datass",rt.RooArgList(x),rt.RooFit.Import(hdatass))
            mcos=rt.RooDataHist("mcos","mcos",rt.RooArgList(x),hmcos)
            mcospdf=rt.RooHistPdf("mcospdf","mcospdf",rt.RooArgSet(x),mcos,0)
            mcss=rt.RooDataHist("mcss","mcss",rt.RooArgList(x),hmcss)
            mcsspdf=rt.RooHistPdf("mcsspdf","mcsspdf",rt.RooArgSet(x),mcss,0)
            Import(w,mcospdf)
            Import(w,mcsspdf)
    
            bgos=w.factory("CMSShape::bgos(x,alphaos[50,40,80],betaos[0.1,0.01,0.25],gammaos[0.05,0.0001,0.2],peak[90])")
            modelos=w.factory("SUM::modelos(fsigos[0.9,0.1,1]*mcospdf,bgos)")

            bgss=w.factory("CMSShape::bgss(x,alphass[50,40,80],betass[0.1,0.01,0.25],gammass[0.05,0.0001,0.2],peak[90])")
            modelss=w.factory("SUM::modelss(fsigss[0.9,0.1,1]*mcsspdf,bgss)")

            cos=rt.TCanvas("cos")
            modelos.fitTo(dataos,rt.RooFit.Range("fit_range"))
            plotos=x.frame(rt.RooFit.Title("os "+binstring))
            dataos.plotOn(plotos)
            modelos.plotOn(plotos)
            modelos.plotOn(plotos,rt.RooFit.Components("bgos"),rt.RooFit.LineStyle(rt.kDashed))
            plotos.Draw()
            cos.SaveAs("fig/CFRate/pt{}to{}_eta{}to{}_OS".format(ptbins[ip],ptbins[ip+1],etabins[ie],etabins[ie+1]).replace(".","p")+".png")

            css=rt.TCanvas("css")
            modelss.fitTo(datass,rt.RooFit.Range("fit_range"))
            plotss=x.frame(rt.RooFit.Title("ss "+binstring)).Clone()
            datass.plotOn(plotss)
            modelss.plotOn(plotss)
            modelss.plotOn(plotss,rt.RooFit.Components("bgss"),rt.RooFit.LineStyle(rt.kDashed))
            plotss.Draw()
            css.SaveAs("fig/CFRate/pt{}to{}_eta{}to{}_SS".format(ptbins[ip],ptbins[ip+1],etabins[ie],etabins[ie+1]).replace(".","p")+".png")
            
            hdataoserr=ctypes.c_double(0.)
            hdataosval=hdataos.IntegralAndError(hdataos.GetXaxis().GetFirst(),hdataos.GetXaxis().GetLast(),hdataoserr)
            hdataoserr=hdataoserr.value
            hdatasserr=ctypes.c_double(0.)
            hdatassval=hdatass.IntegralAndError(hdatass.GetXaxis().GetFirst(),hdatass.GetXaxis().GetLast(),hdatasserr)
            hdatasserr=hdatasserr.value
            fsigosval=w.var("fsigos").getVal()
            fsigoserr=w.var("fsigos").getError()
            fsigssval=w.var("fsigss").getVal()
            fsigsserr=w.var("fsigss").getError()

            hmcoserr=ctypes.c_double(0.)
            hmcosval=hmcos.IntegralAndError(hmcos.GetXaxis().GetFirst(),hmcos.GetXaxis().GetLast(),hmcoserr)
            hmcoserr=hmcoserr.value
            hmcsserr=ctypes.c_double(0.)
            hmcssval=hmcss.IntegralAndError(hmcss.GetXaxis().GetFirst(),hmcss.GetXaxis().GetLast(),hmcsserr)
            hmcsserr=hmcsserr.value
            hmcsscferr=ctypes.c_double(0.)
            hmcsscfval=hmcsscf.IntegralAndError(hmcsscf.GetXaxis().GetFirst(),hmcsscf.GetXaxis().GetLast(),hmcsscferr)
            hmcsscferr=hmcsscferr.value

            this_cfdata,this_cfdata_err=calc_eff(hdatassval*fsigssval,hdataosval*fsigosval,((hdatassval*fsigsserr)**2+(hdatasserr*fsigssval)**2)**0.5,((hdataosval*fsigoserr)**2+(hdataoserr*fsigosval)**2)**0.5)
            this_cfdata*=hmcsscfval/hmcssval
            this_cfdata_err*=hmcsscfval/hmcssval
            
            this_cfmc,this_cfmc_err=calc_eff(hmcssval,hmcosval,hmcsserr,hmcoserr)
            this_cfmc*=hmcsscfval/hmcssval
            this_cfmc_err*=hmcsscfval/hmcssval
            
            if this_cfdata==this_cfdata:
                cfdata.SetBinContent(ie+1,ip+1,this_cfdata)
                cfdata.SetBinError(ie+1,ip+1,this_cfdata_err)
            if this_cfmc==this_cfmc: 
                cfmc.SetBinContent(ie+1,ip+1,this_cfmc)
                cfmc.SetBinError(ie+1,ip+1,this_cfmc_err)
            
            
            hdatamean=hdatass.Clone()
            hdatamean.GetXaxis().SetRangeUser(82,100)
            hmcmean=hmcss.Clone()
            hmcmean.GetXaxis().SetRangeUser(82,100)
            scale=hdatamean.GetMean()/hmcmean.GetMean()
            scaleerr=( (hdatamean.GetMeanError()/hdatamean.GetMean())**2 + (hmcmean.GetMeanError()/hmcmean.GetMean())**2 )**0.5*scale

            scaleerr=2*scale*scaleerr
            scale=scale**2

            print scale,scaleerr
            cfscale.SetBinContent(ie+1,ip+1,scale)
            cfscale.SetBinError(ie+1,ip+1,scaleerr)


            #r=raw_input()


    cfsf=cfdata.Clone("cfsf")
    cfsf.SetTitle("cfsf")
    cfsf.Divide(cfmc)
    cfsf_this=cfsf.Clone("cfsf_this")
    cfsf_this.SetTitle("cfsf_this")
    cfscale_this=cfscale.Clone("cfscale_this")
    cffilepath=os.getenv("SKFlat_WD")+"/data/"+os.getenv("SKFlatV")+"/"+args.era+"/SMP/CFRate.root"
    cfsf_old=None
    cfscale_old=None
    if os.path.exists(cffilepath):
        fold=rt.TFile(cffilepath)
        cfsf_old=fold.Get("cfsf")
        if cfsf_old:
            cfsf_old.SetNameTitle("cfsf_old","cfsf_old")
            cfsf_old.SetDirectory(0)
            h=cfsf_old.Clone()
            for i in range(h.GetNcells()):
                h.SetBinError(i,0)
            cfsf.Multiply(h)
            cfmc.Divide(h)
        cfscale_old=fold.Get("cfscale")
        if cfscale_old:
            cfscale_old.SetNameTitle("cfscale_old","cfscale_old")
            cfscale_old.SetDirectory(0)
            h=cfscale_old.Clone()
            for i in range(h.GetNcells()):
                h.SetBinError(i,0)
            cfscale.Multiply(h)            
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
    cfscale.SetOption("colz text e")
    cfscale.Write("cfscale")
    cfscale_this.SetOption("colz text e")
    cfscale_this.Write("cfscale_this")
    if cfsf_old:
        cfsf_old.SetOption("colz text e")
        cfsf_old.Write("cfsf_old")
    if cfscale_old:
        cfscale_old.SetOption("colz text e")
        cfscale_old.Write("cfscale_old")

    cgcf_eta=rt.TCanvas("cgcf_eta")
    for i in range(len(ptbins)-1):
        data_px=cfdata.ProjectionX("cfdata_px_{}_{}".format(ptbins[i],ptbins[i+1]),i+1,i+1)
        data_px.SetLineColor(i+1)
        data_px.SetLineWidth(2)
        data_px.SetStats(0)
        if i==0: 
            data_px.GetYaxis().SetRangeUser(0,0.06)
            data_px.Draw("hist e")
        else:
            data_px.Draw("same hist e")
        mc_px=cfmc.ProjectionX("cfmc_px_{}_{}".format(ptbins[i],ptbins[i+1]),i+1,i+1)
        mc_px.SetLineColor(i+1)
        mc_px.SetLineStyle(2)
        mc_px.SetLineWidth(2)
        mc_px.SetStats(0)
        mc_px.Draw("same hist e")
    cgcf_eta.Write("cf_eta")

    cgcf_pt=rt.TCanvas("cgcf_pt")
    for i in range(len(etabins)-1):
        data_py=cfdata.ProjectionY("cfdata_py_{}_{}".format(etabins[i],etabins[i+1]),i+1,i+1)
        data_py.SetLineColor(i+1)
        data_py.SetLineWidth(2)
        data_py.SetStats(0)
        if i==0: 
            data_py.GetYaxis().SetRangeUser(0,0.06)
            data_py.Draw("hist e")
        else:
            data_py.Draw("same hist e")
        mc_py=cfmc.ProjectionY("cfmc_py_{}_{}".format(etabins[i],etabins[i+1]),i+1,i+1)
        mc_py.SetLineColor(i+1)
        mc_py.SetLineStyle(2)
        mc_py.SetLineWidth(2)
        mc_py.SetStats(0)
        mc_py.Draw("same hist e")
    cgcf_pt.Write("cf_pt")

    cgcfsf_eta=rt.TCanvas("cgcfsf_eta")
    for i in range(len(ptbins)-1):
        sf_px=cfsf.ProjectionX("cfsf_px_{}_{}".format(ptbins[i],ptbins[i+1]),i+1,i+1)
        sf_px.SetLineColor(i+1)
        sf_px.SetLineWidth(2)
        sf_px.SetStats(0)
        if i==0: 
            sf_px.Draw("hist e")
            sf_px.GetYaxis().SetRangeUser(0.5,2)
        else:
            sf_px.Draw("same hist e")
    cgcfsf_eta.Write("cfsf_eta")

    cgcfsf_pt=rt.TCanvas("cgcfsf_pt")
    for i in range(len(etabins)-1):
        sf_py=cfsf.ProjectionY("cfsf_py_{}_{}".format(etabins[i],etabins[i+1]),i+1,i+1)
        sf_py.SetLineColor(i+1)
        sf_py.SetLineWidth(2)
        sf_py.SetStats(0)
        if i==0: 
            sf_py.Draw("hist e")
            sf_py.GetYaxis().SetRangeUser(0.5,2)
        else:
            sf_py.Draw("same hist e")
    cgcfsf_pt.Write("cfsf_pt")

    f.Close()
    #r=raw_input()
    for i in range(cfsf_this.GetNcells()):
        val=cfsf_this.GetBinContent(i)
        err=cfsf_this.GetBinError(i)
        if val==0: continue
        if val<1-err or val>1+err: return False

    return True


def validate(args):
    plotter=ZpeakPlotter("data ^{}+{}".format(args.dykey,args.bgkey))
    #rt.Verbosity=0

    ptbins=[10,30,40,50,70,100,200]
    etabins=[0.0,1.0,1.4,1.7,2.0,2.5]

    binstring=""
    hdataos=plotter.GetHist(0,"ee"+args.eras+"/dimass","xmin:52 xmax:150 project:z "+binstring)
    hmcos=plotter.GetTH1(plotter.GetHist(1,"ee"+args.eras+"/dimass","xmin:52 xmax:150 project:z "+binstring))
    normsf=hdataos.Integral()/hmcos.Integral()
    plotter.entries[1].weight=normsf
    plotter.DrawPlot("ee"+args.eras+"/dimass","xmin:54 xmax:150 project:z rebin:2 2:widewidey "+binstring)
    plotter.DrawPlot("ee"+args.eras+"/ss_dimass","xmin:54 xmax:150 project:z rebin:4 2:widewidey "+binstring)
    plotter.DrawPlot("ee"+args.eras+"/ss_dimass","suffix:_noCFSF:sim xmin:54 xmax:150 project:z rebin:4 2:widewidey "+binstring)
    plotter.entries[1].weight=1.

    print normsf
    plotter2=ZpeakPlotter("{}*data-{}-{}".format(1/normsf,args.dykey,args.bgkey.replace("+","-")))
    plotter2.DrawPlot("ee"+args.eras+"/ss_dimass","xmin:54 xmax:150 project:z rebin:4 2:widewidey "+binstring)
    r=raw_input()

    for ip in range(len(ptbins)-1):
        binstring="Ymin:{} Ymax:{}".format(ptbins[ip],ptbins[ip+1])
        print binstring
        hdataos=plotter.GetHist(0,"ee"+args.eras+"/dimass","xmin:52 xmax:150 project:z "+binstring)
        hmcos=plotter.GetTH1(plotter.GetHist(1,"ee"+args.eras+"/dimass","xmin:52 xmax:150 project:z "+binstring))
        normsf=hdataos.Integral()/hmcos.Integral()
        plotter.entries[1].weight=normsf
        plotter.DrawPlot("ee"+args.eras+"/dimass","xmin:54 xmax:150 project:z rebin:2 2:widewidey "+binstring)
        plotter.DrawPlot("ee"+args.eras+"/ss_dimass","xmin:54 xmax:150 project:z rebin:4 2:widewidey "+binstring)
        plotter.entries[1].weight=1.
        r=raw_input()

    for ie in range(len(etabins)-1):
        binstring="Xmin:{} Xmax:{} absX".format(etabins[ie],etabins[ie+1])
        print binstring
        hdataos=plotter.GetHist(0,"ee"+args.eras+"/dimass","xmin:52 xmax:150 project:z "+binstring)
        hmcos=plotter.GetTH1(plotter.GetHist(1,"ee"+args.eras+"/dimass","xmin:52 xmax:150 project:z "+binstring))
        normsf=hdataos.Integral()/hmcos.Integral()
        plotter.entries[1].weight=normsf
        plotter.DrawPlot("ee"+args.eras+"/dimass","xmin:54 xmax:150 project:z rebin:2 2:widewidey "+binstring)
        plotter.DrawPlot("ee"+args.eras+"/ss_dimass","xmin:54 xmax:150 project:z rebin:4 2:widewidey "+binstring)
        plotter.entries[1].weight=1.
        r=raw_input()

    # for ip in range(len(ptbins)-1):
    #     for ie in range(len(etabins)-1):
    #         binstring="Xmin:{} Xmax:{} Ymin:{} Ymax:{} absX".format(etabins[ie],etabins[ie+1],ptbins[ip],ptbins[ip+1])
    #         print binstring
    #         hdataos=plotter.GetHist(0,"ee"+args.eras+"/dimass","xmin:52 xmax:150 project:z "+binstring)
    #         hmcos=plotter.GetTH1(plotter.GetHist(1,"ee"+args.eras+"/dimass","xmin:52 xmax:150 project:z "+binstring))
    #         normsf=hdataos.Integral()/hmcos.Integral()
    #         plotter.entries[1].weight=normsf
    #         plotter.DrawPlot("ee"+args.eras+"/dimass","xmin:54 xmax:150 project:z rebin:2 2:widewidey "+binstring)
    #         plotter.DrawPlot("ee"+args.eras+"/ss_dimass","xmin:54 xmax:150 project:z rebin:4 2:widewidey "+binstring)
    #         plotter.entries[1].weight=1.

    #         r=raw_input()

def iterate(args):
    for i in range(args.maxiter):
        print "[CFRate] Iteration {}".format(i)
        if i==0:
            cmd=""
            for sample in args.samples_dy+args.samples_data+args.samples_bg:
                cmd+="SKFlat.py -a ZpeakAnalyzer --skim SkimTree_Dilepton -i {} -e {} -n {} --nmax {} & ".format(sample,args.era,60 if "DY" in sample else 30,args.nmax)
            cmd+="wait;"
            os.system(cmd)
        else:
            cmd=""
            for sample in args.samples_dy:
                cmd+="SKFlat.py -a ZpeakAnalyzer --skim SkimTree_Dilepton -i {} -e {} -n {} & ".format(sample,args.era,args.nmax)
            cmd+="wait;"
            os.system(cmd)

        stop=evaluate(args)

        cffilepath=os.getenv("SKFlat_WD")+"/data/"+os.getenv("SKFlatV")+"/"+args.era+"/SMP/CFRate.root"
        if os.path.exists(cffilepath):
            index=0
            while os.path.exists(cffilepath.replace(".root","_old{}.root".format(index))):
                index+=1
            newpath=cffilepath.replace(".root","_old{}.root".format(index))
            os.system("mv {} {}".format(cffilepath,newpath))
        os.system("mv CFRate.root "+cffilepath)

        if stop:
            print "[CFRate] Stop at iteration {}".format(i)
            break    
        
if __name__=="__main__":
    import argparse

    parser=argparse.ArgumentParser()
    parser.add_argument("action",help="evaluate(eval), validate(val), iterate(iter)")
    parser.add_argument("era",help="2016preVFP(2016a), 2016postVFP(2016b), 2017, 2018")
    parser.add_argument("--maxiter",type=int,default=5,help="maximum iteration")
    parser.add_argument("--nmax",type=int,default=300,help="condor concurrency limit")
    args=parser.parse_args()

    args.samples_data=["DoubleEG"]
    args.samples_dy=["DYJetsToEE_MiNNLO"]
    args.samples_bg=["DYJetsToTauTau_MiNNLO","WW_pythia","WZ_pythia","ZZ_pythia","WJets_MG","TTLL_powheg","SingleTop_tW_top_NoFullyHad","SingleTop_tW_antitop_NoFullyHad","SingleTop_tch_top_Incl","SingleTop_tch_antitop_Incl","SingleTop_sch_Lep","GGToLL"]
    args.dykey="mi"
    args.bgkey="tau_mi+vv+wjets+tt+st+aa"
    if args.era in ["2016preVFP","2016a"]:
        args.era="2016preVFP"
        args.eras="2016a"
    elif args.era in ["2016postVFP","2016b"]:
        args.era="2016postVFP"
        args.eras="2016b"
    elif args.era=="2017":
        args.era="2017"
        args.eras="2017"
    elif args.era=="2018":
        args.era="2018"
        args.eras="2018"
        args.samples_data=["EGamma"]
    else:
        print "unavailable era",args.era
        exit(1)



    if args.action in ["evaludate","eval"]:
        evaluate(args)
    elif args.action in ["validate","val"]:
        validate(args)
    elif args.action in ["iterate","iter"]:
        iterate(args)
    else:
        print "unavailable action",args.action
        exit(1)
