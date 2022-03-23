import os,ctypes,array,copy
import ROOT
#ROOT.gROOT.SetBatch(True)
ROOT.gROOT.LoadMacro("./Plotter/AFBPlotter.cc")
ROOT.TH1.AddDirectory(0)

def GetHist(filename,histname):
    f=ROOT.TFile(filename)
    h=f.Get(histname)
    if not h: print "no {} in {}".format(histname,filename)
    return h

def GetEra(era):
    if era=="2016a": 
        return "2016preVFP"
    elif era=="2016b": 
        return "2016postVFP"
    else: 
        return era

def GetEraShort(era):
    if era=="2016preVFP": 
        return "2016a"
    elif era=="2016postVFP": 
        return "2016b"
    else: 
        return era
    

def Rebin2D(hist,xbins,ybins):
    histnew=ROOT.TH2D(hist.GetName(),hist.GetTitle(),len(xbins)-1,array.array("d",xbins),len(ybins)-1,array.array("d",ybins))
    histnew.SetStats(0)
    for i in range(hist.GetNcells()):
        val=hist.GetBinContent(i)
        err=hist.GetBinError(i)
        binx=ctypes.c_int()
        biny=ctypes.c_int()
        binz=ctypes.c_int()
        hist.GetBinXYZ(i,binx,biny,binz)
        binx=binx.value
        biny=biny.value
        x=hist.GetXaxis().GetBinCenter(binx)
        y=hist.GetYaxis().GetBinCenter(biny)
        
        j=histnew.FindBin(x,y)
        histnew.SetBinContent(j,histnew.GetBinContent(j)+val)
        histnew.SetBinError(j,(histnew.GetBinError(j)**2+err**2)**0.5)
    return histnew

def DrawFakeRate(args):
    c=ROOT.TCanvas()
    c.hists=[]
    c.graphs=[]
    c.leg=ROOT.TLegend(0.2,0.89,0.4,0.75)
    c.leg.SetBorderSize(0)
    for i in range(1,args.hfake.GetNbinsX()+1):
        c.hists+=[args.hfake.ProjectionY("hfake_x{}".format(i),i,i)]
        c.hists[-1].SetLineColor(i)
        c.hists[-1].SetStats(0)
        etamin=args.hfake.GetXaxis().GetBinLowEdge(i)
        etamax=args.hfake.GetXaxis().GetBinLowEdge(i+1)
        c.hists[-1].SetTitle("{}#leq|#eta|<{}".format(etamin,etamax))
        c.leg.AddEntry(c.hists[-1])
            
        if i==1: 
            c.hists[-1].Draw("HIST e")
            if args.channel=="ee":
                c.hists[-1].GetYaxis().SetRangeUser(0,0.3)
                c.hists[-1].GetXaxis().SetTitle("electron p_{T} [GeV]");
            elif args.channel=="mm":
                c.hists[-1].GetYaxis().SetRangeUser(0,5.0)
                c.hists[-1].GetXaxis().SetTitle("muon p_{T} [GeV]");
            c.hists[-1].GetYaxis().SetTitle("transfer factor");
            c.hists[-1].GetXaxis().SetMoreLogLabels()
        else:
            c.hists[-1].Draw("HIST e same")
        c.graphs+=[ROOT.TGraph()]
        for x in range(int(c.hists[-1].GetXaxis().GetBinLowEdge(1)),int(c.hists[-1].GetXaxis().GetBinLowEdge(c.hists[-1].GetXaxis().GetLast()+1))):
            c.graphs[-1].SetPoint(c.graphs[-1].GetN(),x,GetFakeRate(args,args.hfake.GetXaxis().GetBinCenter(i),x))            
        c.graphs[-1].Draw("same")
        c.graphs[-1].SetLineColor(i)
    c.SetLogx()
    c.hists[0].SetTitle(args.era)
    c.leg.Draw()
    c.Update()
    return c

def GetFakeRate(args,eta,pt):
    etamin=args.hfake.GetXaxis().GetBinLowEdge(1)
    etamax=args.hfake.GetXaxis().GetBinUpEdge(args.hfake.GetNbinsX())
    ptmin=args.hfake.GetYaxis().GetBinLowEdge(1)
    ptmax=args.hfake.GetYaxis().GetBinUpEdge(args.hfake.GetNbinsY())
    if eta<etamin: eta=etamin+1e-6
    if eta>etamax: eta=etamax-1e-6
    if pt<ptmin: pt=ptmin+1e-6
    if pt>ptmax: pt=ptmax-1e-6
    eta=args.hfake.GetXaxis().GetBinCenter(args.hfake.GetXaxis().FindBin(eta))
    return args.hfake.Interpolate(eta,pt)

def LL2l(args):
    hl=args.hl.Clone()
    hl.Reset()
    for ih in range(2):
        hal=args.hals[ih]
        for i in range(hal.GetNcells()):
            val=hal.GetBinContent(i)
            err=hal.GetBinError(i)
            binx=ctypes.c_int()
            biny=ctypes.c_int()
            binz=ctypes.c_int()
            hal.GetBinXYZ(i,binx,biny,binz)
            binx=binx.value
            biny=biny.value
            leta=hal.GetXaxis().GetBinCenter(binx)
            lpt=hal.GetYaxis().GetBinCenter(biny)
            if args.channel=="ee":
                if lpt<15: continue
            elif args.channel=="mm":
                if lpt<10: continue
            j=hl.FindBin(leta,lpt)
            hl.SetBinContent(j,hl.GetBinContent(j)+val*GetFakeRate(args,leta,lpt))
            hl.SetBinError(j, ( hl.GetBinError(j)**2 + (err*GetFakeRate(args,leta,lpt))**2 )**0.5 )
    hl.SetName("hl_est")
    return hl

def iterate(args):
    print "[iterate]",args.channel,args.era
    i=0
    while True:
        print "[iterate] iter",i
        hl_estimate=LL2l(args)
        args.hfake.Multiply(args.hl)
        args.hfake.Divide(hl_estimate)
        c=DrawFakeRate(args)
        #raw_input()
        if i==5: 
            raw_input()
            break;
        for j in range(args.hfake.GetNcells()):
            args.hfake.SetBinError(j,0.)
        i+=1

if __name__=="__main__":
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument("-e","--era",default="2016a,2016b,2017,2018",type=str)
    parser.add_argument("-c","--channel",default="ee,mm")
    args_raw=parser.parse_args()

    args_raw.era=args_raw.era.split(",")
    args_raw.channel=args_raw.channel.split(",")

    save=[]
    for era in args_raw.era:
        for channel in args_raw.channel:
            args=copy.deepcopy(args_raw)
            args.plotter=ROOT.AFBPlotter("data-mi-tau_mi-wjets-vv-tttw-aa","FakeAnalyzer")
            args.era=GetEra(era)
            args.erashort=GetEraShort(era)
            if args.era not in ["2016preVFP","2016postVFP","2017","2018"]:
                print "Unknown era {}".format(args.era)
                continue
            args.channel=channel
            args.hals=[]
            if args.channel=="ee":
                for i in range(10):
                    hal=args.plotter.GetHist(0,"EE{}/ss_al{}etapt_noZ".format(args.erashort,i),"noproject")
                    if hal:
                        args.hals+=[hal.Clone()]
                    else:
                        break
                args.hl=args.plotter.GetHist(0,"eE{}/ss_l0etapt_noZ".format(args.erashort),"noproject")
                args.hl=Rebin2D(args.hl,[0,1.5,2.5],[15,25,40,60,100,200,400])
            elif args.channel=="mm":
                for i in range(10):
                    hal=args.plotter.GetHist(0,"MM{}/ss_al{}etapt_noZ".format(args.erashort,i),"noproject")
                    if hal:
                        args.hals+=[hal.Clone()]
                    else:
                        break
                args.hl=args.plotter.GetHist(0,"mM{}/ss_l0etapt_noZ".format(args.erashort),"noproject")
                args.hl=Rebin2D(args.hl,[0,1.0,1.5,2.0,2.5],[10,20,40,60,100,200])
            else:
                print "Unknown channel {}".format(args.channel)
                continue

            args.hfake=args.hl.Clone()
            args.hfake.SetNameTitle("{}{}".format(args.channel,args.era),"{}{}".format(args.channel,args.era))
            args.hfake.Reset()
            for i in range(args.hfake.GetNcells()): args.hfake.SetBinContent(i,0.5)

            #args.plotter.DrawPlot("eE2017/l0etapt_noZ","project:y logx widthweight")
            #args.plotter.DrawPlot("mM2017/l0etapt_noZ","project:y logx widthweight")
            #args.plotter2=ROOT.AFBPlotter("data ^mi+tau_mi+wjets+vv+tttw+aa","FakeAnalyzer")
            #args.plotter2.DrawPlot("mM2017/l0etapt_noZ","project:y logx widthweight")
            #c=ROOT.TCanvas()
            #args.hl.Draw("text colz")
            #raw_input()
            iterate(args)
            save+=[args.hfake.Clone()]
    
    f=ROOT.TFile("FakeRate.root","recreate")
    for h in save:
        h.Write()
    exit(0)
