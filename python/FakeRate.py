import os,ctypes,array,copy
import ROOT
#ROOT.gROOT.SetBatch(True)
ROOT.gROOT.LoadMacro("./Plotter/SKFlatPlotter.cc")
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

def DrawTransferFactor(args,htf):
    c=ROOT.TCanvas()
    c.hists=[]
    c.graphs=[]
    c.leg=ROOT.TLegend(0.2,0.89,0.4,0.75)
    c.leg.SetBorderSize(0)
    for i in range(1,htf.GetNbinsX()+1):
        c.hists+=[htf.ProjectionY("hl0tf_x{}".format(i),i,i)]
        c.hists[-1].SetLineColor(i)
        c.hists[-1].SetStats(0)
        etamin=htf.GetXaxis().GetBinLowEdge(i)
        etamax=htf.GetXaxis().GetBinLowEdge(i+1)
        c.hists[-1].SetTitle("{}#leq|#eta|<{}".format(etamin,etamax))
        c.leg.AddEntry(c.hists[-1])
            
        if i==1: 
            c.hists[-1].Draw("HIST e")
            if args.channel=="ee":
                c.hists[-1].GetYaxis().SetRangeUser(0,1.0)
                c.hists[-1].GetXaxis().SetTitle("electron p_{T} [GeV]");
            elif args.channel=="mm":
                c.hists[-1].GetYaxis().SetRangeUser(0,1.0)
                c.hists[-1].GetXaxis().SetTitle("muon p_{T} [GeV]");
            c.hists[-1].GetYaxis().SetTitle("transfer factor");
            c.hists[-1].GetXaxis().SetMoreLogLabels()
        else:
            c.hists[-1].Draw("HIST e same")
        c.graphs+=[ROOT.TGraph()]
        for x in range(int(c.hists[-1].GetXaxis().GetBinLowEdge(1)),int(c.hists[-1].GetXaxis().GetBinLowEdge(c.hists[-1].GetXaxis().GetLast()+1))):
            c.graphs[-1].SetPoint(c.graphs[-1].GetN(),x,GetTransferFactor(htf,htf.GetXaxis().GetBinCenter(i),x))            
        c.graphs[-1].Draw("same")
        c.graphs[-1].SetLineColor(i)
    c.SetLogx()
    c.hists[0].SetTitle(args.era)
    c.leg.Draw()
    c.Update()
    return c

def GetTransferFactor(htf,eta,pt):
    etamin=htf.GetXaxis().GetBinLowEdge(1)
    etamax=htf.GetXaxis().GetBinUpEdge(htf.GetNbinsX())
    ptmin=htf.GetYaxis().GetBinLowEdge(1)
    ptmax=htf.GetYaxis().GetBinUpEdge(htf.GetNbinsY())
    if eta<etamin: eta=etamin+1e-6
    if eta>etamax: eta=etamax-1e-6
    if pt<ptmin: pt=ptmin+1e-6
    if pt>ptmax: pt=ptmax-1e-6
    eta=htf.GetXaxis().GetBinCenter(htf.GetXaxis().FindBin(eta))
    return htf.Interpolate(eta,pt)

def Apply(htf,hl):
    out=hl.Clone("hlnum_est")
    for i in range(out.GetNcells()):
        val=out.GetBinContent(i)
        err=out.GetBinError(i)
        binx=ctypes.c_int()
        biny=ctypes.c_int()
        binz=ctypes.c_int()
        out.GetBinXYZ(i,binx,biny,binz)
        binx=binx.value
        biny=biny.value
        leta=out.GetXaxis().GetBinCenter(binx)
        lpt=out.GetYaxis().GetBinCenter(biny)
        tf=GetTransferFactor(htf,leta,lpt)
        out.SetBinContent(i,val*tf)
        out.SetBinError(i,err*tf)
    out.SetTitle(out.GetName())
    return out

def evaluate(args):
    print "[evaluate]",args.channel,args.era,args.region
    hl0num=Rebin2D(args.hl0num,args.tfxbins,args.tfybins)
    hl0den=Rebin2D(args.hl0den,args.tfxbins,args.tfybins)
    hl0num.Divide(hl0den)
    hl0num.SetNameTitle(args.hl0tf.GetName(),args.hl0tf.GetTitle())
    args.hl0tf=hl0num

    hl1num=Rebin2D(args.hl1num,args.tfxbins,args.tfybins)
    hl1den=Rebin2D(args.hl1den,args.tfxbins,args.tfybins)
    hl1num.Divide(hl1den)
    hl1num.SetNameTitle(args.hl1tf.GetName(),args.hl1tf.GetTitle())
    args.hl1tf=hl1num
    
    return

def validation(args):
    print "[validation]",args.channel,args.era,args.region
    
    if args.channel=="ee":
        hl0num_val=Rebin2D(args.plotter.GetHist(0,"ee{}/{}/ss_l0etapt".format(args.erashort,args.region),"noproject"),args.tfxbins,args.tfybins)
        hl0den_val=Rebin2D(args.plotter.GetHist(0,"Ee{}/{}/ss_l0etapt".format(args.erashort,args.region),"noproject"),args.tfxbins,args.tfybins)
        hl1num_val=Rebin2D(args.plotter.GetHist(0,"ee{}/{}/ss_l1etapt".format(args.erashort,args.region),"noproject"),args.tfxbins,args.tfybins)
        hl1den_val=Rebin2D(args.plotter.GetHist(0,"eE{}/{}/ss_l1etapt".format(args.erashort,args.region),"noproject"),args.tfxbins,args.tfybins)
    if args.channel=="mm":
        hl0num_val=Rebin2D(args.plotter.GetHist(0,"mm{}/{}/ss_l0etapt".format(args.erashort,args.region),"noproject"),args.tfxbins,args.tfybins)
        hl0den_val=Rebin2D(args.plotter.GetHist(0,"Mm{}/{}/ss_l0etapt".format(args.erashort,args.region),"noproject"),args.tfxbins,args.tfybins)
        hl1num_val=Rebin2D(args.plotter.GetHist(0,"mm{}/{}/ss_l1etapt".format(args.erashort,args.region),"noproject"),args.tfxbins,args.tfybins)
        hl1den_val=Rebin2D(args.plotter.GetHist(0,"mM{}/{}/ss_l1etapt".format(args.erashort,args.region),"noproject"),args.tfxbins,args.tfybins)

    hl0num_val.Divide(hl0den_val)
    cl0=DrawTransferFactor(args,hl0num_val)
    cl0_=DrawTransferFactor(args,args.hl0tf)

    hl1num_val.Divide(hl1den_val)
    cl1=DrawTransferFactor(args,hl1num_val)
    cl1_=DrawTransferFactor(args,args.hl1tf)

    raw_input()
    return
    

def iterate(args):
    print "[iterate]",args.channel,args.era
    i=0
    while True:
        print "[iterate] iter",i
        hl0residue=args.hl0num.Clone("residue")
        hl0residue=Rebin2D(hl0residue,args.tfxbins,args.tfybins)
        hl0residue.Divide(Rebin2D(Apply(args.hl0tf,args.hl0den),args.tfxbins,args.tfybins))
        args.hl0tf.Multiply(hl0residue)

        hl1residue=args.hl1num.Clone("residue")
        hl1residue=Rebin2D(hl1residue,args.tfxbins,args.tfybins)
        hl1residue.Divide(Rebin2D(Apply(args.hl1tf,args.hl1den),args.tfxbins,args.tfybins))
        args.hl1tf.Multiply(hl1residue)

        c=DrawTransferFactor(args)
        #raw_input()
        if i==5: 
            raw_input()
            break;
        for j in range(args.hl0tf.GetNcells()):
            args.hl0tf.SetBinError(j,0.)
        for j in range(args.hl1tf.GetNcells()):
            args.hl1tf.SetBinError(j,0.)
        i+=1

if __name__=="__main__":
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument("-e","--era",default="2016a,2016b,2017,2018",type=str)
    parser.add_argument("-c","--channel",default="ee,mm")
    parser.add_argument("-r","--region",default="cpt/noZ/,noZ/,mpt/noZ/")
    args_raw=parser.parse_args()

    args_raw.era=args_raw.era.split(",")
    args_raw.channel=args_raw.channel.split(",")
    args_raw.region=args_raw.region.split(",")

    save=[]
    for era in args_raw.era:
        for channel in args_raw.channel:
            for region in args_raw.region:
                args=copy.deepcopy(args_raw)
                args.plotter=ROOT.SKFlatPlotter("FakeAnalyzer","data-mi-tau_mi-wjets-vv-tttw-aa")
                args.era=GetEra(era)
                args.erashort=GetEraShort(era)
                if args.era not in ["2016preVFP","2016postVFP","2017","2018"]:
                    print "Unknown era {}".format(args.era)
                    continue
                args.region=region
                if args.region not in ["","noZ/","cpt/","cpt/noZ/","mpt/","mpt/noZ/"]:
                    print "Unknown region {}".format(args.region)
                    continue
                args.channel=channel
                if args.channel=="ee":
                    args.tfxbins=[0.,1.5,2.5]
                    #args.tfxbins=[0,1.0,1.5,2.0,2.5]
                    #args.tfybins=[10.,20.,40.,60.,100.,200.]
                    args.tfybins=[10,12.5,15,17.5,20,22.5,25,27.5,30,35,40,50,100,1000]
                    #args.tfybins=[15.,25.,40.,60.,100.,200.,400.]
                    
                    args.hl0num=args.plotter.GetHist(0,"eE{}/{}/ss_l0etapt".format(args.erashort,region),"noproject")
                    if not args.hl0num: continue
                    args.hl0den=args.plotter.GetHist(0,"EE{}/{}/ss_l0etapt".format(args.erashort,region),"noproject")
                    args.hl0tf=args.hl0num.Clone("{}{}_{}l0tf".format(args.channel,args.era,args.region.replace("/","_")))
                    args.hl0tf=Rebin2D(args.hl0tf,args.tfxbins,args.tfybins)

                    args.hl1num=args.plotter.GetHist(0,"Ee{}/{}/ss_l1etapt".format(args.erashort,region),"noproject")
                    args.hl1den=args.plotter.GetHist(0,"EE{}/{}/ss_l1etapt".format(args.erashort,region),"noproject")
                    args.hl1tf=args.hl1num.Clone("{}{}_{}l1tf".format(args.channel,args.era,args.region.replace("/","_")))
                    args.hl1tf=Rebin2D(args.hl1tf,args.tfxbins,args.tfybins)

                elif args.channel=="mm":
                    #args.tfxbins=[0,1.0,1.5,2.0,2.5]
                    args.tfxbins=[0,1.5,2.5]
                    #args.tfybins=[10.,20.,40.,60.,100.,200.]
                    args.tfybins=[10,12.5,15,17.5,20,22.5,25,27.5,30,35,40,50,100,1000]
                    
                    args.hl0num=args.plotter.GetHist(0,"mM{}/{}/ss_l0etapt".format(args.erashort,args.region),"noproject")
                    if not args.hl0num: continue
                    args.hl0den=args.plotter.GetHist(0,"MM{}/{}/ss_l0etapt".format(args.erashort,args.region),"noproject")
                    args.hl0tf=args.hl0num.Clone("{}{}_{}l0tf".format(args.channel,args.era,args.region.replace("/","_")))
                    args.hl0tf=Rebin2D(args.hl0tf,args.tfxbins,args.tfybins)
                    
                    args.hl1num=args.plotter.GetHist(0,"Mm{}/{}/ss_l1etapt".format(args.erashort,args.region),"noproject")
                    args.hl1den=args.plotter.GetHist(0,"MM{}/{}/ss_l1etapt".format(args.erashort,args.region),"noproject")
                    args.hl1tf=args.hl1num.Clone("{}{}_{}l1tf".format(args.channel,args.era,args.region.replace("/","_")))
                    args.hl1tf=Rebin2D(args.hl1tf,args.tfxbins,args.tfybins)
                else:
                    print "Unknown channel {}".format(args.channel)
                    continue
                    
                args.hl0tf.SetTitle(args.hl0tf.GetName())
                args.hl1tf.SetTitle(args.hl1tf.GetName())
                args.hl0tf.Reset()
                args.hl1tf.Reset()
                for i in range(args.hl0tf.GetNcells()): args.hl0tf.SetBinContent(i,0.5)
                for i in range(args.hl1tf.GetNcells()): args.hl1tf.SetBinContent(i,0.5)
                
                #args.plotter.DrawPlot("eE2017/l0etapt_noZ","project:y logx widthweight")
                #args.plotter.DrawPlot("mM2017/l0etapt_noZ","project:y logx widthweight")
                #args.plotter2=ROOT.AFBPlotter("data ^mi+tau_mi+wjets+vv+tttw+aa","FakeAnalyzer")
                #args.plotter2.DrawPlot("mM2017/l0etapt_noZ","project:y logx widthweight")
                #c=ROOT.TCanvas()
                #args.hl.Draw("text colz")
                #raw_input()
                
                #iterate(args)
                evaluate(args)
                #validation(args)
                #c=DrawTransferFactor(args,0)
                #c2=DrawTransferFactor(args,1)
                #raw_input()
                
                save+=[args.hl0tf.Clone()]
                save+=[args.hl1tf.Clone()]
    
    f=ROOT.TFile("FakeTF.root","recreate")
    for h in save:
        h.Write()
    exit(0)
