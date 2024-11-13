import os,ctypes,array,copy
import ROOT
#ROOT.gROOT.SetBatch(True)
ROOT.gROOT.LoadMacro("./Plotter/FakePlotter2.cc")
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


def MakeTF(args,region="",option=""):
    l0tf=None
    l1tf=None
    if args.channel=="mm":
        n=4
        events=[[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                events[i][j]=args.plotter.GetHist(0,"mm{}/wp{}{}/{}ss_dimass".format(args.erashort,i,j,region),"rebin:{20,70} "+option).GetBinContent(1)
            print events[i]

        if "s2m0" in option:
            num=args.plotter.GetHist(0,"mm{}/wp12/{}ss_l0etapt".format(args.erashort,region),"noproject "+option)
            hnum=Rebin2D(num,args.tfxbins,args.tfybins)
            den=args.plotter.GetHist(0,"mm{}/wp22/{}ss_l0etapt".format(args.erashort,region),"noproject "+option)
            hden=Rebin2D(den,args.tfxbins,args.tfybins)
            hnum.Divide(hden)
            hnum.Scale((events[1][2]+events[0][2])/events[1][2])
            l0tf=hnum.Clone()

            num=args.plotter.GetHist(0,"mm{}/wp21/{}ss_l1etapt".format(args.erashort,region),"noproject "+option)
            hnum=Rebin2D(num,args.tfxbins,args.tfybins)
            den=args.plotter.GetHist(0,"mm{}/wp22/{}ss_l1etapt".format(args.erashort,region),"noproject "+option)
            hden=Rebin2D(den,args.tfxbins,args.tfybins)
            hnum.Divide(hden)
            hnum.Scale((events[2][1]+events[2][0])/events[2][1])
            l1tf=hnum.Clone()

        else:
            num=args.plotter.GetHist(0,"mm{}/wp[01]2/{}ss_l0etapt".format(args.erashort,region),"noproject "+option)
            hnum=Rebin2D(num,args.tfxbins,args.tfybins)
            den=args.plotter.GetHist(0,"mm{}/wp22/{}ss_l0etapt".format(args.erashort,region),"noproject "+option)
            hden=Rebin2D(den,args.tfxbins,args.tfybins)
            hnum.Divide(hden)
            l0tf=hnum.Clone()

            num=args.plotter.GetHist(0,"mm{}/wp2[01]/{}ss_l1etapt".format(args.erashort,region),"noproject "+option)
            hnum=Rebin2D(num,args.tfxbins,args.tfybins)
            den=args.plotter.GetHist(0,"mm{}/wp22/{}ss_l1etapt".format(args.erashort,region),"noproject "+option)
            hden=Rebin2D(den,args.tfxbins,args.tfybins)
            hnum.Divide(hden)
            l1tf=hnum.Clone()

    if args.channel=="ee":
        n=6
        events=[[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                events[i][j]=args.plotter.GetHist(0,"ee{}/wp{}{}/{}ss_dimass".format(args.erashort,i,j,region),"rebin:{20,70} "+option).GetBinContent(1)
            print events[i]

            
        if "s2m0" in option:
            num=args.plotter.GetHist(0,"ee{}/wp14/{}ss_l0etapt".format(args.erashort,region),"noproject "+option)
            hnum=Rebin2D(num,args.tfxbins,args.tfybins)
            den=args.plotter.GetHist(0,"ee{}/wp54/{}ss_l0etapt".format(args.erashort,region),"noproject "+option)
            hden=Rebin2D(den,args.tfxbins,args.tfybins)
            hnum.Divide(hden)
            hnum.Scale((events[0][4]+events[1][4])/(events[1][4]))
            l0tf=hnum.Clone()

            num=args.plotter.GetHist(0,"ee{}/wp41/{}ss_l1etapt".format(args.erashort,region),"noproject "+option)
            hnum=Rebin2D(num,args.tfxbins,args.tfybins)
            den=args.plotter.GetHist(0,"ee{}/wp45/{}ss_l1etapt".format(args.erashort,region),"noproject "+option)
            hden=Rebin2D(den,args.tfxbins,args.tfybins)
            hnum.Divide(hden)
            hnum.Scale((events[4][0]+events[4][1])/(events[4][1]))
            l1tf=hnum.Clone()

        else:
            num=args.plotter.GetHist(0,"ee{}/wp[01]4/{}ss_l0etapt".format(args.erashort,region),"noproject "+option)
            hnum=Rebin2D(num,args.tfxbins,args.tfybins)
            den=args.plotter.GetHist(0,"ee{}/wp54/{}ss_l0etapt".format(args.erashort,region),"noproject "+option)
            hden=Rebin2D(den,args.tfxbins,args.tfybins)
            hnum.Divide(hden)
            l0tf=hnum.Clone()

            num=args.plotter.GetHist(0,"ee{}/wp4[01]/{}ss_l1etapt".format(args.erashort,region),"noproject "+option)
            hnum=Rebin2D(num,args.tfxbins,args.tfybins)
            den=args.plotter.GetHist(0,"ee{}/wp45/{}ss_l1etapt".format(args.erashort,region),"noproject "+option)
            hden=Rebin2D(den,args.tfxbins,args.tfybins)
            hnum.Divide(hden)
            l1tf=hnum.Clone()

    return l0tf,l1tf

def MakeSet(hists,name):
    print "[MakeSet]", name
    rt=[]
    if not hists[0][0]: 
        print "no hist... pass"
        return rt
    hist=hists[0][0].Clone(name)
    hist.SetTitle(name)
    for i in range(1,len(hists)):
        for ibin in range(hist.GetNcells()):
            val=hist.GetBinContent(ibin)
            err=hist.GetBinError(ibin)
            maxdiff=0
            for j in range(len(hists[i])):
                diff=abs(hists[i][j].GetBinContent(ibin)-val)
                if diff>maxdiff:
                    maxdiff=diff
            hist.SetBinError(ibin,(err**2+maxdiff**2)**0.5)
    rt+=[hist]
    for i in range(len(hists)):
        for j in range(len(hists[i])):
            hists[i][j].SetName(name+"_s{}m{}".format(i,j))
            hists[i][j].SetTitle(name+"_s{}m{}".format(i,j))
            rt+=[hists[i][j].Clone()]
    return rt

if __name__=="__main__":
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument("-e","--era",default="2016a,2016b,2017,2018",type=str)
    #parser.add_argument("-e","--era",default="2018",type=str)
    #parser.add_argument("-c","--channel",default="ee,mm,ej,mj")
    parser.add_argument("-c","--channel",default="ee,mm")
    #parser.add_argument("-r","--region",default=",noZ/,0bjet/noZ/,nbjet/noZ/,cpt/noZ/,cpt/0bjet/noZ/,cpt/nbjet/noZ/,mpt/noZ/,mpt/0bjet/noZ/,mpt/nbjet/noZ/")
    #parser.add_argument("-r","--region",default="metcut/,0bjet/metcut/,nbjet/metcut/,noZ/metcut/,0bjet/noZ/metcut/,nbjet/noZ/metcut/")
    parser.add_argument("-r","--region",default="metcut/,nbjet/metcut/,noZ/metcut/,nbjet/noZ/metcut/")
    args_raw=parser.parse_args()

    args_raw.era=args_raw.era.split(",")
    args_raw.channel=args_raw.channel.split(",")
    args_raw.region=args_raw.region.split(",")

    save=[]
    plotter=ROOT.FakePlotter2("data-mi-tau_mi-wjets-vv-tt-st-aa-mg10")
    for era in args_raw.era:
        for channel in args_raw.channel:
            for region in args_raw.region:
                args=copy.deepcopy(args_raw)
                args.plotter=plotter
                args.era=GetEra(era)
                args.erashort=GetEraShort(era)
                if args.era not in ["2016preVFP","2016postVFP","2017","2018"]:
                    print "Unknown era {}".format(args.era)
                    continue
                args.region=region
                args.channel=channel
                args.tfxbins=[0.,1.5,2.5]
                if args.channel=="ee":
                    args.tfybins=[10,12.5,15,17.5,20,22.5,25,27.5,30,35,40,50,100,200]
                
                elif args.channel=="mm":
                    args.tfybins=[10,12.5,15,17.5,20,22.5,25,27.5,30,35,40,50,60,100,200]

                else:
                    print "Unknown channel {}".format(args.channel)
                    continue

                hl0tf_s0m0,hl1tf_s0m0=MakeTF(args,region)
                hl0tf_s1m0,hl1tf_s1m0=MakeTF(args,region,"scale:1.1:data")
                hl0tf_s1m1,hl1tf_s1m1=MakeTF(args,region,"scale:0.9:data")
                # hl0tf_s2m0,hl1tf_s2m0=MakeTF(args,region,"s2m0")
                # save+=MakeSet([[hl0tf_s0m0],[hl0tf_s1m0,hl0tf_s1m1],[hl0tf_s2m0]],"{}{}_{}l0tf".format(args.channel,args.era,args.region.replace("/","_")))
                # save+=MakeSet([[hl1tf_s0m0],[hl1tf_s1m0,hl1tf_s1m1],[hl1tf_s2m0]],"{}{}_{}l1tf".format(args.channel,args.era,args.region.replace("/","_")))
                save+=MakeSet([[hl0tf_s0m0],[hl0tf_s1m0,hl0tf_s1m1]],"{}{}_{}l0tf".format(args.channel,args.era,args.region.replace("/","_")))
                save+=MakeSet([[hl1tf_s0m0],[hl1tf_s1m0,hl1tf_s1m1]],"{}{}_{}l1tf".format(args.channel,args.era,args.region.replace("/","_")))
                        
    f=ROOT.TFile("FakeTF.root","recreate")
    for h in save:
        # for i in range(h.GetNcells()):
        #     val=h.GetBinContent(i)
        #     if val<0:
        #         h.SetBinContent(i,0.)
        #         h.SetBinError(i,0.)
        h.Write()
    exit(0)
