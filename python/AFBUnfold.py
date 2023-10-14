#!/usr/bin/env python

import array,os,sys,re
import ROOT,ctypes

ROOT.TH1.AddDirectory(0)
ROOT.TH1.SetDefaultSumw2(1)

ROOT.gROOT.ProcessLine('#include"AFBPlotter.cc"')
SystematicSuffixes=dict(ROOT.AFBPlotter("").GetSystematicSuffixes("totalsys"))
DEBUG=0

class Config(object):
    def __init__(self,histname,data=None,sim=None):
        self.histname=histname[:histname.find("_dressed")+8]
        self.suffix=histname[histname.find("_dressed")+8:]
        self.option=""
        if self.suffix!="":
            self.option=SystematicSuffixes[self.suffix]

        words=histname.split("/")
        self.channel=words[0][:2]
        self.era=words[0][2:]
        self.region=words[1]
        self.bins=[
            array.array("d",[52.,3000.]),
            array.array("d",[-2.4,2.4]),
            array.array("d",[0.,650.]),
            array.array("d",[-1.,0.,1.]),
        ]
        if "dimass" in histname:
            self.merge_axis=None
            self.merge_bins=array.array("d",[0,1])
            self.primary_axis=0
            self.primary_bins_skflat=array.array("d",ROOT.AFBAnalyzer.afb_mbin)
            if self.region=="0bjet":
                self.primary_bins=array.array("d",ROOT.AFBAnalyzer.afb_mbin)
                self.bins[self.primary_axis]=self.primary_bins[::2]
            else:
                self.primary_bins=array.array("d",[52,60,65,70,77,90,106,120,140,175,200,240,280,340,400,600,3000])
                self.bins[self.primary_axis]=array.array("d",[52,65,77,106,140,200,280,400,3000])
            self.matrixnames=[re.sub(r"genfid_dimass[CSRecoil]*_dressed","response_afbm",self.histname)]
        elif "dirap" in histname:
            self.merge_axis=0
            self.merge_bins=array.array("d",[52.,77.,106.,280.,3000.])
            self.bins[self.merge_axis]=self.merge_bins[:]
            self.primary_axis=1
            self.primary_bins_skflat=array.array("d",ROOT.AFBAnalyzer.afb_ybin)
            if self.region=="0bjet":
                self.primary_bins=array.array("d",[-2.4,-2.0,-1.6,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.6,2.0,2.4])
                self.bins[self.primary_axis]=array.array("d",[0.0,0.4,0.8,1.2,1.6,2.0,2.4])
            else:
                self.primary_bins=array.array("d",[-2.4,-1.6,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.6,2.4])
                self.bins[self.primary_axis]=array.array("d",[0.0,0.4,0.8,1.2,1.6,2.4])
            self.matrixnames=[re.sub(r"genfid_dirap[CSRecoil]*_dressed","response_afby",self.histname)+"_m"+str(i) for i in range(len(self.merge_bins)-1)]
        elif "dipt" in histname:
            self.merge_axis=0
            self.merge_bins=array.array("d",[52.,77.,106.,280.,3000.])
            self.bins[self.merge_axis]=self.merge_bins[:]
            self.primary_axis=2
            self.primary_bins_skflat=array.array("d",ROOT.AFBAnalyzer.afb_ptbin)
            if self.region=="0bjet":
                self.primary_bins=array.array("d",ROOT.AFBAnalyzer.afb_ptbin)
                self.bins[self.primary_axis]=self.primary_bins[::2]
            else:
                self.primary_bins=array.array("d",[0,2,10,20,28,40,50,60,70,80,90,100,120,140,190,650])
                self.bins[self.primary_axis]=array.array("d",[0,2,20,40,60,80,100,140,650])
            self.matrixnames=[re.sub(r"genfid_dipt[CSRecoil]*_dressed","response_afbpt",self.histname)+"_m"+str(i) for i in range(len(self.merge_bins)-1)]
        self.primary_nbin=len(self.primary_bins)-1
        self.merge_nbin=len(self.merge_bins)-1

        if data==None:
            if self.region=="0bjet":
                data="data-tau_mi-vv-wjets-tt-st-qcdss-aa"
            elif self.region=="nbjet":
                data="data-mi-tau_mi-vv-wjets-ttlj-st-qcdss-aa"
            else:
                print "[Config::Init] Unknown region"
                exit(1)
        if sim==None:
            if self.region=="0bjet":
                sim="mi"
            elif self.region=="nbjet":
                sim="ttll"
            else:
                print "[Config::Init] Unknown region"
                exit(1)
        self.plotter=ROOT.AFBPlotter(data+" "+sim)

        self.outfilename=os.environ["SKFlat_WD"]+"/AFBResult/result_"+data.split("-",1)[0]+"_"+sim+"_"+histname.replace("/","_")+self.suffix+".root"

    def PrimaryAxisStr(self):
        axes=["x","y","z","u"]
        return axes[self.primary_axis]
    def MergeAxisStr(self):
        axes=["x","y","z","u"]
        return axes[self.merge_axis]

def GetUnfoldBinCenter(bins,i):
    if type(i) is ctypes.c_int:
        i=i.value
    x=None
    cost=None
    if i<1: 
        x=bins[0]-(bins[1]-bins[0])/2.
        cost=0.
    elif i>2*(len(bins)-1):
        x=bins[-1]+(bins[-1]-bins[-2])/2.
        cost=0.
    else:
        if i>len(bins)-1:
            cost=0.5
            i-=len(bins)-1
        else:
            cost=-0.5
        x=(bins[i]+bins[i-1])/2.
    return x,cost

def RebinResponseMatrix(matrix,xbins_old,ybins_old,xbins_new,ybins_new):
    rt=ROOT.TH2D(matrix.GetName(),matrix.GetTitle(),2*(len(xbins_new)-1),1,1+2*(len(xbins_new)-1),2*(len(ybins_new)-1),1,1+2*(len(ybins_new)-1))
    for i in range(matrix.GetNcells()):
        ix=ctypes.c_int(-1)
        iy=ctypes.c_int(-1)
        iz=ctypes.c_int(-1)
        matrix.GetBinXYZ(i,ix,iy,iz)
        x,costx=GetUnfoldBinCenter(xbins_old,ix)
        y,costy=GetUnfoldBinCenter(ybins_old,iy)
        if xbins_new[0]==0 and xbins_old[0]!=0:
            x=abs(x)
        if ybins_new[0]==0 and ybins_old[0]!=0:
            y=abs(y)
        ix_new=ROOT.AFBAnalyzer.GetUnfoldBin(len(xbins_new)-1,xbins_new,x,costx)
        iy_new=ROOT.AFBAnalyzer.GetUnfoldBin(len(ybins_new)-1,ybins_new,y,costy)        

        val0=matrix.GetBinContent(i)
        err0=matrix.GetBinError(i)

        val1=rt.GetBinContent(ix_new,iy_new)
        err1=rt.GetBinError(ix_new,iy_new)
        
        rt.SetBinContent(ix_new,iy_new,val0+val1)
        rt.SetBinError(ix_new,iy_new,(err0**2+err1**2)**0.5)
        #print i,x,costx,y,costy,ix,iy,ix_new,iy_new,val0,val1,rt.Integral()

    print matrix.Integral(),rt.Integral()
    return rt
        
        
def GetUnfoldedAFBHists(config):
    config.unfolded_afb_hists=[]
    for h in config.unfolded_hists:
        print h.GetName(),h.GetTitle()
        histname=str(h.GetName()).replace("unfolded","unfoldedafb")
        afb_hist=ROOT.TH1D(histname,histname,len(config.bins[config.primary_axis])-1,config.bins[config.primary_axis])
        if hasattr(h,"cov"):
            afb_hist.cov=ROOT.TH2D(histname+"_cov",histname+"_cov",len(config.bins[config.primary_axis])-1,config.bins[config.primary_axis],len(config.bins[config.primary_axis])-1,config.bins[config.primary_axis])
        n=afb_hist.GetNbinsX()
        for i in range(1,n+1):
            nf=h.GetBinContent(i+n)
            nb=h.GetBinContent(i)
            afb_hist.SetBinContent(i,(nf-nb)/(nf+nb))
            if hasattr(h,"cov"):
                ef2=h.cov.GetBinContent(i+n,i+n)
                eb2=h.cov.GetBinContent(i,i)
                efeb=h.cov.GetBinContent(i+n,i)
                for j in range(i,n+1):
                    nfj=h.GetBinContent(j+n)
                    nbj=h.GetBinContent(j)
                    efefj=h.cov.GetBinContent(i+n,j+n)
                    efebj=h.cov.GetBinContent(i+n,j)
                    ebefj=h.cov.GetBinContent(i,j+n)
                    ebebj=h.cov.GetBinContent(i,j)
                    afb_hist.cov.SetBinContent(i,j,4./(nf+nb)**2/(nfj+nbj)**2*(nb*nbj*efefj-nb*nfj*efebj-nf*nbj*ebefj+nf*nfj*ebebj))
                    if i!=j:
                        afb_hist.cov.SetBinContent(j,i,4./(nf+nb)**2/(nfj+nbj)**2*(nb*nbj*efefj-nb*nfj*efebj-nf*nbj*ebefj+nf*nfj*ebebj))
            else:
                ef2=h.GetBinError(i+n)**2
                eb2=h.GetBinError(i)**2
                efeb=0.
            afb_hist.SetBinError(i,2./(nf+nb)**2*(ef2*nb**2+eb2*nf**2-2*nf*nb*efeb)**0.5)
                
        config.unfolded_afb_hists+=[afb_hist]
    return config.unfolded_afb_hists

def ClosureTest(config):
    hist=config.plotter.GetHist(1,"mm2017/nbjet/genfid_myptcostRecoil_dressed")
    hist_check=config.plotter.GetHist(1,"mm2017/nbjet/genfid_myptcost_dressed_check")
    hist_check.SetLineColor(4)
    hist.Draw()
    hist_check.Draw("same")
    raw_input()
    for im in range(config.merge_nbin):
        histname=config.histname.replace("genfid_","").replace("_dressed","")
        print histname
        str_project=" project:{} ".format(config.PrimaryAxisStr())
        str_merge_bin=" {merge_axis}min:{low} {merge_axis}max:{high} ".format(merge_axis=config.MergeAxisStr().upper(),low=config.merge_bins[im],high=config.merge_bins[im+1])
        str_rebin= " rebin:{"+",".join(map(str,config.primary_bins))+"} "
        forward=config.plotter.GetHist(1,histname," Umin:0 Umax:1 "+str_project+str_merge_bin+str_rebin)
        backward=config.plotter.GetHist(1,histname," Umin:-1 Umax:0 "+str_project+str_merge_bin+str_rebin)
        hist_input=ROOT.TH1D(forward.GetName(),forward.GetTitle(),2*config.primary_nbin,1,2*config.primary_nbin+1)
        for i in range(1,forward.GetNbinsX()+1):
            ibin=ROOT.AFBAnalyzer.GetUnfoldBin(config.primary_nbin,config.primary_bins,forward.GetBinCenter(i),0.5)
            hist_input.SetBinContent(ibin,forward.GetBinContent(i)+hist_input.GetBinContent(ibin))
            hist_input.SetBinError(ibin,(forward.GetBinError(i)**2+hist_input.GetBinError(ibin)**2)**0.5)
            ibin=ROOT.AFBAnalyzer.GetUnfoldBin(config.primary_nbin,config.primary_bins,backward.GetBinCenter(i),-0.5)
            hist_input.SetBinContent(ibin,backward.GetBinContent(i)+hist_input.GetBinContent(ibin))
            hist_input.SetBinError(ibin,(backward.GetBinError(i)**2+hist_input.GetBinError(ibin)**2)**0.5)

        histname=config.histname
        print histname
        str_project=" project:{} ".format(config.PrimaryAxisStr())
        str_merge_bin=" {merge_axis}min:{low} {merge_axis}max:{high} ".format(merge_axis=config.MergeAxisStr().upper(),low=config.merge_bins[im],high=config.merge_bins[im+1]) 
        str_rebin= " rebin:{"+",".join(map(str,config.primary_bins))+"} "
        forward=config.plotter.GetHist(1,histname," Umin:0 Umax:1 "+str_project+str_merge_bin+str_rebin)
        backward=config.plotter.GetHist(1,histname," Umin:-1 Umax:0 "+str_project+str_merge_bin+str_rebin)
        hist_true=ROOT.TH1D(forward.GetName(),forward.GetTitle(),2*config.primary_nbin,1,2*config.primary_nbin+1)
        for i in range(1,forward.GetNbinsX()+1):
            ibin=ROOT.AFBAnalyzer.GetUnfoldBin(config.primary_nbin,config.primary_bins,forward.GetBinCenter(i),0.5)
            hist_true.SetBinContent(ibin,forward.GetBinContent(i)+hist_true.GetBinContent(ibin))
            hist_true.SetBinError(ibin,(forward.GetBinError(i)**2+hist_true.GetBinError(ibin)**2)**0.5)
            ibin=ROOT.AFBAnalyzer.GetUnfoldBin(config.primary_nbin,config.primary_bins,backward.GetBinCenter(i),-0.5)
            hist_true.SetBinContent(ibin,backward.GetBinContent(i)+hist_true.GetBinContent(ibin))
            hist_true.SetBinError(ibin,(backward.GetBinError(i)**2+hist_true.GetBinError(ibin)**2)**0.5)
        
        matrixname=config.matrixnames[im]
        print matrixname
        matrix=config.plotter.GetHist(1,matrixname,"noproject")

        c1=ROOT.TCanvas()
        hist_input.Draw()
        hist_input.SetLineColor(4)
        hy=matrix.ProjectionY("hy")
        hy.Draw("same")
        print "hist_input:",hist_input.Integral(),"hy:",hy.Integral()
        
        c2=ROOT.TCanvas()
        hist_true.Draw()
        hist_true.SetLineColor(4)
        hx=matrix.ProjectionX("hx")
        hx.Draw("same")
        print "hist_true:",hist_true.Integral(),"hx:",hx.Integral()

        raw_input()

def Unfold(matrix_orig,hist_orig,savecov=False):
    matrix=matrix_orig.Clone()
    hist=hist_orig.Clone()
    matrix.SetBinContent(0,0)
    matrix.SetBinError(0,0)
    hsub=matrix.ProjectionY("sub",0,0)
    hist.Add(hsub,-1)
    for i in range(matrix.GetYaxis().GetNbins()+2):
        matrix.SetBinContent(0,i,0)
        matrix.SetBinError(0,i,0)

    if DEBUG:
        hy=matrix.ProjectionY()
        hy.SetMarkerSize(0.4)
        hy.SetMarkerStyle(20)
        for i in range(hist.GetNcells()):
            #if i%30==0:
            #    hist.SetBinContent(i,hist.GetBinContent(i)/2)
            print i,hist.GetBinContent(i),hist.GetBinError(i),hist.GetBinContent(i)-hy.GetBinContent(i)
        hist.Draw()
        hist.SetLineColor(4)
        hist.SetLineWidth(2)
        hy.Draw("same")
        raw_input()
        #matrix.Draw("colz")
        #raw_input()
    unfold=ROOT.TUnfoldDensity(matrix,ROOT.TUnfold.kHistMapOutputHoriz,ROOT.TUnfold.kRegModeNone,ROOT.TUnfold.kEConstraintArea)
    unfold.SetInput(hist);
    unfold.DoUnfold(0)
    unfolded=unfold.GetOutput(hist.GetName()+"_unfolded")
    if savecov:
        unfolded.cov=unfold.GetEmatrixTotal(unfolded.GetName()+"_cov")
    if DEBUG:
        prob=unfold.GetProbabilityMatrix("prob")
        prob.Draw("colz")
        raw_input()
        unfolded.Draw()
        unfolded.SetLineColor(4)
        unfolded.SetLineWidth(2)
        hx=matrix.ProjectionX()
        hx.SetMarkerSize(0.4)
        hx.SetMarkerStyle(20)
        hx.Draw("same")
        print hx.Integral(1,30),hx.GetBinContent(0)
        print unfolded.Integral(1,30),unfolded.GetBinContent(0)
        raw_input()
        hist.Draw()
        folded=unfold.GetFoldedOutput("folded")
        folded.Draw("same")
        raw_input()
    return unfolded

def GetUnfoldedHists(config):
    config.unfolded_hists=[]
    for im in range(config.merge_nbin):
        histname=config.histname.replace("genfid_","").replace("_dressed","")
        str_project=" project:{} ".format(config.PrimaryAxisStr())
        str_merge_bin=" "
        if config.merge_axis is not None:
            str_merge_bin=" {merge_axis}min:{low} {merge_axis}max:{high} ".format(merge_axis=config.MergeAxisStr().upper(),low=config.merge_bins[im],high=config.merge_bins[im+1])
        str_rebin= " rebin:{"+",".join(map(str,config.primary_bins))+"} "
        print histname, config.option, str_project,str_merge_bin,str_rebin
        forward=config.plotter.GetHist(0,histname,config.option+" Umin:0 Umax:1 "+str_project+str_merge_bin+str_rebin)
        backward=config.plotter.GetHist(0,histname,config.option+" Umin:-1 Umax:0 "+str_project+str_merge_bin+str_rebin)
        hist_input=ROOT.TH1D(forward.GetName(),forward.GetTitle(),2*config.primary_nbin,1,2*config.primary_nbin+1)
        for i in range(1,forward.GetNbinsX()+1):
            ibin=ROOT.AFBAnalyzer.GetUnfoldBin(config.primary_nbin,config.primary_bins,forward.GetBinCenter(i),0.5)
            hist_input.SetBinContent(ibin,forward.GetBinContent(i)+hist_input.GetBinContent(ibin))
            hist_input.SetBinError(ibin,(forward.GetBinError(i)**2+hist_input.GetBinError(ibin)**2)**0.5)
            ibin=ROOT.AFBAnalyzer.GetUnfoldBin(config.primary_nbin,config.primary_bins,backward.GetBinCenter(i),-0.5)
            hist_input.SetBinContent(ibin,backward.GetBinContent(i)+hist_input.GetBinContent(ibin))
            hist_input.SetBinError(ibin,(backward.GetBinError(i)**2+hist_input.GetBinError(ibin)**2)**0.5)
        
        matrixname=config.matrixnames[im]
        print matrixname, config.option+" noproject"
        matrix=config.plotter.GetHist(1,matrixname,config.option+" noproject")
        matrix=RebinResponseMatrix(matrix,config.primary_bins_skflat,config.primary_bins_skflat,config.bins[config.primary_axis],config.primary_bins)
        hist_unfolded=Unfold(matrix,hist_input,config.suffix=="")
        hist_unfolded.SetName(matrixname.replace("response","unfolded")+config.suffix)
        hist_unfolded.SetTitle(matrixname.replace("response","unfolded")+config.suffix)
        if hasattr(hist_unfolded,"cov"):
            hist_unfolded.cov.SetName(hist_unfolded.GetName()+"_cov")
            hist_unfolded.cov.SetTitle(hist_unfolded.GetTitle()+"_cov")
        config.unfolded_hists+=[hist_unfolded]
        config.plotter.pdir=ROOT.TDirectory("plotdir","plotdir")
    return config.unfolded_hists

def Fluctuate(hist,name=""):
    rt=hist.Clone(name)
    rt.Reset()
    for i in range(hist.GetNcells()):
        rt.SetBinContent(i,ROOT.gRandom.Gaus(hist.GetBinContent(i),hist.GetBinError(i)*1.414))
        rt.SetBinError(i,hist.GetBinError(i))
    return rt

def GetBias(config):
    config.bias={}
    for im in range(config.merge_nbin):
        matrixname=config.matrixnames[im]
        print matrixname
        matrix=config.plotter.GetHist(1,matrixname,"noproject")
        matrix.Rebin2D(2,1) ## FIXME: define rebin function with bin array
        for i in range(100):
            print "toy",i
            toy=Fluctuate(matrix,matrixname+"_toy{}".format(i))
            hist_input=toy.ProjectionY(matrixname+"_toy{}".format(i)+"_py")
            hist_unfolded=Unfold(matrix,hist_input)
            hist_true=toy.ProjectionX(matrixname+"_toy{}".format(i)+"_px")

            for ibin in range(1,hist_true.GetNbinsX()+1):
                biasname=matrixname.replace("response_","bias_")+"_bin{}".format(ibin)
                if biasname not in config.bias:
                    config.bias[biasname]=ROOT.TH1D(biasname,biasname,40,-10,10)
                if hist_unfolded.GetBinError(ibin):
                    config.bias[biasname].Fill((hist_unfolded.GetBinContent(ibin)-hist_true.GetBinContent(ibin))/hist_unfolded.GetBinError(ibin))
                else:
                    print "zero error",biasname,hist_unfolded.GetBinContent(ibin),hist_unfolded.GetBinError(ibin),hist_true.GetBinContent(ibin)
                
    return config.bias

def GetHist4D(config):
    histname=config.histname+config.suffix
    print histname
    hist=ROOT.TH4D(histname,histname,
                   len(config.bins[0])-1,config.bins[0],
                   len(config.bins[1])-1,config.bins[1],
                   len(config.bins[2])-1,config.bins[2],
                   len(config.bins[3])-1,config.bins[3],)
    xs=[1500.,0.,300.]
    for im in range(len(config.merge_bins)-1):
        if config.merge_axis:
            xs[config.merge_axis]=(config.merge_bins[im]+config.merge_bins[im+1])/2
        n=len(config.bins[config.primary_axis])-1
        for i in range(n):
            xs[config.primary_axis]=(config.bins[config.primary_axis][i]+config.bins[config.primary_axis][i+1])/2
            ibin=hist.FindBin(xs[0],xs[1],xs[2],0.5)
            hist.SetBinContent(ibin,config.unfolded_hists[im].GetBinContent(1+i+n))
            hist.SetBinError(ibin,config.unfolded_hists[im].GetBinError(1+i+n))
            ibin=hist.FindBin(xs[0],xs[1],xs[2],-0.5)
            hist.SetBinContent(ibin,config.unfolded_hists[im].GetBinContent(1+i))
            hist.SetBinError(ibin,config.unfolded_hists[im].GetBinError(1+i))
                        
    config.hist4d=hist
    return hist


def WriteHist(hist):
    histname=hist.GetName()
    dirname=os.path.dirname(histname)
    basename=os.path.basename(histname)
    if not ROOT.gFile.Get(dirname):
        ROOT.gFile.mkdir(dirname)
    ROOT.gFile.cd(dirname)
    hist.Write(basename)
    if hasattr(hist,"cov"):
        WriteHist(hist.cov)

def Run(config):
    GetUnfoldedHists(config)
    GetUnfoldedAFBHists(config)
    GetHist4D(config)
    if not os.path.exists(os.path.dirname(config.outfilename)):
        os.makedirs(os.path.dirname(config.outfilename))
    f=ROOT.TFile(config.outfilename,"recreate")
    for h in config.unfolded_hists+config.unfolded_afb_hists+[config.hist4d]:
        print h.GetName()
        WriteHist(h)

    # GetBias(config)
    # hist_bias_mean=ROOT.TH1D(config.matrixname.replace("response_","bias_")+"_mean",config.matrixname.replace("response_","bias_")+"_mean",100,-1,1)
    # hist_bias_std=ROOT.TH1D(config.matrixname.replace("response_","bias_")+"_std",config.matrixname.replace("response_","bias_")+"_std",100,0,2)
    # for key in sorted(list(config.bias.keys())):
    #     WriteHist(config.bias[key])
    #     hist_bias_mean.Fill(config.bias[key].GetMean())
    #     hist_bias_std.Fill(config.bias[key].GetStdDev())
    # WriteHist(hist_bias_mean)
    # WriteHist(hist_bias_std)
         
def RunCondor(arg):
    os.system('condor_submit $SKFlat_WD/AFBResult/condor.jds -a arguments={}'.format(arg))
    
def Merge():
    files=os.listdir(os.environ["SKFlat_WD"]+"/AFBResult")
    files=filter(lambda x:x.startswith("result_"),files)
    n=900
    for i in range(len(files)/n+1):
        os.system("cd $SKFlat_WD/AFBResult; hadd -f final_{}.root {}".format(i," ".join(files[i*n:(i+1)*n])))
    os.system("hadd -f $SKFlat_WD/AFBResult/final.root $SKFlat_WD/AFBResult/final_*.root")
    os.system("rm $SKFlat_WD/AFBResult/final_*.root")
    return 

if __name__=="__main__":
    #ClosureTest(Config("mm2017/nbjet/genfid_dirapRecoil_dressed"))

    if len(sys.argv)>1:
        if sys.argv[1]=="merge":
            Merge()
            exit()
        else:
            c=Config(sys.argv[1])
            Run(c)
            exit()

    suffixes=list(SystematicSuffixes.keys())
    for channel in ["ee","mm"]:
        for era in ["2016a","2016b","2017","2018"]:
            #for suffix in [""]:
            for suffix in [""]+suffixes:
                print channel,era,suffix
                RunCondor(channel+era+"/0bjet/genfid_dimassCS_dressed"+suffix)
                RunCondor(channel+era+"/0bjet/genfid_diptCS_dressed"+suffix)
                RunCondor(channel+era+"/0bjet/genfid_dirapCS_dressed"+suffix)
                RunCondor(channel+era+"/nbjet/genfid_dimassRecoil_dressed"+suffix)
                RunCondor(channel+era+"/nbjet/genfid_diptRecoil_dressed"+suffix)
                RunCondor(channel+era+"/nbjet/genfid_dirapRecoil_dressed"+suffix)
    os.system("condor_wait $SKFlat_WD/AFBResult/condor.log")
    Merge()    
    exit()

    outputpath=os.environ["SKFlatOutputDir"]+os.environ["SKFlatV"]+"/AFBAnalyzer/"
