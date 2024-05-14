#!/usr/bin/env python

import array,os,sys,re
import ROOT,ctypes

ROOT.TH1.AddDirectory(0)
ROOT.TH1.SetDefaultSumw2(1)

ROOT.gROOT.ProcessLine('#include"AFBSystPlotter.cc"')
SystematicSuffixes=dict(ROOT.AFBSystPlotter("").GetSystematicSuffixes("totalsys"))
DEBUG=0
grid_mbin=ROOT.AFBAnalyzer.grid_mbin

class Config(object):
    def __init__(self,histname,suffix="",data=None,sim=None):
        self.histname=histname
        self.suffix=suffix
        self.option=""
        if self.suffix!="":
            self.option=SystematicSuffixes[self.suffix]

        words=histname.split("/")
        self.channel=words[0][:2]
        self.era=words[0][2:]
        self.region=words[1]
        if "dimass" in histname:
            self.matrixname=re.sub(r"/dimass","/response_afbm",self.histname)
            self.bins_matrix=array.array("d",ROOT.AFBAnalyzer.afb_mbin)
            if self.region=="0bjet":
                self.bins_reco=array.array("d",ROOT.AFBAnalyzer.unfold_0bjet_mbin_reco)
                self.bins_gen=array.array("d",ROOT.AFBAnalyzer.unfold_0bjet_mbin_gen)
            else:
                self.bins_reco=array.array("d",ROOT.AFBAnalyzer.unfold_nbjet_mbin_reco)
                self.bins_gen=array.array("d",ROOT.AFBAnalyzer.unfold_nbjet_mbin_gen)
        elif "dirap" in histname:
            self.matrixname=re.sub(r"/dirap","/response_afby",self.histname)
            self.bins_matrix=array.array("d",ROOT.AFBAnalyzer.afb_ybin)
            if self.region=="0bjet":
                self.bins_reco=array.array("d",ROOT.AFBAnalyzer.unfold_0bjet_ybin_reco)
                self.bins_gen=array.array("d",ROOT.AFBAnalyzer.unfold_0bjet_ybin_gen)
            else:
                self.bins_reco=array.array("d",ROOT.AFBAnalyzer.unfold_nbjet_ybin_reco)
                self.bins_gen=array.array("d",ROOT.AFBAnalyzer.unfold_nbjet_ybin_gen)
        elif "dipt" in histname:
            self.matrixname=re.sub(r"/dipt","/response_afbpt",self.histname)
            self.bins_matrix=array.array("d",ROOT.AFBAnalyzer.afb_ptbin)
            if self.region=="0bjet":
                self.bins_reco=array.array("d",ROOT.AFBAnalyzer.unfold_0bjet_ptbin_reco)
                self.bins_gen=array.array("d",ROOT.AFBAnalyzer.unfold_0bjet_ptbin_gen)
            else:
                self.bins_reco=array.array("d",ROOT.AFBAnalyzer.unfold_nbjet_ptbin_reco)
                self.bins_gen=array.array("d",ROOT.AFBAnalyzer.unfold_nbjet_ptbin_gen)
        self.nbin_reco=len(self.bins_reco)-1
        self.nbin_gen=len(self.bins_gen)-1

        if data==None:
            if self.region=="0bjet":
                data="data-tau_mi-vv-wjets-tt-st-qcd-aa"
            elif self.region=="nbjet":
                data="data-mi-tau_mi-vv-wjets-ttlj-st-qcd-aa"
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
        self.plotter=ROOT.AFBSystPlotter(data+" "+sim)
        self.hists=[]

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

    #print matrix.Integral(),rt.Integral()
    return rt
        
        
def GetUnfoldedAFBHist(config):
    h=config.unfolded_hist
    histname=str(h.GetName()).replace("unfolded","unfoldedafb")
    afb_hist=ROOT.TH1D(histname,histname,len(config.bins_gen)-1,config.bins_gen)
    if hasattr(h,"cov"):
        afb_hist.cov=ROOT.TH2D(histname+"_cov",histname+"_cov",len(config.bins_gen)-1,config.bins_gen,len(config.bins_gen)-1,config.bins_gen)
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
                
    config.unfolded_afb_hist=afb_hist
    config.hists+=[afb_hist]
    return config.unfolded_afb_hist

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

def GetInputHist(config):
    histname=config.histname
    str_project=" project:x "
    if "dirap" in histname: str_project+=" absx "
    str_rebin= " rebin:{"+",".join(map(str,config.bins_reco))+"} "
    print "[GetInputHist]",histname, config.option, str_project,str_rebin
    forward=config.plotter.GetHist(0,histname,config.option+" Ymin:0 Ymax:1 "+str_project+str_rebin)
    backward=config.plotter.GetHist(0,histname,config.option+" Ymin:-1 Ymax:0 "+str_project+str_rebin)
    input_hist=ROOT.TH1D(forward.GetName(),forward.GetTitle(),2*config.nbin_reco,1,2*config.nbin_reco+1)
    for i in range(1,forward.GetNbinsX()+1):
        ibin=ROOT.AFBAnalyzer.GetUnfoldBin(config.nbin_reco,config.bins_reco,forward.GetBinCenter(i),0.5)
        input_hist.SetBinContent(ibin,forward.GetBinContent(i)+input_hist.GetBinContent(ibin))
        input_hist.SetBinError(ibin,(forward.GetBinError(i)**2+input_hist.GetBinError(ibin)**2)**0.5)
        ibin=ROOT.AFBAnalyzer.GetUnfoldBin(config.nbin_reco,config.bins_reco,backward.GetBinCenter(i),-0.5)
        input_hist.SetBinContent(ibin,backward.GetBinContent(i)+input_hist.GetBinContent(ibin))
        input_hist.SetBinError(ibin,(backward.GetBinError(i)**2+input_hist.GetBinError(ibin)**2)**0.5)
    config.input_hist=input_hist
    return input_hist

def GetMatrix(config):
    matrixname=config.matrixname
    print matrixname, config.option+" noproject"
    matrix=config.plotter.GetHist(1,matrixname,config.option+" noproject")
    matrix=RebinResponseMatrix(matrix,config.bins_matrix,config.bins_matrix,config.bins_gen,config.bins_reco)
    config.matrix=matrix
    return matrix

def GetUnfoldedHist(config):
    matrix=config.matrix
    input_hist=config.input_hist
    histname=config.matrixname.replace("response","unfolded")
    unfolded_hist=Unfold(matrix,input_hist,config.suffix=="")
    unfolded_hist.SetName(histname+config.suffix)
    unfolded_hist.SetTitle(histname+config.suffix)
    if hasattr(unfolded_hist,"cov"):
        unfolded_hist.cov.SetName(unfolded_hist.GetName()+"_cov")
        unfolded_hist.cov.SetTitle(unfolded_hist.GetTitle()+"_cov")
    config.unfolded_hist=unfolded_hist
    config.plotter.pdir=ROOT.TDirectory("plotdir","plotdir")
    config.hists+=[config.unfolded_hist]
    return config.unfolded_hist

def GetGenAFBHist(config):
    matrixX=config.matrix.ProjectionX()
    histname=config.matrixname.replace("response_","gen_")+config.suffix
    gen_afb_hist=ROOT.TH1D(histname,histname,len(config.bins_gen)-1,config.bins_gen)
    n=gen_afb_hist.GetNbinsX()
    for i in range(1,n+1):
        nf=matrixX.GetBinContent(i+n)
        nb=matrixX.GetBinContent(i)
        gen_afb_hist.SetBinContent(i,(nf-nb)/(nf+nb))
        ef2=matrixX.GetBinError(i+n)**2
        eb2=matrixX.GetBinError(i)**2
        efeb=0.
        gen_afb_hist.SetBinError(i,2./(nf+nb)**2*(ef2*nb**2+eb2*nf**2-2*nf*nb*efeb)**0.5)
                
    config.gen_afb_hist=gen_afb_hist
    config.hists+=[gen_afb_hist]
    return gen_afb_hist

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
    dirname=os.path.dirname(histname).replace("201[678][ab]?","Run2")
    basename=os.path.basename(histname)
    if not ROOT.gFile.Get(dirname):
        ROOT.gFile.mkdir(dirname)
    ROOT.gFile.cd(dirname)
    hist.Write(basename)
    if hasattr(hist,"cov"):
        WriteHist(hist.cov)

def WriteHists(config,outfilename):
    f=ROOT.TFile(outfilename,"update")
    for h in config.hists:
        print h.GetName()
        WriteHist(h)
    f.Close()

def Run(suffix):
    if suffix=="nominal": suffix=""
    outfilename=os.environ["SKFlat_WD"]+"/AFBResult/result"+suffix+".root"
    if os.path.exists(outfilename):
        os.remove(outfilename)

    histnames=[]
    for channel in ["ee","mm","ll"]:
        for era in ["2016a","2016b","2017","2018","201[678][ab]?"]:
            #if channel=="ll" and era!="201[678][ab]?": continue
            for region in ["/0bjet","/nbjet"]:
                histnames+=[channel+era+region+"/dimass"]
                histnames+=[channel+era+region+"/dirap_m0"]
                histnames+=[channel+era+region+"/dirap_m1"]
                histnames+=[channel+era+region+"/dirap_m2"]
                histnames+=[channel+era+region+"/dirap_m3"]
                histnames+=[channel+era+region+"/dipt_m0"]
                histnames+=[channel+era+region+"/dipt_m1"]
                histnames+=[channel+era+region+"/dipt_m2"]
                histnames+=[channel+era+region+"/dipt_m3"]

    for histname in histnames:
        config=Config(histname,suffix)
        GetInputHist(config)
        GetMatrix(config)
        GetUnfoldedHist(config)
        GetUnfoldedAFBHist(config)
        GetGenAFBHist(config)
        #GetHist4D(config)
        if not os.path.exists(os.path.dirname(outfilename)):
            os.makedirs(os.path.dirname(outfilename))
        WriteHists(config,outfilename)

    # GetBias(config)
    # hist_bias_mean=ROOT.TH1D(config.matrixname.replace("response_","bias_")+"_mean",config.matrixname.replace("response_","bias_")+"_mean",100,-1,1)
    # hist_bias_std=ROOT.TH1D(config.matrixname.replace("response_","bias_")+"_std",config.matrixname.replace("response_","bias_")+"_std",100,0,2)
    # for key in sorted(list(config.bias.keys())):
    #     WriteHist(config.bias[key])
    #     hist_bias_mean.Fill(config.bias[key].GetMean())
    #     hist_bias_std.Fill(config.bias[key].GetStdDev())
    # WriteHist(hist_bias_mean)
    # WriteHist(hist_bias_std)

def AssertSame(hist1,hist2,check_error=True):
    for i in range(1,hist1.GetNbinsX()+1):
        val1=hist1.GetBinContent(i)
        err1=hist1.GetBinError(i)
        val2=hist2.GetBinContent(i)
        err2=hist2.GetBinError(i)
        if err1==0 or err2==0:
            print "[Warning] bin",i,"zero error"
            exit(1)
        if abs((val1-val2)/val2)>1e-6:
            print "[Warning] bin",i,"different value",val1,val2
            exit(1)
        if check_error:
            if abs((err1-err2)/err2)>1e-6:
                print "[Warning] bin",i,"different error",err1,err2
                exit(1)
    return

def TestClosure(suffix):
    if suffix=="nominal": suffix=""
    outfilename=os.environ["SKFlat_WD"]+"/AFBResult/result"+suffix+".root"
    if os.path.exists(outfilename):
        os.remove(outfilename)

    histnames=[]
    for channel in ["ee","mm"]:
        for era in ["2016a","2016b","2017","2018","201[678][ab]?"]:
            for region in ["/0bjet","/nbjet"]:
                histnames+=[channel+era+region+"/dimass"]
                histnames+=[channel+era+region+"/dirap_m0"]
                histnames+=[channel+era+region+"/dirap_m1"]
                histnames+=[channel+era+region+"/dirap_m2"]
                histnames+=[channel+era+region+"/dirap_m3"]
                histnames+=[channel+era+region+"/dipt_m0"]
                histnames+=[channel+era+region+"/dipt_m1"]
                histnames+=[channel+era+region+"/dipt_m2"]
                histnames+=[channel+era+region+"/dipt_m3"]

    for histname in histnames:
        toydata="mi" if "0bjet" in histname else "ttll"
        config=Config(histname,suffix,data=toydata)
        GetInputHist(config)
        GetMatrix(config)
        GetUnfoldedHist(config)
        GetUnfoldedAFBHist(config)
        GetGenAFBHist(config)
        if not os.path.exists(os.path.dirname(outfilename)):
            os.makedirs(os.path.dirname(outfilename))
        WriteHists(config,outfilename)

        ## input test
        matrixY=config.matrix.ProjectionY()
        # matrixY.SetOption("e2")
        # matrixY.SetLineColor(2)
        # matrixY.SetFillStyle(3002)
        # matrixY.SetFillColor(3)
        # matrixY.SetDirectory(0)
        #hists=[config.input_hist,matrixY]
        #c=config.plotter.DrawPlot(hists)
        #raw_input()
        AssertSame(config.input_hist,matrixY)

        ## Closure test
        # config.gen_afb_hist.SetOption("e2")
        # config.gen_afb_hist.SetLineColor(2)
        # config.gen_afb_hist.SetFillStyle(3002)
        # config.gen_afb_hist.SetFillColor(3)
        # config.gen_afb_hist.SetDirectory(0)
        # hists=[config.unfolded_afb_hist,config.gen_afb_hist]
        # c=config.plotter.DrawPlot(hists)
        # raw_input()
        AssertSame(config.unfolded_afb_hist,config.gen_afb_hist,check_error=False)

def RunCondor(arg):
    os.system('condor_submit $SKFlat_WD/AFBResult/condor.jds -a arguments={} > /dev/null'.format(arg))

def GetConditionNumber(histname):
    config=Config(histname)
    matrixname=config.matrixname
    matrix=config.plotter.GetHist(1,matrixname,config.option+" noproject")
    matrix=RebinResponseMatrix(matrix,config.bins_matrix,config.bins_matrix,config.bins_gen,config.bins_reco)
    proj=matrix.ProjectionX("proj",1,matrix.GetNbinsY())    
    response=ROOT.TMatrixD(1,matrix.GetNbinsX(),1,matrix.GetNbinsY())
    for i in range(1,matrix.GetNbinsX()+1):
        for j in range(1,matrix.GetNbinsY()+1):
            if proj.GetBinContent(i):
                response[i][j]=matrix.GetBinContent(i,j)/proj.GetBinContent(i)
            else:
                print histname,i
                response[i][j]=0.
    response.T()
    svd=ROOT.TDecompSVD(response)
    return svd.Condition()

def PrintConditionNumberAll():
    histnames=[]
    for channel in ["ee","mm"]:
        for era in ["2016a","2016b","2017","2018"]:
            for region in ["/0bjet","/nbjet"]:
                histnames+=[channel+era+region+"/dimass"]
                histnames+=[channel+era+region+"/dirap_m0"]
                histnames+=[channel+era+region+"/dirap_m1"]
                histnames+=[channel+era+region+"/dirap_m2"]
                histnames+=[channel+era+region+"/dirap_m3"]
                histnames+=[channel+era+region+"/dipt_m0"]
                histnames+=[channel+era+region+"/dipt_m1"]
                histnames+=[channel+era+region+"/dipt_m2"]
                histnames+=[channel+era+region+"/dipt_m3"]

    for histname in histnames:
        print histname,GetConditionNumber(histname)
    
def Merge():
    files=os.listdir(os.environ["SKFlat_WD"]+"/AFBResult")
    files=sorted(filter(lambda x:x.startswith("result"),files))
    n=900
    for i in range(len(files)/n+1):
        os.system("cd $SKFlat_WD/AFBResult; hadd -f final_{}.root {}".format(i," ".join(files[i*n:(i+1)*n])))
    os.system("hadd -f $SKFlat_WD/AFBResult/final.root $SKFlat_WD/AFBResult/final_*.root")
    os.system("rm $SKFlat_WD/AFBResult/final_*.root")
    return 

def SaveResponseAll(path):
    histnames=[]
    for channel in ["ee","mm"]:
        for era in ["2016a","2016b","2017","2018","201[678][ab]?"]:
        #for era in ["201[678][ab]?"]:
            for region in ["/0bjet","/nbjet"]:
                histnames+=[channel+era+region+"/dimass"]
                histnames+=[channel+era+region+"/dirap_m0"]
                histnames+=[channel+era+region+"/dirap_m1"]
                histnames+=[channel+era+region+"/dirap_m2"]
                histnames+=[channel+era+region+"/dirap_m3"]
                histnames+=[channel+era+region+"/dipt_m0"]
                histnames+=[channel+era+region+"/dipt_m1"]
                histnames+=[channel+era+region+"/dipt_m2"]
                histnames+=[channel+era+region+"/dipt_m3"]

    for histname in histnames:
        config=Config(histname)
        matrixname=config.matrixname
        matrix=config.plotter.GetHist(1,matrixname,config.option+" noproject")
        matrix=RebinResponseMatrix(matrix,config.bins_matrix,config.bins_matrix,config.bins_gen,config.bins_reco)
        proj=matrix.ProjectionX("proj",1,matrix.GetNbinsY())    
        response=ROOT.TMatrixD(1,matrix.GetNbinsX(),1,matrix.GetNbinsY())
        for i in range(1,matrix.GetNbinsX()+1):
            for j in range(1,matrix.GetNbinsY()+1):
                if proj.GetBinContent(i):
                    response[i][j]=matrix.GetBinContent(i,j)/proj.GetBinContent(i)
                else:
                    print histname,i
                    response[i][j]=0.
        c=ROOT.gROOT.MakeDefCanvas()
        h=ROOT.TH2D(response)
        h.SetStats(0)
        dilepton="#mu#mu" if "mm" in histname else "ee"
        variable="m("+dilepton+")"
        if "dirap" in histname: variable="y("+dilepton+")"
        elif "dipt" in histname: variable="p_{T}("+dilepton+")"
        h.GetYaxis().SetTitle("GEN "+variable+" bin index")
        h.GetXaxis().SetTitle("RECO "+variable+" bin index")
        h.SetMinimum(-1e-2)
        h.Draw("colz")
        if "2016a" in histname: era="2016preVFP"
        elif "2016b" in histname: era="2016postVFP"
        elif "2017" in histname: era="2017"
        elif "2018" in histname: era="2018"
        elif "201[678][ab]?" in histname: era="Run2"
        config.plotter.DrawPreliminary(c,era,"","nolumi")
        latex=ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextColor(ROOT.kGray)
        if "0bjet" in histname:
            latex.DrawLatex(0.17,0.8,"DY")
        else:
            latex.DrawLatex(0.17,0.8,"t#bar{t}")
        cut=""
        if "_m0" in histname:
            latex.DrawLatex(0.17,0.75,"{} #leq m < {} GeV".format(grid_mbin[0],grid_mbin[1]))
        elif "_m1" in histname:
            latex.DrawLatex(0.17,0.75,"{} #leq m < {} GeV".format(grid_mbin[1],grid_mbin[2]))
        elif "_m2" in histname:
            latex.DrawLatex(0.17,0.75,"{} #leq m < {} GeV".format(grid_mbin[2],grid_mbin[3]))
        elif "_m3" in histname:
            latex.DrawLatex(0.17,0.75,"{} #leq m < {} GeV".format(grid_mbin[3],grid_mbin[4]))
            
        print histname
        #raw_input()
        c.SaveAs(path+"/"+histname.replace("/","_").replace("201[678][ab]?","Run2")+".png")
        c.SaveAs(path+"/"+histname.replace("/","_").replace("201[678][ab]?","Run2")+".pdf")
    

if __name__=="__main__":
    # SaveResponseAll("")
    # exit()
    # TestClosure("nominal")

    if len(sys.argv)>1:
        if sys.argv[1]=="merge":
            Merge()
            exit()
        else:
            Run(sys.argv[1])
            exit()


    condor_jds=os.environ["SKFlat_WD"]+"/AFBResult/condor.jds"
    if not os.path.exists(condor_jds):
        if not os.path.exists(os.path.dirname(condor_jds)):
            os.makedirs(os.path.dirname(condor_jds))
        with open(condor_jds,"w") as f:
            f.write(
'''
executable = $ENV(SKFlat_WD)/python/AFBUnfold.py
log = $ENV(SKFlat_WD)/AFBResult/condor.log
output = $ENV(SKFlat_WD)/AFBResult/condor.out
error = $ENV(SKFlat_WD)/AFBResult/condor.err
getenv = true
should_transfer_files = yes
jobbatchname = AFBResult
queue 1
'''
            )

    suffixes=sorted(list(SystematicSuffixes.keys()))
    for suffix in ["nominal"]+suffixes:
        RunCondor(suffix)
    os.system("condor_wait $SKFlat_WD/AFBResult/condor.log")
    Merge()    
    exit()

    outputpath=os.environ["SKFlatOutputDir"]+os.environ["SKFlatV"]+"/AFBAnalyzer/"
