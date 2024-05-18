import os,array,copy
import numpy as np
import ROOT
ROOT.TH1.AddDirectory(0)
ROOT.TH1.SetDefaultSumw2(1)
ROOT.gROOT.ProcessLine('#include"AFBSystPlotter.cc"')

_plotter=ROOT.AFBSystPlotter("mi ttll")
SYSNAME="totalsys"
#SYSNAME="dypdf"
SystematicSuffixes=dict(_plotter.GetSystematicSuffixes(SYSNAME)).keys()
Systematics=dict(_plotter.systematics)
grid_mbin=ROOT.AFBAnalyzer.grid_mbin

def Variation2Suffix(variation):
    return str(ROOT.Plotter.Variation2Suffix(variation))

class AFBMeasurement:
    def __init__(self,value,stat):
        self.value=value
        self.stat=stat
        self.syst={}
    
    def GetValue(self):
        return self.value

    def GetStatError(self):
        return self.stat
    
    def SetSystError(self,key,error,isValue=False):
        if isValue:
            error-=self.value
        self.syst[key]=error
        return

    def GetSystError(self,key):
        return self.syst[key]    

    def __mul__(self,other):
        rt=copy.deepcopy(self)
        rt.value*=other
        rt.stat*=other
        for key in rt.syst:
            rt.syst[key]*=other
        return rt
    __rmul__=__mul__

    def __add__(self,other):
        rt=AFBMeasurement(self.value+other.value,(self.stat**2+other.stat**2)**0.5)
        keys=set(self.syst.keys()+other.syst.keys())
        for key in keys:
            rt.syst[key]=self.syst[key]+other.syst[key]
        return rt

    def __sub__(self,other):
        return self.__add__(-1*other)


class AFBMeasurements:
    def __init__(self,filename,histname,combine="blue"):
        self.combine=combine
        self.default_syst=SYSNAME
        self.histname=histname
        f=ROOT.TFile(filename)
        h=f.Get(histname)
        print histname
        self.measurements=[]
        ncells=h.GetNcells()
        for i in range(ncells):
            self.measurements+=[AFBMeasurement(h.GetBinContent(i),h.GetBinError(i))]
        self.bins=array.array("d",[h.GetBinLowEdge(i) for i in range(1,h.GetNbinsX()+2)])
        hcov=f.Get(histname+"_cov")
        if hcov:
            self.cov_stat=np.zeros((ncells,ncells))
            for i in range(ncells):
                for j in range(ncells):
                    self.cov_stat[i][j]=hcov.GetBinContent(i,j)
        else:
            self.cov_stat=np.zeros((ncells,ncells))
            for i in range(ncells):
                self.cov_stat[i][i]=h.GetBinError(i)**2
            
        for suffix in SystematicSuffixes:
            h=f.Get(histname+suffix)
            if not h:
                print histname+suffix
            for i in range(ncells):
                self.measurements[i].SetSystError(suffix,h.GetBinContent(i),isValue=True)
        self.syst={}
        self.dist=None
        if "unfoldedafb" in histname:
            self.dist=AFBMeasurements(filename,histname.replace("unfoldedafb","unfolded"),combine=combine)
        return

    def __len__(self):
        return len(self.measurements)

    def __str__(self,syst=[]):
        if len(syst)==0:
            syst=[self.default_syst]
        out=[]
        out+=[ "\t".join( ["bin","val","stat"]+syst ) ]
        for i in range(len(self.measurements)):
            out+=[ "\t".join( map(str,[i,self.GetValue(i),self.GetStatError(i)]+[self.GetSystError(i,key) for key in syst]) ) ]
        return "\n".join(out)

    def Dist2AFB(self):
        if not self.dist:
            print "No distribution"
            return
        ncells=len(self.measurements)
        n=ncells-2
        self.measurements=[AFBMeasurement(0,0) for i in range(ncells)]
        for i in range(1,n+1):
            nf=self.dist.GetValue(i+n)
            nb=self.dist.GetValue(i)
            afb=(nf-nb)/(nf+nb)
            self.measurements[i].value=afb
            ef2=self.dist.cov_stat[i+n,i+n]
            eb2=self.dist.cov_stat[i,i]
            efeb=self.dist.cov_stat[i+n,i]
            for j in range(i,n+1):
                nfj=self.dist.GetValue(j+n)
                nbj=self.dist.GetValue(j)
                efefj=self.dist.cov_stat[i+n,j+n]
                efebj=self.dist.cov_stat[i+n,j]
                ebefj=self.dist.cov_stat[i,j+n]
                ebebj=self.dist.cov_stat[i,j]
                self.cov_stat[i,j]=4./(nf+nb)**2/(nfj+nbj)**2*(nb*nbj*efefj-nb*nfj*efebj-nf*nbj*ebefj+nf*nfj*ebebj)
                if i!=j:
                    self.cov_stat[j,i]=4./(nf+nb)**2/(nfj+nbj)**2*(nb*nbj*efefj-nb*nfj*efebj-nf*nbj*ebefj+nf*nfj*ebebj)
                else:
                    self.measurements[i].error=self.cov_stat[i,i]**0.5
        for suffix in SystematicSuffixes:
            for i in range(ncells):
                if i==0 or i==ncells-1:
                    afb=0
                else:
                    nf=self.dist.measurements[i+n].value+self.dist.measurements[i+n].syst[suffix]
                    nb=self.dist.measurements[i].value+self.dist.measurements[i].syst[suffix]
                    afb=(nf-nb)/(nf+nb)
                self.measurements[i].SetSystError(suffix,afb,isValue=True)
        self.syst={}
        return

    def EvalSystError(self,key):
        key=str(key)
        syst=Systematics[key]
        if syst.type==ROOT.Systematic.Type.MULTI:
            subkeys=map(str,syst.keys)
            for subkey in subkeys:
                self.EvalSystError(subkey)
            self.syst[key]=sum([self.syst[subkey] for subkey in subkeys])
        else:
            self.syst[key]=AFBMeasurements.CalcCov(key,self.measurements)
        return

    @staticmethod
    def CalcCov(key,measurements):
        n=len(measurements)
        rt=np.zeros((n,n))
        key=str(key)
        syst=Systematics[key]
        if syst.type==ROOT.Systematic.Type.MULTI:
            subkeys=map(str,syst.keys)
            for subkey in subkeys:
                rt+=AFBMeasurements.CalcCov(subkey,measurements)
        elif syst.type==ROOT.Systematic.Type.ENVELOPE:
            subkeys=map(Variation2Suffix,syst.variations)
            vec=max([[x.GetSystError(subkey) for x in measurements] for subkey in subkeys],key=sum)
            rt=np.outer(vec,vec)
        elif syst.type==ROOT.Systematic.Type.GAUSSIAN:
            subkeys=map(Variation2Suffix,syst.variations)
            vecs=[[x.GetSystError(subkey) for x in measurements] for subkey in subkeys]
            covs=map(lambda x:np.outer(x,x),vecs)
            rt=sum(covs)/len(covs)
        elif syst.type==ROOT.Systematic.Type.HESSIAN:
            subkeys=map(Variation2Suffix,syst.variations)
            vecs=[[x.GetSystError(subkey) for x in measurements] for subkey in subkeys]
            covs=map(lambda x:np.outer(x,x),vecs)
            rt=sum(covs)
        elif syst.type==ROOT.Systematic.Type.CORRELATED:
            subkeys=map(Variation2Suffix,syst.variations)
            vec=[sum([x.GetSystError(subkey) for subkey in subkeys]) for x in measurements]
            rt=np.outer(vec,vec)
        else:
            print "[Error] Unknown systematic type",syst.type,syst.title
            exit(1)
        return rt

    def GetSystError(self,i,key=""):
        if key=="":
            key=self.default_syst
        if key not in self.syst:
            self.EvalSystError(key)
        return self.syst[key][i][i]**0.5

    def GetValue(self,i):
        return self.measurements[i].GetValue()

    def GetStatError(self,i):
        return self.cov_stat[i][i]**0.5

    def GetTotalError(self,i):
        return (self.GetStatError(i)**2+self.GetSystError(i,self.default_syst)**2)**0.5

    def GetSystCov(self,key,i,j):
        if key not in self.syst:
            self.EvalSystError(key)
        return self.syst[key][i][j]
        
    def GetCov(self,i,j):
        return self.cov_stat[i][j]+self.GetSystCov(self.default_syst,i,j)

    def GetHists(self):
        hist_stat=ROOT.TH1D(self.histname,self.histname,len(self.bins)-1,self.bins)
        hist_stat.SetDirectory(0)
        hist_total=ROOT.TH1D(self.histname,self.histname,len(self.bins)-1,self.bins)
        hist_total.SetDirectory(0)
        for i in range(hist_stat.GetNcells()):
            hist_stat.SetBinContent(i,self.GetValue(i))
            hist_stat.SetBinError(i,self.GetStatError(i))
            hist_total.SetBinContent(i,self.GetValue(i))
            hist_total.SetBinError(i,self.GetTotalError(i))
        hists=ROOT.Hists()
        hists.push_back(hist_stat)
        hists.push_back(hist_total)
        hists.save=[hist_stat,hist_total]
        return hists

    def __sub__(self,other):
        rt=copy.deepcopy(self)
        rt.measurements=[rt.measurements[i]-other.measurements[i] for i in range(len(rt.measurements))]
        rt.cov_stat+=other.cov_stat
        rt.syst={}
        return rt

    def __add__(self,other):
        rt=copy.deepcopy(self)
        rt.measurements=[rt.measurements[i]+other.measurements[i] for i in range(len(rt.measurements))]
        rt.cov_stat+=other.cov_stat
        rt.syst={}
        return rt

    def __mod__(self,other):
        rt=copy.deepcopy(self)
        n=len(rt.measurements)
        if n!=len(other.measurements):
            print "Inconsistent number of bins",n,len(other.measurements)
        cov_stat0=copy.deepcopy(self.cov_stat)
        cov_stat1=copy.deepcopy(other.cov_stat)
        for i in range(n):
            if self.GetStatError(i)==0:
                if other.GetStatError(i)==0:
                    ws=[0.5,0.5]
                else:
                    ws=[0.,1.]
            else:
                if other.GetStatError(i)==0:
                    ws=[1.,0.]
                else:
                    if self.combine=="simple":
                        ws=[1/self.cov_stat[i][i],1/other.cov_stat[i][i]]
                    elif self.combine=="blue":
                        cov=AFBMeasurements.CalcCov(self.default_syst,[self.measurements[i],other.measurements[i]])
                        cov[0][0]+=self.cov_stat[i][i]
                        cov[1][1]+=other.cov_stat[i][i]
                        try:
                            covI=np.linalg.inv(cov)
                        except:
                            print cov
                            print self.GetValue(i), self.GetStatError(i), other.GetValue(i), other.GetStatError(i)
                            exit()
                        ws=np.matmul(covI,np.ones((2,1)))
                        ws=ws.T[0]
                    else:
                        print "Unknown combine method",self.combine
                        exit(1)
                    ws/=sum(ws)
            rt.measurements[i]=self.measurements[i]*ws[0]+other.measurements[i]*ws[1]
            cov_stat0[i,:]*=ws[0]
            cov_stat0[:,i]*=ws[0]
            cov_stat1[i,:]*=ws[1]
            cov_stat1[:,i]*=ws[1]
        rt.cov_stat=cov_stat0+cov_stat1
        rt.syst={}
        if rt.dist:
            if other.dist:
                rt.dist+=other.dist
            else:
                print "[Error] [AFBMeasurements::Combined] there is self.dist but no other.dist"
                exit(1)
        return rt
        
    def Compare(self,other):
        diff=self-other
        diff.EvalSystError(self.default_syst)
        cov=diff.syst[self.default_syst]+diff.cov_stat
        try:
            covI=np.linalg.inv(cov[1:-1,1:-1])
        except:
            return {}
        d=np.array([[diff.GetValue(i) for i in range(1,len(diff.measurements)-1)]])
        chi2=np.matmul(np.matmul(d,covI),d.T)[0][0]
        ndf=len(d[0])
        return {"chi2":chi2, "ndf":ndf, "pvalue":ROOT.TMath.Prob(chi2,ndf)}

    def Save(self,filename=""):
        systs=["efficiencySF","bcharge","electronenergy","muonmomentum","dytheory","tttheory"]
        out=["\\begin{tabular} { "+" ".join(["r"]*(6+len(systs)))+" }"]
        out+=["bin & min & max & val & stat & syst & "+" & ".join(systs)+" \\\\"]
        for i in range(1,len(self.measurements)-1):
            #line="{} & {} & {} & {} & {} & {} & ".format(i,self.bins[i-1],self.bins[i],round(self.GetValue(i),4),round(self.GetStatError(i),4),round(self.GetSystError(i),4)) \
            line="{} & {} & {} & {} & {} & {} & ".format(i,self.bins[i-1],self.bins[i],"blind",round(self.GetStatError(i),4),round(self.GetSystError(i),4)) \
            +" & ".join([format(self.GetSystError(i,syst),".4f") for syst in systs]) \
            +" \\\\"
            out+=[line]
        out+=["\\end{tabular}"]
        out="\n".join(out)
        if filename!="":
            dirname=os.path.dirname(filename)
            if not os.path.exists(dirname):
                os.makedirs(dirname)
            with open(filename,"w") as f:
                f.write(out)
        else:
            print out
        return

    def SaveCov(self,filename=""):
        out=["\\begin{tabular} { "+" ".join(["r" for i in range(1,len(self.measurements)-1+4)])+" }"]
        out+=["bin & min & max & val & "+" & ".join([str(i) for i in range(1,len(self.measurements)-1)])+" \\\\"]
        for i in range(1,len(self.measurements)-1):
            #line="{} & {} & {} & {} & ".format(i,self.bins[i-1],self.bins[i],round(self.GetValue(i),4)) \
            line="{} & {} & {} & {} & ".format(i,self.bins[i-1],self.bins[i],"blind") \
            +" & ".join([format(self.GetCov(i,j),".8f") for j in range(1,len(self.measurements)-1)]) \
            +" \\\\"
            out+=[line]
        out+=["\\end{tabular}"]
        # out="\n".join(map(lambda x:x.replace("\\","\\\\"),out))
        out="\n".join(out)
        if filename!="":
            dirname=os.path.dirname(filename)
            if not os.path.exists(dirname):
                os.makedirs(dirname)
            with open(filename,"w") as f:
                f.write(out)
        else:
            print out
        return
        
        
def SaveCompareEraAll(inputpath,outputpath):
    colors=[ROOT.kBlack,ROOT.kRed,ROOT.kGreen+1,ROOT.kBlue,ROOT.kYellow+1,ROOT.kMagenta,ROOT.kCyan,ROOT.kGray,ROOT.kPink+1,ROOT.kSpring+1,ROOT.kAzure+1,ROOT.kOrange+1,ROOT.kViolet+1,ROOT.kTeal+1,ROOT.kWhite,ROOT.kGray+1,ROOT.kGray+3]
    _plotter.plotdir=outputpath
    for channel in ["mm","ee"]:
        for region in ["0bjet","nbjet"]:
            for histname in ["unfoldedafb_afbm","unfoldedafb_afby_m0","unfoldedafb_afby_m1","unfoldedafb_afby_m2","unfoldedafb_afby_m3","unfoldedafb_afbpt_m0","unfoldedafb_afbpt_m1","unfoldedafb_afbpt_m2","unfoldedafb_afbpt_m3"]:
                ms=[AFBMeasurements(inputpath,channel+era+"/"+region+"/"+histname) for era in ["2016a","2016b","2017","2018"]]
                ms=[ms[0]%ms[1]%ms[2]%ms[3]]+ms
                hists=map(lambda x:x.GetHists(),ms)
                for i,hh in enumerate(hists):
                    for h in hh:
                        h.SetLineColor(colors[i])
                        h.SetOption("e1")
                        h.SetName(h.GetName().split("/")[0])
                        if i==0:
                            h.SetName("Run2")
                            h.SetLineWidth(2)
                p=ROOT.Plot()
                p.SetOption("TLleg")
                if "afbm" in histname:
                    p.SetOption("logx")
                    p.SetOption("xtitle:'m(ll) [GeV]'")
                elif "afby" in histname:
                    p.SetOption("xtitle:'y(ll)'")
                elif "afbpt" in histname:
                    p.SetOption("logx")
                    p.SetOption("xtitle:'p_{T}(ll) [GeV]'")
                    
                p.hists=ROOT.vector("Hists")(hists)
                c=ROOT.TCanvas()
                _plotter.DrawDiff(p)
                #raw_input()
                _plotter.SaveCanvas(c,"diff/"+channel+"_"+region+"_"+histname+".png")

def SaveUnfoldedPlotAll(inputpath,outputpath):
    _plotter.plotdir=outputpath
    for region in ["0bjet","nbjet"]:
        for histname in ["unfoldedafb_afbm","unfoldedafb_afby_m0","unfoldedafb_afby_m1","unfoldedafb_afby_m2","unfoldedafb_afby_m3","unfoldedafb_afbpt_m0","unfoldedafb_afbpt_m1","unfoldedafb_afbpt_m2","unfoldedafb_afbpt_m3"]:
            print region,histname
            mdatas=[]
            msims=[]
            for channel in ["ee","mm"]:
                for era in ["2016a","2016b","2017","2018"]:
                #for era in ["Run2"]:
                    mdatas+=[AFBMeasurements(inputpath,channel+era+"/"+region+"/"+histname)]
                    msims+=[AFBMeasurements(inputpath,channel+era+"/"+region+"/"+histname.replace("unfoldedafb_","gen_"))]
            mdata=reduce(lambda x,y:x%y,mdatas)
            msim=reduce(lambda x,y:x%y,msims)
            #m=ms[-1]
            hdata=mdata.GetHists()
            hsim=msim.GetHists()

            if region=="0bjet":
                plot_option="ytitle:'A_{FB}^{CS}'"
                plot_option+=" ymin:-0.09 ymax:0.49"
                for i in range(len(hsim)):
                    hsim[i].SetOption("hist e1")
                    hsim[i].SetLineColor(2)
                    hsim[i].SetName("DY POWHEG MiNNLO_{PS}+Pythia8+PHOTOS")                    
            else:
                plot_option="ytitle:'A_{FB}^{Recoil}'"
                plot_option+=" ymin:0.0 ymax:1.49"
                for i in range(len(hsim)):
                    hsim[i].SetOption("hist e1")
                    hsim[i].SetLineColor(6)
                    hsim[i].SetName("t\bar{t} POWHEG+Pythia8")

            if "afbm" in histname:
                plot_option+=" xtitle:m(ll) logx"
            elif "afby" in histname:
                plot_option+=" xtitle:y(ll)"
            elif "afbpt" in histname:
                plot_option+=" xtitle:p_{T}(ll) logx"


            for i in range(len(hdata)):
                h=hdata[i]
                h.SetName("data (blind)")
                h.SetOption("e1")
                h.SetMarkerStyle(20)
                h.SetMarkerSize(0.7)
                h.SetLineColor(1)
                h.SetMarkerColor(1)
                for ib in range(h.GetNcells()):
                   h.SetBinContent(ib,hsim[0].GetBinContent(ib))

            p=ROOT.Plot()
            p.SetOption(plot_option)
            p.hists=ROOT.vector("Hists")([hdata,hsim])
            c=ROOT.TCanvas()
            _plotter.DrawCompare(p)
            _plotter.DrawPreliminary(c,"Run2")
            savename=histname.replace("unfoldedafb_afb","afb_")
            savename=savename.replace("y_m","y").replace("pt_m","pt")
            _plotter.SaveCanvas(c,region+"_"+savename+".png",False)
            _plotter.SaveCanvas(c,region+"_"+savename+".pdf")
            #raw_input()

def SaveDeltaPlotAll(inputpath,outputpath):
    _plotter.plotdir=outputpath
    for region in ["0bjet","nbjet"]:
        for histname in ["unfoldedafb_afbm","unfoldedafb_afby_m0","unfoldedafb_afby_m1","unfoldedafb_afby_m2","unfoldedafb_afby_m3","unfoldedafb_afbpt_m0","unfoldedafb_afbpt_m1","unfoldedafb_afbpt_m2","unfoldedafb_afbpt_m3"]:
            print region,histname
            mdatas=[]
            msims=[]
            for channel in ["ee","mm"]:
                mdatass=[]
                msimss=[]
                for era in ["2016a","2016b","2017","2018"]:
                #for era in ["Run2"]:
                    mdatass+=[AFBMeasurements(inputpath,channel+era+"/"+region+"/"+histname)]
                    msimss+=[AFBMeasurements(inputpath,channel+era+"/"+region+"/"+histname.replace("unfoldedafb_","gen_"))]
                mdatas+=[reduce(lambda x,y:x%y,mdatass)]
                msims+=[reduce(lambda x,y:x%y,msimss)]
            mdata=mdatas[1]-mdatas[0]
            msim=msims[1]-msims[0]
            hdata=mdata.GetHists()
            hsim=msim.GetHists()

            if region=="0bjet":
                plot_option="ytitle:'#Delta A_{FB}^{CS}'"
                for i in range(len(hsim)):
                    hsim[i].SetOption("hist e1")
            else:
                plot_option="ytitle:'#Delta A_{FB}^{Recoil}'"
                for i in range(len(hsim)):
                    hsim[i].SetOption("hist e1")

            if "afbm" in histname:
                plot_option+=" xtitle:m(ll) logx"
            elif "afby" in histname:
                plot_option+=" xtitle:y(ll)"
            elif "afbpt" in histname:
                plot_option+=" xtitle:p_{T}(ll) logx"

            for i in range(len(hdata)):
                h=hdata[i]
                h.SetName("data (blind)")
                h.SetOption("e1")
                h.SetMarkerStyle(20)
                h.SetMarkerSize(0.7)
                h.SetLineColor(1)
                h.SetMarkerColor(1)
                for ib in range(h.GetNcells()):
                    h.SetBinContent(ib,0.)

            p=ROOT.Plot()
            p.SetOption(plot_option)
            #p.hists=ROOT.vector("Hists")([hdata,hsim])
            p.hists=ROOT.vector("Hists")([hdata])
            c=ROOT.TCanvas()
            _plotter.DrawSig(p)
            _plotter.DrawPreliminary(c,"Run2")
            savename=histname.replace("unfoldedafb_afb","dafb_")
            savename=savename.replace("y_m","y").replace("pt_m","pt")
            _plotter.SaveCanvas(c,region+"_"+savename+".png",False)
            _plotter.SaveCanvas(c,region+"_"+savename+".pdf")
            #raw_input()

def SaveTableAll(inputpath,outputpath):
    for region in ["0bjet","nbjet"]:
        for histname in ["unfoldedafb_afbm","unfoldedafb_afby_m0","unfoldedafb_afby_m1","unfoldedafb_afby_m2","unfoldedafb_afby_m3","unfoldedafb_afbpt_m0","unfoldedafb_afbpt_m1","unfoldedafb_afbpt_m2","unfoldedafb_afbpt_m3"]:
            ms=[]
            for channel in ["ee","mm"]:
                #for era in ["Run2"]:
                for era in ["2016a","2016b","2017","2018"]:
                    ms+=[AFBMeasurements(inputpath,channel+era+"/"+region+"/"+histname)]
            m=reduce(lambda x,y:x%y,ms)
            m.Save(outputpath+"/"+region+"/"+histname+".tex")

def SaveTableCovAll(inputpath,outputpath):
    for region in ["0bjet","nbjet"]:
        for histname in ["unfoldedafb_afbm","unfoldedafb_afby_m0","unfoldedafb_afby_m1","unfoldedafb_afby_m2","unfoldedafb_afby_m3","unfoldedafb_afbpt_m0","unfoldedafb_afbpt_m1","unfoldedafb_afbpt_m2","unfoldedafb_afbpt_m3"]:
            ms=[]
            for channel in ["ee","mm"]:
                #for era in ["Run2"]:
                for era in ["2016a","2016b","2017","2018"]:
                    ms+=[AFBMeasurements(inputpath,channel+era+"/"+region+"/"+histname)]
            m=reduce(lambda x,y:x%y,ms)
            m.SaveCov(outputpath+"/"+region+"/"+histname+"_cov.tex")

def SaveCorrelationAll(inputpath,outputpath):
    _plotter.plotdir=outputpath
    for region in ["0bjet","nbjet"]:
        for histname in ["unfoldedafb_afbm","unfoldedafb_afby_m0","unfoldedafb_afby_m1","unfoldedafb_afby_m2","unfoldedafb_afby_m3","unfoldedafb_afbpt_m0","unfoldedafb_afbpt_m1","unfoldedafb_afbpt_m2","unfoldedafb_afbpt_m3"]:
            ms=[]
            for channel in ["ee","mm"]:
                #for era in ["Run2"]:
                for era in ["2016a","2016b","2017","2018"]:
                    ms+=[AFBMeasurements(inputpath,channel+era+"/"+region+"/"+histname)]
            m=reduce(lambda x,y:x%y,ms)
            cov=[[m.GetCov(i,j) for i in range(len(m))] for j in range(len(m))]
            h=ROOT.TH2D("correlation","correlation",len(m)-2,1,len(m)-1,len(m)-2,1,len(m)-1)
            for i in range(1,len(m)-1):
                for j in range(1,len(m)-1):
                    h.SetBinContent(i,j,m.GetCov(i,j)/(m.GetCov(i,i)*m.GetCov(j,j))**0.5)
            c=ROOT.gROOT.MakeDefCanvas()
            h.SetStats(0)
            variable="m(ll)"
            if "dirap" in histname: variable="y(ll)"
            elif "dipt" in histname: variable="p_{T}(ll)"
            h.GetYaxis().SetTitle(variable+" bin index")
            h.GetXaxis().SetTitle(variable+" bin index")
            h.SetMaximum(1.01)
            h.SetMinimum(-1.01)
            h.Draw("colz")
            _plotter.DrawPreliminary(c,"Run2")
            latex=ROOT.TLatex()
            latex.SetNDC()
            latex.SetTextColor(ROOT.kWhite)
            if "0bjet" in histname:
                latex.DrawLatex(0.17,0.8,"A_{FB}^{CS}")
            else:
                latex.DrawLatex(0.17,0.8,"A_{FB}^{Recoil}")
            cut=""
            if "_m0" in histname:
                latex.DrawLatex(0.17,0.75,"{} #leq m < {} GeV".format(grid_mbin[0],grid_mbin[1]))
            elif "_m1" in histname:
                latex.DrawLatex(0.17,0.75,"{} #leq m < {} GeV".format(grid_mbin[1],grid_mbin[2]))
            elif "_m2" in histname:
                latex.DrawLatex(0.17,0.75,"{} #leq m < {} GeV".format(grid_mbin[2],grid_mbin[3]))
            elif "_m3" in histname:
                latex.DrawLatex(0.17,0.75,"{} #leq m < {} GeV".format(grid_mbin[3],grid_mbin[4]))
            savename=histname.replace("unfoldedafb_afb","afb_")
            savename=savename.replace("y_m","y").replace("pt_m","pt")
            savename+="_correlation"
            _plotter.SaveCanvas(c,region+"_"+savename+".png",False)
            _plotter.SaveCanvas(c,region+"_"+savename+".pdf")
    return

def Compare(mss,base=0):
    N=len(mss)
    gss=[]
    for ms in mss:
        n=len(ms)
        gs=[ROOT.TGraphErrors() for _ in range(n*2)]
        colors=[1,2,4,6,7,8,9]
        for i in range(n):
            gs[i].SetLineColor(colors[i])
            gs[i+n].SetLineColor(colors[i])
        for i in range(1,len(ms[0])-1):
            vals=np.array([m.GetValue(i) for m in ms])
            errs=np.array([m.GetTotalError(i) for m in ms])
            errs_stat=np.array([m.GetStatError(i) for m in ms])
            if errs[1]==0: continue
            vals,errs,errs_stat=(vals-vals[base])/errs[base],errs/errs[base],errs_stat/errs[base]
            width=0.5/(n-1) if n>0 else 0
            start=width/2
            for j in range(n):
                gs[j].SetPoint(i-1,i-start+width*j,vals[j])
                gs[j].SetPointError(i-1,0,errs[j])
                gs[n+j].SetPoint(i-1,i-start+width*j,vals[j])
                gs[n+j].SetPointError(i-1,0,errs_stat[j])
        gss+=[gs]
    c=ROOT.TCanvas()
    c.hists=[]
    c.Divide(1,N)
    top=0.9
    margin=0.03
    hight=(0.8-3*margin)/N
    for i in range(N):
        c.cd(i+1)
        ROOT.gPad.SetPad(0,top,1,top-hight-margin)
        ROOT.gPad.SetTopMargin(0)
        ROOT.gPad.SetBottomMargin(margin/(hight+margin))
        ROOT.gPad.SetFillStyle(0)
        n=gss[i][0].GetN()
        sigma=ROOT.TH1D("sigma","",100,0.1,n+0.9)
        for ib in range(sigma.GetNcells()):
            sigma.SetBinError(ib,1)
        sigma.SetFillStyle(3001)
        sigma.SetFillColor(ROOT.kGreen)
        sigma.GetYaxis().SetTickLength(0.003)
        sigma.GetYaxis().SetRangeUser(-2.99,2.99)
        sigma.GetYaxis().SetNdivisions(205)
        sigma.GetYaxis().SetLabelSize(0.2)
        sigma.GetXaxis().SetLabelSize(0.2)
        sigma.Draw("e2")
        sigma.SetStats(0)
        sigma.SetDirectory(0)
        n_next=-1 
        if i+1<N: 
            n_next=gss[i+1][0].GetN()
        if n!=n_next:
            top-=hight+margin
        else:
            top-=hight
            sigma.GetXaxis().SetLabelSize(0)
        c.hists+=[sigma]
        for g in gss[i]:
            g.Draw("p same")
    c.gss=gss
    latex=ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.03)
    latex.SetTextFont(42)
    c.cd()
    latex.DrawLatex(0.8,0.05,"Bin index")
    return c

def SaveCombineSchemeAll(inputpath,outputpath,combine="blue",dist2afb=False):
    histnames=["unfoldedafb_afbm","unfoldedafb_afby_m0","unfoldedafb_afby_m1","unfoldedafb_afby_m2","unfoldedafb_afby_m3","unfoldedafb_afbpt_m0","unfoldedafb_afbpt_m1","unfoldedafb_afbpt_m2","unfoldedafb_afbpt_m3"]
    for region in ["0bjet","nbjet"]:
        mss=[]
        for histname in histnames:
        #for histname in ["unfoldedafb_afbm"]:
            ms=[]
            for channel in ["ee","mm"]:
                for era in ["2016a","2016b","2017","2018"]:
                    ms+=[AFBMeasurements(inputpath,channel+era+"/"+region+"/"+histname,combine=combine)]
            m0=reduce(lambda x,y:x%y,ms)
            
            ms=[]
            for channel in ["ll"]:
                for era in ["2016a","2016b","2017","2018"]:
                    ms+=[AFBMeasurements(inputpath,channel+era+"/"+region+"/"+histname,combine=combine)]
            m1=reduce(lambda x,y:x%y,ms)

            ms=[]
            for channel in ["ee","mm"]:
                for era in ["Run2"]:
                    ms+=[AFBMeasurements(inputpath,channel+era+"/"+region+"/"+histname,combine=combine)]
            m2=reduce(lambda x,y:x%y,ms)
    
            m3=AFBMeasurements(inputpath,"llRun2/"+region+"/"+histname,combine=combine)

            ms=[m0,m1,m2,m3]
            if dist2afb:
                for m in ms:
                    m.Dist2AFB()
            mss+=[ms]
        c=Compare(mss,base=0)
        ROOT.Plotter.DrawPreliminary(c,"Run2")
        latex=ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextFont(42)
        for i in range(len(mss)):
            c.cd(i+1)
            mass_range=""
            if "_m0" in histnames[i]:
                mass_range=";52#leq m_{ll}<77 GeV"
            if "_m1" in histnames[i]:
                mass_range=";77#leq m_{ll}<106 GeV"
            if "_m2" in histnames[i]:
                mass_range=";106#leq m_{ll}<200 GeV"
            if "_m3" in histnames[i]:
                mass_range=";200#leq m_{ll}<3000 GeV"

            variable="m_{ll}"
            if "afby" in histnames[i]:
                variable="|y_{ll}|"
            if "afbpt" in histnames[i]:
                variable="p_{T}^{ll}"                
            latex.SetTextAlign(13)
            latex.SetTextSize(0.12)
            latex.DrawLatex(0.16,0.99,"A_{FB}("+variable+mass_range+")")
        c.cd()
        latex.SetTextSize(0.03)
        latex.SetTextAlign(11)
        latex.SetTextAngle(90)
        latex.DrawLatex(0.05,0.5,"#frac{A_{FB}-A_{FB}^{nominal}}{#sigma(A_{FB}^{nominal})}")
        leg=ROOT.TLegend(0.13,0.03,0.77,0.1)
        leg.SetBorderSize(0)
        leg.SetNColumns(2)
        leg.AddEntry(c.gss[0][0],"era/channel-dependent unfolding (nominal)","l")
        leg.AddEntry(c.gss[0][1],"era-dependent unfolding","l")
        leg.AddEntry(c.gss[0][2],"channel-dependent unfolding","l")
        leg.AddEntry(c.gss[0][3],"unfolding at once","l")
        leg.Draw()
        c.leg=leg
        #raw_input()
        _plotter.plotdir=outputpath
        _plotter.SaveCanvas(c,region+"_combine_scheme.png",False)
        _plotter.SaveCanvas(c,region+"_combine_scheme.pdf")
    return

def SaveCombineMethodAll(inputpath,outputpath,channels=["ee","mm"],eras=["2016a","2016b","2017","2018"]):
    histnames=["unfoldedafb_afbm","unfoldedafb_afby_m0","unfoldedafb_afby_m1","unfoldedafb_afby_m2","unfoldedafb_afby_m3","unfoldedafb_afbpt_m0","unfoldedafb_afbpt_m1","unfoldedafb_afbpt_m2","unfoldedafb_afbpt_m3"]
    for region in ["0bjet","nbjet"]:
        mss=[]
        for histname in histnames:
        #for histname in ["unfoldedafb_afbm"]:
            ms=[]
            for combine,dist2afb in [("simple",False),("blue",False),("blue",True)]:
                m=[]
                for channel in channels:
                    for era in eras:
                        m+=[AFBMeasurements(inputpath,channel+era+"/"+region+"/"+histname,combine=combine)]
                m=reduce(lambda x,y:x%y,m)
                if dist2afb:
                    m.Dist2AFB()
                ms+=[m]
            mss+=[ms]
        c=Compare(mss,base=1)
        ROOT.Plotter.DrawPreliminary(c,"Run2")
        latex=ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextFont(42)
        for i in range(len(mss)):
            c.cd(i+1)
            mass_range=""
            if "_m0" in histnames[i]:
                mass_range=";52#leq m_{ll}<77 GeV"
            if "_m1" in histnames[i]:
                mass_range=";77#leq m_{ll}<106 GeV"
            if "_m2" in histnames[i]:
                mass_range=";106#leq m_{ll}<200 GeV"
            if "_m3" in histnames[i]:
                mass_range=";200#leq m_{ll}<3000 GeV"

            variable="m_{ll}"
            if "afby" in histnames[i]:
                variable="|y_{ll}|"
            if "afbpt" in histnames[i]:
                variable="p_{T}^{ll}"                
            latex.SetTextAlign(13)
            latex.SetTextSize(0.12)
            latex.DrawLatex(0.16,0.99,"A_{FB}("+variable+mass_range+")")
        c.cd()
        latex.SetTextSize(0.03)
        latex.SetTextAlign(11)
        latex.SetTextAngle(90)
        latex.DrawLatex(0.05,0.5,"#frac{A_{FB}-A_{FB}^{nominal}}{#sigma(A_{FB}^{nominal})}")
        leg=ROOT.TLegend(0.13,0.03,0.77,0.1)
        leg.SetBorderSize(0)
        leg.SetNColumns(2)
        leg.AddEntry(c.gss[0][0],"stat. weighted average","l")
        leg.AddEntry(c.gss[0][1],"BLUE method (nominal)","l")
        leg.AddEntry(c.gss[0][2],"stacked","l")
        leg.Draw()
        c.leg=leg
        #raw_input()
        _plotter.plotdir=outputpath
        _plotter.SaveCanvas(c,region+"_combine_method.png",False)
        _plotter.SaveCanvas(c,region+"_combine_method.pdf")
    return

def CompareGen(inputpath,outputpath):
    histnames=["gen_afbm","gen_afby_m0","gen_afby_m1","gen_afby_m2","gen_afby_m3","gen_afbpt_m0","gen_afbpt_m1","gen_afbpt_m2","gen_afbpt_m3"]
    for region in ["0bjet","nbjet"]:
        mss=[]
        for histname in histnames:
            ms=[AFBMeasurements(inputpath,"eeRun2/"+region+"/"+histname),AFBMeasurements(inputpath,"mmRun2/"+region+"/"+histname)]
            mss+=[ms]
        c=Compare(mss)
        ROOT.Plotter.DrawPreliminary(c,"","","nolumi noera")
        latex=ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextFont(42)
        for i in range(len(mss)):
            c.cd(i+1)
            mass_range=""
            if "_m0" in histnames[i]:
                mass_range=";52#leq m_{ll}<77 GeV"
            if "_m1" in histnames[i]:
                mass_range=";77#leq m_{ll}<106 GeV"
            if "_m2" in histnames[i]:
                mass_range=";106#leq m_{ll}<200 GeV"
            if "_m3" in histnames[i]:
                mass_range=";200#leq m_{ll}<3000 GeV"

            variable="m_{ll}"
            if "afby" in histnames[i]:
                variable="|y_{ll}|"
            if "afbpt" in histnames[i]:
                variable="p_{T}^{ll}"                
            latex.SetTextAlign(13)
            latex.SetTextSize(0.12)
            latex.DrawLatex(0.16,0.99,"A_{FB}("+variable+mass_range+")")
        c.cd()
        latex.SetTextSize(0.03)
        latex.SetTextAlign(11)
        latex.SetTextAngle(90)
        latex.DrawLatex(0.05,0.5,"#frac{A_{FB}-A_{FB}^{nominal}}{#sigma(A_{FB}^{nominal})}")
        leg=ROOT.TLegend(0.13,0.03,0.77,0.1)
        leg.SetBorderSize(0)
        leg.SetNColumns(2)
        leg.AddEntry(c.gss[0][0],"gen ee","l")
        leg.AddEntry(c.gss[0][1],"gen #mu#mu","l")
        leg.Draw()
        c.leg=leg
        #raw_input()
        _plotter.plotdir=outputpath
        _plotter.SaveCanvas(c,region+"_compare_gen.png",False)
        _plotter.SaveCanvas(c,region+"_compare_gen.pdf")
    return

if __name__=="__main__":
    # a=AFBMeasurements("AFBResult/final.root","mm2018/0bjet/unfoldedafb_afbm")
    # print a
    # a.Dist2AFB()
    # print a

    #SaveCompareEraAll("AFBResult/final.root","fig/AFBMeasurements/diff")
    #SaveUnfoldedPlotAll("AFBResult/final.root","fig/AFBMeasurements")
    # SaveCombineSchemeAll("AFBResult/final.root","test_simple0",combine="simple",dist2afb=False)
    # SaveCombineSchemeAll("AFBResult/final.root","test_simple1",combine="simple",dist2afb=True)
    # SaveCombineSchemeAll("AFBResult/final.root","test_blue0",combine="blue",dist2afb=False)
    # SaveCombineSchemeAll("AFBResult/final.root","test_blue1",combine="blue",dist2afb=True)
    #SaveCombineMethodAll("AFBResult/final.root","test_method")
    #SaveCombineMethodAll("AFBResult/final.root","test_method_eradep",eras=["2016a","2016b","2017","2018"])
    CompareGen("AFBResult/final.root","fig/AFBMeasurements")
    #SaveCorrelationAll("AFBResult/final.root","fig/AFBMeasurements")
    # syst="triggerSF_mode1 triggerSF_interpolation".split()
    # m=AFBMeasurements("AFBResult/final.root","ll2018/0bjet/unfoldedafb_afby_m0")
    # print m.__str__(syst)
    # m=AFBMeasurements("AFBResult/final.root","mmRun2/nbjet/unfoldedafb_afby_m1")
    # print m.__str__(syst)
    pass

    # hists=ROOT.vector("Hists")()
    # h2017=m2017.GetHists()
    # h2018=m.GetHists()
    # hists.push_back(h2017)
    # hists.push_back(h2018)
    # p=ROOT.Plot()
    # p.hists=hists
    # print hists.size()
    # print hists.at(0).size()
    # for i in range(hists.size()):
    #     for j in range(hists.at(i).size()):
    #         print hists.at(i).at(j).GetName()
    # _plotter.DrawCompare(p)
    # raw_input()
    # print m
