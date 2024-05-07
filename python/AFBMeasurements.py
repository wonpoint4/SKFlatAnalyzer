import os,array,copy
import numpy as np
import ROOT
ROOT.TH1.AddDirectory(0)
ROOT.TH1.SetDefaultSumw2(1)
ROOT.gROOT.ProcessLine('#include"AFBSystPlotter.cc"')

_plotter=ROOT.AFBSystPlotter("mi ttll")
SystematicSuffixes=dict(_plotter.GetSystematicSuffixes("totalsys")).keys()
Systematics=dict(_plotter.systematics)

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
    def __init__(self,filename,histname):
        self.default_syst="totalsys"
        #self.default_syst="z0weight"
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
                    ws/=sum(ws)
            rt.measurements[i]=self.measurements[i]*ws[0]+other.measurements[i]*ws[1]
            cov_stat0[i,:]*=ws[0]
            cov_stat0[:,i]*=ws[0]
            cov_stat1[i,:]*=ws[1]
            cov_stat1[:,i]*=ws[1]
        rt.cov_stat=cov_stat0+cov_stat1
        rt.syst={}
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
                #for era in ["2016a","2016b","2017","2018"]:
                for era in ["Run2"]:
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
                #for era in ["2016a","2016b","2017","2018"]:
                for era in ["Run2"]:
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
                for era in ["2016a","2016b","2017","2018"]:
                    ms+=[AFBMeasurements(inputpath,channel+era+"/"+region+"/"+histname)]
            m=reduce(lambda x,y:x%y,ms)
            m.Save(outputpath+"/"+region+"/"+histname+".tex")

def SaveTableCovAll(inputpath,outputpath):
    for region in ["0bjet","nbjet"]:
        for histname in ["unfoldedafb_afbm","unfoldedafb_afby_m0","unfoldedafb_afby_m1","unfoldedafb_afby_m2","unfoldedafb_afby_m3","unfoldedafb_afbpt_m0","unfoldedafb_afbpt_m1","unfoldedafb_afbpt_m2","unfoldedafb_afbpt_m3"]:
            ms=[]
            for channel in ["ee","mm"]:
                for era in ["2016a","2016b","2017","2018"]:
                    ms+=[AFBMeasurements(inputpath,channel+era+"/"+region+"/"+histname)]
            m=reduce(lambda x,y:x%y,ms)
            m.SaveCov(outputpath+"/"+region+"/"+histname+"_cov.tex")

def SaveCombineScheme(inputpath,outputpath):
    region="nbjet"
    gss=[]
    for histname in ["unfoldedafb_afbm","unfoldedafb_afby_m0","unfoldedafb_afby_m1","unfoldedafb_afby_m2","unfoldedafb_afby_m3","unfoldedafb_afbpt_m0","unfoldedafb_afbpt_m1","unfoldedafb_afbpt_m2","unfoldedafb_afbpt_m3"]:
    #for histname in ["unfoldedafb_afbm"]:
        ms=[]
        for channel in ["ee","mm"]:
            for era in ["2016a","2016b","2017","2018"]:
                ms+=[AFBMeasurements(inputpath,channel+era+"/"+region+"/"+histname)]
        m0=reduce(lambda x,y:x%y,ms)

        ms=[]
        for channel in ["ee","mm"]:
            for era in ["Run2"]:
                ms+=[AFBMeasurements(inputpath,channel+era+"/"+region+"/"+histname)]
        m1=reduce(lambda x,y:x%y,ms)

        ms=[]
        for channel in ["ll"]:
            for era in ["2016a","2016b","2017","2018"]:
                ms+=[AFBMeasurements(inputpath,channel+era+"/"+region+"/"+histname)]
        m2=reduce(lambda x,y:x%y,ms)
    
        m3=AFBMeasurements(inputpath,"llRun2/"+region+"/"+histname)

        gs=[ROOT.TGraphErrors(),ROOT.TGraphErrors(),ROOT.TGraphErrors(),ROOT.TGraphErrors(),ROOT.TGraphErrors(),ROOT.TGraphErrors(),ROOT.TGraphErrors(),ROOT.TGraphErrors()]
        gs[0].SetLineColor(2)
        gs[1].SetLineColor(1)
        gs[2].SetLineColor(4)
        gs[3].SetLineColor(6)
        gs[4+0].SetLineColor(2)
        gs[4+1].SetLineColor(1)
        gs[4+2].SetLineColor(4)
        gs[4+3].SetLineColor(6)
        for i in range(1,len(m0)-1):
            vals=np.array([m0.GetValue(i),m1.GetValue(i),m2.GetValue(i),m3.GetValue(i)])
            errs=np.array([m0.GetTotalError(i),m1.GetTotalError(i),m2.GetTotalError(i),m3.GetTotalError(i)])
            errs_stat=np.array([m0.GetStatError(i),m1.GetStatError(i),m2.GetStatError(i),m3.GetStatError(i)])
            if errs[1]==0: continue
            vals,errs,errs_stat=(vals-vals[1])/errs[1],errs/errs[1],errs_stat/errs[1]
            for j in range(4):
                gs[j].SetPoint(i-1,i-0.35+0.1*j,vals[j])
                gs[j].SetPointError(i-1,0,errs[j])
                gs[4+j].SetPoint(i-1,i-0.35+0.1*j,vals[j])
                gs[4+j].SetPointError(i-1,0,errs_stat[j])
        gss+=[gs]
    c=ROOT.TCanvas()
    c.hists=[]
    c.Divide(1,9)
    for i in range(len(gss)):
        c.cd(i+1)
        ROOT.gPad.SetTopMargin(0)
        ROOT.gPad.SetBottomMargin(0)
        n=gss[i][0].GetN()
        sigma=ROOT.TH1D("sigma","sigma",n,0,n)
        for ib in range(sigma.GetNcells()):
            sigma.SetBinError(ib,1)
        sigma.SetFillStyle(3001)
        sigma.SetFillColor(ROOT.kGreen)
        sigma.GetYaxis().SetTickLength(0)
        sigma.GetYaxis().SetRangeUser(-1.99,1.99)
        sigma.Draw("e2")
        sigma.SetStats(0)
        sigma.SetDirectory(0)
        c.hists+=[sigma]
        for g in gss[i]:
            g.Draw("p same")
    raw_input()

if __name__=="__main__":
    #SaveCompareEraAll("AFBResult/final.root","fig/AFBMeasurements/diff")
    #SaveUnfoldedPlotAll("AFBResult/final.root","fig/AFBMeasurements")
    SaveCombineScheme("AFBResult/final.root","fig/AFBMeasurements")
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
