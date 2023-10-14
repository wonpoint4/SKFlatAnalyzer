import os,array,copy
import numpy as np
import ROOT
ROOT.TH1.AddDirectory(0)
ROOT.TH1.SetDefaultSumw2(1)
ROOT.gROOT.ProcessLine('#include"AFBPlotter.cc"')

_plotter=ROOT.AFBPlotter("mi ttll")
SystematicSuffixes=dict(_plotter.GetSystematicSuffixes("totalsys")).keys()
Systematics=dict(_plotter.systematics)

def Variation2Suffix(variation):
    variation=str(variation)
    words=variation.split(":")
    Type=words[0]
    tag=words[2]
    if Type=="replace":
        return "_"+tag+words[1].split("->")[1]
    elif Type=="scale":
        return "_"+tag+"_scale"
    print "[Error] [Variation2Suffix] Cannot convert "+variation
    return None

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
        for suffix in SystematicSuffixes:
            h=f.Get(histname+suffix)
            for i in range(ncells):
                self.measurements[i].SetSystError(suffix,h.GetBinContent(i),isValue=True)
        self.syst={}
        return

    def __str__(self,syst=[]):
        if len(syst)==0:
            syst=[self.default_syst]
        out=[]
        out+=[ "\t".join( ["bin","val","stat"]+syst ) ]
        for i in range(len(self.measurements)):
            out+=[ "\t".join( map(str,[i,self.GetValue(i),self.GetStatError(i)]+[self.GetSystError(key,i) for key in syst]) ) ]
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
        else:
            print "[Error] Unknown systematic type",syst.type,syst.title
            exit(1)
        return rt

    def GetSystError(self,key,i):
        if key not in self.syst:
            self.EvalSystError(key)
        return self.syst[key][i][i]**0.5

    def GetValue(self,i):
        return self.measurements[i].GetValue()

    def GetStatError(self,i):
        return self.cov_stat[i][i]**0.5

    def GetTotalError(self,i):
        return (self.GetStatError(i)**2+self.GetSystError(self.default_syst,i)**2)**0.5

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

    def Save(self,filename):
        out=["\\begin{tabular} { "+" ".join(["r" for i in range(1,len(self.measurements)-1+4)])+" }"]
        for i in range(1,len(self.measurements)-1):
            line="{} & {} & {} & {} & ".format(i,self.bins[i-1],self.bins[i],round(self.GetValue(i),4)) \
            +" & ".join([str(round(self.GetCov(i,j),4)) for j in range(1,len(self.measurements)-1)]) \
            +" \\\\"
            out+=[line]
        out+=["\\end{tabular}"]
        # out="\n".join(map(lambda x:x.replace("\\","\\\\"),out))
        out="\n".join(out)
        print out
        
def CompareEraAll():
    colors=[ROOT.kBlack,ROOT.kRed,ROOT.kGreen+1,ROOT.kBlue,ROOT.kYellow+1,ROOT.kMagenta,ROOT.kCyan,ROOT.kGray,ROOT.kPink+1,ROOT.kSpring+1,ROOT.kAzure+1,ROOT.kOrange+1,ROOT.kViolet+1,ROOT.kTeal+1,ROOT.kWhite,ROOT.kGray+1,ROOT.kGray+3]
    _plotter.plotdir="fig/AFBMeasurements/"
    for channel in ["mm","ee"]:
        for region in ["0bjet","nbjet"]:
            for histname in ["unfoldedafb_afbm","unfoldedafb_afby_m0","unfoldedafb_afby_m1","unfoldedafb_afby_m2","unfoldedafb_afby_m3","unfoldedafb_afbpt_m0","unfoldedafb_afbpt_m1","unfoldedafb_afbpt_m2","unfoldedafb_afbpt_m3"]:
                ms=[AFBMeasurements("AFBResult/final.root",channel+era+"/"+region+"/"+histname) for era in ["2016a","2016b","2017","2018"]]
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

if __name__=="__main__":
    # inputpath=os.environ["SKFlatOutputDir"]+os.environ["SKFlatV"]+"/AFBAnalyzer/"

    # m=AFBMeasurements("AFBResult/final.root","mm2018/nbjet/unfoldedafb_afbm")
    # m.Save("")
    # exit()

    CompareEraAll()

    # m2017=AFBMeasurements("AFBResult/final.root","mm2017/nbjet/unfoldedafb_afbm")    
    # print m.Compare(m2017)
    _plotter.plotdir="fig/AFBMeasurements/"
    for region in ["0bjet","nbjet"]:
        for histname in ["unfoldedafb_afbm","unfoldedafb_afby_m0","unfoldedafb_afby_m1","unfoldedafb_afby_m2","unfoldedafb_afby_m3","unfoldedafb_afbpt_m0","unfoldedafb_afbpt_m1","unfoldedafb_afbpt_m2","unfoldedafb_afbpt_m3"]:
            print region,histname
            ms=[]
            for channel in ["ee","mm"]:
                for era in ["2016a","2016b","2017","2018"]:
                    ms+=[AFBMeasurements("AFBResult/final.root",channel+era+"/"+region+"/"+histname)]
            m=reduce(lambda x,y:x%y,ms)
            #m=ms[-1]
            hdata=m.GetHists()
            str_bins=[str(m.GetHists().at(0).GetBinLowEdge(i)) for i in range(1,m.GetHists().at(0).GetNbinsX()+2)]
            str_bins="{"+",".join(str_bins)+"}"
            

            sim_histname=""
            sim_histoption="AFB "
            plot_option="ytitle:'A_{FB}'"
            if region=="0bjet":
                sim_index=0
                afbtype="CS"
                #sim_histoption+=" sysname:dytheory"
                plot_option+=" ymin:-0.09 ymax:0.49"
            else:
                sim_index=1
                afbtype="Recoil"
                #sim_histoption+=" sysname:tttheory" 
                plot_option+=" ymin:0.0 ymax:1.49"
            if "afbm" in histname:
                sim_histname="dimass"
                sim_histoption+=" project:x"
                sim_histoption+=" rebinX:"+str_bins
                plot_option+=" xtitle:m(ll) logx"
            elif "afby" in histname:
                sim_histname="dirap"
                sim_histoption+=" project:y"
                sim_histoption+=" rebinY:"+str_bins
                plot_option+=" xtitle:y(ll)"
            elif "afbpt" in histname:
                sim_histname="dipt"
                sim_histoption+=" project:z"
                sim_histoption+=" rebinZ:"+str_bins
                plot_option+=" xtitle:p_{T}(ll) logx"
            sim_histname="[em][em]201[678][ab]?/"+region+"/genfid_"+sim_histname+afbtype+"_dressed"
            #sim_histname="mm2018/"+region+"/genfid_"+sim_histname+afbtype+"_dressed"
            if "_m0" in histname:
                sim_histoption+=" Xmin:52 Xmax:77"
            elif "_m1" in histname:
                sim_histoption+=" Xmin:77 Xmax:106"
            elif "_m2" in histname:
                sim_histoption+=" Xmin:106 Xmax:280"
            elif "_m3" in histname:
                sim_histoption+=" Xmin:280 Xmax:3000"
            hsim=_plotter.GetHistSys(sim_index,sim_histname,sim_histoption)
            for h in hsim:
                h.SetOption("hist e1")

            for h in hdata:
                h.SetName("RunII (blind)")
                h.SetOption("e1")
                h.SetMarkerStyle(20)
                h.SetMarkerSize(0.7)
                h.SetLineColor(1)
                h.SetMarkerColor(1)
                for i in range(h.GetNcells()):
                    h.SetBinContent(i,hsim[0].GetBinContent(i))
            print sim_histname,sim_histoption
            p=ROOT.Plot()
            p.SetOption(plot_option)
            p.hists=ROOT.vector("Hists")([hdata,hsim])
            c=ROOT.TCanvas()
            _plotter.DrawCompare(p)
            _plotter.SaveCanvas(c,region+"_"+histname+".png")
            #raw_input()


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

 
