#!/usr/bin/env python

import sys,os,math
import numpy as np
from array import array
import ROOT as rt
import pickle

class AFBResult:
    def __init__(self):
        self.channel=""
        self.year=0
        self.m_range=(0.,0.)
        self.y_range=(0.,0.)
        self.pt_range=(0.,0.)
        self.absy=False
        self.data_val=0.
        self.data_err_stat=0.
        self.sim_val=0.
        self.sim_err_stat=0.
        self.sim_err_syst=0.
        self.sim_syst_raw={}
        self.sim_syst={}
        self.sigma=0
        self.entries=""
        if not hasattr(rt,"AFBPlotter"):
            rt.gInterpreter.ProcessLine(".L "+os.getenv("SKFlat_WD")+"/Plotter/AFBPlotter.cc")
            rt.Verbosity=2
    
    def Init(self,channel,year,m_range,y_range,pt_range,absy=True,entries=""):
        if entries=="": plotter=rt.AFBPlotter()
        else: plotter=rt.AFBPlotter(entries.replace(","," "))
        if channel!="ee" and channel!="mm":
            print("Wrong channel "+str(channel))
            exit(1)
        if year!=2016 and year!=2017 and year!=2018:
            print("Wrong year "+str(year))
            exit(1)

        plot=plotter.MakePlot(channel+str(year)+"/m[{:.0f},{:.0f}]/y[{:.1f},{:.1f}]/pt[{:.0f},{:.0f}]/costhetaCS".format(m_range[0],m_range[1],y_range[0],y_range[1],pt_range[0],pt_range[1]),"project:u")
        if absy: plot+="absy"

        self.data_val,self.data_err_stat=plotter.GetAFBAndError(plotter.GetHist(plotter.entries[0],plot))
        self.sim_val,self.sim_err_stat=plotter.GetAFBAndError(plotter.GetHist(plotter.entries[1],plot))

        ## SYS
        suffixes=["_PUweight_up","_PUweight_down",
                  "_prefireweight_up","_prefireweight_down",
                  "_nozptweight",
                  "_RECOSF_up","_RECOSF_down",
                  "_IDSF_up","_IDSF_down",
                  "_ISOSF_up","_ISOSF_down",
                  "_triggerSF_up","_triggerSF_down"]
        for suffix in suffixes:
            self.sim_syst_raw[suffix]=plotter.GetAFBAndError(plotter.GetHist(plotter.entries[1],plot+("suffix:"+suffix+" varibit:12")))[0]

        self.EvalSyst("PUweight",["_PUweight_up","_PUweight_down"])
        self.EvalSyst("prefireweight",["_prefireweight_up","_prefireweight_down"])
        self.EvalSyst("zptweight",["_nozptweight"])

        self.EvalSyst("RECOSF",["_RECOSF_up","_RECOSF_down"])
        self.EvalSyst("IDSF",["_IDSF_up","_IDSF_down"])
        self.EvalSyst("ISOSF",["_ISOSF_up","_ISOSF_down"])
        self.EvalSyst("triggerSF",["_triggerSF_up","_triggerSF_down"])
        self.CombineSyst("efficiencySF",["RECOSF","IDSF","ISOSF","triggerSF"])

        ## PDFSYS
        suffixes=["_scalevariation0","_scalevariation1","_scalevariation2","_scalevariation3","_scalevariation4","_scalevariation6","_scalevariation8",
                  "_alphaS_up","_alphaS_down"]
        suffixes+=["_pdf"+str(i) for i in range(100)]
        for suffix in suffixes:
            self.sim_syst_raw[suffix]=plotter.GetAFBAndError(plotter.GetHist(plotter.entries[1],plot+("suffix:"+suffix+" varibit:4")))[0]

        self.EvalSyst("scalevariation",["_scalevariation0","_scalevariation1","_scalevariation2","_scalevariation3","_scalevariation4","_scalevariation6","_scalevariation8"])
        self.EvalSyst("alphaS",["_alphaS_up","_alphaS_down"])
        if year==2016: self.EvalSyst("pdf",["_pdf"+str(i) for i in range(100)],"gaussian")
        else: self.EvalSyst("pdf",["_pdf"+str(i) for i in range(100)],"hessian")
        self.CombineSyst("PDF",["pdf","alphaS"])

        ## Xsec SYS
        plotter.Setup(plotter.mode.replace("+tt","+1.06*tt"))
        self.sim_syst_raw["ttXsec"]=plotter.GetAFBAndError(plotter.GetHist(plotter.entries[1],plot))[0]
        self.EvalSyst("ttXsec",["ttXsec"])
        

        self.sim_err_syst=0
        for sys in self.sim_syst.values():
            self.sim_err_syst+=sys**2
        self.sim_err_syst=math.sqrt(self.sim_err_syst)

        if self.data_err_stat or self.sim_err_stat:
            self.sigma=(self.data_val-self.sim_val)/math.sqrt(self.data_err_stat**2+self.sim_err_stat**2+self.sim_err_syst**2)

        self.channel=channel
        self.year=year
        self.m_range=m_range
        self.y_range=y_range
        self.pt_range=pt_range
        self.entries=entries
        self.Print(3)

    def EvalSyst(self,name,suffixes,method="envelope"):
        if method=="envelope":
            max_diff=0
            for suffix in suffixes:
                if abs(self.sim_syst_raw[suffix]-self.sim_val)>max_diff:
                    max_diff=abs(self.sim_syst_raw[suffix]-self.sim_val)
            self.sim_syst[name]=max_diff
        elif method=="gaussian":
            square=0
            for suffix in suffixes:
                square+=pow(self.sim_syst_raw[suffix]-self.sim_val,2)
            ms=square/len(suffixes)
            self.sim_syst[name]=math.sqrt(ms)
        elif method=="hessian":
            square=0
            for suffix in suffixes:
                square+=pow(self.sim_syst_raw[suffix]-self.sim_val,2)
            self.sim_syst[name]=math.sqrt(square)
        else:
            print("unkown method "+method)

    def CombineSyst(self,name,systnames):
        self.sim_syst[name]=math.sqrt(sum([self.sim_syst[x]**2 for x in systnames]))
        for x in systnames: self.sim_syst.pop(x,None)

    def Print(self,verbosity=1):
        if verbosity>0: 
            out='''{} {} m[{:.0f},{:.0f}] y[{:.1f},{:.1f}] pt[{:.0f},{:.0f}] : {:.2f} sigma
'''.format(self.channel,self.year,self.m_range[0],self.m_range[1],self.y_range[0],self.y_range[1],self.pt_range[0],self.pt_range[1],self.sigma)
        if verbosity>1: 
            out+='''  data: {:.5f} +- {:.5f}
  sim:  {:.5f} +- {:.5f} +- {:.5f}
'''.format(self.data_val,self.data_err_stat,self.sim_val,self.sim_err_stat,self.sim_err_syst)
        if verbosity>2:
            for name,val in self.sim_syst.items():
                out+='''    {:<15}: +- {:.5f}
'''.format(name,val)
        if verbosity>3:
            for name,val in self.sim_syst_raw.items():
                out+='''      {:<17}: {:.5f}
'''.format(name,val)
        sys.stdout.write(out)

    def GetName(self,previous=None,order=["channel","year","m","y","pt"]):
        names=[]
        for key in order:
            if key=="channel":
                if not previous or self.channel!=previous.channel: names+=[self.channel]
            elif key=="year":
                if not previous or self.year!=previous.year: names+=[str(self.year)]
            elif key=="m":
                if not previous or self.m_range!=previous.m_range: names+=["m[{:.0f},{:.0f}]".format(self.m_range[0],self.m_range[1])]
            elif key=="y":
                if not previous or self.y_range!=previous.y_range: names+=["y[{:.1f},{:.1f}]".format(self.y_range[0],self.y_range[1])]
            elif key=="pt":
                if not previous or self.pt_range!=previous.pt_range: names+=["pt[{:.0f},{:.0f}]".format(self.pt_range[0],self.pt_range[1])]
            else:
                print("unknown key: "+key)
                exit(1)                
        return " ".join(names)
        
class AFBReport:
    def __init__(self,condor=False):
        self.results=[]
        self.condor=condor
        self.order=["channel","year","m","y","pt"]

    def Copy(self):
        a=AFBReport(self.condor)
        a.results=self.results
        a.order=self.order
        return a

    def Init(self,channels=["ee"],years=[2016],mbins=[60,80],ybins=[0.0,1.2],ptbins=[0,20],entries=""):
        self.results=[]
        if self.condor:
            workingdir=os.popen("mktemp -d -p .").read().strip()
            print("Submitting jobs")
            for channel in channels:
                for year in years:
                    for im in range(len(mbins)-1):
                        for iy in range(len(ybins)-1):
                            for ip in range(len(ptbins)-1):
                                runname="run_{}_{}_m{}_y{}_p{}".format(channel,year,im,iy,ip)
                                arguments="create {}/{}.pkl --channels {} --years {} --mbins {},{} --ybins {},{} --ptbins {},{}".format(workingdir,runname,channel,year,mbins[im],mbins[im+1],ybins[iy],ybins[iy+1],ptbins[ip],ptbins[ip+1])
                                if entries!="": arguments+=" --entries "+entries
                                os.system(
'''condor_submit -batch-name AFBReport &> /dev/null <<EOF
executable = {0}
output = {1}/{2}.out
error = {1}/{2}.err
log = {1}/condor.log
getenv = true
arguments = {3}
queue
EOF
'''.format(os.popen("which AFBReport.py").read().strip(),workingdir,runname,arguments))
                                sys.stdout.write(".")
                                sys.stdout.flush()
            print("Waiting jobs")
            os.system("condor_wait "+workingdir+"/condor.log")

            if os.popen("egrep -rIi 'error|warning' "+workingdir).read()!="":
                print("Error during condor run. Check logs in "+workingdir)
                print(os.popen("egrep -rIi 'error|warning' "+workingdir).read())
                exit(1)

            for channel in channels:
                for year in years:
                    for im in range(len(mbins)-1):
                        for iy in range(len(ybins)-1):
                            for ip in range(len(ptbins)-1):
                                runname="run_{}_{}_m{}_y{}_p{}".format(channel,year,im,iy,ip)
                                this_report=AFBReport()
                                this_report.Load(workingdir+"/"+runname+".pkl")
                                self.Add(this_report)
            os.system("rm -r "+workingdir)
        else:
            for channel in channels:
                for year in years:
                    for im in range(len(mbins)-1):
                        for iy in range(len(ybins)-1):
                            for ip in range(len(ptbins)-1):
                                result=AFBResult()
                                result.Init(channel,year,(mbins[im],mbins[im+1]),(ybins[iy],ybins[iy+1]),(ptbins[ip],ptbins[ip+1]),entries=entries)
                                self.results+=[result]


    def Add(self,b):
        self.results+=b.results

    def Load(self,path):
        with open(path) as f:
            self.results+=pickle.load(f).results

    def Save(self,path):
        with open(path,"w") as f:
            pickle.dump(self,f)

    def Print(self,verbosity=1):
        for result in self.results:
            result.Print(verbosity)

    def PlotVal(self):
        results=self.results
        years=sorted(list(set([x.year for x in results])))
        channels=sorted(list(set([x.channel for x in results])))
        mbins=sorted(list(set([x.m_range[0] for x in results]+[x.m_range[1] for x in results])))
        ybins=sorted(list(set([x.y_range[0] for x in results]+[x.y_range[1] for x in results])))
        ptbins=sorted(list(set([x.pt_range[0] for x in results]+[x.pt_range[1] for x in results])))
        nbins=len(results)
        if len(years)==1 and len(channels)==1:
            numericbin=True
            if len(mbins)>2 and len(ybins)<3 and len(ptbins)<3: bins=mbins
            elif len(ybins)>2 and len(mbins)<3 and len(ptbins)<3: bins=ybins
            elif len(ptbins)>2 and len(mbins)<3 and len(ybins)<3: bins=ptbins
            else: 
                numericbin=False
                bins=range(nbins+1)
        else:
            numericbin=False
            bins=range(nbins+1)

        c=rt.TCanvas()
        c.cd()
        if not numericbin: c.SetBottomMargin(0.4)
        data=rt.TH1D("data","data",nbins,array('d',bins))
        data.SetMarkerStyle(20)
        sim=rt.TH1D("sim","sim",nbins,array('d',bins))
        sim.SetLineColor(2)
        sim.SetFillStyle(3001)
        sim.SetFillColor(2)
        sim_total=rt.TH1D("sim_total","sim_total",nbins,array('d',bins))
        sim_total.SetLineColor(2)
        sim_total.SetFillStyle(3001)
        sim_total.SetFillColor(3)
        for i in range(len(results)):
            result=results[i]
            if math.isnan(result.data_val): continue
            if math.isnan(result.sim_val): continue
            data.SetBinContent(i+1,result.data_val)
            data.SetBinError(i+1,result.data_err_stat)
            sim.SetBinContent(i+1,result.sim_val)
            sim.SetBinError(i+1,result.sim_err_stat)
            sim_total.SetBinContent(i+1,result.sim_val)
            sim_total.SetBinError(i+1,math.sqrt(result.sim_err_stat**2+result.sim_err_syst**2))
            if not numericbin:
                sim_total.GetXaxis().SetBinLabel(i+1,result.GetName(results[i-1],order=self.order))
                sim_total.GetXaxis().LabelsOption("v")

        sim_total.SetStats(0)
        sim_total.Draw("e2")
        sim_total.GetYaxis().SetTitle("val")
        sim.Draw("same e2")
        data.Draw("same e")
        return c, sim_total, sim, data

    def PlotSyst(self):
        results=self.results
        years=sorted(list(set([x.year for x in results])))
        channels=sorted(list(set([x.channel for x in results])))
        mbins=sorted(list(set([x.m_range[0] for x in results]+[x.m_range[1] for x in results])))
        ybins=sorted(list(set([x.y_range[0] for x in results]+[x.y_range[1] for x in results])))
        ptbins=sorted(list(set([x.pt_range[0] for x in results]+[x.pt_range[1] for x in results])))
        nbins=len(results)
        if len(years)==1 and len(channels)==1:
            numericbin=True
            if len(mbins)>2 and len(ybins)<3 and len(ptbins)<3: bins=mbins
            elif len(ybins)>2 and len(mbins)<3 and len(ptbins)<3: bins=ybins
            elif len(ptbins)>2 and len(mbins)<3 and len(ybins)<3: bins=ptbins
            else: 
                numericbin=False
                bins=range(nbins+1)
        else:
            numericbin=False
            bins=range(nbins+1)

        c=rt.TCanvas()
        c.cd()
        c.SetLogy()
        if not numericbin: c.SetBottomMargin(0.4)
        hists={}
        for i,syst in enumerate(results[0].sim_syst.keys()):
            hists[syst]=rt.TH1D(syst,syst,nbins,array('d',bins))
            hists[syst].SetLineColor(i+1)
            if i==4: hists[syst].SetLineColor(rt.kYellow+1)
            hists[syst].SetLineWidth(2)
        for i in range(len(results)):
            result=results[i]
            if math.isnan(result.data_val): continue
            if math.isnan(result.sim_val): continue
            for name,hist in hists.items():
                #hist.SetBinContent(i+1,abs(result.sim_syst[name]/result.sim_err_syst))
                hist.SetBinContent(i+1,abs(result.sim_syst[name]))
                hist.SetBinError(i+1,0)
                if not numericbin:
                    hist.GetXaxis().SetBinLabel(i+1,result.GetName(results[i-1],order=self.order))
                    hist.GetXaxis().LabelsOption("v")

        legend=rt.TLegend(0.1,0.9,0.35,0.65)
        for i,hist in enumerate(hists.values()):
            hist.SetStats(0)
            if i==0: 
                hist.Draw("hist")
                hist.GetYaxis().SetRangeUser(0.00001,0.1)
            else: hist.Draw("hist same")
            hist.GetYaxis().SetTitle("Error")
            legend.AddEntry(hist)
        legend.Draw()
        return [c]+hists.values()+[legend]

    def PlotDiff(self):
        nbins=len(self.results)
        print(nbins)
        c=rt.TCanvas()
        c.cd()
        c.SetBottomMargin(0.4)
        data=rt.TH1D("data","data",nbins,0,nbins)
        data.SetMarkerStyle(20)
        sim=rt.TH1D("sim","sim",nbins,0,nbins)
        sim.SetLineColor(2)
        sim.SetFillStyle(3001)
        sim.SetFillColor(2)
        sim_total=rt.TH1D("sim_total","sim_total",nbins,0,nbins)
        sim_total.SetLineColor(2)
        sim_total.SetFillStyle(3001)
        sim_total.SetFillColor(3)
        for i in range(len(self.results)):
            result=self.results[i]
            if math.isnan(result.data_val): continue
            if math.isnan(result.sim_val): continue
            data.SetBinContent(i+1,result.data_val-result.sim_val)
            data.SetBinError(i+1,result.data_err_stat)
            sim.SetBinContent(i+1,result.sim_val-result.sim_val)
            sim.SetBinError(i+1,result.sim_err_stat)
            sim_total.SetBinContent(i+1,result.sim_val-result.sim_val)
            sim_total.SetBinError(i+1,math.sqrt(result.sim_err_stat**2+result.sim_err_syst**2))
            sim_total.GetXaxis().SetBinLabel(i+1,result.GetName(self.results[i-1],order=self.order))
        sim_total.GetXaxis().LabelsOption("v")
        sim_total.SetStats(0)
        sim_total.Draw("e2")
        sim_total.GetYaxis().SetTitle("diff")
        sim.Draw("same e2")
        data.Draw("same e")
        return c, sim_total, sim, data

    def PlotSig(self):
        nbins=len(self.results)
        print(nbins)
        c=rt.TCanvas()
        c.cd()
        c.SetBottomMargin(0.4)
        data=rt.TH1D("data","data",nbins,0,nbins)
        data.SetMarkerStyle(20)
        sim=rt.TH1D("sim","sim",nbins,0,nbins)
        sim.SetLineColor(2)
        sim.SetFillStyle(3001)
        sim.SetFillColor(2)
        sim_total=rt.TH1D("sim_total","sim_total",nbins,0,nbins)
        sim_total.SetLineColor(2)
        sim_total.SetFillStyle(3001)
        sim_total.SetFillColor(3)
        for i in range(len(self.results)):
            result=self.results[i]
            if math.isnan(result.data_val): continue
            if math.isnan(result.sim_val): continue
            scale=math.sqrt(result.data_err_stat**2+result.sim_err_stat**2+result.sim_err_syst**2)
            if math.isnan(scale): continue
            data.SetBinContent(i+1,(result.data_val-result.sim_val)/scale)
            data.SetBinError(i+1,result.data_err_stat/scale)
            sim.SetBinContent(i+1,(result.sim_val-result.sim_val)/scale)
            sim.SetBinError(i+1,result.sim_err_stat/scale)
            sim_total.SetBinContent(i+1,(result.sim_val-result.sim_val)/scale)
            sim_total.SetBinError(i+1,math.sqrt(result.sim_err_stat**2+result.sim_err_syst**2)/scale)
            sim_total.GetXaxis().SetBinLabel(i+1,result.GetName(self.results[i-1],order=self.order))
        sim_total.GetXaxis().LabelsOption("v")
        sim_total.GetYaxis().SetRangeUser(-5,5)
        sim_total.SetStats(0)
        sim_total.Draw("e2")
        sim_total.GetYaxis().SetTitle("#sigma")
        sim.Draw("same e2")
        data.Draw("same e")
        c.SetGridy()
        return c, sim_total, sim, data

    def Sort(self,keys):
        for key in reversed(keys):
            if key=="channel": 
                self.results=sorted(self.results,key=lambda x: x.channel)
            elif key=="year": 
                self.results=sorted(self.results,key=lambda x: x.year)
            elif key=="m": 
                self.results=sorted(self.results,key=lambda x: x.m_range)
            elif key=="y": 
                self.results=sorted(self.results,key=lambda x: x.y_range)
            elif key=="pt": 
                self.results=sorted(self.results,key=lambda x: x.pt_range)
            else:
                print("unknown key: "+key)
                exit(1)
        old_order=self.order
        self.order=keys
        for key in old_order:
            if key not in self.order: self.order+=[key]

    def Sorted(self,keys):
        a=self.Copy()
        a.Sort(keys)
        return a

    def Select(self,channels=None,years=None,m=None,y=None,pt=None):
        newresults=[]
        for result in self.results:
            if channels is not None:
                if not result.channel in channels: continue
            if years is not None:
                if not result.year in years: continue
            if m is not None:
                if result.m_range[0]<m[0]: continue
                if result.m_range[1]>m[-1]: continue
            if y is not None:
                if result.y_range[0]<y[0]: continue
                if result.y_range[1]>y[-1]: continue
            if pt is not None:
                if result.pt_range[0]<pt[0]: continue
                if result.pt_range[1]>pt[-1]: continue
            newresults+=[result]
        self.results=newresults

    def Selected(self,channels=None,years=None,m=None,y=None,pt=None):
        a=self.Copy()
        a.Select(channels,years,m,y,pt)
        return a

    def Combine(self,keys=[]):
        result_group=[]
        for result in self.results:
            done=False
            for i in range(len(result_group)):
                example=result_group[i][0]
                if "channel" not in keys and example.channel!=result.channel: continue
                if "year" not in keys and example.year!=result.year: continue
                if example.m_range!=result.m_range: continue
                if example.y_range!=result.y_range: continue
                if example.pt_range!=result.pt_range: continue
                result_group[i]+=[result]
                done=True
                break
            if not done:
                result_group+=[[result]]
        new_results=[]
        for results in result_group:
            nresult=len(results)
            unit=np.mat(np.ones(shape=(nresult,1)))
            a=AFBResult()
            
            data_cov=self.GetCovarianceMatrix([x.data_val for x in results],[x.data_err_stat for x in results])
            data_w=data_cov.I*unit/(unit.T*data_cov.I*unit)
            a.data_val=float(np.mat([x.data_val for x in results])*data_w)
            a.data_err_stat=math.sqrt(data_w.T*data_cov*data_w)

            sim_cov=self.GetCovarianceMatrix([x.sim_val for x in results],[x.sim_err_stat for x in results],[x.sim_syst for x in results])
            sim_w=sim_cov.I*unit/(unit.T*sim_cov.I*unit)
            a.sim_val=float(np.mat([x.sim_val for x in results])*sim_w)
            a.sim_err_stat=math.sqrt(sum([sim_w[i]**2*results[i].sim_err_stat**2 for i in range(nresult)]))
            a.sim_err_syst=math.sqrt(sim_w.T*sim_cov*sim_w-a.sim_err_stat**2)
            for key in results[0].sim_syst.keys():
                a.sim_syst[key]=float(np.mat([x.sim_syst[key] for x in results])*sim_w)

            a.channel=",".join(set([x.channel for x in results]))
            a.year=",".join(sorted(list(set([str(x.year) for x in results]))))
            a.m_range=results[0].m_range
            a.y_range=results[0].y_range
            a.pt_range=results[0].pt_range
            new_results+=[a]
            if a.data_err_stat or a.sim_err_stat:
                a.sigma=(a.data_val-a.sim_val)/math.sqrt(a.data_err_stat**2+a.sim_err_stat**2+a.sim_err_syst**2)
        self.results=new_results            

    def Combined(self,keys=[]):
        a=self.Copy()
        a.Combine(keys)
        return a

    def GetCovarianceMatrix(self,vals=[],errs_stat=[],systs=[]):
        nval=len(vals)
        cov=np.mat(np.zeros(shape=(nval,nval)))
        for i in range(nval):
            cov[i,i]+=errs_stat[i]**2
        syst_names=set([])
        for syst in systs:
            syst_names|=set(syst.keys())
        for name in syst_names:
            sigma=np.mat(np.zeros(shape=(nval,1)))
            for i in range(nval):
                syst=systs[i]
                if name in syst:
                    sigma[i]=syst[name]
            cov+=sigma*sigma.T
        return cov

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("action",help="create, print, plot, select, merge")
    parser.add_argument("path",nargs="+",help="filename for action")
    parser.add_argument("--condor",action="store_true",help="[create] Use condor")
    parser.add_argument("--entries",default="",type=str,help="[create] ")
    parser.add_argument("--channels",type=str,help="[all] ex) ee,mm")
    parser.add_argument("--years",type=str,help="[all] ex) 2016,2017,2018")
    parser.add_argument("--mbins",type=str,help="[all] ex) 60,80,100,120,200,500,3000")
    parser.add_argument("--ybins",type=str,help="[all] ex) 0.0,1.2,2.4")
    parser.add_argument("--ptbins",type=str,help="[all] ex) 0,20,50,100,650")
    parser.add_argument("--printlevel",default=1,type=int,help="[print] ")
    parser.add_argument("--sort",type=str,help="[print,plot] ex) channel,year,m,y,pt")
    parser.add_argument("--combine",type=str,help="[print,plot] ex) channel,year,m,y,pt")
    parser.add_argument("--plots",type=str,help="[plot] ")
    args=parser.parse_args()

    if args.channels is not None: args.channels=[x for x in args.channels.split(",")]
    if args.years is not None: args.years=[int(x) for x in args.years.split(",")]
    if args.mbins is not None: args.mbins=[float(x) for x in args.mbins.split(",")]
    if args.ybins is not None: args.ybins=[float(x) for x in args.ybins.split(",")]
    if args.ptbins is not None: args.ptbins=[float(x) for x in args.ptbins.split(",")]
    if args.sort is not None: args.sort=[x for x in args.sort.split(",")]
    if args.combine is not None: args.combine=[x for x in args.combine.split(",")]

    if args.action=="create":
        if len(args.path)!=1:
            print("'create' action needs one path (target) ")
            exit(1)

        if args.channels is None: args.channels=["ee","mm"]
        if args.years is None: args.years=[2016,2017,2018]
        if args.mbins is None: args.mbins=[60,80,100,120,200,500,3000]
        if args.ybins is None: args.ybins=[0.0,1.2,2.4]
        if args.ptbins is None: args.ptbins=[0,20,50,100,650]
        if args.entries=="m": 
            args.entries="data,amcM+tau_amcM+vv+wjets+tt+tw"
        elif args.entries=="pt":
            args.entries="data,amcPt+tau_amcPt+vv+wjets+tt+tw"
        
        a=AFBReport(args.condor)
        a.Init(args.channels,args.years,args.mbins,args.ybins,args.ptbins,entries=args.entries)
        a.Print(2)
        a.Save(args.path[0])

    elif args.action=="print":
        if len(args.path)!=1:
            print("'print' action needs one path (source) ")
            exit(1)

        with open(args.path[0]) as f:
            a=pickle.load(f)

        a.Select(channels=args.channels,years=args.years,m=args.mbins,y=args.ybins,pt=args.ptbins)
        if args.sort is not None: a.Sort(args.sort)
        if args.combine is not None: a.Combine(args.combine)

        a.Print(args.printlevel)


    elif args.action=="plot":
        if len(args.path)<1:
            print("'plot' action needs at least one path (source [target]) ")
            exit(1)

        with open(args.path[0]) as f:
            a=pickle.load(f)
        
        a.Select(channels=args.channels,years=args.years,m=args.mbins,y=args.ybins,pt=args.ptbins)
        if args.sort is not None: a.Sort(args.sort)
        if args.combine is not None: a.Combine(args.combine)

        c0=a.PlotVal()
        c1=a.PlotDiff()
        c2=a.PlotSig()
        c3=a.PlotSyst()
        if len(args.path)>1:
            if len(args.path[1].split(":"))==2:
                filepath=args.path[1].split(":")[0]
                tdirpath=args.path[1].split(":")[1]
            else:
                filepath=args.path[1]
                tdirpath=""
            f=rt.TFile(filepath,"update")
            tdir=f.GetDirectory(tdirpath)
            if tdir==None:
                f.mkdir(tdirpath)
            f.cd(tdirpath)
            for obj in c0: obj.Write()
            for obj in c3[:-1]: obj.Write()
        else:
            b=raw_input("press any key: ")

    elif args.action=="select":
        if len(args.path)!=2:
            print("'select' action needs two paths (source target) ")
            exit(1)

        with open(args.path[0]) as f:
            a=pickle.load(f)
        
        a.Select(channels=args.channels,years=args.years,m=args.mbins,y=args.ybins,pt=args.ptbins)
        if args.sort is not None: a.Sort(args.sort)

        a.Print(1)
        a.Save(args.path[1])

    elif args.action=="merge":
        if len(args.path)<2:
            print("'merge' action needs at least two paths (sources... target) ")
            exit(1)

        with open(args.path[0]) as f:
            a=pickle.load(f)

        for path in args.path[1:-1]:
            with open(path) as f:
                b=pickle.load(f)
                a.results+=b.results
        
        a.Select(channels=args.channels,years=args.years,m=args.mbins,y=args.ybins,pt=args.ptbins)
        if args.sort is not None: a.Sort(args.sort)

        a.Print(1)
        a.Save(args.path[-1])

    elif args.action=="smpv":
        with open("results.pkl") as f:
            a=pickle.load(f)
        a.Sort(["m","pt","y"])
        a.Select(m=[60,120])
        mbins=sorted(list(set([x.m_range[0] for x in a.results]+[x.m_range[1] for x in a.results])))
        ybins=sorted(list(set([x.y_range[0] for x in a.results]+[x.y_range[1] for x in a.results])))
        ptbins=sorted(list(set([x.pt_range[0] for x in a.results]+[x.pt_range[1] for x in a.results])))
        indiv=[]
        comb=[]
        for im in range(len(mbins)-1):
            for iy in range(len(ybins)-1):
                indiv+=[a.Selected(m=[mbins[im],mbins[im+1]],y=[ybins[iy],ybins[iy+1]]).PlotVal()]
                comb+=[a.Combined(["year","channel"]).Selected(m=[mbins[im],mbins[im+1]],y=[ybins[iy],ybins[iy+1]]).PlotVal()]
                comb[-1][1].GetXaxis().GetLabels().Delete()
                comb[-1][1].SetBins(4,np.array(ptbins))
                comb[-1][2].SetBins(4,np.array(ptbins))
                comb[-1][3].SetBins(4,np.array(ptbins))
                
        b=raw_input("press any key: ")

    elif args.action=="test":
        with open(args.path[0]) as f:
            a=pickle.load(f)
        
        a.Select(channels=args.channels,years=args.years,m=args.mbins,y=args.ybins,pt=args.ptbins)
        if args.sort is not None: a.Sort(args.sort)
        
        a.Print(3)
        for result in a.results:
            result.CombineSyst("efficiency SF",["RECOSF","IDSF","ISOSF","triggerSF"])
        a.Print(3)

    else:
        print("unknown action")
        exit(1)
