import math
import numpy as np
import ROOT
ROOT.gROOT.ProcessLine('#include"BBPlotter.cc"')
ROOT.Plotter.SetupStyle()

p=ROOT.BBPlotter("eff")

def calc(a0,a1,a2,a3):
    N=a0+a1+a2+a3
    bm=(2*a0+a1+a2)/N
    bp=(a1+a2+2*a3)/N
    return (bm+(bm**2-4*a0/N)**0.5)/2,(bp+(bp**2-4*a3/N)**0.5)/2

def calcWithCov(a0,e0,a1,e1,a2,e2,a3,e3):
    cov=np.zeros((2,2))
    nominal=np.array(calc(a0,a1,a2,a3))
    stat0=np.array(calc(a0+e0,a1,a2,a3))
    stat1=np.array(calc(a0,a1+e1,a2,a3))
    stat2=np.array(calc(a0,a1,a2+e2,a3))
    stat3=np.array(calc(a0,a1,a2,a3+e3))
    cov+=np.outer(stat0-nominal,stat0-nominal)
    cov+=np.outer(stat1-nominal,stat1-nominal)
    cov+=np.outer(stat2-nominal,stat2-nominal)
    cov+=np.outer(stat3-nominal,stat3-nominal)
    return nominal,cov

def GetAccuracy(ientry,channel,option=""):
    h=p.GetHist(ientry,channel+"/charge",option)
    return calcWithCov(h.GetBinContent(1),h.GetBinError(1),h.GetBinContent(2),h.GetBinError(2),h.GetBinContent(3),h.GetBinError(3),h.GetBinContent(4),h.GetBinError(4))

def calc_eff(valp,valf,errp=None,errf=None):
    if valp+valf==0: return 0.,0.
    eff=valp/(valp+valf)
    if errp is None and errf is None:
        return eff
    err = 1/(valp+valf)**2*math.sqrt(errp*errp*valf*valf+errf*errf*valp*valp)
    return eff,err

def GetTrueAccuracy(channel,option=""):
    hm=p.GetHist(1,channel+"/b[01]mcorrect",option)
    am,em=calc_eff(hm.GetBinContent(2),hm.GetBinContent(1),hm.GetBinError(2),hm.GetBinError(1))
    hp=p.GetHist(1,channel+"/b[01]pcorrect",option)
    ap,ep=calc_eff(hp.GetBinContent(2),hp.GetBinContent(1),hp.GetBinError(2),hp.GetBinError(1))
    return (am,ap),(em,ep)
    
def DrawAccuracy(channels):
    c=ROOT.gROOT.MakeDefCanvas()
    c.SetLeftMargin(0.25)

    gdata=ROOT.TGraphErrors()
    gsim=ROOT.TGraphErrors()
    gtrue=ROOT.TGraphErrors()
    for i in range(len(channels)):
        channel=channels[i]
        g=gdata
        value,cov=GetAccuracy(0,channel)
        g.SetPoint(2*i,value[0],2*i+0.5)
        g.SetPointError(2*i,cov[0][0]**0.5,0)
        g.SetPoint(2*i+1,value[1],2*i+1+0.5)
        g.SetPointError(2*i+1,cov[1][1]**0.5,0)
    
        g=gsim
        value,cov=GetAccuracy(1,channel)
        g.SetPoint(2*i,value[0],2*i+0.5+0.1)
        g.SetPointError(2*i,cov[0][0]**0.5,0)
        g.SetPoint(2*i+1,value[1],2*i+1+0.5+0.1)
        g.SetPointError(2*i+1,cov[1][1]**0.5,0)

        value,cov=GetAccuracy(1,channel,"suffix:_FSR_up:ttll")
        g.SetPoint(2*i+2*len(channel),value[0],2*i+0.5+0.12)
        g.SetPointError(2*i+2*len(channel),cov[0][0]**0.5,0)
        g.SetPoint(2*i+1+2*len(channel),value[1],2*i+1+0.5+0.12)
        g.SetPointError(2*i+1+2*len(channel),cov[1][1]**0.5,0)

        value,cov=GetAccuracy(1,channel,"suffix:_FSR_down:ttll")
        g.SetPoint(2*i+4*len(channel),value[0],2*i+0.5+0.14)
        g.SetPointError(2*i+4*len(channel),cov[0][0]**0.5,0)
        g.SetPoint(2*i+1+4*len(channel),value[1],2*i+1+0.5+0.14)
        g.SetPointError(2*i+1+4*len(channel),cov[1][1]**0.5,0)
    
        g=gtrue
        value,error=GetTrueAccuracy(channel)
        g.SetPoint(2*i,value[0],2*i+0.5+0.2)
        g.SetPointError(2*i,error[0],0)
        g.SetPoint(2*i+1,value[1],2*i+1+0.5+0.2)
        g.SetPointError(2*i+1,error[1],0)

        value,error=GetTrueAccuracy(channel,"suffix:_FSR_up:ttll")
        g.SetPoint(2*i+2*len(channel),value[0],2*i+0.5+0.22)
        g.SetPointError(2*i+2*len(channel),cov[0][0]**0.5,0)
        g.SetPoint(2*i+1+2*len(channel),value[1],2*i+1+0.5+0.22)
        g.SetPointError(2*i+1+2*len(channel),cov[1][1]**0.5,0)

        value,error=GetTrueAccuracy(channel,"suffix:_FSR_down:ttll")
        g.SetPoint(2*i+4*len(channel),value[0],2*i+0.5+0.24)
        g.SetPointError(2*i+4*len(channel),cov[0][0]**0.5,0)
        g.SetPoint(2*i+1+4*len(channel),value[1],2*i+1+0.5+0.24)
        g.SetPointError(2*i+1+4*len(channel),cov[1][1]**0.5,0)
    
    hframe=ROOT.TH2D("hframe","",100,0.60,0.65,len(channels)*2,0,len(channels)*2);
    for i in range(len(channels)):
        hframe.GetYaxis().SetBinLabel(i*2+1,"#alpha^{-} ("+channels[i]+")")
        hframe.GetYaxis().SetBinLabel(i*2+2,"#alpha^{+} ("+channels[i]+")")
    hframe.SetStats(0)
    hframe.Draw()
    gdata.SetMarkerStyle(20)
    gdata.SetMarkerSize(0.6)
    gdata.SetMarkerColor(1)
    gdata.SetLineColor(1)
    gdata.Draw("same p")
    
    gsim.SetMarkerStyle(20)
    gsim.SetMarkerSize(0.6)
    gsim.SetMarkerColor(2)
    gsim.SetLineColor(2)
    gsim.Draw("same p")
    
    gtrue.SetMarkerStyle(20)
    gtrue.SetMarkerSize(0.6)
    gtrue.SetMarkerColor(4)
    gtrue.SetLineColor(4)
    gtrue.Draw("same p")
    
    raw_input()

DrawAccuracy(["me201[678][ab]?","mm201[678][ab]?","ee201[678][ab]?"])
#DrawAccuracy([channel+era for channel in ["me","mm","ee"] for era in ["2016a","2016b","2017","2018","201[678][ab]?"]])
#DrawAccuracy([channel+era for channel in ["[me][me]"] for era in ["2016a","2016b","2017","2018","201[678][ab]?"]])
