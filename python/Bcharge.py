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
    print(nominal)
    print("cov",cov)
    values,vectors=np.linalg.eig(cov)
    vector0=np.array([vectors[0][0],vectors[1][0]])*values[0]**0.5
    vector1=np.array([vectors[0][1],vectors[1][1]])*values[1]**0.5
    print(vector0)
    print(vector1)
    print(np.outer(vector0,vector0)+np.outer(vector1,vector1))
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

def GetTrueAccuracyByType(channel,type,option=""):
    hm=p.GetHist(1,channel+"/b[01]mcorrect_type{}".format(type),option)
    #hm=p.GetHist(1,channel+"/b1mcorrect_type{}".format(type),option)
    am,em=calc_eff(hm.GetBinContent(2),hm.GetBinContent(1),hm.GetBinError(2),hm.GetBinError(1))
    hp=p.GetHist(1,channel+"/b[01]pcorrect_type{}".format(type),option)
    #hp=p.GetHist(1,channel+"/b1pcorrect_type{}".format(type),option)
    ap,ep=calc_eff(hp.GetBinContent(2),hp.GetBinContent(1),hp.GetBinError(2),hp.GetBinError(1))
    return (am,ap),(em,ep)
    
def DrawAccuracy(channels):
    channels.reverse()
    c=ROOT.gROOT.MakeDefCanvas()
    c.SetLeftMargin(0.2)

    gdata=ROOT.TGraphErrors()
    gsim=ROOT.TGraphErrors()
    gtrue=ROOT.TGraphErrors()
    gdata=[ROOT.TGraphErrors(),ROOT.TGraphErrors()]
    gsim=[ROOT.TGraphErrors(),ROOT.TGraphErrors()]
    gtrue=[ROOT.TGraphErrors(),ROOT.TGraphErrors()]
    for i in range(len(channels)):
        channel=channels[i]
        g=gdata
        value,cov=GetAccuracy(0,channel)
        for j in range(2):
            g[j].SetPoint(i,value[j],2*i+j+0.5)
            g[j].SetPointError(i,cov[j][j]**0.5,0)
    
        g=gsim
        value,cov=GetAccuracy(1,channel)
        for j in range(2):
            g[j].SetPoint(i,value[j],2*i+j+0.4)
            g[j].SetPointError(i,cov[j][j]**0.5,0)

        # value,cov=GetAccuracy(1,channel,"suffix:_FSR_up:ttll")
        # g.SetPoint(2*i+2*len(channel),value[0],2*i+0.5+0.12)
        # g.SetPointError(2*i+2*len(channel),cov[0][0]**0.5,0)
        # g.SetPoint(2*i+1+2*len(channel),value[1],2*i+1+0.5+0.12)
        # g.SetPointError(2*i+1+2*len(channel),cov[1][1]**0.5,0)

        # value,cov=GetAccuracy(1,channel,"suffix:_FSR_down:ttll")
        # g.SetPoint(2*i+4*len(channel),value[0],2*i+0.5+0.14)
        # g.SetPointError(2*i+4*len(channel),cov[0][0]**0.5,0)
        # g.SetPoint(2*i+1+4*len(channel),value[1],2*i+1+0.5+0.14)
        # g.SetPointError(2*i+1+4*len(channel),cov[1][1]**0.5,0)
    
        g=gtrue
        value,error=GetTrueAccuracy(channel)
        for j in range(2):
            g[j].SetPoint(i,value[j],2*i+j+0.3)
            g[j].SetPointError(i,cov[j][j]**0.5,0)

        # value,error=GetTrueAccuracy(channel,"suffix:_FSR_up:ttll")
        # g.SetPoint(2*i+2*len(channel),value[0],2*i+0.5+0.22)
        # g.SetPointError(2*i+2*len(channel),cov[0][0]**0.5,0)
        # g.SetPoint(2*i+1+2*len(channel),value[1],2*i+1+0.5+0.22)
        # g.SetPointError(2*i+1+2*len(channel),cov[1][1]**0.5,0)

        # value,error=GetTrueAccuracy(channel,"suffix:_FSR_down:ttll")
        # g.SetPoint(2*i+4*len(channel),value[0],2*i+0.5+0.24)
        # g.SetPointError(2*i+4*len(channel),cov[0][0]**0.5,0)
        # g.SetPoint(2*i+1+4*len(channel),value[1],2*i+1+0.5+0.24)
        # g.SetPointError(2*i+1+4*len(channel),cov[1][1]**0.5,0)
    
    hframe=ROOT.TH2D("hframe","",100,0.60,0.65,len(channels)*2+1,0,len(channels)*2+1);
    for i in range(len(channels)):
        title=channels[i]
        title=title.replace("201[678][ab]?"," Run2")
        title=title.replace("m","#mu")
        hframe.GetYaxis().SetBinLabel(i*2+1,"#alpha^{#minus} ("+title+")")
        hframe.GetYaxis().SetBinLabel(i*2+2,"#alpha^{#plus} ("+title+")")
    hframe.SetStats(0)
    hframe.Draw()
    hframe.GetXaxis().SetTitle("Accuracy")
    hframe.GetYaxis().SetLabelSize(hframe.GetYaxis().GetLabelSize()*1.5)

    leg=ROOT.TLegend(c.GetLeftMargin()+0.01,0.79,0.99-c.GetRightMargin(),0.89)
    leg.AddEntry(gdata[0],"data")
    leg.AddEntry(gsim[0],"simulation")
    leg.AddEntry(gtrue[0],"simulation (gen info)")
    leg.SetBorderSize(0)
    leg.Draw()

    gdata[0].SetMarkerStyle(20)
    gdata[0].SetMarkerSize(0.8)
    gdata[0].SetMarkerColor(1)
    gdata[0].SetLineColor(1)
    gdata[0].Draw("same p")

    gdata[1].SetMarkerStyle(24)
    gdata[1].SetMarkerSize(0.8)
    gdata[1].SetMarkerColor(1)
    gdata[1].SetLineColor(1)
    gdata[1].Draw("same p")
    
    gsim[0].SetMarkerStyle(20)
    gsim[0].SetMarkerSize(0.6)
    gsim[0].SetMarkerColor(2)
    gsim[0].SetLineColor(2)
    gsim[0].Draw("same p")

    gsim[1].SetMarkerStyle(24)
    gsim[1].SetMarkerSize(0.6)
    gsim[1].SetMarkerColor(2)
    gsim[1].SetLineColor(2)
    gsim[1].Draw("same p")
    
    gtrue[0].SetMarkerStyle(20)
    gtrue[0].SetMarkerSize(0.6)
    gtrue[0].SetMarkerColor(4)
    gtrue[0].SetLineColor(4)
    gtrue[0].Draw("same p")

    gtrue[1].SetMarkerStyle(24)
    gtrue[1].SetMarkerSize(0.6)
    gtrue[1].SetMarkerColor(4)
    gtrue[1].SetLineColor(4)
    gtrue[1].Draw("same p")

    #raw_input()
    c.hists=[gdata,gsim,gtrue,hframe,leg]
    return c

def DrawAccuracyByType(channels,types):
    c=ROOT.gROOT.MakeDefCanvas()
    c.SetLeftMargin(0.25)

    gtrue=ROOT.TGraphErrors()
    for i in range(len(types)*len(channels)):
        ic=i//len(types)
        it=i%len(types)        
        channel=channels[ic]
        g=gtrue
        value,error=GetTrueAccuracyByType(channel,types[it])
        g.SetPoint(2*i,value[0],2*i+0.5+0.2)
        g.SetPointError(2*i,error[0],0)
        g.SetPoint(2*i+1,value[1],2*i+1+0.5+0.2)
        g.SetPointError(2*i+1,error[1],0)
    
    hframe=ROOT.TH2D("hframe","",100,0.60,0.70,len(channels)*len(types)*2,0,len(channels)*len(types)*2);
    for i in range(len(types)*len(channels)):
        ic=i//len(types)
        it=i%len(types)
        hframe.GetYaxis().SetBinLabel(i*2+1,"#alpha^{-} ("+channels[ic]+", type "+str(types[it])+")")
        hframe.GetYaxis().SetBinLabel(i*2+2,"#alpha^{+} ("+channels[ic]+", type "+str(types[it])+")")
    hframe.SetStats(0)
    hframe.Draw()
    gtrue.SetMarkerStyle(20)
    gtrue.SetMarkerSize(0.6)
    gtrue.SetMarkerColor(4)
    gtrue.SetLineColor(4)
    gtrue.Draw("same p")
    
    raw_input()

if __name__=="__main__":
    #value,cov=GetAccuracy(0,"me201[678][ab]?")
    #value,cov=GetAccuracy(1,"me201[678][ab]?")
    DrawAccuracy(["mn201[678][ab]?","en201[678][ab]?","me201[678][ab]?","mm201[678][ab]?","ee201[678][ab]?"])
    #DrawAccuracy(["ee201[678][ab]?","mm201[678][ab]?","me201[678][ab]?"])
    #DrawAccuracy([channel+era for channel in ["mn","en","me","mm","ee"] for era in ["2016a","2016b","2017","2018","201[678][ab]?"]])
    #DrawAccuracy([channel+era for channel in ["[me]n","[me][me]"] for era in ["2016a","2016b","2017","2018","201[678][ab]?"]])
    
    #DrawAccuracyByType(["mn201[678][ab]?","en201[678][ab]?","me201[678][ab]?","mm201[678][ab]?","ee201[678][ab]?"],[0,1,2])
    #DrawAccuracyByType([channel+era for channel in ["mn","en"] for era in ["2016a","2016b","2017","2018","201[678][ab]?"]],[0,1,2])
