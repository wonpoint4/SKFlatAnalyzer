import array,sys
import ROOT
from scipy import stats
from scipy.optimize import minimize
from scipy.optimize import basinhopping,shgo,brute
import numpy as np

ROOT.gROOT.ProcessLine('#include"ZpeakPlotter.cc"')
p=ROOT.ZpeakPlotter("data-qcdss mi+tau_mi+vv+wjets+tt+st+aa")

def Apply(hist,pars):
    scale,centralRes,leftFrac,leftRes,rightFrac,rightRes=pars
    #print "[Apply]",hist.GetName(),scale,centralRes,leftFrac,leftRes,rightFrac,rightRes
    bins=[hist.GetBinLowEdge(i) for i in range(hist.GetNcells())]
    bins+=[np.inf]
    bins[0]=-np.inf

    precision=0.1
    h=hist.Clone()
    for i in range(1,hist.GetNcells()-1):
        center=hist.GetBinCenter(i)
        val=hist.Interpolate(center*(1-scale))/(1+scale)
        err=(val/hist.GetBinContent(i))**0.5*hist.GetBinError(i)
        h.SetBinContent(i,val)
        h.SetBinError(i,err)

    ws=[0. for _ in range(h.GetNcells())]
    w2s=[0. for _ in range(h.GetNcells())]
    for i in range(1,h.GetNcells()-1):
        center=h.GetBinCenter(i)
        centralFrac=1-leftFrac-rightFrac
        if abs(centralRes)<1e-10:
            ws[i]+=h.GetBinContent(i)*centralFrac
            w2s[i]+=h.GetBinError(i)**2*centralFrac
        else:
            frac=centralFrac*(stats.norm.cdf((bins[i+1]/center-1)/centralRes)-stats.norm.cdf((bins[i]/center-1)/centralRes))
            ws[i]+=h.GetBinContent(i)*frac
            w2s[i]+=h.GetBinError(i)**2*frac
            for j in range(i-1,-1,-1):
                frac=centralFrac*(stats.norm.cdf((bins[j+1]/center-1)/centralRes)-stats.norm.cdf((bins[j]/center-1)/centralRes))
                ws[j]+=h.GetBinContent(i)*frac
                w2s[j]+=h.GetBinError(i)**2*frac
                if h.GetBinContent(i)*frac < precision : break
            for j in range(i+1,h.GetNcells()):
                frac=centralFrac*(stats.norm.cdf((bins[j+1]/center-1)/centralRes)-stats.norm.cdf((bins[j]/center-1)/centralRes))
                ws[j]+=h.GetBinContent(i)*frac
                w2s[j]+=h.GetBinError(i)**2*frac
                if h.GetBinContent(i)*frac < precision : break

        if abs(leftRes)<1e-10:
            ws[i]+=h.GetBinContent(i)*leftFrac
            w2s[i]+=h.GetBinError(i)**2*leftFrac
        else:
            frac=2*leftFrac*(0.5-stats.norm.cdf((bins[i]/center-1)/leftRes))
            ws[i]+=h.GetBinContent(i)*frac
            w2s[i]+=h.GetBinError(i)**2*frac
            for j in range(i-1,-1,-1):
                frac=2*leftFrac*(stats.norm.cdf((bins[j+1]/center-1)/leftRes)-stats.norm.cdf((bins[j]/center-1)/leftRes))
                ws[j]+=h.GetBinContent(i)*frac
                w2s[j]+=h.GetBinError(i)**2*frac
                if h.GetBinContent(i)*frac < precision : break

        if abs(rightRes)<1e-10:
            ws[i]+=h.GetBinContent(i)*rightFrac
            w2s[i]+=h.GetBinError(i)**2*rightFrac
        else:
            frac=2*rightFrac*(stats.norm.cdf((bins[i+1]/center-1)/rightRes)-0.5)
            ws[i]+=h.GetBinContent(i)*frac
            w2s[i]+=h.GetBinError(i)**2*frac
            for j in range(i+1,h.GetNcells()):
                frac=2*rightFrac*(stats.norm.cdf((bins[j+1]/center-1)/rightRes)-stats.norm.cdf((bins[j]/center-1)/rightRes))
                ws[j]+=h.GetBinContent(i)*frac
                w2s[j]+=h.GetBinError(i)**2*frac
                if h.GetBinContent(i)*frac < precision : break

    for i in range(1,len(ws)-1):
        h.SetBinContent(i,ws[i])
        h.SetBinError(i,w2s[i]**0.5)
        
    return h

def ApplyBoth(hdata,hsim,pars):
    scale,centralRes,leftFrac,leftRes,rightFrac,rightRes=pars
    scaleData,scaleSim=-scale,0.
    centralResData,centralResSim=(0.,centralRes) if centralRes>0 else (abs(centralRes),0.)
    leftFracData,leftFracSim=(0.,leftFrac) if leftRes>0 else (leftFrac,0.)
    leftResData,leftResSim=(0.,leftRes) if leftRes>0 else (abs(leftRes),0.)
    rightFracData,rightFracSim=(0.,rightFrac) if rightRes>0 else (rightFrac,0.)
    rightResData,rightResSim=(0.,rightRes) if rightRes>0 else (abs(rightRes),0.)

    hdata=Apply(hdata,[scaleData,centralResData,leftFracData,leftResData,rightFracData,rightResData])
    hsim=Apply(hsim,[scaleSim,centralResSim,leftFracSim,leftResSim,rightFracSim,rightResSim])
    return hdata,hsim

def GetChi2(hdata,hsim,pars,bins=[5800,6700,7400,8000,8600,9200,10000,11200]):
    hdata,hsim=ApplyBoth(hdata,hsim,pars)
    chi2=0
    ndf=0
    #newbins=array.array("d",[5800,6300,6700,7100,7400,7700,8000,8300,8600,8900,9200,9600,10000,10600,11200])
    newbins=array.array("d",bins)
    xmin,xmax=newbins[0],newbins[-1]
    hdata=hdata.Rebin(len(newbins)-1,"hdata",newbins)
    hsim=hsim.Rebin(len(newbins)-1,"hdata",newbins)
    hdata.GetXaxis().SetRangeUser(xmin,xmax)
    hsim.GetXaxis().SetRangeUser(xmin,xmax)
    hsim.Scale(hdata.Integral()/hsim.Integral())
    for i in range(hdata.FindBin(xmin),hdata.FindBin(xmax)):
        chi2+=(hdata.GetBinContent(i)-hsim.GetBinContent(i))**2/(hdata.GetBinError(i)**2+hsim.GetBinError(i)**2)
        ndf+=1
    return chi2,ndf

def Merge(hists):
    hist=None
    for h in hists:
        if not h: continue
        if not hist:
            hist=h.Clone()
        else:
            hist.Add(h)
    return hist

def Test(channel="mm2016a",option=""):
    etabins=array.array("d",[0.0,0.4,0.8,1.2,1.6,2.0,2.5])    
    hdatas=[[p.GetHist(0,channel+"/lmetalpetam2",option+" project:z xmin:5800 xmax:11200 Xmin:{} Xmax:{} Ymin:{} Ymax:{}".format(etabins[ix],etabins[ix+1],etabins[iy],etabins[iy+1])) for iy in range(len(etabins)-1)] for ix in range(len(etabins)-1)]
    hsims=[[p.GetHist(1,channel+"/lmetalpetam2",option+" project:z xmin:5800 xmax:11200 Xmin:{} Xmax:{} Ymin:{} Ymax:{}".format(etabins[ix],etabins[ix+1],etabins[iy],etabins[iy+1])) for iy in range(len(etabins)-1)] for ix in range(len(etabins)-1)]

    target=5
    hdata=Merge([hdatas[ix][iy] for ix in range(len(etabins)-1) for iy in range(len(etabins)-1) if target in [iy]])
    hsim=Merge([hsims[ix][iy] for ix in range(len(etabins)-1) for iy in range(len(etabins)-1) if target in [iy]])

    bounds=((-0.1,0.1),(-0.1,0.1),(0,0.1),(-1,1),(0,0.1),(-1,1))
    bounds_central=((-0.001,0.001),(-0.01,0.01))
    #result=minimize(lambda x:GetChi2(hdata,hsim,x)[0],[0,-0.003,0.003,0.3,0.002,0.3],bounds=bounds,options={'disp':True})
    #result=minimize(lambda x:GetChi2(hdata,hsim,list(x)+[0,0,0,0],bins=[7400,7700,8000,8300,8600,8900,9200])[0],[0.0,0.0],bounds=bounds_central,options={'disp':True})
    result=shgo(lambda x:GetChi2(hdata,hsim,list(x)+[0,0,0,0],bins=[7400,7700,8000,8300,8600,8900,9200])[0],sampling_method='sobol',n=20,bounds=bounds_central,options={'disp':True},minimizer_kwargs={'disp':True,'options':{'disp':True}})
    #result=brute(lambda x:GetChi2(hdata,hsim,x)[0],bounds,Ns=3,disp=True)
    print result
    pars=list(result.x.copy())+[0,0,0,0]

    print pars
    print GetChi2(hdata,hsim,[0,0,0,0,0,0])
    print GetChi2(hdata,hsim,pars)

    this_hdata,this_hsim=ApplyBoth(hdata,hsim,pars)

    hdata.GetXaxis().SetRangeUser(5800,11200)
    hsim.GetXaxis().SetRangeUser(5800,11200)
    hsim.Scale(hdata.Integral()/hsim.Integral())
    this_hdata.GetXaxis().SetRangeUser(5800,11200)
    this_hsim.GetXaxis().SetRangeUser(5800,11200)
    this_hsim.Scale(this_hdata.Integral()/this_hsim.Integral())
    p.DrawPlot([hdata,hsim],"norm chi2 2:finey")
    p.DrawPlot([this_hdata,this_hsim],"norm chi2 2:finey")
    raw_input()

def Test2(channel="mm2016a",option=""):
    etabins=array.array("d",[0.0,0.4,0.8,1.2,1.6,2.0,2.5])
    target=2
    hdata=Merge([p.GetHist(0,channel+"/lmetalpetam2",option+" project:z xmin:5800 xmax:11200 Ymin:{} Ymax:{}".format(etabins[target],etabins[target+1])),p.GetHist(0,channel+"/lmetalpetam2",option+" project:z xmin:5800 xmax:11200 Xmin:{} Xmax:{}".format(etabins[target],etabins[target+1]))])
    hsim=Merge([p.GetHist(1,channel+"/lmetalpetam2",option+" project:z xmin:5800 xmax:11200 Ymin:{} Ymax:{}".format(etabins[target],etabins[target+1])),p.GetHist(1,channel+"/lmetalpetam2",option+" project:z xmin:5800 xmax:11200 Xmin:{} Xmax:{}".format(etabins[target],etabins[target+1]))])

    #mm2016a 5
    #this_hdata,this_hsim=ApplyBoth(hdata,hsim,[ 0.00026291, -0.00924583,  0.010691,    0.11306612,  0.00677987,  0.1711329 ])
    #this_hdata,this_hsim=ApplyBoth(hdata,hsim,[ 0.00026291, -0.00924583,  0,    0,  0,  0 ])
    
    this_hdata,this_hsim=ApplyBoth(hdata,hsim,[0.00023113001816805237, 0.00650435953211277, 0.002952356612697262, -0.1254816377070748, 0.0023439877811243944, -0.1378167760282065])
    
    #ee2016b
    #this_hdata,this_hsim=ApplyBoth(hdata,hsim,[1.2081525587878775e-05, 0.001976545272106665, 0.02, -0.010980567520955857, 0.001000000025953911, -0.19724580640465905])
    #this_hdata,this_hsim=ApplyBoth(hdata,hsim,[-3.2081525587878775e-04, 0.005976545272106665, 0.02, -0.030980567520955857, 0.001000000025953911, -0.19724580640465905])

    hdata.GetXaxis().SetRangeUser(5800,11200)
    hsim.GetXaxis().SetRangeUser(5800,11200)
    hsim.Scale(hdata.Integral()/hsim.Integral())
    this_hdata.GetXaxis().SetRangeUser(5800,11200)
    this_hsim.GetXaxis().SetRangeUser(5800,11200)
    this_hsim.Scale(this_hdata.Integral()/this_hsim.Integral())
    p.DrawPlot([hdata,hsim],"norm chi2 2:finey")
    p.DrawPlot([this_hdata,this_hsim],"norm chi2 2:finey")
    raw_input()


def Eval(channel,option=""):
    etabins=array.array("d",[0.0,0.4,0.8,1.2,1.6,2.0,2.5])
    hscale=ROOT.TH1D(channel+"_scale",channel+"_scale",len(etabins)-1,etabins)
    hcentralRes=ROOT.TH1D(channel+"_centralRes",channel+"_centralRes",len(etabins)-1,etabins)
    hleftFrac=ROOT.TH1D(channel+"_leftFrac",channel+"_leftFrac",len(etabins)-1,etabins)
    hleftRes=ROOT.TH1D(channel+"_leftRes",channel+"_leftRes",len(etabins)-1,etabins)
    hrightFrac=ROOT.TH1D(channel+"_rightFrac",channel+"_rightFrac",len(etabins)-1,etabins)
    hrightRes=ROOT.TH1D(channel+"_rightRes",channel+"_rightRes",len(etabins)-1,etabins)

    hdatas=[[p.GetHist(0,channel+"/lmetalpetam2",option+" project:z xmin:5800 xmax:11200 Xmin:{} Xmax:{} Ymin:{} Ymax:{}".format(etabins[ix],etabins[ix+1],etabins[iy],etabins[iy+1])) for iy in range(len(etabins)-1)] for ix in range(len(etabins)-1)]
    hsims=[[p.GetHist(1,channel+"/lmetalpetam2",option+" project:z xmin:5800 xmax:11200 Xmin:{} Xmax:{} Ymin:{} Ymax:{}".format(etabins[ix],etabins[ix+1],etabins[iy],etabins[iy+1])) for iy in range(len(etabins)-1)] for ix in range(len(etabins)-1)]
    pars=[[0.0,0.0,0.0,0.0,0.0,0.0] for i in range(len(etabins)-1)]

    niter=2
    for it in range(niter):
        print "[Iteration]",it
        this_hmdatas=[[None for ix in range(len(etabins)-1)] for iy in range(len(etabins)-1)]
        this_hmsims=[[None for ix in range(len(etabins)-1)] for iy in range(len(etabins)-1)]
        this_hpdatas=[[None for ix in range(len(etabins)-1)] for iy in range(len(etabins)-1)]
        this_hpsims=[[None for ix in range(len(etabins)-1)] for iy in range(len(etabins)-1)]
        print "Prepare"
        for ix in range(len(etabins)-1):
            for iy in range(len(etabins)-1):
                this_hmdatas[ix][iy],this_hmsims[ix][iy]=ApplyBoth(hdatas[ix][iy],hsims[ix][iy],pars[ix])
                this_hpdatas[ix][iy],this_hpsims[ix][iy]=ApplyBoth(hdatas[ix][iy],hsims[ix][iy],pars[iy])
        for i in range(len(etabins)-1):
            hdata=Merge([this_hmdatas[ix][iy] for ix in range(len(etabins)-1) for iy in range(len(etabins)-1) if i==iy]+[this_hpdatas[ix][iy] for ix in range(len(etabins)-1) for iy in range(len(etabins)-1) if i==ix])
            hsim=Merge([this_hmsims[ix][iy] for ix in range(len(etabins)-1) for iy in range(len(etabins)-1) if i==iy]+[this_hpsims[ix][iy] for ix in range(len(etabins)-1) for iy in range(len(etabins)-1) if i==ix])

            bounds=((-0.001,0.001),(-0.01,0.01),(0.001,0.02),(-0.3,0.3),(0.001,0.02),(-0.3,0.3))
            print "Fit central",i
            if it==0:
                result=shgo(lambda x:GetChi2(hdata,hsim,list(x)+pars[i][2:],bins=[7700,7800,7900,8000,8100,8200,8300,8400,8500,8600,8700,8800,8900])[0],bounds[:2],sampling_method='sobol',n=20,options={'disp':True})
            elif it==niter-1:
                result=shgo(lambda x:GetChi2(hdata,hsim,list(x)+pars[i][2:],bins=[7700,7800,7900,8000,8100,8200,8300,8400,8500,8600,8700,8800,8900])[0],bounds[:2],sampling_method='sobol',n=20,options={'disp':True})
                #result=minimize(lambda x:GetChi2(hdata,hsim,x)[0],pars[i],method='TNC',bounds=bounds,options={'disp':True})
            else:
                result=shgo(lambda x:GetChi2(hdata,hsim,list(x)+pars[i][2:],bins=[7700,7800,7900,8000,8100,8200,8300,8400,8500,8600,8700,8800,8900])[0],bounds[:2],sampling_method='sobol',n=20,options={'disp':True})
            pars[i][0],pars[i][1]=result.x

            print "Fit left",i
            if it==0:
                result=shgo(lambda x:GetChi2(hdata,hsim,pars[i][:2]+list(x)+pars[i][4:],bins=[5800,6700,7400,8000,8600,9200])[0],bounds[2:4],sampling_method='sobol',n=20,options={'disp':True})
            elif it==niter-1:
                result=shgo(lambda x:GetChi2(hdata,hsim,pars[i][:2]+list(x)+pars[i][4:],bins=[5800,6700,7400,8000,8600,9200])[0],bounds[2:4],sampling_method='sobol',n=20,options={'disp':True})
                #result=minimize(lambda x:GetChi2(hdata,hsim,x)[0],pars[i],method='TNC',bounds=bounds,options={'disp':True})
            else:
                result=shgo(lambda x:GetChi2(hdata,hsim,pars[i][:2]+list(x)+pars[i][4:],bins=[5800,6700,7400,8000,8600,9200])[0],bounds[2:4],sampling_method='sobol',n=20,options={'disp':True})
            pars[i][2],pars[i][3]=result.x

            print "Fit right",i
            if it==0:
                result=shgo(lambda x:GetChi2(hdata,hsim,pars[i][:4]+list(x),bins=[7400,8000,8600,9200,10000,11200])[0],bounds[4:],sampling_method='sobol',n=20,options={'disp':True})
            elif it==niter-1:
                result=shgo(lambda x:GetChi2(hdata,hsim,pars[i][:4]+list(x),bins=[7400,8000,8600,9200,10000,11200])[0],bounds[4:],sampling_method='sobol',n=20,options={'disp':True})
                #result=minimize(lambda x:GetChi2(hdata,hsim,x)[0],pars[i],method='TNC',bounds=bounds,options={'disp':True})
            else:
                result=shgo(lambda x:GetChi2(hdata,hsim,pars[i][:4]+list(x),bins=[7400,8000,8600,9200,10000,11200])[0],bounds[4:],sampling_method='sobol',n=20,options={'disp':True})
            pars[i][4],pars[i][5]=result.x

            print "before",GetChi2(hdata,hsim,[0,0,0,0,0,0])
            print "after",GetChi2(hdata,hsim,pars[i])
            print "pars",pars[i]

        for par in pars:
            print par

    for i in range(len(etabins)-1):
        hscale.SetBinContent(i+1,pars[i][0])
        hcentralRes.SetBinContent(i+1,pars[i][1])
        hleftFrac.SetBinContent(i+1,pars[i][2])
        hleftRes.SetBinContent(i+1,pars[i][3])
        hrightFrac.SetBinContent(i+1,pars[i][4])
        hrightRes.SetBinContent(i+1,pars[i][5])

    f=ROOT.TFile(channel+".root","recreate")
    hscale.Write()
    hcentralRes.Write()
    hleftFrac.Write()
    hleftRes.Write()
    hrightFrac.Write()
    hrightRes.Write()

def Validation(channel,option=""):
    etabins=array.array("d",[0.0,0.4,0.8,1.2,1.6,2.0,2.5])
    f=ROOT.TFile(channel+".root")
    hscale=f.Get(channel+"_scale")
    hcentralRes=f.Get(channel+"_centralRes")
    hleftFrac=f.Get(channel+"_leftFrac")
    hleftRes=f.Get(channel+"_leftRes")
    hrightFrac=f.Get(channel+"_rightFrac")
    hrightRes=f.Get(channel+"_rightRes")
    hpars=[hscale,hcentralRes,hleftFrac,hleftRes,hrightFrac,hrightRes]

    hdatas=[[p.GetHist(0,channel+"/lmetalpetam2",option+" project:z xmin:5800 xmax:11200 Xmin:{} Xmax:{} Ymin:{} Ymax:{}".format(etabins[ix],etabins[ix+1],etabins[iy],etabins[iy+1])) for iy in range(len(etabins)-1)] for ix in range(len(etabins)-1)]
    hsims=[[p.GetHist(1,channel+"/lmetalpetam2",option+" project:z xmin:5800 xmax:11200 Xmin:{} Xmax:{} Ymin:{} Ymax:{}".format(etabins[ix],etabins[ix+1],etabins[iy],etabins[iy+1])) for iy in range(len(etabins)-1)] for ix in range(len(etabins)-1)]
    pars=[[hpars[ip].GetBinContent(i+1) for ip in range(len(hpars))] for i in range(len(etabins)-1)]
    
    this_hdatas=[[None for ix in range(len(etabins)-1)] for iy in range(len(etabins)-1)]
    this_hsims=[[None for ix in range(len(etabins)-1)] for iy in range(len(etabins)-1)]
    for ix in range(len(etabins)-1):
        for iy in range(len(etabins)-1):
            this_hdatas[ix][iy],this_hsims[ix][iy]=ApplyBoth(hdatas[ix][iy],hsims[ix][iy],pars[ix])
            this_hdatas[ix][iy],this_hsims[ix][iy]=ApplyBoth(this_hdatas[ix][iy],this_hsims[ix][iy],pars[iy])
    for i in range(len(etabins)-1):
        hdata=Merge([hdatas[ix][iy] for ix in range(len(etabins)-1) for iy in range(len(etabins)-1) if i in [ix,iy]])
        hsim=Merge([hsims[ix][iy] for ix in range(len(etabins)-1) for iy in range(len(etabins)-1) if i in [ix,iy]])
        this_hdata=Merge([this_hdatas[ix][iy] for ix in range(len(etabins)-1) for iy in range(len(etabins)-1) if i in [ix,iy]])
        this_hsim=Merge([this_hsims[ix][iy] for ix in range(len(etabins)-1) for iy in range(len(etabins)-1) if i in [ix,iy]])
        
        hdata.GetXaxis().SetRangeUser(5800,11200)
        hsim.GetXaxis().SetRangeUser(5800,11200)
        hsim.Scale(hdata.Integral()/hsim.Integral())
        this_hdata.GetXaxis().SetRangeUser(5800,11200)
        this_hsim.GetXaxis().SetRangeUser(5800,11200)
        this_hsim.Scale(this_hdata.Integral()/this_hsim.Integral())
        p.DrawPlot([hdata,hsim],"norm chi2 2:finey")
        p.DrawPlot([this_hdata,this_hsim],"norm chi2 2:finey")
        p.DrawPlot(channel+"/dimass","norm chi2 2:finey project:z rebin:{{54,66,76,82,86,89.5,92.5,96,100,106,116,150}} widthweight absX Xmin:{} Xmax:{}".format(etabins[i],etabins[i+1]))
        p.DrawPlot(channel+"/dimass_roccor_residual","norm chi2 2:finey project:z rebin:{{54,66,76,82,86,89.5,92.5,96,100,106,116,150}} widthweight absX Xmin:{} Xmax:{}".format(etabins[i],etabins[i+1]))

    raw_input()
    
    
if __name__=="__main__":    
    if sys.argv[1]=="eval":
        Eval(sys.argv[2])
        exit()
    elif sys.argv[1]=="validation":
        Validation(sys.argv[2])
        exit()
    elif sys.argv[1]=="test":
        Test(sys.argv[2])
        exit()
    elif sys.argv[1]=="test2":
        Test2(sys.argv[2])
        exit()
    exit()
    #Eval("ee2016b")
    #Validation("mm2016a")
    #Validation("ee2016b")
    Test()
    #Test2()
    exit()
    hists=[]
    #for channel in ["ee2016a","ee2016b","ee2017","ee2018","mm2016a","mm2016b","mm2017","mm2018"]:
    #for channel in ["ee2016a","mm2016a"]:
    for channel in ["ee2016b","mm2016b"]:
        hists+=Eval(channel)
    f=ROOT.TFile("test.root","recreate")
    for h in hists:
        h.Write()
    f.Close()
    exit()

cs=[]

cs+=[test("mm2016a")]
# cs+=[test("mm2016a","Xmin:20 Xmax:30")]
# cs+=[test("mm2016a","Xmin:30 Xmax:40")]
# cs+=[test("mm2016a","Xmin:40 Xmax:50")]
# cs+=[test("mm2016a","Xmin:50 Xmax:1000")]
cs+=[test("mm2016a","suffix:_muonmomentum_residual")]
# cs+=[test("mm2016a","Xmin:20 Xmax:30 suffix:_muonmomentum_residual")]
# cs+=[test("mm2016a","Xmin:30 Xmax:40 suffix:_muonmomentum_residual")]
# cs+=[test("mm2016a","Xmin:40 Xmax:50 suffix:_muonmomentum_residual")]
# cs+=[test("mm2016a","Xmin:50 Xmax:1000 suffix:_muonmomentum_residual")]
raw_input()

# cs+=[test("mm2016a")]
# cs+=[test("mm2016a","suffix:_muonmomentum_residual")]
# cs+=[test("mm2016a","suffix:_residual:sim")]
# cs+=[test("ee2016a")]
# cs+=[test("ee2016a","suffix:_electronenergy_residual")]
# cs+=[test("ee2016a","suffix:_residual:sim")]
# raw_input()
# exit()

for channel in ["ee2016a","ee2016b","ee2017","ee2018","mm2016a","mm2016b","mm2017","mm2018"]:
#for channel in ["mm2016a"]:
    cs+=[test(channel)]
raw_input()
f=ROOT.TFile("test.root","recreate")
for c in cs:
    for h in c.hists:
        h.Write()
f.Close()
