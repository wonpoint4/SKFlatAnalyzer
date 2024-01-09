## eff residual systematic set with iteration
import os,sys,re
import ROOT
ROOT.gROOT.LoadMacro("./Plotter/AFBPlotter.cc")
ROOT.TH1.AddDirectory(0)
ROOT.TH1.SetDefaultSumw2(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

def Make2D(h4d):
    h2d=ROOT.TH2D(h4d.GetName(),h4d.GetTitle(),h4d.GetXaxis().GetNbins(),h4d.GetXaxis().GetXbins().GetArray(),h4d.GetYaxis().GetNbins(),h4d.GetYaxis().GetXbins().GetArray())
    for ix in range(h4d.GetXaxis().GetNbins()+2):
        for iy in range(h4d.GetYaxis().GetNbins()+2):
            for iz in range(h4d.GetZaxis().GetNbins()+2):
                for iu in range(h4d.GetUaxis().GetNbins()+2):
                    val=h4d.GetBinContent(ix,iy,iz,iu)
                    err=h4d.GetBinError(ix,iy,iz,iu)
                    h2d.SetBinContent(ix,iy,h2d.GetBinContent(ix,iy)+val)
                    h2d.SetBinError(ix,iy,(h2d.GetBinError(ix,iy)**2+err**2)**0.5)
                    h2d.SetBinContent(iz,iu,h2d.GetBinContent(iz,iu)+val)
                    h2d.SetBinError(iz,iu,(h2d.GetBinError(iz,iu)**2+err**2)**0.5)
    return h2d

def Apply(h4d,h2d):
    for ix in range(h4d.GetXaxis().GetNbins()+2):
        for iy in range(h4d.GetYaxis().GetNbins()+2):
            for iz in range(h4d.GetZaxis().GetNbins()+2):
                for iu in range(h4d.GetUaxis().GetNbins()+2):
                    sf=1.
                    if h2d.GetBinContent(ix,iy): sf*=h2d.GetBinContent(ix,iy)
                    if h2d.GetBinContent(iz,iu): sf*=h2d.GetBinContent(iz,iu)
                    h4d.SetBinContent(ix,iy,iz,iu,h4d.GetBinContent(ix,iy,iz,iu)*sf)
                    h4d.SetBinError(ix,iy,iz,iu,h4d.GetBinError(ix,iy,iz,iu)*sf)
    return

def GetChi2(h1,h2):
    chi2=0.
    ndf=0
    for i in range(h1.GetNcells()):
        val1=h1.GetBinContent(i)
        err1=h1.GetBinError(i)
        val2=h2.GetBinContent(i)
        err2=h2.GetBinError(i)
        #print val1,err1,val2,err2
        if err1 or err2:
            chi2+=(val1-val2)**2/(err1**2+err2**2)
            ndf+=1
    return chi2,ndf,ROOT.TMath.Prob(chi2,ndf)

def MoveOverflow(h4d):
    nx=h4d.GetXaxis().GetNbins()
    ny=h4d.GetXaxis().GetNbins()
    nz=h4d.GetXaxis().GetNbins()
    nu=h4d.GetXaxis().GetNbins()
    for iy,iz,iu in [(iy,iz,iu) for iy in range(ny+2) for iz in range(nz+2) for iu in range(nu+2)]:
        val0=h4d.GetBinContent(0,iy,iz,iu)
        err0=h4d.GetBinError(0,iy,iz,iu)
        val1=h4d.GetBinContent(1,iy,iz,iu)
        err1=h4d.GetBinError(1,iy,iz,iu)
        h4d.SetBinContent(1,iy,iz,iu,val0+val1)
        h4d.SetBinError(1,iy,iz,iu,(err0**2+err1**2)**0.5)
        val0=h4d.GetBinContent(nx+1,iy,iz,iu)
        err0=h4d.GetBinError(nx+1,iy,iz,iu)
        val1=h4d.GetBinContent(nx,iy,iz,iu)
        err1=h4d.GetBinError(nx,iy,iz,iu)
        h4d.SetBinContent(nx,iy,iz,iu,val0+val1)
        h4d.SetBinError(nx,iy,iz,iu,(err0**2+err1**2)**0.5)
    for ix,iz,iu in [(ix,iz,iu) for ix in range(nx+2) for iz in range(nz+2) for iu in range(nu+2)]:
        val0=h4d.GetBinContent(ix,0,iz,iu)
        err0=h4d.GetBinError(ix,0,iz,iu)
        val1=h4d.GetBinContent(ix,1,iz,iu)
        err1=h4d.GetBinError(ix,1,iz,iu)
        h4d.SetBinContent(ix,1,iz,iu,val0+val1)
        h4d.SetBinError(ix,1,iz,iu,(err0**2+err1**2)**0.5)
        val0=h4d.GetBinContent(ix,ny+1,iz,iu)
        err0=h4d.GetBinError(ix,ny+1,iz,iu)
        val1=h4d.GetBinContent(ix,ny,iz,iu)
        err1=h4d.GetBinError(ix,ny,iz,iu)
        h4d.SetBinContent(ix,ny,iz,iu,val0+val1)
        h4d.SetBinError(ix,ny,iz,iu,(err0**2+err1**2)**0.5)
    for ix,iy,iu in [(ix,iy,iu) for ix in range(nx+2) for iy in range(ny+2) for iu in range(nu+2)]:
        val0=h4d.GetBinContent(ix,iy,0,iu)
        err0=h4d.GetBinError(ix,iy,0,iu)
        val1=h4d.GetBinContent(ix,iy,1,iu)
        err1=h4d.GetBinError(ix,iy,1,iu)
        h4d.SetBinContent(ix,iy,1,iu,val0+val1)
        h4d.SetBinError(ix,iy,1,iu,(err0**2+err1**2)**0.5)
        val0=h4d.GetBinContent(ix,iy,nz+1,iu)
        err0=h4d.GetBinError(ix,iy,nz+1,iu)
        val1=h4d.GetBinContent(ix,iy,nz,iu)
        err1=h4d.GetBinError(ix,iy,nz,iu)
        h4d.SetBinContent(ix,iy,nz,iu,val0+val1)
        h4d.SetBinError(ix,iy,nz,iu,(err0**2+err1**2)**0.5)
    for ix,iy,iz in [(ix,iy,iz) for ix in range(nx+2) for iy in range(ny+2) for iz in range(nz+2)]:
        val0=h4d.GetBinContent(ix,iy,iz,0)
        err0=h4d.GetBinError(ix,iy,iz,0)
        val1=h4d.GetBinContent(ix,iy,iz,1)
        err1=h4d.GetBinError(ix,iy,iz,1)
        h4d.SetBinContent(ix,iy,iz,1,val0+val1)
        h4d.SetBinError(ix,iy,iz,1,(err0**2+err1**2)**0.5)
        val0=h4d.GetBinContent(ix,iy,iz,nu+1)
        err0=h4d.GetBinError(ix,iy,iz,nu+1)
        val1=h4d.GetBinContent(ix,iy,iz,nu)
        err1=h4d.GetBinError(ix,iy,iz,nu)
        h4d.SetBinContent(ix,iy,iz,nu,val0+val1)
        h4d.SetBinError(ix,iy,iz,nu,(err0**2+err1**2)**0.5)

def addResidual(infilename):
    print infilename
    outfilename=infilename.replace(".root","_residual.root")

    if "2016preVFP" in infilename:
        era="2016a"
    elif "2016postVFP" in infilename:
        era="2016b"
    elif "2017" in infilename:
        era="2017"
    elif "2018" in infilename:
        era="2018"
    
    if "Electron" in infilename:
        channel="el"
    elif "Muon" in infilename:
        channel="mu"

    os.system("cp {} {}".format(infilename,outfilename));
    f=ROOT.TFile(outfilename,"update")
    hsim=f.Get("sim")
    hsf_origin=f.Get("sf")

    plotter=ROOT.AFBPlotter("data mi+tau_mi+vv+wjets+tt+st+qcdss+aa","EfficiencyValidation")
    hdata4d=plotter.GetHist(0,channel+era+"/m80to100/lpetaptlmetapt","noproject")
    #MoveOverflow(hdata4d)
    hsim4d=plotter.GetHist(1,channel+era+"/m80to100/lpetaptlmetapt","noproject")
    #MoveOverflow(hsim4d)
    scale=hdata4d.Integral(0,-1,0,-1,0,-1,0,-1)/hsim4d.Integral(0,-1,0,-1,0,-1,0,-1)
    hsim4d.Scale(scale)

    hsf=hsf_origin.Clone("hsf")
    chi2_old=1e6
    for i in range(10):
        hdata2d=Make2D(hdata4d)
        hsim2d=Make2D(hsim4d)
        chi2,ndf,prob=GetChi2(hdata2d,hsim2d)
        print i, chi2,ndf,prob
        this_hsf=hdata2d.Clone("this_hsf")
        this_hsf.Divide(hsim2d)
        for j in range(this_hsf.GetNcells()):
            #print this_hsf.GetBinContent(j), this_hsf.GetBinError(j)
            this_hsf.SetBinError(j,0)
        Apply(hsim4d,this_hsf)
        hsf.Multiply(this_hsf)
        if (chi2_old-chi2)/chi2_old<0.05:
            break
        chi2_old=chi2

    ## fluctuataion
#    hdata2d=Make2D(hdata4d)
#    hsim2d=Make2D(hsim4d)
#    chi2,ndf,prob=GetChi2(hdata2d,hsim2d)
#    print "before fluctuation", chi2, ndf, prob
#    errscale=(1-(chi2/ndf))**0.5
#    this_hsf=hdata2d.Clone("this_hsf")
#    this_hsf.Divide(hsim2d)
#    for j in range(this_hsf.GetNcells()):
#        err=this_hsf.GetBinError(j)
#        this_hsf.SetBinContent(j,ROOT.gRandom.Gaus(this_hsf.GetBinContent(j),err*errscale))
#        this_hsf.SetBinError(j,0)
#    Apply(hsim4d,this_hsf)
#    hsf.Multiply(this_hsf)

    print "final", GetChi2(Make2D(hdata4d),Make2D(hsim4d))

    hdata=hsim.Clone("data")
    hdata.Multiply(hsf)
    f.cd()
    for h in [hdata, hsim, hsf]:
        for i in range(h.GetNcells()):
            h.SetBinError(i,0)
    iset=max([int(re.match("sf_s([0-9]+)",key.GetName()).group(1)) for key in f.GetListOfKeys() if re.match("sf_s([0-9]+)",key.GetName())])+1
    hdata.SetNameTitle("data_s{}m0".format(iset),"residual")
    hsim.SetNameTitle("sim_s{}m0".format(iset),"residual")
    hsf.SetNameTitle("sf_s{}m0".format(iset),"residual")
    hdata.Write()
    hsim.Write()
    hsf.Write()
    f.Close()

    
def GetCurrentEffFileName(era,channel):
    if channel=="Electron":
        key="Muon_MediumID_trkIsoLoose"
    elif channel=="Muon":
        key="Electron_MediumID"
    else:
        print "Unknown channal ",channel
        exit(1)

    filename=os.popen("cat $SKFlat_WD/data/$SKFlatV/"+era+"/ID/eff.conf|egrep '"+key+"[^_]'|awk '{print $3}'").read().strip()
    filename=filename.replace("_residual","")
    return os.environ["SKFlat_WD"]+"/data/"+os.environ["SKFlatV"]+"/"+era+"/ID/"+filename

if __name__=="__main__":
    if sys.argv[1]=="all":
        for era in ["2016preVFP","2016postVFP","2017","2018"]:
            addResidual(GetCurrentEffFileName(era,"Electron"))
            addResidual(GetCurrentEffFileName(era,"Muon"))
    else:
        addResidual(sys.argv[1])
