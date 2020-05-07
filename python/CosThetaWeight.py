import sys,os
import ROOT

def GetCosThetaWeight(hden,hnum):
    if hden.ClassName() is not "TH4D":
        print("Source is not TH4D... return NULL")
        return None
    if hnum.ClassName() is not "TH4D":
        print("Target is not TH4D... return NULL")
        return None
    if not hden.CheckConsistency(hnum):
        print("Source and Target are not consistent... return NULL")
        return None
    
    hout=hden.Clone("CosThetaWeight")
    hout.Reset()
    nx=hout.GetXaxis().GetNbinsX()
    ny=hout.GetYaxis().GetNbinsY()
    nz=hout.GetZaxis().GetNbinsZ()
    for ix in range(nx+2):
        for iy in range(ny+2):
            for iz in range(nz+2):
                this_hden=hden.ProjectionU("d_{}_{}_{}".format(ix,iy,iz),ix,ix,iy,iy,iz,iz)
                this_hnum=hnum.ProjectionU("n_{}_{}_{}".format(ix,iy,iz),ix,ix,iy,iy,iz,iz)
                this_hden.Scale(1/this_hden.Integral())
                this_hnum.Scale(1/this_hnum.Integral())
                this_hnum.Divide(this_hden)
                
                
