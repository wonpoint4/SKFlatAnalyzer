import math
import numpy as np
import ROOT
ROOT.gROOT.ProcessLine('#include"dybPlotter.cc"')
ROOT.Plotter.SetupStyle()

dyb = ROOT.dybPlotter()
AFB_data_nominal = ROOT.TH1D()
AFB_mc_nominal = ROOT.TH1D()
sin2w_values = [0.23151, 0.23154, 0.23157, 0.2230, 0.2300, 0.2305, 0.2310, 0.2315, 0.2320, 0.2325, 0.2330]
sin2w_indice = [3, 4, 5, 6, 7, 0, 1, 2, 8, 9, 10]
chargeBins = ["y[0,5]/", "y[0,0.1]/", "y[0.1,0.2]/", "y[0.2,0.6]/", "y[0.6,1]/", "y[1,3]/", "y[3,5]/"]
NmassBins = 30 # 52 ~ 200 GeV instead of 52 ~ 3000 GeV. See afb_mbin[afb_mbinnum+1] in dybAnalyzer.h

def getdAFB(AFB_nominal, AFB):
    dAFB = []
    for iBin in range(NmassBins):
        dAFB.append(AFB_nominal.GetBinContent(iBin + 1) - AFB.GetBinContent(iBin + 1))

    return np.array(dAFB).reshape(len(dAFB), 1)

def calChi2WithCov(channel, syst=""):
    #chargeBins = ["y[0,5]/"]
    chi2s = np.zeros((len(chargeBins), len(sin2w_values)))

    for ch in range(len(chargeBins)):
        chargeBin = chargeBins[ch]
        AFB_data_nominal = dyb.GetHist(0, channel+chargeBin+"AFBrecoil"+syst+"(x)", "")
        AFB_mc_nominal = dyb.GetHist(1, channel+chargeBin+"AFBrecoil"+syst+"(x)", "")
        cov_stat = np.zeros((NmassBins, NmassBins))
        for iBin in range(NmassBins):
            cov_stat[iBin][iBin] += AFB_data_nominal.GetBinError(iBin + 1)**2
            cov_stat[iBin][iBin] += AFB_mc_nominal.GetBinError(iBin + 1)**2

        print("channel+chargeBin : "+channel+chargeBin+", cov_stat : ")
        print(cov_stat)
        inv_cov_stat = np.linalg.inv(cov_stat)
        print(inv_cov_stat)

        for sin in range(len(sin2w_indice)):
            iSin = sin2w_indice[sin]
            AFB_mc_sin2w = dyb.GetHist(1, channel+chargeBin+"AFBrecoil"+syst+"(x)", "suffix:_sthw2_%i:dy" % iSin);
            dAFB = getdAFB(AFB_mc_nominal, AFB_mc_sin2w)
            print("sin2w_variation : ", iSin, "sin2w = ",sin2w_values[iSin], ", dAFB = ", dAFB)
            print("shape of dAFB : ", np.shape(dAFB), "shape of dAFB.transpose() : ", np.shape(dAFB.transpose()))
            chi2 = np.linalg.multi_dot([dAFB.transpose(), inv_cov_stat, dAFB])
            print("chi2 = ", chi2[0][0])
            chi2s[ch][sin] = chi2[0][0]

    print(chi2s)

if __name__=="__main__":
    calChi2WithCov("[em][em]201[678][ab]?/")
