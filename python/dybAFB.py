import math
import numpy as np
import matplotlib.pyplot as plt
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
allsysts = [
    ["_jet_scale_up", "_jet_scale_down"], ["_jet_smear_up", "_jet_smear_down"],
    ["_prefireweight_up", "_prefireweight_down"], ["_PUweight_up", "_PUweight_down"], ["_PUjetSF_up", "_PUjetSF_down"],
    ["_btagSF_hup", "_btagSF_hdown"], ["_btagSF_lup", "_btagSF_ldown"],
    ["_btagSF_hcorr"], ["_btagSF_huncorr2016a"], ["_btagSF_huncorr2016b"], ["_btagSF_huncorr2017"], ["_btagSF_huncorr2018"],
    ["_btagSF_lcorr"], ["_btagSF_luncorr2016a"], ["_btagSF_luncorr2016b"], ["_btagSF_luncorr2017"], ["_btagSF_luncorr2018"],
    ["_bChargeSF1_up"], ["_bChargeSF1_down"],
]
systnames = ["JES up, down", "JER up, down", "Prefiring up, down", "PUreweight up, down", "PUjetIDSF up, down",
    "btagSF hup, down", "btagSF lup, down", "btagSF hcorr", "btagSF huncorr2016a", "btagSF huncorr2016b", "btagSF huncorr2017", "btagSF huncorr2018",
    "btagSF lcorr", "btagSF luncorr2016a", "btagSF luncorr2016b", "btagSF luncorr2017", "btagSF luncorr2018", "bChargeSF1_up", "bChargeSF1_down",
]

def calPrecision(chi2s, nametag):
    if chi2s.shape[0] == 1 or "1D" in nametag: chi2s = chi2s[0]
    else: chi2s = chi2s.sum(axis=0) - chi2s[0]

    sin2ws = []
    for i in sin2w_indice:
        sin2ws.append(sin2w_values[i])

    fit = np.polyfit(sin2ws, chi2s, 2)
    pol2 = np.poly1d(fit)
    x = np.linspace(0.22300, 0.23800, 15000)
    y = pol2(x)
    miny = 100
    minx = 0
    x1 = 0
    x2 = 0
    for i in range(len(y)):
        if y[i] < miny:
            miny = y[i] # Find the minimum chi2
            minx = x[i] # and sintheta2

    for i in range(len(y)):
        yi = y[i]
        diff1 = 0.01
        diff2 = 0.01
        if i < len(x)/2 and abs(yi - (miny + 1)) < diff1:
            x1 = x[i]
            diff1 = abs(yi - (miny + 1))
        if i > len(x)/2 and abs(yi - (miny + 1)) < diff2:
            x2 = x[i]
            diff2 = abs(yi - (miny + 1))

    print("sin2w range : ", x1, x2)
    print("sin2w central : ", (x1 + x2) / 2, "1 sigma : ",(x2 - x1) / 2)
    plt.plot(sin2ws, chi2s, 'o', color='black')
    plt.plot(x, y, color='blue')
    z=np.full(len(y), miny+1)
    plt.plot(x, z, color='red')

    plt.xlabel("sin$^{2}\\theta$")
    plt.ylabel("$\\chi^2$")
    plt.title(r"$\chi^2$ fitting", fontsize=15)
    plt.legend(loc='upper right')
    plt.grid()
    plt.savefig("./precision_"+nametag+".png", dpi=300, bbox_inches='tight')
    plt.close()

def getdAFB(AFB_nominal, AFB):
    dAFB = []
    for iBin in range(NmassBins):
        dAFB.append(AFB_nominal.GetBinContent(iBin + 1) - AFB.GetBinContent(iBin + 1))

    return np.array(dAFB).reshape(len(dAFB), 1)

def calChi2WithCov(channel, option=""):
    if "1D" in option: chargeBins = ["y[0,5]/"]
    else: chargeBins = ["y[0,5]/", "y[0,0.1]/", "y[0.1,0.2]/", "y[0.2,0.6]/", "y[0.6,1]/", "y[1,3]/", "y[3,5]/"]
    chi2s = np.zeros((len(chargeBins), len(sin2w_values)))

    for ch in range(len(chargeBins)):
        chargeBin = chargeBins[ch]
        AFB_data_nominal = dyb.GetHist(0, channel+chargeBin+"AFBrecoil(x)", "")
        AFB_mc_nominal = dyb.GetHist(1, channel+chargeBin+"AFBrecoil(x)", "")
        cov_stat = np.zeros((NmassBins, NmassBins))
        cov  = np.zeros((NmassBins, NmassBins))
        for iBin in range(NmassBins):
            cov_stat[iBin][iBin] += AFB_data_nominal.GetBinError(iBin + 1)**2

        print("")
        print("")
        print("channel+chargeBin : "+channel+chargeBin+", cov_stat_data : ")
        print(cov_stat)
        print("Sqrt of Trace(Cstat_data) = ", np.trace(cov_stat)**0.5, "Det(Cstat_data) = ", np.linalg.det(cov_stat), "Condition(Cstat_data) = ", np.linalg.cond(cov_stat))
        for iBin in range(NmassBins):
            cov_stat[iBin][iBin] += AFB_mc_nominal.GetBinError(iBin + 1)**2

        print("channel+chargeBin : "+channel+chargeBin+", cov_stat_total : ")
        print(cov_stat)
        print("Sqrt of Trace(Cstat_total) = ", np.trace(cov_stat)**0.5, "Det(Cstat_total) = ", np.linalg.det(cov_stat), "Condition(Cstat_total) = ", np.linalg.cond(cov_stat))
        cov += cov_stat

        if "syst" in option:
            print("\n Systematics")
            for systs in range(len(allsysts)):
                print("\n -"+systnames[systs])
                cov_syst_bigger = np.zeros((NmassBins, NmassBins))
                for syst in range(len(allsysts[systs])):
                    AFB_mc_syst = dyb.GetHist(1, channel+chargeBin+"AFBrecoil"+allsysts[systs][syst]+"(x)", "");
                    dAFB = getdAFB(AFB_mc_nominal, AFB_mc_syst)
                    cov_syst = np.outer(dAFB, dAFB)
                    print("cov syst", systs, syst)
                    print(cov_syst)
                    if np.trace(cov_syst) > np.trace(cov_syst_bigger): cov_syst_bigger = cov_syst ## Choose the cov_systematic with the larger trace
                    cov += cov_syst_bigger
                    print(systnames[systs]+", Sqrt of Trace(Csyst) = ", np.trace(cov_syst_bigger)**0.5)

        inv_cov = np.linalg.inv(cov)
        print("\n channel+chargeBin : "+channel+chargeBin+", cov : ")
        print(cov)
        print("\n inv_cov : ")
        print(inv_cov)
        print("")
        print("Sqrt of Trace(cov) = ", np.trace(cov)**0.5, ", Sqrt of Trace(inv_cov) = ", np.trace(inv_cov)**0.5, "Det(inv_cov) = ", np.linalg.det(inv_cov), "Condition(inv_cov) = ", np.linalg.cond(inv_cov))
        print(np.dot(cov, inv_cov))

        for sin in range(len(sin2w_indice)):
            iSin = sin2w_indice[sin]
            AFB_mc_sin2w = dyb.GetHist(1, channel+chargeBin+"AFBrecoil(x)", "suffix:_sthw2_%i:dy" % iSin);
            dAFB = getdAFB(AFB_mc_nominal, AFB_mc_sin2w)
            print("")
            print(" sin2w_variation : ", iSin, "sin2w = ",sin2w_values[iSin], ", dAFB = ", dAFB)
            chi2 = np.linalg.multi_dot([dAFB.transpose(), inv_cov, dAFB])
            print("chi2 = ", chi2[0][0])
            chi2s[ch][sin] = chi2[0][0]

    print(chi2s)
    return chi2s

if __name__=="__main__":
    chi2s_stat = calChi2WithCov("[em][em]201[678][ab]?/", "") # "1D"
    chi2s_total = calChi2WithCov("[em][em]201[678][ab]?/", "syst") # "1D syst"

    calPrecision(chi2s_stat, "1D_stat")
    calPrecision(chi2s_stat, "2D_stat")
    calPrecision(chi2s_total, "1D_total")
    calPrecision(chi2s_total, "2D_total")
