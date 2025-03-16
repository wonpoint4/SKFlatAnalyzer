import math
import numpy as np
import matplotlib.pyplot as plt
import ROOT
ROOT.gROOT.ProcessLine('#include"dybPlotter.cc"')
ROOT.Plotter.SetupStyle()

dyb = ROOT.dybPlotter()
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
    "btagSF lcorr", "btagSF luncorr2016a", "btagSF luncorr2016b", "btagSF luncorr2017", "btagSF luncorr2018", "bChargeSF eig0", "bChargeSF1 eig1",
]

def calPrecision(chi2s, nametag):
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

    print(nametag, "sin2w central : ", (x1 + x2) / 2, "1 sigma : ",(x2 - x1) / 2, "sin2w range : ", x1, x2)
    plt.plot(sin2ws, chi2s, 'o', color='black')
    plt.plot(x, y, color='blue')
    z=np.full(len(y), miny+1)
    plt.plot(x, z, color='red')

    plt.title(r"$\chi^{2}$ Fitting ("+nametag+")", fontsize=15)
    plt.xlabel("$sin^{2}\\theta^{l}_{eff}$")
    plt.ylabel("$\chi^{2}$")
    plt.text(0.224, 0.1, "$\chi^{2}_{max}$ = %.3f" % chi2s[0])
    plt.text(0.231, chi2s[0] * 0.79, "$\sin^{2}\\theta^{l}_{eff}$ = %.5f $\\pm$ %.5f" % ((x1 + x2) / 2, (x2 - x1) / 2))

    plt.legend(loc='upper right')
    plt.grid()
    plt.savefig("./precision_"+nametag+".png", dpi=300, bbox_inches='tight')
    plt.close()

def getdAFB(AFB_nominal, AFB):
    dAFB = []
    for iBin in range(NmassBins):
        dAFB.append(AFB_nominal.GetBinContent(iBin + 1) - AFB.GetBinContent(iBin + 1))

    return np.array(dAFB).reshape(len(dAFB), 1)

def caldAFBs(channel):
    dAFBs = [[] for i in range(len(chargeBins))] # dAFBs[chargeBins][sin2w scenarios]

    for ch in range(len(chargeBins)):
        dAFB = np.zeros(NmassBins)
        for sin in range(len(sin2w_indice)):
            iSin = sin2w_indice[sin]
            AFB_mc_nominal = dyb.GetHist(1, channel+chargeBins[ch]+"AFBrecoil(x)", "")
            AFB_mc_sin2w = dyb.GetHist(1, channel+chargeBins[ch]+"AFBrecoil(x)", "suffix:_sthw2_%i:dy" % iSin);
            dAFB = getdAFB(AFB_mc_nominal, AFB_mc_sin2w)
            dAFBs[ch].append(dAFB) # dAFB per each sin2w scenario
            print("\n dAFBs of "+channel+chargeBins[ch]+", sin2w_variation : %d" % iSin)
            #print("sin2w = ",sin2w_values[iSin], ", dAFB = ", dAFBs[ch][sin])

    return dAFBs

def calCovs(channel, option=""):
    if "1D" in option: chargeBins = ["y[0,5]/"]
    else: chargeBins = ["y[0,5]/", "y[0,0.1]/", "y[0.1,0.2]/", "y[0.2,0.6]/", "y[0.6,1]/", "y[1,3]/", "y[3,5]/"]
    covs = [[] for i in range(len(chargeBins))] # covs[chargeBins][each cov terms]

    for ch in range(len(chargeBins)):
        AFB_data_nominal = dyb.GetHist(0, channel+chargeBins[ch]+"AFBrecoil(x)", "")
        AFB_mc_nominal = dyb.GetHist(1, channel+chargeBins[ch]+"AFBrecoil(x)", "")

        cov_stat_data = np.zeros((NmassBins, NmassBins))
        cov_stat_mc = np.zeros((NmassBins, NmassBins))
        for iBin in range(NmassBins):
            cov_stat_data[iBin][iBin] = AFB_data_nominal.GetBinError(iBin + 1)**2
            cov_stat_mc[iBin][iBin] = AFB_mc_nominal.GetBinError(iBin + 1)**2
        covs[ch].append(cov_stat_data) # cov_stat_data
        covs[ch].append(cov_stat_mc) # cov_stat_mc
        print("\n channel+chargeBin : "+channel+chargeBins[ch])
        print("Sqrt of Trace(Cstat_data) = ", np.trace(covs[ch][0])**0.5, "Det(Cstat_data) = ", np.linalg.det(covs[ch][0]), "Condition(Cstat_data) = ", np.linalg.cond(covs[ch][0]))
        print("Sqrt of Trace(Cstat_mc) = ", np.trace(covs[ch][1])**0.5, "Det(Cstat_mc) = ", np.linalg.det(covs[ch][1]), "Condition(Cstat_mc) = ", np.linalg.cond(covs[ch][1]))

        if "syst" in option:
            print("\n Systematics")
            for systs in range(len(allsysts)):
                print("\n -"+systnames[systs])
                cov_syst_bigger = np.zeros((NmassBins, NmassBins))
                for syst in range(len(allsysts[systs])):
                    AFB_mc_syst = dyb.GetHist(1, channel+chargeBins[ch]+"AFBrecoil"+allsysts[systs][syst]+"(x)", "");
                    dAFB = getdAFB(AFB_mc_nominal, AFB_mc_syst)
                    cov_syst = np.outer(dAFB, dAFB)
                    print("cov syst", systs, syst)
                    #print(cov_syst)
                    if np.trace(cov_syst) > np.trace(cov_syst_bigger): cov_syst_bigger = cov_syst ## Choose the cov_systematic with the larger trace
                covs[ch].append(cov_syst_bigger) # each cov_syst
                print(systnames[systs]+", Sqrt of Trace(Csyst) = ", np.trace(cov_syst_bigger)**0.5)

    return covs

def calChi2sWithCov(covs, dAFBs, mode=-1):
    cov = np.zeros((NmassBins, NmassBins))
    for i in range(len(covs)):
        if i == mode: continue # For (N-1) syst uncertainties
        else: cov += covs[i]

    inv_cov = np.linalg.inv(cov)
    print("\n cov : ")
    #print(cov)
    print("\n inv_cov : ")
    #print(inv_cov)
    print("")
    print("Sqrt of Trace(cov) = ", np.trace(cov)**0.5, ", Sqrt of Trace(inv_cov) = ", np.trace(inv_cov)**0.5, "Det(inv_cov) = ", np.linalg.det(inv_cov), "Condition(inv_cov) = ", np.linalg.cond(inv_cov))
    #print(np.dot(cov, inv_cov))

    chi2s = []
    for i in range(len(sin2w_indice)):
        chi2 = np.linalg.multi_dot([dAFBs[i].transpose(), inv_cov, dAFBs[i]])
        chi2s.append(chi2[0][0])

    print("chi2s = ", chi2s)
    return chi2s

if __name__=="__main__":
    channel = "[em][em]201[678][ab]?/"
    dAFBs = caldAFBs(channel)

    covs_stat = calCovs(channel, "") # "" or "1D"
    chi2s_1D_stat = calChi2sWithCov(covs_stat[0], dAFBs[0])
    chi2s_2D_stat = [0] * len(sin2w_indice)
    for i in range(1, len(chargeBins)):
        chi2s = calChi2sWithCov(covs_stat[i], dAFBs[i])
        chi2s_2D_stat = [chi2s_2D_stat[j] + chi2s[j] for j in range(len(chi2s_2D_stat))]

    calPrecision(chi2s_1D_stat, "1D_stat")
    calPrecision(chi2s_2D_stat, "2D_stat")

    covs_total = calCovs(channel, "syst") # "syst" or "1D syst"
    chi2s_1D_total = calChi2sWithCov(covs_total[0], dAFBs[0])
    chi2s_2D_total = [0] * len(sin2w_indice)
    for i in range(1, len(chargeBins)):
        chi2s = calChi2sWithCov(covs_total[i], dAFBs[i])
        chi2s_2D_total = [chi2s_2D_total[j] + chi2s[j] for j in range(len(chi2s_2D_total))]

    calPrecision(chi2s_1D_total, "1D_total")
    calPrecision(chi2s_2D_total, "2D_total")

    # (N-1) stat + syst uncertainties
    for n in range(len(covs_total[0])):
        chi2s_1D_total_N_1 = calChi2sWithCov(covs_total[0], dAFBs[0], n)
        chi2s_2D_total_N_1 = [0] * len(sin2w_indice)
        for i in range(1, len(chargeBins)):
            chi2s = calChi2sWithCov(covs_total[i], dAFBs[i], n)
            chi2s_2D_total_N_1 = [chi2s_2D_total_N_1[j] + chi2s[j] for j in range(len(chi2s_2D_total_N_1))]

        nosource = "no"
        if n == 0: nosource += "DataStat"
        elif n == 1: nosource += "MCStat"
        else: nosource += systnames[n-2].replace(", ", "").replace(" ", "_")
        calPrecision(chi2s_1D_total_N_1, ("1D_total_%i_" % n)+nosource)
        calPrecision(chi2s_2D_total_N_1, ("2D_total_%i_" % n)+nosource)
