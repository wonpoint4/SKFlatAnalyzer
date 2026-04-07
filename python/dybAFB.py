import os
import math
import numpy as np
import matplotlib.pyplot as plt
import ROOT
ROOT.gROOT.ProcessLine('#include"dybPlotter.cc"')
ROOT.Plotter.SetupStyle()

dyb = ROOT.dybPlotter("data ^dyb_mi+dyB_mi+dyall+ttall+ewkall", "dybAnalyzer_backup")
sin2w_values = [0.23151, 0.23154, 0.23157, 0.2230, 0.2300, 0.2305, 0.2310, 0.2315, 0.2320, 0.2325, 0.2330]
sin2w_indice = [3, 4, 5, 6, 7, 0, 1, 2, 8, 9, 10]
sin2w_values = [0.22654, 0.22854, 0.23054, 0.23104, 0.23154, 0.23204, 0.23254, 0.23454, 0.23654]
sin2w_indice = [1, 2, 3, 4, 5, 6, 7, 8, 9]
chargeBins = ["y[0,5]/", "y[0,0.1]/", "y[0.1,0.2]/", "y[0.2,0.6]/", "y[0.6,1]/", "y[1,3]/", "y[3,5]/"]
nMassbins = 12 # 52 ~ 500 GeV. See afb_mbin[afb_mbinnum+1] in dybAnalyzer.h
xsec_unc = {
    #"dy" : [1.7, -1.8],
    "wjets" : [3.8, -3.8],
    "ttll" : [4.8, -6.1],
    "ttlj" : [4.8, -6.1],
    "tw" : [5.4, -5.4],
    "stt" : [4.2, -3.6],
    "sts" : [3.9, -3.5],
    "ww" : [2.5, -2.2],
    "wz" : [6.1, -6.1],
    "zz" : [4.9, -4.9],
    "aa" : [30, -30],
    "qcd" : [30, -30],
}
lumi_unc = {
    "2016"    : [0.985, 0.985, 0.0,   0.0],
    "2017"    : [0.0,   0.0,   0.378, 0.0],
    "2018"    : [0.0,   0.0,   0.0,   0.439],
    "17and18" : [0.0,   0.0,   0.626, 0.582],
    "161718"  : [0.742, 0.742, 0.369, 0.414],
}
systematics = {
    # Stat
    "stat_Data" : [["stat_Data"]],
    "stat_MC" :   [["stat_MC"]],
    # SYS
    "JES" :       [["_jet_scale"+updown+era for updown in ["_up", "_down"]] for era in [":2016preVFP", ":2016postVFP", ":2017", ":2018"]],
    "JER" :       [["_jet_smear"+updown+era for updown in ["_up", "_down"]] for era in [":2016preVFP", ":2016postVFP", ":2017", ":2018"]],
    "Prefiring" : [["_prefireweight_up", "_prefireweight_down"]],
    "PU" :        [["_PUweight_up", "_PUweight_down"]],
    "PUIDSF" :    [["_PUjetSF_up", "_PUjetSF_down"]],
    "btagSF" :    [["_btagSF"+flavor+corr] for flavor in ["_h", "_l"] for corr in ["corr", "uncorr:2016preVFP", "uncorr:2016postVFP", "uncorr:2017", "uncorr:2018"]],
    "bChargeSF" : [["_bChargeSF1"+updown+bCh] for updown in ["_up", "_down"] for bCh in ["0", "1", "2", "3", "4", "5"]],
    "Lumi" :      [["lumi_"+eras+updown for updown in ["_up", "_down"]] for eras in lumi_unc.keys()],
    "CFSF" :      [["_CFSF_up", "_CFSF_down"]],
    # PDFSYS
    "Scales" :    [["_scalevariation%d:%s" % (i, j) for i in [0, 1, 2, 3, 4, 6, 8]] for j in ["dy", "tt"]], # 0=(1, 1), 5=(2, 0.5), and 7=(0.5, 2)
    "AlphaS" :    [["_alphaS_up", "_alphaS_down"]],
    "ISR" :       [["_ISR_up", "_ISR_down"]],
    "FSR" :       [["_FSR_up", "_FSR_down"]],
    #"PDF" :       [["_pdf%d" % i] for i in range(100)],
    "PDF0" :      [["_pdf%d" % i] for i in range(20)],
    "PDF1" :      [["_pdf%d" % i] for i in range(20,40)],
    "PDF2" :      [["_pdf%d" % i] for i in range(40,60)],
    "PDF3" :      [["_pdf%d" % i] for i in range(60,80)],
    "PDF4" :      [["_pdf%d" % i] for i in range(80,100)],
    "Bkgs" :      [["norm_"+bkgs+updown for updown in ["_up", "_down"]] for bkgs in xsec_unc.keys()],
    "Toppt" :     [["_noToppt"]],
    "Zpt" :       [["_Zpt"]],
    "Weak" :      [["_noWeak"]],
    # LEPSYS
    "MuTracking" : [["_muonTrackingeffSF_s%dm0" % i] for i in [1, 2, 3, 4, 5, 6, 9, 10, 14, 15, 16]] + [["_muonTrackingeffSF_s%dm0" % i, "_muonTrackingeffSF_s%dm1" % i] for i in [7, 8, 11, 12, 13]],
    "MuRECO" :     [["_muonRECOeffSF_s%dm0" % i] for i in [1, 2, 3, 4, 5, 6, 8, 9, 10, 14, 15, 16]] + [["_muonRECOeffSF_s%dm0" % i, "_muonRECOeffSF_s%dm1" % i] for i in [7, 11, 12, 13]],
    "MuID" :       [["_muonIDeffSF_s%dm0" % i] for i in [1, 2, 3, 4, 5, 6, 9, 10, 14, 15, 16, 17]] + [["_muonIDeffSF_s%dm0" % i, "_muonIDeffSF_s%dm1" % i] for i in [7, 8, 11, 12, 13]],
    "MuTrigger" :  [["_muonTriggereffSF_s%dm0" % i] for i in [1, 2, 3, 4, 5, 6, 9, 10, 14, 15]] + [["_muonTriggereffSF_s%dm0" % i, "_muonTriggereffSF_s%dm1" % i] for i in [7, 8, 11, 12, 13]],
    "ElRECO" :     [["_electronRECOeffSF_s%dm0" % i] for i in [1, 2, 3, 4, 5, 6, 9, 10, 14, 15, 16]] + [["_electronRECOeffSF_s%dm0" % i, "_electronRECOeffSF_s%dm1" % i] for i in [7, 8, 11, 12, 13]],
    "ElID" :       [["_electronIDeffSF_s%dm0" % i] for i in [1, 2, 3, 4, 5, 6, 9, 10, 14, 15, 16, 17, 18]] + [["_electronIDeffSF_s%dm0" % i, "_electronIDeffSF_s%dm1" % i] for i in [7, 8, 11, 12, 13]],
    "ElTrigger" :  [["_electronTriggereffSF_s%dm0" % i] for i in [1, 2, 3, 4, 5, 6, 9, 10, 14, 15, 16]] + [["_electronTriggereffSF_s%dm0" % i, "_electronTriggereffSF_s%dm1" % i] for i in [7, 8, 11, 12, 13]],
    "EffStatReplica" :    [["_lepeffSF_stat%d" % i] for i in range(20)],
    "RoccoRStatReplica" : [["_MuonMomentum_s1m%d" % i] for i in range(40)],
    "RoccoR" :            [["_MuonMomentum_s%dm0" % i for i in range(2, 6)]],
    "AepcorStatReplica" : [["_ElectronEnergy_s1m%d" % i] for i in range(40)],
    "Aepcor" :            [["_ElectronEnergy_s%dm0" % i for i in range(2, 9)]],
}

def findIntersection(x, y, sigma=1):
    miny = 999
    minx = 0
    x1 = 0
    x2 = 0
    for i in range(len(y)):
        if y[i] < miny:
            miny = y[i] # Find the minimum chi2
            minx = x[i] # and sintheta2

    diff1 = 999
    diff2 = 999
    for i in range(len(y)):
        if x[i] < minx and abs(y[i] - (miny + sigma**2)) < diff1:
            x1 = x[i]
            diff1 = abs(y[i] - (miny + sigma**2))
        if x[i] > minx and abs(y[i] - (miny + sigma**2)) < diff2:
            x2 = x[i]
            diff2 = abs(y[i] - (miny + sigma**2))

    return x1, x2, minx, miny

def calPrecision(chi2s_stat, chi2s_total, nametag=""):
    sin2ws = sin2w_values#[]
    #for i in sin2w_indice:
    #    sin2ws.append(sin2w_values[i])

    sigma = 1
    xmin = 0.22200
    xmax = 0.24200
    nxbin = int((xmax - xmin) / 0.00001 * 100)
    xicenter = int((0.23154 - xmin) / 0.00001 * 100)
    x = np.linspace(xmin, xmax, nxbin)

    # stat-only
    fit_stat = np.polyfit(sin2ws, chi2s_stat, 2)
    pol2_stat = np.poly1d(fit_stat)
    y_stat = pol2_stat(x)
    x1_stat, x2_stat, minx_stat, miny_stat = findIntersection(x, y_stat, sigma)
    y_stat -= miny_stat
    chi2s_stat -= miny_stat

    # stat+syst
    fit_total = np.polyfit(sin2ws, chi2s_total, 2)
    pol2_total = np.poly1d(fit_total)
    y_total = pol2_total(x)
    x1_total, x2_total, minx_total, miny_total = findIntersection(x, y_total, sigma)
    y_total -= miny_total
    chi2s_total -= miny_total

    central_stat = minx_stat
    unc_stat = (x2_stat - x1_stat) / 2
    central_total = minx_total
    unc_total = (x2_total - x1_total) / 2
    unc_syst = (unc_total**2 - unc_stat**2)**0.5

    plt.plot(sin2ws, chi2s_stat, 'o', color='black')
    plt.plot(sin2ws, chi2s_total, 'o', color='black')
    pol2_stat, = plt.plot(x, y_stat, color='darkviolet', label="stat-only")
    pol2_total, = plt.plot(x, y_total, color='blue', label="stat+syst")
    z = np.full(len(y_stat), sigma)
    plt.plot(x, z, color='red')

    ymax = max(y_stat[-1], y_stat[0], y_total[-1], y_total[0])
    #plt.title(r"$\chi^{2}$ Fitting ("+nametag+")", fontsize=15)
    plt.xlabel("$sin^{2}\\theta^{l}_{eff}$")
    plt.ylabel("$\Delta\chi^{2}$")
    plt.text(0.23514 + (0.23514 - x[0]) * -1.0, 0.3, "$\chi^{2}_{min}$ = %.3f, %.3f" % (miny_stat, miny_total))
    plt.text(0.23514 + (0.23514 - x[0]) * -0.84, ymax * 0.8, "$\sin^{2}\\theta^{l}_{eff}$ = %.5f $\\pm$ %.5f (stat-only)" % (central_stat, unc_stat))
    plt.text(0.23514 + (0.23514 - x[0]) * -0.84, ymax * 0.7, "$\sin^{2}\\theta^{l}_{eff}$ = %.5f $\\pm$ %.5f (stat) $\\pm$ %.5f (syst)" % (central_total, unc_stat, unc_syst))
    plt.text(0.23514 + (0.23514 - x[0]) * -0.65, ymax * 0.65, "= %.5f $\\pm$ %.5f (total)" % (central_total, unc_total))

    plt.legend(handles=[pol2_stat, pol2_total], loc='upper right')
    plt.grid()
    # CMS Style
    yRange = y_stat[0] - y_stat[xicenter]
    plt.text(0.23514 + (0.23514 - x[0]) * -1.07, ymax * 1.07, r"$\bf{CMS}$ Preliminary", fontsize=14)
    plt.text(0.23514 + (0.23514 - x[0]) * -0.42, ymax * 1.07, r"$\it{Working\ in\ progress}$", fontsize=10, color='red')
    plt.text(0.23514 + (0.23514 - x[0]) * 0.45, ymax * 1.12, r"$\bf{Run\ II}$", fontsize=12)
    plt.text(0.23514 + (0.23514 - x[0]) * 0.14, ymax * 1.07, r"138 fb$^{-1}$ (13 TeV)", fontsize=12)
    #plt.text(0.23514 + (0.23514 - x[0]) * 0.32, ymax * 1.12, r"$\bf{Run\ 2018}$", fontsize=12)
    #plt.text(0.23514 + (0.23514 - x[0]) * 0.14, ymax * 1.07, r"59.6 fb$^{-1}$ (13 TeV)", fontsize=12) # 19.5, 16.8, 42.1, 59.6

    if "no" not in nametag: plt.savefig("./precision_"+nametag+".pdf", dpi=300, bbox_inches='tight')
    plt.close()

    print(nametag, "sin2w central : ", central_total, "1 sigma (stat): ", unc_stat, "1 sigma (total): ", unc_total, "sin2w range : ", x1_total, x2_total)
    return unc_stat, unc_total

def getdAFB(AFB_ref, AFB):
    dAFB = []
    for iBin in range(nMassbins):
        dAFB.append(AFB_ref.GetBinContent(iBin + 1) - AFB.GetBinContent(iBin + 1))
        #print(iBin+1, "th Bin : nominal AFB = ", AFB_ref.GetBinContent(iBin + 1), "syst AFB = ", AFB.GetBinContent(iBin + 1), "diff = ", AFB_ref.GetBinContent(iBin + 1) - AFB.GetBinContent(iBin + 1))

    return np.array(dAFB)

## dAFBs[sin2w scenarios][chargeBins] - numpy 1D array with dimass bins
def getdAFBs_sin2w(channel):
    dAFBs = [[] for i in range(len(sin2w_indice))]
    dAFBs_full = []

    for sin in range(len(sin2w_indice)):
        #iSin = sin2w_indice[sin]
        iSin = sin2w_values[sin]
        dAFB_full = np.array([])
        for ch in range(len(chargeBins)):
            AFB_data = dyb.GetHist(0, channel+chargeBins[ch]+"AFBrecoil(x)", "")
            AFB_mc = dyb.GetHist(1, channel+chargeBins[ch]+"AFBrecoil(x)", "")
            #AFB_mc_sin2w = dyb.GetHist(1, channel+chargeBins[ch]+"AFBrecoil(x)", "suffix:_sthw2_%i:dy" % iSin)
            AFB_mc_sin2w = dyb.GetHist(1, channel+chargeBins[ch]+"AFBrecoil(x)", "suffix:_Recoil_weakNLOHO_s2eff_%i:dy" % (iSin * 1e5))
            for i in range(1, AFB_data.GetNbinsX() + 1):
                if AFB_data.GetBinCenter(i) < 120: AFB_data.SetBinContent(i, AFB_mc.GetBinContent(i)) # Blinded under 120 GeV
            dAFB = getdAFB(AFB_data, AFB_mc_sin2w)
            #dAFB = getdAFB(AFB_mc, AFB_mc_sin2w)
            dAFBs[sin].append(dAFB)
            if ch != 0: dAFB_full = np.append(dAFB_full, dAFB)
            print("\n dAFBs of "+channel+chargeBins[ch], (", sin2w_variation : %d, trace(dAFB) of " % iSin), (dAFB * dAFB).sum(), ", len(dAFB) = ", len(dAFB))
            print(dAFB)
        dAFBs_full.append(dAFB_full)
        print("\n dAFBs of "+channel, (", sin2w_variation : %d, trace(dAFB_full) of " % iSin), (dAFB_full * dAFB_full).sum(), ", len(dAFB_full) = ", len(dAFB_full))
        print(dAFB_full)

    return np.array(dAFBs), np.array(dAFBs_full)

## dAFBs{systkey}[term][syst][chargeBins] - numpy 1D array with dimass bins
def getdAFBs_syst(channel, dAFBs, dAFBs_full, missingSyst=""):
    for systkey, list_syst in systematics.items():
        if missingSyst != "" and systkey != missingSyst: continue
        dAFBs[systkey] = []
        dAFBs_full[systkey] = []
        for term in range(len(list_syst)):
            dAFBs[systkey].append([])
            dAFBs_full[systkey].append([])
            for syst in range(len(list_syst[term])):
                dAFBs[systkey][term].append([])
                dAFBs_full[systkey][term].append([])
                dAFB_full = np.array([])
                for ch in range(len(chargeBins)):
                    print("\n channel+chargeBin : "+channel+chargeBins[ch])
                    dAFB = []
                    if "stat_Data" in systkey:
                        AFB_data_nominal = dyb.GetHist(0, channel+chargeBins[ch]+"AFBrecoil(x)", "")
                        for iBin in range(nMassbins):
                            dAFB.append(AFB_data_nominal.GetBinError(iBin + 1))
                        dAFBs[systkey][term][syst].append(np.array(dAFB))
                        print("trace(dAFB_stat_Data) = ", (dAFBs["stat_Data"][term][syst][ch] * dAFBs["stat_Data"][term][syst][ch]).sum(), ", len(dAFB_stat_Data) = ", len(dAFBs["stat_Data"][term][syst][ch]))
                        if ch != 0: dAFB_full = np.append(dAFB_full, dAFB)
                    elif "stat_MC" in systkey:
                        AFB_mc_nominal = dyb.GetHist(1, channel+chargeBins[ch]+"AFBrecoil(x)", "")
                        for iBin in range(nMassbins):
                            dAFB.append(AFB_mc_nominal.GetBinError(iBin + 1))
                        dAFBs[systkey][term][syst].append(np.array(dAFB))
                        print("trace(dAFB_stat_MC) = ", (dAFBs["stat_MC"][term][syst][ch] * dAFBs["stat_MC"][term][syst][ch]).sum(), ", len(dAFB_stat_MC) = ", len(dAFBs["stat_MC"][term][syst][ch]))
                        if ch != 0: dAFB_full = np.append(dAFB_full, dAFB)
                    else: # Systematics
                        AFB_mc_nominal = dyb.GetHist(1, channel+chargeBins[ch]+"AFBrecoil(x)", "")
                        syststr = "suffix:"+list_syst[term][syst]
                        if "norm" in list_syst[term][syst]:
                            for process, uncs in xsec_unc.items():
                                if process in list_syst[term][syst]: syststr = "scale:%.3f:%s" % (1 + uncs[(0 if "up" in list_syst[term][syst] else 1)] * 0.01, process)
                        elif "lumi" in list_syst[term][syst]:
                            for era, uncs in lumi_unc.items():
                                if era in list_syst[term][syst]:
                                    unc = [1 + a * 0.01 * (1 if "up" in list_syst[term][syst] else -1) for a in uncs]
                                    syststr = "scale:%.3f:2016preVFP scale:%.3f:2016postVFP scale:%.3f:2017 scale:%.3f:2018" % (unc[0], unc[1], unc[2], unc[3])
                        print("syststr = "+syststr)
                        AFB_mc_syst = dyb.GetHist(1, channel+chargeBins[ch]+"AFBrecoil(x)", syststr)
                        dAFB_syst = getdAFB(AFB_mc_nominal, AFB_mc_syst)
                        dAFBs[systkey][term][syst].append(dAFB_syst)
                        print(list_syst[term][syst]+", trace(dAFB_syst) = ", (dAFBs[systkey][term][syst][ch] * dAFBs[systkey][term][syst][ch]).sum(), ", len(dAFB_syst) = ", len(dAFBs[systkey][term][syst][ch]))
                        if ch != 0: dAFB_full = np.append(dAFB_full, dAFB_syst)
                    #print(dAFBs[systkey][term][syst][ch])

                dAFBs_full[systkey][term][syst].append(dAFB_full)
                print(list_syst[term][syst]+", trace(dAFB_syst_full) = ", (dAFBs_full[systkey][term][syst][0] * dAFBs_full[systkey][term][syst][0]).sum(), ", len(dAFB_syst_full) = ", len(dAFBs_full[systkey][term][syst][0]))
                #print(dAFBs_full[systkey][term][syst])

    return dAFBs, dAFBs_full

def calChi2sWithCov(dAFBs_sin2w, dAFBs_syst, chargeBin=0, statOnly=False, N_1=""):
    dim = dAFBs_sin2w.shape[1] # 30 or 180
    cov = np.zeros((dim, dim))
    for systkey, list_syst in dAFBs_syst.items():
        if systkey == N_1: continue # For (N-1) syst uncertainties
        elif statOnly and "stat_" not in systkey: continue # For Stat-only uncertainties
        else:
            cov_term = np.zeros((dim, dim))
            for term in range(len(list_syst)):
                dAFBs_trace = -1
                cov_syst_bigger = np.zeros((dim, dim))
                for syst in range(len(list_syst[term])):
                    trace = (dAFBs_syst[systkey][term][syst][chargeBin] * dAFBs_syst[systkey][term][syst][chargeBin]).sum()
                    #print(syst, "trace(cov_syst_cand) = ", trace)
                    if trace > dAFBs_trace:
                        dAFBs_trace = trace
                        cov_syst_bigger = np.outer(dAFBs_syst[systkey][term][syst][chargeBin], dAFBs_syst[systkey][term][syst][chargeBin])
                        if "stat_" in systkey: cov_syst_bigger = np.diag(np.diag(cov_syst_bigger))
                #print(cov_syst_bigger)
                #print(term, "trace(cov_syst_bigger) = ", np.trace(cov_syst_bigger))
                cov_term += cov_syst_bigger
            print("trace(cov_term) of "+systkey+" = ", np.trace(cov_term), "len(Csyst) = ", len(cov_term))
            if "Replica" in systkey:
                cov_term = cov_term / len(list_syst)
                print("cov_term divided by ", len(list_syst), ", and then trace(cov_term) of "+systkey+" = ", np.trace(cov_term), "len(Csyst) = ", len(cov_term))
            cov += cov_term
            #for i in range(dim):
            #    print("cov_term : ", i, "th bin = ", cov_term[i][i]**0.5)
            #for i in range(dim):
            #    print("cov : ", i, "th bin = ", cov[i][i]**0.5)

    inv_cov = np.linalg.inv(cov)
    print("\n cov : ")
    #print(cov)
    print("\n inv_cov : ")
    #print(inv_cov)
    print("")
    print("trace(cov) = ", np.trace(cov), ", trace(inv_cov) = ", np.trace(inv_cov), "len(inv_cov) = ", len(inv_cov))
    #print(np.dot(cov, inv_cov))

    chi2s = []
    for i in range(len(sin2w_indice)):
        dAFB = dAFBs_sin2w[i].reshape(dim, 1)
        chi2 = np.linalg.multi_dot([dAFB.transpose(), inv_cov, dAFB])
        chi2s.append(chi2[0][0])

    print("chi2s = ", chi2s)
    return chi2s

if __name__=="__main__":
    channel = "[em][em]201[678][ab]?/"
    #channel = "mm201[678][ab]?/"
    npz_files_tag = "_12bins"

    ## dAFBs_sin2w[sin2w scenarios][chargeBins]
    channel_path = channel.replace("?", "").replace("/", "").replace("[", "").replace("]", "")
    dAFBs_sin2w_npz = channel_path+"_dAFBs_sin2w"+npz_files_tag+".npz"
    if not os.path.exists(dAFBs_sin2w_npz):
        dAFBs_sin2w, dAFBs_sin2w_full = getdAFBs_sin2w(channel)
        np.savez(dAFBs_sin2w_npz, X = dAFBs_sin2w, Y = dAFBs_sin2w_full)
        print("New "+dAFBs_sin2w_npz+" is saved")
    else:
        npz_sin2w = np.load(dAFBs_sin2w_npz)
        print(dAFBs_sin2w_npz+" is loaded")
        dAFBs_sin2w = npz_sin2w['X']
        dAFBs_sin2w_full = npz_sin2w['Y']

    ## dAFBs_syst{systematics}[chargeBins]
    dAFBs_syst = {}
    dAFBs_syst_full = {}
    dAFBs_syst_npz = channel_path+"_dAFBs_syst"+npz_files_tag+"_hadded.npz"
    if not os.path.exists(dAFBs_syst_npz):
        dAFBs_syst, dAFBs_syst_full = getdAFBs_syst(channel, dAFBs_syst, dAFBs_syst_full)
        np.savez(dAFBs_syst_npz, X = dAFBs_syst, Y = dAFBs_syst_full)
        print("New "+dAFBs_syst_npz+" is saved")
    else:
        npz_syst = np.load(dAFBs_syst_npz, allow_pickle=True)
        print(dAFBs_syst_npz+" is loaded")
        dAFBs_syst = npz_syst['X'][()]
        dAFBs_syst_full = npz_syst['Y'][()]
        for syst in systematics.keys():
            if syst not in dAFBs_syst:
                print(syst+" is missing in "+dAFBs_syst_npz)
                dAFBs_syst, dAFBs_syst_full = getdAFBs_syst(channel, dAFBs_syst, dAFBs_syst_full, syst)
                np.savez(dAFBs_syst_npz, X = dAFBs_syst, Y = dAFBs_syst_full)
                print(syst+" is added in "+dAFBs_syst_npz+", "+dAFBs_syst_npz+" is updated")

    ## Chi2s (1D, 2D and full2D)
    chi2s_1D_stat = calChi2sWithCov(dAFBs_sin2w[:, 0], dAFBs_syst, 0, True)
    chi2s_1D_total = calChi2sWithCov(dAFBs_sin2w[:, 0], dAFBs_syst, 0)
    chi2s_2D_stat = [0] * len(sin2w_indice)
    chi2s_2D_total = [0] * len(sin2w_indice)
    for i in range(1, len(chargeBins)):
        chi2s_stat = calChi2sWithCov(dAFBs_sin2w[:, i], dAFBs_syst, i, True)
        chi2s_2D_stat = [chi2s_2D_stat[j] + chi2s_stat[j] for j in range(len(chi2s_2D_stat))]
        chi2s_total = calChi2sWithCov(dAFBs_sin2w[:, i], dAFBs_syst, i)
        chi2s_2D_total = [chi2s_2D_total[j] + chi2s_total[j] for j in range(len(chi2s_2D_total))]
    chi2s_full2D_stat = calChi2sWithCov(dAFBs_sin2w_full, dAFBs_syst_full, 0, True)
    chi2s_full2D_total = calChi2sWithCov(dAFBs_sin2w_full, dAFBs_syst_full, 0)

    unc_1D_stat, unc_1D_total = calPrecision(chi2s_1D_stat, chi2s_1D_total, "1D")
    unc_2D_stat, unc_2D_total = calPrecision(chi2s_2D_stat, chi2s_2D_total, "2D")
    unc_full2D_stat, unc_full2D_total = calPrecision(chi2s_full2D_stat, chi2s_full2D_total, "full2D")

    print("1D precision = %.5f (stat) pm %.5f (syst) = %.5f (total)" % (unc_1D_stat, (unc_1D_total**2 - unc_1D_stat**2)**0.5, unc_1D_total))
    print("2D precision = %.5f (stat) pm %.5f (syst) = %.5f (total)" % (unc_2D_stat, (unc_2D_total**2 - unc_2D_stat**2)**0.5, unc_2D_total))
    print("full2D precision = %.5f (stat) pm %.5f (syst) = %.5f (total)" % (unc_full2D_stat, (unc_full2D_total**2 - unc_full2D_stat**2)**0.5, unc_full2D_total))

    ## (N-1) stat + syst uncertainties
    n = 0
    for syst in dAFBs_syst.keys():
        if "stat_" in syst: continue
        chi2s_1D_total_N_1 = calChi2sWithCov(dAFBs_sin2w[:, 0], dAFBs_syst, 0, False, syst)
        chi2s_2D_total_N_1 = [0] * len(sin2w_indice)
        for i in range(1, len(chargeBins)):
            chi2s = calChi2sWithCov(dAFBs_sin2w[:, i], dAFBs_syst, i, False, syst)
            chi2s_2D_total_N_1 = [chi2s_2D_total_N_1[j] + chi2s[j] for j in range(len(chi2s_2D_total_N_1))]
        chi2s_full2D_total_N_1 = calChi2sWithCov(dAFBs_sin2w_full, dAFBs_syst_full, 0, False, syst)

        unc_1D_stat, unc_1D_total_N_1 = calPrecision(chi2s_1D_stat, chi2s_1D_total_N_1, ("1D_%i_" % n)+"no"+syst)
        unc_2D_stat, unc_2D_total_N_1 = calPrecision(chi2s_2D_stat, chi2s_2D_total_N_1, ("2D_%i_" % n)+"no"+syst)
        unc_full2D_stat, unc_full2D_total_N_1 = calPrecision(chi2s_full2D_stat, chi2s_full2D_total_N_1, ("full2D_%i_" % n)+"no"+syst)

        print("Impact of "+syst+" on 1D precision = %.5f " % (unc_1D_total**2 - unc_1D_total_N_1**2)**0.5)
        print("Impact of "+syst+" on 2D precision = %.5f " % (unc_2D_total**2 - unc_2D_total_N_1**2)**0.5)
        print("Impact of "+syst+" on full2D precision = %.5f " % (unc_full2D_total**2 - unc_full2D_total_N_1**2)**0.5)
        n += 1
