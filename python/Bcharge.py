import math
import numpy as np
import ROOT
ROOT.gROOT.ProcessLine('#include"BBPlotter.cc"')
ROOT.gROOT.ProcessLine('#include"ttljPlotter.cc"')
ROOT.Plotter.SetupStyle()

p=ROOT.BBPlotter("eff")
ttlj = ROOT.ttljPlotter("correct_ttlj wrong_ttlj unmatched_ttlj")
ttlj_gen = ROOT.ttljPlotter("ttlj_2b")
forNorm = ROOT.ttljPlotter("data mc")

xsec_unc = {
    "dy" : [1.7, -1.8],
    "wjets" : [3.8, -3.8],
    "ttll" : [4.8, -6.1],
    "ttlj" : [4.8, -6.1],
    "ttjj" : [4.8, -6.1],
    "tw" : [5.4, -5.4],
    "stt" : [4.2, -3.6],
    "sts" : [3.9, -3.5],
    "ww" : [2.5, -2.2],
    "wz" : [6.1, -6.1],
    "zz" : [4.9, -4.9],
    #"aa" : [30, -30],
    "qcd" : [30, -30],
}
systematics = [
    ["_donorm"],
    ["_noSelQ"],
    ["_lumi_up", "_lumi_down"],
    ["_jetpt25", "_jetpt55"],
    ["_jeteta5", "_jeteta1p5"],
    ["_jet_scale_up", "_jet_scale_down"], ["_jet_smear_up", "_jet_smear_down"],
    #["_prefireweight_up", "_prefireweight_down"], ["_PUweight_up", "_PUweight_down"], ["_PUjetSF_up", "_PUjetSF_down"],
    #["_btagSF_hup", "_btagSF_hdown"], ["_btagSF_lup", "_btagSF_ldown"],
    #["_btagSF_hcorr"], ["_btagSF_huncorr2016a"], ["_btagSF_huncorr2016b"], ["_btagSF_huncorr2017"], ["_btagSF_huncorr2018"],
    #["_btagSF_lcorr"], ["_btagSF_luncorr2016a"], ["_btagSF_luncorr2016b"], ["_btagSF_luncorr2017"], ["_btagSF_luncorr2018"],
    ["_scalevariation%d" % i for i in [0, 1, 2, 3, 4, 6, 8]], # 0=(1, 1), 5=(2, 0.5), and 7=(0.5, 2)
    ["_FSR_up", "_FSR_down"],
] + [["norm_"+bkgs+updown for updown in ["_up", "_down"]] for bkgs in xsec_unc.keys()]

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

### For Liklihood ratio method
def getfc(channel, chargeBin="", syst=""):
    suffix = "" if "norm" in syst or "lumi" in syst else "suffix:"+syst

    correct = ttlj.GetHist(0, channel+"/reco[bB]Charge"+chargeBin, suffix).Integral()
    wrong = ttlj.GetHist(1, channel+"/reco[bB]Charge"+chargeBin, suffix).Integral()
    fc = correct / (correct + wrong)
    fc_e = fc * (1 - fc) / (correct + wrong)

    data = forNorm.GetHist(0, channel+"/reco[bB]Charge"+chargeBin, suffix).Integral()
    mc = forNorm.GetHist(1, channel+"/reco[bB]Charge"+chargeBin, suffix).Integral()
    print("getfc function : fc = ", fc, ", fc_e = ", fc_e, ", data = ", data, ", mc = ", mc, ", norm = ", data / mc)

    return fc, fc_e, (data / mc)

def calc_withLR(fc, fp, fm):
    ap = (fc * (fp + fc -1) + (1 - fc) * (fm + fc -1)) / (2 * fc - 1)
    am = ((1 - fc) * (fp + fc -1) + fc * (fm + fc -1)) / (2 * fc - 1)

    return ap, am

def calcWithCov_withLR(fcs, fps, fms, fc_stat, fp_stat, fm_stat):
    cov_stat = np.zeros((2, 2))
    cov = np.zeros((2, 2))
    print("fc", fcs)
    fc_nominal = fcs[-1][0]
    fp_nominal = fps[-1][0]
    fm_nominal = fms[-1][0]
    nominal = np.array(calc_withLR(fc_nominal, fp_nominal, fm_nominal))
    stat0 = np.array(calc_withLR(fc_nominal + fc_stat, fp_nominal, fm_nominal))
    stat1 = np.array(calc_withLR(fc_nominal, fp_nominal + fp_stat, fm_nominal))
    stat2 = np.array(calc_withLR(fc_nominal, fp_nominal, fm_nominal + fm_stat))
    cov_stat += np.outer(stat0 - nominal, stat0 - nominal)
    cov_stat += np.outer(stat1 - nominal, stat1 - nominal)
    cov_stat += np.outer(stat2 - nominal, stat2 - nominal)
    print("nomianl : ", nominal, ", stat unc : ", [math.sqrt(cov_stat[0][0]), math.sqrt(cov_stat[1][1])])
    print("cov", cov_stat)

    values, vectors = np.linalg.eig(cov_stat)
    vector0 = np.array([vectors[0][0], vectors[1][0]]) * values[0] ** 0.5
    vector1 = np.array([vectors[0][1], vectors[1][1]]) * values[1] ** 0.5
    print("stat vect0", vector0)
    print("stat vect1", vector1)
    print(np.outer(vector0, vector0) + np.outer(vector1, vector1))

    cov += cov_stat
    if len(fcs) == 1: return nominal, cov_stat, cov

    print("\nSystematics")
    for systs in range(len(fcs) - 1):
        cov_syst_bigger = np.zeros((2, 2))
        for syst in range(len(fcs[systs])):
            dsyst = np.array(calc_withLR(fcs[systs][syst], fps[systs][syst], fms[systs][syst])) - nominal
            cov_syst = np.outer(dsyst, dsyst)
            print(systs, syst, ", syst : "+systematics[systs][syst]+", alpha : ", dsyst + nominal, ", d(alpha) : ", dsyst, ", trace(cov) : ", np.trace(cov_syst))
            print(cov_syst)
            if np.trace(cov_syst) > np.trace(cov_syst_bigger): cov_syst_bigger = cov_syst ## Choose the cov_systematic with the larger trace
        cov += cov_syst_bigger

    print("nomianl : ", nominal, ", stat+syst unc : ", [math.sqrt(cov[0][0]), math.sqrt(cov[1][1])])
    print("cov", cov)

    values, vectors = np.linalg.eig(cov)
    vector0 = np.array([vectors[0][0], vectors[1][0]]) * values[0] ** 0.5
    vector1 = np.array([vectors[0][1], vectors[1][1]]) * values[1] ** 0.5
    print("stat+syst vect0", vector0)
    print("stat+syst vect1", vector1)
    print(np.outer(vector0, vector0) + np.outer(vector1, vector1))

    return nominal, cov_stat, cov

def GetAccuracy_withLR(ientry, channel, chargeBin="", option=""):
    allsysts = [[""]]
    if "syst" in option: allsysts = systematics + allsysts # Nominal last

    print("chargeBin : "+chargeBin+", ientry = ", ientry, ", option = "+option+"\n")
    fc, fp, fm = 0, 0, 0
    fc_e, fp_e, fm_e = -1, -1, -1
    fcs, fps, fms = [[]], [[]], [[]]

    for systs in range(len(allsysts)):
        for syst in allsysts[systs]:
            fc, fc_stat, norm = getfc(channel, chargeBin, syst)
            if "donorm" not in syst: norm = 1.
            if "lumi_up" in syst: norm *= (100 + 1.616477652180815) / 100
            elif "lumi_down" in syst: norm *= (100 - 1.616477652180815) / 100
            a = ROOT.ttljPlotter("data_sub ttlj", norm)

            suffix = "" if "norm" in syst or "lumi" in syst else "suffix:"+syst
            if "norm" in syst:
                for process, uncs in xsec_unc.items():
                    if process in syst: suffix = "scale:%.3f:%s" % (1 + uncs[(0 if "up" in syst else 1)] * 0.01, process)

            hp = a.GetHist(ientry, channel+"/recoBCharge"+chargeBin, suffix)
            hm = a.GetHist(ientry, channel+"/recobCharge"+chargeBin, suffix)
            print("syst = "+syst+", suffix = "+suffix+", hp = ", hp.Integral(), ", fpxhp = ", hp.GetBinContent(2), ", hm = ", hm.Integral(), ", fmxhm = ", hm.GetBinContent(1), ", fc = ", fc, ", norm = ", norm)
            hp.Scale(1. / hp.Integral())
            hm.Scale(1. / hm.Integral())
            fcs[systs].append(fc)
            fps[systs].append(hp.GetBinContent(2))
            fp_stat = hp.GetBinError(2)
            fms[systs].append(hm.GetBinContent(1))
            fm_stat = hm.GetBinError(1)
        if systs < len(allsysts) - 1:
            fcs.append([])
            fps.append([])
            fms.append([])

    nominal, cov_stat, cov = calcWithCov_withLR(fcs, fps, fms, fc_stat, fp_stat, fm_stat)
    return nominal, cov_stat, cov

def GetTrueAccuracy_withLR(channel, chargeBin="", option=""):
    allsysts = [[""]]
    if "syst" in option: allsysts = allsysts + systematics # Nominal first

    ap, am, stat_ep, stat_em, ep, em = 0, 0, -1, -1, 0, 0
    for systs in range(len(allsysts)):
        syst_ep_bigger, syst_em_bigger = 0, 0
        for syst in allsysts[systs]:
            suffix = "" if "norm" in syst or "lumi" in syst else "suffix:"+syst
            hp = ttlj_gen.GetHist(0, channel+"/genBCharge"+chargeBin+"_L[pm]", suffix)
            hm = ttlj_gen.GetHist(0, channel+"/genbCharge"+chargeBin+"_L[pm]", suffix)
            hp.Scale(1. / hp.Integral())
            hm.Scale(1. / hm.Integral())

            if syst == "":
                ap = hp.GetBinContent(2)
                stat_ep = hp.GetBinError(2)
                am = hm.GetBinContent(1)
                stat_em = hm.GetBinError(1)
            else:
                dep = hp.GetBinContent(2) - ap
                dem = hm.GetBinContent(1) - am
                if (dep**2 + dem**2) > (syst_ep_bigger**2 + syst_em_bigger**2): ## Choose the dep, dem with the larger square sum (~ trace)
                    syst_ep_bigger = dep
                    syst_em_bigger = dem

        ep = (ep**2 + syst_ep_bigger**2)**0.5
        em = (em**2 + syst_em_bigger**2)**0.5

    ep = (ep**2 + stat_ep**2)**0.5
    em = (em**2 + stat_em**2)**0.5

    print("True Accuracy nominal = ", (ap, am), "stat error = ", (stat_ep, stat_em), "stat+syst error = ", (ep, em))
    return (ap, am), (stat_ep, stat_em), (ep, em)

def DrawAccuracy_withLR(channels, Xrange=(0.6, 0.67), tag="", option=""):
    gdata = [ROOT.TGraphErrors(), ROOT.TGraphErrors()]
    gsim = [ROOT.TGraphErrors(), ROOT.TGraphErrors()]
    gtrue = [ROOT.TGraphErrors(), ROOT.TGraphErrors()]
    gdata_tot_unc = [ROOT.TGraphErrors(), ROOT.TGraphErrors()]
    gsim_tot_unc = [ROOT.TGraphErrors(), ROOT.TGraphErrors()]
    gtrue_tot_unc = [ROOT.TGraphErrors(), ROOT.TGraphErrors()]

    chBins = ["_[0-3]", "_0", "_1", "_2", "_3", "_4", "_5"]
    ptBins = ["_pt0", "_pt1", "_pt2", "_pt3", "_pt4"]
    for i in range(len(channels)):
        channel = channels[i]
        print("\n\n@@@ Channel : "+channel+" started @@@\n")

        chargeBin = ""
        for Bin in chBins:
            if Bin in channel:
                channel = channel.replace(Bin, "")
                chargeBin = Bin.replace("_", "")
                break
        for Bin in ptBins:
            if Bin in channel:
                channel = channel.replace(Bin, "")
                chargeBin += Bin
                break

        alphas = []
        betas = []
        alphas_error = []
        betas_error = []
        value, cov_stat, cov = GetAccuracy_withLR(0, channel, chargeBin, option)
        for j in range(2):
            gdata[j].SetPoint(i, value[j], 2 * i + 1 - j + 0.5)
            gdata[j].SetPointError(i, cov_stat[j][j]**0.5, 0)
            gdata_tot_unc[j].SetPoint(i, value[j], 2 * i + 1 - j + 0.5)
            gdata_tot_unc[j].SetPointError(i, cov[j][j]**0.5, 0)
            alphas.append(value[j])
            betas.append(1 - value[j])
            alphas_error.append(cov_stat[j][j]**0.5)
            betas_error.append(cov_stat[j][j]**0.5)

        value, cov_stat, cov = GetAccuracy_withLR(1, channel, chargeBin, option)
        for j in range(2):
            gsim[j].SetPoint(i, value[j], 2 * i + 1 - j + 0.4)
            gsim[j].SetPointError(i, cov_stat[j][j]**0.5, 0)
            gsim_tot_unc[j].SetPoint(i, value[j], 2 * i + 1 - j + 0.4)
            gsim_tot_unc[j].SetPointError(i, cov[j][j]**0.5, 0)
            alphas[j] *= 1. / value[j]
            betas[j] *= 1. / (1 - value[j])
            alphas_error[j] *= 1. / (value[j])
            betas_error[j] *= 1. / (1 - value[j])

        value, stat_e, e = GetTrueAccuracy_withLR(channel, chargeBin, option)
        for j in range(2):
            gtrue[j].SetPoint(i, value[j], 2 * i + 1 - j + 0.3)
            gtrue[j].SetPointError(i, stat_e[j], 0)
            gtrue_tot_unc[j].SetPoint(i, value[j], 2 * i + 1 - j + 0.3)
            gtrue_tot_unc[j].SetPointError(i, e[j], 0)

        print("alpha ratio : ", alphas, "+-", alphas_error, ", (1-alpha) ratio : ", betas, "+-", betas_error)

    c = ROOT.gROOT.MakeDefCanvas()
    c.SetLeftMargin(0.2)
    c.SetRightMargin(0.05)
    c.SetBottomMargin(0.1)

    hframe = ROOT.TH2D("hframe", "", 100, Xrange[0], Xrange[1], len(channels) * 2 + 2, 0, len(channels) * 2 + 2);
    is_era_compared = True if [a for a in channels if "201[678][ab]?" in a] and [a for a in channels if "2018" in a] else False
    is_lep_compared = True if [a for a in channels if "[em]" in a] and [a for a in channels if "e201" in a] else False
    #is_nPV_compared = True if [a for a in channels if "Run2s" in a] and [a for a in channels if "2018" in a] else False
    #is_pt_compared = False
    is_LeftMargin_wide = False
    for i in range(len(channels)):
        title = channels[i]
        title = title.replace("201[678][ab]?", "Run2")
        if is_era_compared:
            title = title.replace("2016a", "2016a")
            title = title.replace("2016b", "2016b")
            title = title.replace("2016[ab]?", "2016")
            title = title.replace("2017", "2017")
            title = title.replace("2018", "2018")
        else:
            title = title.replace("2016a", "")
            title = title.replace("2016b", "")
            title = title.replace("2016[ab]?", "")
            title = title.replace("2017", "")
            title = title.replace("2018", "")

        if is_lep_compared:
            title = title.replace("[em]", "l")
            title = title.replace("e", "e")
            title = title.replace("m", "#mu")
        else: title = title.replace("[em]", "")

        title = title.replace("Run2F", "nPV #leq 10")
        title = title.replace("Run2S", "nPV#in(10,20]")
        title = title.replace("Run2L", "nPV#in(20,30]")
        title = title.replace("Run2M", "nPV#in(30,40]")
        title = title.replace("Run2H", "nPV#in(40,50]")
        title = title.replace("Run2V", "50 < nPV")

        title = title.replace("Run2_pt0", "pt(b)#in[25,35)")
        title = title.replace("Run2_pt1", "pt(b)#in[35,50)")
        title = title.replace("Run2_pt2", "pt(b)#in[50,80)")
        title = title.replace("Run2_pt3", "pt(b)#in[80,120)")
        title = title.replace("Run2_pt4", "120 #leq pt(b)")

        if not is_era_compared: title = title.replace("Run2", "")

        chargeBin = ""
        for Bin in chBins:
            if Bin in title:
                title = title.replace(Bin, "")
                chargeBin = Bin.replace("_", "")
                break

        if title != "":
            is_LeftMargin_wide = True
            title = " ("+title+")"
        alpha = "#alpha"
        if chargeBin != "":
            if "[0-3]" in chargeBin: alpha = "#alpha_{ j}"
            elif "4" in chargeBin: alpha = "#alpha_{ #mu}"
            elif "5" in chargeBin: alpha = "#alpha_{ e}"
            else: alpha = "#alpha_{ j,"+chargeBin+"}"
        hframe.GetYaxis().SetBinLabel(i * 2 + 1, alpha+"^{ #minus}"+title)
        hframe.GetYaxis().SetBinLabel(i * 2 + 2, alpha+"^{ #plus}"+title)

    if is_LeftMargin_wide == True: c.SetLeftMargin(0.08)
    hframe.SetStats(0)
    hframe.Draw()
    hframe.GetXaxis().SetTitle("Accuracy")
    hframe.GetYaxis().SetLabelSize(hframe.GetYaxis().GetLabelSize() * 1.5)

    leg = ROOT.TLegend(c.GetLeftMargin() + 0.01, 0.79, 0.99 - c.GetRightMargin(), 0.89)
    leg.AddEntry(gdata[0], "data")
    leg.AddEntry(gsim[0], "simulation")
    leg.AddEntry(gtrue[0], "simulation (truth)")
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

    # For stat+syst uncertainties
    gdata_tot_unc[0].SetMarkerStyle(20)
    gdata_tot_unc[0].SetMarkerSize(0.8)
    gdata_tot_unc[0].SetMarkerColor(1)
    gdata_tot_unc[0].SetLineColor(1)
    gdata_tot_unc[0].Draw("same p")

    gdata_tot_unc[1].SetMarkerStyle(24)
    gdata_tot_unc[1].SetMarkerSize(0.8)
    gdata_tot_unc[1].SetMarkerColor(1)
    gdata_tot_unc[1].SetLineColor(1)
    gdata_tot_unc[1].Draw("same p")

    gsim_tot_unc[0].SetMarkerStyle(20)
    gsim_tot_unc[0].SetMarkerSize(0.6)
    gsim_tot_unc[0].SetMarkerColor(2)
    gsim_tot_unc[0].SetLineColor(2)
    gsim_tot_unc[0].Draw("same p")

    gsim_tot_unc[1].SetMarkerStyle(24)
    gsim_tot_unc[1].SetMarkerSize(0.6)
    gsim_tot_unc[1].SetMarkerColor(2)
    gsim_tot_unc[1].SetLineColor(2)
    gsim_tot_unc[1].Draw("same p")

    gtrue_tot_unc[0].SetMarkerStyle(20)
    gtrue_tot_unc[0].SetMarkerSize(0.6)
    gtrue_tot_unc[0].SetMarkerColor(4)
    gtrue_tot_unc[0].SetLineColor(4)
    gtrue_tot_unc[0].Draw("same p")

    gtrue_tot_unc[1].SetMarkerStyle(24)
    gtrue_tot_unc[1].SetMarkerSize(0.6)
    gtrue_tot_unc[1].SetMarkerColor(4)
    gtrue_tot_unc[1].SetLineColor(4)
    gtrue_tot_unc[1].Draw("same p")

    # From DrawPreliminary() in Plotter.cc
    latex = ROOT.TLatex()
    latex.SetTextSize(0.04);
    latex.SetNDC();
    latex.SetTextAlign(11);
    leftmargin = c.GetLeftMargin();
    rightmargin = c.GetRightMargin();
    topmargin = c.GetTopMargin();
    latex.DrawLatex(0.01 + leftmargin, 1.01 - topmargin, "CMS #bf{#it{Preliminary}}")

    latex.SetTextSize(0.035);
    latex.SetTextAlign(31);
    if [a for a in channels if "201[678][ab]?" in a]:
      latex.DrawLatex(1 - rightmargin, 1.01 - topmargin, "138 fb^{-1} (13 TeV)")
      latex.DrawLatex(1 - rightmargin, 1.05 - topmargin, "Run II")
    elif [a for a in channels if "2018" in a]:
      latex.DrawLatex(1 - rightmargin, 1.01 - topmargin, "59.8 fb^{-1} (13 TeV)")
      latex.DrawLatex(1 - rightmargin, 1.05 - topmargin, "Run 2018")
    elif [a for a in channels if "2017" in a]:
      latex.DrawLatex(1 - rightmargin, 1.01 - topmargin, "41.5 fb^{-1} (13 TeV)")
      latex.DrawLatex(1 - rightmargin, 1.05 - topmargin, "Run 2017")
    elif [a for a in channels if "2016[ab]?" in a]:
      latex.DrawLatex(1 - rightmargin, 1.01 - topmargin, "36.3 fb^{-1} (13 TeV)")
      latex.DrawLatex(1 - rightmargin, 1.05 - topmargin, "Run 2016")
    elif [a for a in channels if "2016b" in a]:
      latex.DrawLatex(1 - rightmargin, 1.01 - topmargin, "16.8 fb^{-1} (13 TeV)")
      latex.DrawLatex(1 - rightmargin, 1.05 - topmargin, "Run 2016postVFP")
    elif [a for a in channels if "2016a" in a]:
      latex.DrawLatex(1 - rightmargin, 1.01 - topmargin, "19.5 fb^{-1} (13 TeV)")
      latex.DrawLatex(1 - rightmargin, 1.05 - topmargin, "Run 2016preVFP")

    latex.SetTextSize(0.035);
    latex.SetTextColor(2)
    latex.SetTextAlign(21)
    latex.DrawLatex(0.6 * (leftmargin - rightmargin) + 0.5, 1.01 - topmargin, "#it{Working in progress}")
    c.Update()

    #raw_input()
    c.hists = [gdata, gsim, gtrue, gdata_tot_unc, gsim_tot_unc, gtrue_tot_unc, hframe, leg, latex]
    nametag = ""
    if tag != "": nametag = nametag+"_"+tag
    if option != "": nametag = nametag+"_"+option
    c.SaveAs("Accuracies"+nametag+".png")
    c.SaveAs("Accuracies"+nametag+".pdf")

    return c

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
    
    hframe=ROOT.TH2D("hframe","",100,0.60,0.67,len(channels)*2+1,0,len(channels)*2+1);
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

    #DrawAccuracy(["mn201[678][ab]?","en201[678][ab]?","me201[678][ab]?","mm201[678][ab]?","ee201[678][ab]?"])
    #DrawAccuracy(["ee201[678][ab]?","mm201[678][ab]?","me201[678][ab]?"])
    #DrawAccuracy([channel+era for channel in ["mn","en","me","mm","ee"] for era in ["2016a","2016b","2017","2018","201[678][ab]?"]])
    #DrawAccuracy([channel+era for channel in ["[me]n","[me][me]"] for era in ["2016a","2016b","2017","2018","201[678][ab]?"]])
    
    #DrawAccuracyByType(["mn201[678][ab]?","en201[678][ab]?","me201[678][ab]?","mm201[678][ab]?","ee201[678][ab]?"],[0,1,2])
    #DrawAccuracyByType([channel+era for channel in ["mn","en"] for era in ["2016a","2016b","2017","2018","201[678][ab]?"]],[0,1,2])

    leps = ["[em]", "e", "m"]
    eras = ["2016a", "2016b", "2017", "2018", "201[678][ab]?"]
    chargeBins = ["", "_4", "_5", "_[0-3]", "_0", "_1", "_2", "_3"]
    ptBins = ["", "_pt0", "_pt1", "_pt2", "_pt3", "_pt4"]
    option = "syst" # "syst" or ""

    #DrawAccuracy_withLR([leps[0]+eras[0]+chargeBin for chargeBin in chargeBins], (0.51, 0.84), eras[0]+"s", option)
    #DrawAccuracy_withLR([leps[0]+eras[1]+chargeBin for chargeBin in chargeBins], (0.51, 0.84), eras[1]+"s", option)
    #DrawAccuracy_withLR([leps[0]+eras[2]+chargeBin for chargeBin in chargeBins], (0.51, 0.84), eras[2]+"s", option)
    #DrawAccuracy_withLR([leps[0]+eras[3]+chargeBin for chargeBin in chargeBins], (0.51, 0.84), eras[3]+"s", option)
    DrawAccuracy_withLR([leps[0]+eras[4]+chargeBin for chargeBin in chargeBins], (0.51, 0.84), "Run2ss", option)

    #DrawAccuracy_withLR([leps[0]+eras[4]+chargeBins[0]+ptBin for ptBin in ptBins], (0.51, 0.84), "pt_chargeBin", option)
    #DrawAccuracy_withLR([leps[0]+eras[4]+chargeBins[1]+ptBin for ptBin in ptBins], (0.51, 0.84), "pt_chargeBin4", option)
    #DrawAccuracy_withLR([leps[0]+eras[4]+chargeBins[2]+ptBin for ptBin in ptBins], (0.51, 0.84), "pt_chargeBin5", option)
    #DrawAccuracy_withLR([leps[0]+eras[4]+chargeBins[3]+ptBin for ptBin in ptBins], (0.51, 0.84), "pt_chargeBin0123", option)
    #DrawAccuracy_withLR([leps[0]+eras[4]+chargeBins[4]+ptBin for ptBin in ptBins], (0.51, 0.84), "pt_chargeBin0", option)
    #DrawAccuracy_withLR([leps[0]+eras[4]+chargeBins[5]+ptBin for ptBin in ptBins], (0.51, 0.84), "pt_chargeBin1", option)
    #DrawAccuracy_withLR([leps[0]+eras[4]+chargeBins[6]+ptBin for ptBin in ptBins], (0.51, 0.84), "pt_chargeBin2", option)
    #DrawAccuracy_withLR([leps[0]+eras[4]+chargeBins[7]+ptBin for ptBin in ptBins], (0.51, 0.84), "pt_chargeBin3", option)

    #DrawAccuracy_withLR([leps[0]+eras[4]+"F"+chargeBin for chargeBin in chargeBins], (0.505, 0.84), "Run2F", option)
    #DrawAccuracy_withLR([leps[0]+eras[4]+"S"+chargeBin for chargeBin in chargeBins], (0.505, 0.84), "Run2S", option)
    #DrawAccuracy_withLR([leps[0]+eras[4]+"L"+chargeBin for chargeBin in chargeBins], (0.505, 0.84), "Run2L", option)
    #DrawAccuracy_withLR([leps[0]+eras[4]+"M"+chargeBin for chargeBin in chargeBins], (0.505, 0.84), "Run2M", option)
    #DrawAccuracy_withLR([leps[0]+eras[4]+"H"+chargeBin for chargeBin in chargeBins], (0.505, 0.84), "Run2H", option)
    #DrawAccuracy_withLR([leps[0]+eras[4]+"V"+chargeBin for chargeBin in chargeBins], (0.505, 0.84), "Run2V", option)

    #DrawAccuracy_withLR([leps[0]+eras[4]+nPV+chargeBins[0] for nPV in ["", "F", "S", "L", "M", "H", "V"]], (0.59, 0.66), "nPV_chargeBin", option)
    #DrawAccuracy_withLR([leps[0]+eras[4]+nPV+chargeBins[1] for nPV in ["", "F", "S", "L", "M", "H", "V"]], (0.70, 0.79), "nPV_chargeBin4", option)
    #DrawAccuracy_withLR([leps[0]+eras[4]+nPV+chargeBins[2] for nPV in ["", "F", "S", "L", "M", "H", "V"]], (0.63, 0.79), "nPV_chargeBin5", option)
    #DrawAccuracy_withLR([leps[0]+eras[4]+nPV+chargeBins[3] for nPV in ["", "F", "S", "L", "M", "H", "V"]], (0.57, 0.645), "nPV_chargeBin0123", option)
    #DrawAccuracy_withLR([leps[0]+eras[4]+nPV+chargeBins[4] for nPV in ["", "F", "S", "L", "M", "H", "V"]], (0.505, 0.55), "nPV_chargeBin0", option)
    #DrawAccuracy_withLR([leps[0]+eras[4]+nPV+chargeBins[5] for nPV in ["", "F", "S", "L", "M", "H", "V"]], (0.54, 0.62), "nPV_chargeBin1", option)
    #DrawAccuracy_withLR([leps[0]+eras[4]+nPV+chargeBins[6] for nPV in ["", "F", "S", "L", "M", "H", "V"]], (0.61, 0.72), "nPV_chargeBin2", option)
    #DrawAccuracy_withLR([leps[0]+eras[4]+nPV+chargeBins[7] for nPV in ["", "F", "S", "L", "M", "H", "V"]], (0.62, 0.84), "nPV_chargeBin3", option)

    DrawAccuracy_withLR([leps[0]+era for era in eras], (0.605, 0.66), "Eras", option)
    #DrawAccuracy_withLR([lep+era for era in eras for lep in leps], (0.605, 0.66), "Eras_Leps", option)
    #chargeBins = ["", "_[0-3]", "_4", "_5"]
    #DrawAccuracy_withLR([leps[0]+era+chargeBin for era in eras for chargeBin in chargeBins], (0.59, 0.78), "Eras_Bins", option)
