import math
import numpy as np
import ROOT
ROOT.gROOT.ProcessLine('#include"BBPlotter.cc"')
ROOT.gROOT.ProcessLine('#include"ttljPlotter.cc"')
ROOT.Plotter.SetupStyle()

p=ROOT.BBPlotter("eff")
ttlj = ROOT.ttljPlotter("correct_ttlj wrong_ttlj unmatched_ttlj")
ttlj_gen = ROOT.ttljPlotter("ttlj_kin")
forNorm = ROOT.ttljPlotter("data mc")

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
    correct = ttlj.GetHist(0, channel+"/lepbjetCharge"+chargeBin+"Easy_Lm"+syst).Integral()
    wrong = ttlj.GetHist(1, channel+"/lepbjetCharge"+chargeBin+"Easy_Lm"+syst).Integral()
    fc = correct / (correct + wrong)
    fc_e = fc * (1 - fc) / (correct + wrong)

    data = forNorm.GetHist(0, channel+"/lepbjetCharge"+chargeBin+"Easy_Lm"+("" if "jet_scale" not in syst else syst)).Integral()
    mc = forNorm.GetHist(1, channel+"/lepbjetCharge"+chargeBin+"Easy_Lm"+syst).Integral()

    return fc, fc_e, data / mc

def calc_withLR(fc, flep, fhad):
    ap = (fc * (flep + fc -1) + (1 - fc) * (fhad + fc -1)) / (2 * fc - 1)
    am = ((1 - fc) * (flep + fc -1) + fc * (fhad + fc -1)) / (2 * fc - 1)

    return ap, am

def calcWithCov_withLR(fcs, fleps, fhads, fc_stat, flep_stat, fhad_stat):
    cov_stat = np.zeros((2, 2))
    cov = np.zeros((2, 2))
    print("fc", fcs)
    fc_nominal = fcs[-1][0]
    flep_nominal = fleps[-1][0]
    fhad_nominal = fhads[-1][0]
    nominal = np.array(calc_withLR(fc_nominal, flep_nominal, fhad_nominal))
    stat0 = np.array(calc_withLR(fc_nominal + fc_stat, flep_nominal, fhad_nominal))
    stat1 = np.array(calc_withLR(fc_nominal, flep_nominal + flep_stat, fhad_nominal))
    stat2 = np.array(calc_withLR(fc_nominal, flep_nominal, fhad_nominal + fhad_stat))
    cov_stat += np.outer(stat0 - nominal, stat0 - nominal)
    cov_stat += np.outer(stat1 - nominal, stat1 - nominal)
    cov_stat += np.outer(stat2 - nominal, stat2 - nominal)
    print("nomianl", nominal, "stat unc", [math.sqrt(cov_stat[0][0]), math.sqrt(cov_stat[1][1])])
    print("cov", cov_stat)

    values, vectors= np.linalg.eig(cov_stat)
    vector0 = np.array([vectors[0][0], vectors[1][0]]) * values[0] ** 0.5
    vector1 = np.array([vectors[0][1], vectors[1][1]]) * values[1] ** 0.5
    print("stat vect0", vector0)
    print("stat vect1", vector1)
    print(np.outer(vector0, vector0) + np.outer(vector1, vector1))

    cov += cov_stat
    if len(fcs) == 1: return nominal, cov_stat, cov

    print("\nSystematics")
    for systs in range(len(fcs) - 1):
        if systs == 0: print(" -JES up, down")
        elif systs == 1: print(" -JER up, down")
        elif systs == 2: print(" -Prefiring up, down")
        elif systs == 3: print(" -PU reweight up, down")
        elif systs == 4: print(" -PUjetID SF up, down")
        elif systs == 5: print(" -btagSF h, l up, down")
        elif systs == 7: print(" -btagSF h, l corr, uncorr")
        cov_syst_bigger = np.zeros((2,2))
        for syst in range(len(fcs[systs])):
            dsyst = np.array(calc_withLR(fcs[systs][syst], fleps[systs][syst], fhads[systs][syst])) - nominal
            cov_syst = np.outer(dsyst, dsyst)
            print("cov syst", systs, syst)
            print(cov_syst)
            if np.trace(cov_syst) > np.trace(cov_syst_bigger): cov_syst_bigger = cov_syst ## Choose the cov_systematic with the larger trace
        cov += cov_syst_bigger

    print("nomianl", nominal, "stat+syst unc", [math.sqrt(cov[0][0]), math.sqrt(cov[1][1])])
    print("cov", cov)
    values, vectors= np.linalg.eig(cov)
    vector0 = np.array([vectors[0][0], vectors[1][0]]) * values[0] ** 0.5
    vector1 = np.array([vectors[0][1], vectors[1][1]]) * values[1] ** 0.5
    print("stat+syst vect0", vector0)
    print("stat+syst vect1", vector1)
    print(np.outer(vector0, vector0) + np.outer(vector1, vector1))

    return nominal, cov_stat, cov

def GetAccuracy_withLR(ientry, channel, chargeBin="", option=""):
    allsysts = [[""]]
    if "syst" in option:
        allsysts = [
            ["_jet_scale_up", "_jet_scale_down"], ["_jet_smear_up", "_jet_smear_down"],
            ["_prefireweight_up", "_prefireweight_down"], ["_PUweight_up", "_PUweight_down"], ["_PUjetSF_up", "_PUjetSF_down"],
            ["_btagSF_hup", "_btagSF_hdown"], ["_btagSF_lup", "_btagSF_ldown"],
            ["_btagSF_hcorr"], ["_btagSF_huncorr2016a"], ["_btagSF_huncorr2016b"], ["_btagSF_huncorr2017"], ["_btagSF_huncorr2018"],
            ["_btagSF_lcorr"], ["_btagSF_luncorr2016a"], ["_btagSF_luncorr2016b"], ["_btagSF_luncorr2017"], ["_btagSF_luncorr2018"],
            [""]
        ]

    print("chargeBin : "+chargeBin+", option = "+option)
    fc, flep, fhad = 0, 0, 0
    fc_e, flep_e, fhad_e = -1, -1, -1
    fcs, fleps, fhads = [[]], [[]], [[]]

    for systs in range(len(allsysts)):
        for syst in allsysts[systs]:
            fc, fc_stat, norm = getfc(channel, chargeBin, syst)
            a = ROOT.ttljPlotter("data_sub ttlj", norm)
            hlep = a.GetHist(ientry, channel+"/lepbjetCharge"+chargeBin+"Easy_Lm"+("" if "jet_scale" not in syst and ientry == 0 else syst))
            hhad = a.GetHist(ientry, channel+"/hadbjetCharge"+chargeBin+"Easy_Lm"+("" if "jet_scale" not in syst and ientry == 0 else syst))
            hlep.Scale(1. / hlep.Integral())
            hhad.Scale(1. / hhad.Integral())
            fcs[systs].append(fc)
            fleps[systs].append(hlep.GetBinContent(2))
            flep_stat = hlep.GetBinError(2)
            fhads[systs].append(hhad.GetBinContent(1))
            fhad_stat = hhad.GetBinError(1)
        if systs < len(allsysts) - 1:
            fcs.append([])
            fleps.append([])
            fhads.append([])

    nominal, cov_stat, cov = calcWithCov_withLR(fcs, fleps, fhads, fc_stat, flep_stat, fhad_stat)
    return nominal, cov_stat, cov

def GetTrueAccuracy_withLR(channel, chargeBin="", option=""):
    allsysts = [[""]]
    if "syst" in option:
        allsysts.extend([
            ["_jet_scale_up", "_jet_scale_down"], ["_jet_smear_up", "_jet_smear_down"],
            ["_prefireweight_up", "_prefireweight_down"], ["_PUweight_up", "_PUweight_down"], ["_PUjetSF_up", "_PUjetSF_down"],
            ["_btagSF_hup", "_btagSF_hdown"], ["_btagSF_lup", "_btagSF_ldown"],
            ["_btagSF_hcorr"], ["_btagSF_huncorr2016a"], ["_btagSF_huncorr2016b"], ["_btagSF_huncorr2017"], ["_btagSF_huncorr2018"],
            ["_btagSF_lcorr"], ["_btagSF_luncorr2016a"], ["_btagSF_luncorr2016b"], ["_btagSF_luncorr2017"], ["_btagSF_luncorr2018"],
        ])

    ap, am, stat_ep, stat_em, ep, em = 0, 0, -1, -1, 0, 0
    for systs in range(len(allsysts)):
        syst_ep_bigger, syst_em_bigger = 0, 0
        for syst in allsysts[systs]:
            hp = ttlj_gen.GetHist(0, channel+"/genbbarjetCharge"+chargeBin+"Easy_Lm"+syst, option)
            hm = ttlj_gen.GetHist(0, channel+"/genbjetCharge"+chargeBin+"Easy_Lm"+syst, option)
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

    Bins = ["_[0-3]", "_0", "_1", "_2", "_3", "_4", "_5"]
    for i in range(len(channels)):
        channel = channels[i]
        print("\n\n@@@ Channel : "+channel+" started @@@\n")

        chargeBin = ""
        for Bin in Bins:
            if Bin in channel:
                channel = channel.replace(Bin, "")
                chargeBin = Bin.replace("_", "")
                break

        value, cov_stat, cov = GetAccuracy_withLR(0, channel, chargeBin, option)
        for j in range(2):
            gdata[j].SetPoint(i, value[j], 2 * i + 1 - j + 0.5)
            gdata[j].SetPointError(i, cov_stat[j][j]**0.5, 0)
            gdata_tot_unc[j].SetPoint(i, value[j], 2 * i + 1 - j + 0.5)
            gdata_tot_unc[j].SetPointError(i, cov[j][j]**0.5, 0)

        value, cov_stat, cov = GetAccuracy_withLR(1, channel, chargeBin, option)
        for j in range(2):
            gsim[j].SetPoint(i, value[j], 2 * i + 1 - j + 0.4)
            gsim[j].SetPointError(i, cov_stat[j][j]**0.5, 0)
            gsim_tot_unc[j].SetPoint(i, value[j], 2 * i + 1 - j + 0.4)
            gsim_tot_unc[j].SetPointError(i, cov[j][j]**0.5, 0)

        value, stat_e, e = GetTrueAccuracy_withLR(channel, chargeBin, option)
        for j in range(2):
            gtrue[j].SetPoint(i, value[j], 2 * i + 1 - j + 0.3)
            gtrue[j].SetPointError(i, stat_e[j], 0)
            gtrue_tot_unc[j].SetPoint(i, value[j], 2 * i + 1 - j + 0.3)
            gtrue_tot_unc[j].SetPointError(i, e[j], 0)

    c = ROOT.gROOT.MakeDefCanvas()
    c.SetLeftMargin(0.2)

    hframe = ROOT.TH2D("hframe", "", 100, Xrange[0], Xrange[1], len(channels) * 2 + 2, 0, len(channels) * 2 + 2);
    for i in range(len(channels)):
        title = channels[i]
        title = title.replace("201[678][ab]?"," Run2")
        title = title.replace("E","e")
        title = title.replace("m","#mu")
        title = title.replace("[e#mu] Run2F", "nPV <= 10")
        title = title.replace("[e#mu] Run2S", "nPV(10,20]")
        title = title.replace("[e#mu] Run2L", "nPV(20,30]")
        title = title.replace("[e#mu] Run2M", "nPV(30,40]")
        title = title.replace("[e#mu] Run2H", "nPV(40,50]")
        title = title.replace("[e#mu] Run2V", "nPV > 50")
        chargeBin = ""
        for Bin in Bins:
            if Bin in title:
                title = title.replace(Bin, "")
                chargeBin = Bin.replace("_", "")
                break
        alpha = "#alpha"
        if chargeBin != "":
            if chargeBin == "[0-3]": alpha = "#alpha_{j}"
            elif chargeBin == "4": alpha = "#alpha_{#mu}"
            elif chargeBin == "5": alpha = "#alpha_{e}"
            else: alpha = "#alpha_{j,"+chargeBin+"}"
        hframe.GetYaxis().SetBinLabel(i * 2 + 1, alpha+"^{#minus} ("+title+")")
        hframe.GetYaxis().SetBinLabel(i * 2 + 2, alpha+"^{#plus} ("+title+")")

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

    #raw_input()
    c.hists = [gdata, gsim, gtrue, gdata_tot_unc, gsim_tot_unc, gtrue_tot_unc, hframe, leg]
    nametag = ""
    if tag != "": nametag = nametag+"_"+tag
    if option != "": nametag = nametag+"_"+option
    c.SaveAs("Accuracies"+nametag+".png")

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

    leps = ["[Em]", "E", "m"]
    eras = ["2016a", "2016b", "2017", "2018", "201[678][ab]?"]
    chargeBins = ["", "_4", "_5", "_[0-3]", "_0", "_1", "_2", "_3"]
    option = "syst" # "syst" or ""

    DrawAccuracy_withLR([leps[0]+eras[0]+chargeBin for chargeBin in chargeBins], (0.51, 0.84), eras[0]+"s", option)
    DrawAccuracy_withLR([leps[0]+eras[1]+chargeBin for chargeBin in chargeBins], (0.51, 0.84), eras[1]+"s", option)
    DrawAccuracy_withLR([leps[0]+eras[2]+chargeBin for chargeBin in chargeBins], (0.51, 0.84), eras[2]+"s", option)
    DrawAccuracy_withLR([leps[0]+eras[3]+chargeBin for chargeBin in chargeBins], (0.51, 0.84), eras[3]+"s", option)
    DrawAccuracy_withLR([leps[0]+eras[4]+chargeBin for chargeBin in chargeBins], (0.51, 0.84), "Run2s", option)

    DrawAccuracy_withLR([leps[0]+eras[4]+"F"+chargeBin for chargeBin in chargeBins], (0.48, 0.86), "Run2F", option)
    DrawAccuracy_withLR([leps[0]+eras[4]+"S"+chargeBin for chargeBin in chargeBins], (0.48, 0.86), "Run2S", option)
    DrawAccuracy_withLR([leps[0]+eras[4]+"L"+chargeBin for chargeBin in chargeBins], (0.48, 0.86), "Run2L", option)
    DrawAccuracy_withLR([leps[0]+eras[4]+"M"+chargeBin for chargeBin in chargeBins], (0.48, 0.86), "Run2M", option)
    DrawAccuracy_withLR([leps[0]+eras[4]+"H"+chargeBin for chargeBin in chargeBins], (0.48, 0.86), "Run2H", option)
    DrawAccuracy_withLR([leps[0]+eras[4]+"V"+chargeBin for chargeBin in chargeBins], (0.48, 0.86), "Run2V", option)

    DrawAccuracy_withLR([leps[0]+eras[4]+nPV+chargeBins[0] for nPV in ["", "F", "S", "L", "M", "H", "V"]], (0.58, 0.66), "chargeBin", option)
    DrawAccuracy_withLR([leps[0]+eras[4]+nPV+chargeBins[1] for nPV in ["", "F", "S", "L", "M", "H", "V"]], (0.66, 0.83), "chargeBin4", option)
    DrawAccuracy_withLR([leps[0]+eras[4]+nPV+chargeBins[2] for nPV in ["", "F", "S", "L", "M", "H", "V"]], (0.65, 0.81), "chargeBin5", option)
    DrawAccuracy_withLR([leps[0]+eras[4]+nPV+chargeBins[3] for nPV in ["", "F", "S", "L", "M", "H", "V"]], (0.56, 0.65), "chargeBin0123", option)
    DrawAccuracy_withLR([leps[0]+eras[4]+nPV+chargeBins[4] for nPV in ["", "F", "S", "L", "M", "H", "V"]], (0.48, 0.56), "chargeBin0", option)
    DrawAccuracy_withLR([leps[0]+eras[4]+nPV+chargeBins[5] for nPV in ["", "F", "S", "L", "M", "H", "V"]], (0.52, 0.63), "chargeBin1", option)
    DrawAccuracy_withLR([leps[0]+eras[4]+nPV+chargeBins[6] for nPV in ["", "F", "S", "L", "M", "H", "V"]], (0.60, 0.73), "chargeBin2", option)
    DrawAccuracy_withLR([leps[0]+eras[4]+nPV+chargeBins[7] for nPV in ["", "F", "S", "L", "M", "H", "V"]], (0.55, 0.85), "chargeBin3", option)

    DrawAccuracy_withLR([leps[0]+era for era in eras], (0.6, 0.67), "Eras", option)
    DrawAccuracy_withLR([lep+era for era in eras for lep in leps], (0.6, 0.67), "Eras_Leps", option)
    chargeBins = ["", "_[0-3]", "_4", "_5"]
    DrawAccuracy_withLR([leps[0]+era+chargeBin for era in eras for chargeBin in chargeBins], (0.59, 0.79), "Eras_Bins", option)
