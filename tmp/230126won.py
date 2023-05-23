import ROOT
ROOT.gROOT.ProcessLine("#include\"EfficiencyPlotter.cc\"")
ROOT.TH1.AddDirectory(0)
ROOT.gROOT.SetBatch(True)

if __name__=="__main__":
    p=ROOT.EfficiencyPlotter()
    p.plotdir="fig/230126won"
    for channel in ["ee","mm","el","mu"]:
        for era in ["2016a","2016b","2017","2018","201[678][ab]?"]:
            common="preliminary "
            lepton=channel[0].replace("m","#mu")
            channel_=lepton+lepton
            era_=era.replace("201[678][ab]?","Run2")
            p.SavePlot("{}{}/m52to150/dimass".format(channel,era),common+"xtitle:'m({}) [GeV]' save:{}{}_dimass.png,pdf".format(channel_,channel,era_))
            p.SavePlot("{}{}/m52to150/dimass".format(channel,era),common+"xtitle:'m({}) [GeV]' save:{}{}_dimass_notriggerSF.png,pdf suffix:_notriggerSF".format(channel_,channel,era_))

            p.SavePlot("{}{}/m80to100/l0pt".format(channel,era),common+"xtitle:'Leading {} p_{{T}} [GeV]' xmax:100 save:{}{}_l0pt.png,pdf".format(lepton,channel,era_))
            p.SavePlot("{}{}/m80to100/l0pt".format(channel,era),common+"xtitle:'Leading {} p_{{T}} [GeV]' xmax:100 save:{}{}_l0pt_notriggerSF.png,pdf suffix:_notriggerSF".format(lepton,channel,era_))
            
            p.SavePlot("{}{}/m80to100/l0eta".format(channel,era),common+"xtitle:'Leading {} #eta' xmin:-2.4 xmax:2.4 rebin:2 BMleg save:{}{}_l0eta.png,pdf".format(channel_[0],channel,era_))
            p.SavePlot("{}{}/m80to100/l0eta".format(channel,era),common+"xtitle:'Leading {} #eta' xmin:-2.4 xmax:2.4 rebin:2 BMleg save:{}{}_l0eta_notriggerSF.png,pdf suffix:_notriggerSF".format(channel_[0],channel,era_))
