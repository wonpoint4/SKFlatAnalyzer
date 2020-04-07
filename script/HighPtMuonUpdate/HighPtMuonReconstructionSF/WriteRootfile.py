import os,ROOT
from array import array

Years = [
"2016",
"2017",
"2018",
]

PBins = [50, 100, 150, 200, 300, 400, 600, 1500, 3500]
EtaBins = [0, 1.6, 2.4]

for Year in Years:

  print '@@@@ Year : '+Year
  print 'Eta\tEta\tP\tP\tSF\tErr'

  outroot = ROOT.TFile('HighPtMuonRecoSF_Run'+Year+'.root','RECREATE')

  histSF = ROOT.TH2D('HighPtMuonRecoSF', '', len(PBins)-1, array('d',PBins), len(EtaBins)-1, array('d',EtaBins))

  lines = open(Year+'.txt').readlines()
  for i in range(0,len(lines)):

    Eta_l = EtaBins[i]
    Eta_r = EtaBins[i+1]

    line = lines[i].strip('\n')
    words = line.split('\t')

    for j in range(0,len(words)):

      P_l = PBins[j]
      P_r = PBins[j+1]

      word = words[j]
      ValandErr = word.replace('(','').replace(')','').split()

      Value = ""
      Err = ""
      if '-' in word:
        Value = "1"
        Err = "0"
      else:
        Value = ValandErr[0]
        Err = ValandErr[1]

      histSF.SetBinContent(j+1, i+1, float(Value))
      histSF.SetBinError(j+1, i+1, float(Err))

      print str(Eta_l)+'\t'+str(Eta_r)+'\t'+str(P_l)+'\t'+str(P_r)+'\t'+Value+'\t'+Err

  histSF.Write()

  outroot.Close()
