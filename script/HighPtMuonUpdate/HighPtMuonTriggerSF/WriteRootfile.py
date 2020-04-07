import os,ROOT,math
from array import array

Years = [
"2016",
"2017",
"2018",
]

PBins = [52, 55, 60, 120, 200, 300, 500]
EtaBins = [0, 0.9, 1.2, 2.1, 2.4]

originalfile = ROOT.TFile('OutFile-v20190510-Combined-Run2016BtoH_Run2017BtoF_Run2018AtoD-M120to10000.root')

histnames = [
'Data',
'DY',
]

for Year in Years:

  print '@@@@ Year : '+Year
  print 'Eta\tEta\tP\tP\tSF\tErr'

  outroot = ROOT.TFile('HighPtMuonTriggerSF_Run'+Year+'.root','RECREATE')

  #### Data

  h_val_Data = originalfile.Get('Eff_'+Year+'_Data_var')
  h_err_Data_Lower = originalfile.Get('Eff_'+Year+'_Data_errorLower')
  h_err_Data_Upper = originalfile.Get('Eff_'+Year+'_Data_errorUpper')
  h_Data = h_val_Data.Clone()

  for ix in range(0, h_err_Data_Lower.GetXaxis().GetNbins()):
    for iy in range(0, h_err_Data_Lower.GetYaxis().GetNbins()):

      err_Lower = h_err_Data_Lower.GetBinContent(ix+1,iy+1)
      err_Upper = h_err_Data_Upper.GetBinContent(ix+1,iy+1)

      newerr = math.sqrt( (err_Lower*err_Lower+err_Upper*err_Upper)/2. )

      h_Data.SetBinError(ix+1, iy+1, newerr)

      #print '%1.5f\t%1.5f\t%1.5f' % ( err_Lower, err_Upper, newerr )

  h_Data.SetName('Eff_Data')
  h_Data.Write()

  #### MC

  h_val_MC = originalfile.Get('Eff_'+Year+'_DY_var')
  h_err_MC_Lower = originalfile.Get('Eff_'+Year+'_DY_errorLower')
  h_err_MC_Upper = originalfile.Get('Eff_'+Year+'_DY_errorUpper')
  h_MC = h_val_MC.Clone()

  for ix in range(0, h_err_MC_Lower.GetXaxis().GetNbins()):
    for iy in range(0, h_err_MC_Lower.GetYaxis().GetNbins()):

      err_Lower = h_err_MC_Lower.GetBinContent(ix+1,iy+1)
      err_Upper = h_err_MC_Upper.GetBinContent(ix+1,iy+1)

      newerr = math.sqrt( (err_Lower*err_Lower+err_Upper*err_Upper)/2. )

      h_MC.SetBinError(ix+1, iy+1, newerr)

      #print '%1.5f\t%1.5f\t%1.5f' % ( err_Lower, err_Upper, newerr )

  h_MC.SetName('Eff_MC')
  h_MC.Write()

  #### SF

  h_val_SF = originalfile.Get('SF_'+Year+'_var')
  h_err_SF_Lower = originalfile.Get('SF_'+Year+'_errorLower')
  h_err_SF_Upper = originalfile.Get('SF_'+Year+'_errorUpper')
  h_SF = h_val_SF.Clone()

  for ix in range(0, h_err_SF_Lower.GetXaxis().GetNbins()):
    for iy in range(0, h_err_SF_Lower.GetYaxis().GetNbins()):

      err_Lower = h_err_SF_Lower.GetBinContent(ix+1,iy+1)
      err_Upper = h_err_SF_Upper.GetBinContent(ix+1,iy+1)

      newerr = math.sqrt( (err_Lower*err_Lower+err_Upper*err_Upper)/2. )

      h_SF.SetBinError(ix+1, iy+1, newerr)

      #print '%1.5f\t%1.5f\t%1.5f' % ( err_Lower, err_Upper, newerr )

  h_SF.SetName('SF')
  h_SF.Write()



  outroot.Close()
