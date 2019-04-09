import ROOT,os

def CheckBadFile(filepath):
  f = ROOT.TFile(filepath)

  IsZombie = f.IsZombie()
  IsRecovered = f.TestBit(1024)

  return (IsZombie or IsRecovered)
