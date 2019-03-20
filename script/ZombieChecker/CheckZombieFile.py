import ROOT,os

def CheckZombieFile(filepath):
  f = ROOT.TFile(filepath)
  return f.IsZombie()
