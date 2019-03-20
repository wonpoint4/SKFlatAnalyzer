import os
from CheckZombieFile import *

SKFlat_WD = os.environ['SKFlat_WD']
SKFlatV = os.environ['SKFlatV']
DATA_DIR = SKFlat_WD+'/data/'+SKFlatV+'/'

USER = 'jhchoi'
Year = "2016"

PDs = open('samples_'+Year+'.txt').readlines()

for PD in PDs:

  if "#" in PD:
    continue

  PD = PD.strip('\n')
  PD = PD.replace(':','_')
  Files = open(DATA_DIR+'/'+Year+'/Sample/ForKISTI/'+PD+'.txt')

  print "###################"
  print PD
  print "###################"
  for File in Files:
    File = File.replace('root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/'+USER+'/SKFlat/','/xrootd_user/'+USER+'/xrootd/SKFlat/').strip('\n')
    if CheckZombieFile(File):
      print '@@ '+File
      #print '@@ --> removed'
      #os.system('rm '+File)
