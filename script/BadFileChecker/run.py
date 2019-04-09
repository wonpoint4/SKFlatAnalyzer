import os
from CheckBadFile import *

SKFlat_WD = os.environ['SKFlat_WD']
SKFlatV = os.environ['SKFlatV']
DATA_DIR = SKFlat_WD+'/data/'+SKFlatV+'/'

USER = 'jalmond'
Year = "2018"

PDs = open('samples_'+Year+'.txt').readlines()

out = open('remove_command.sh','w')

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

    filepath_removable = File.replace('root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/'+USER+'/SKFlat/','/xrootd_user/'+USER+'/xrootd/SKFlat/').strip('\n')
    filepath_access = File.replace('root://cms-xrdr.sdfarm.kr:1094///xrd','/xrootd').strip('\n')

    if CheckBadFile(filepath_access):
      out.write('rm '+filepath_removable+'\n')
      #print '@@ '+filepath_removable
      #print '@@ --> removed'
      #os.system('rm '+filepath_removable)

out.close()
