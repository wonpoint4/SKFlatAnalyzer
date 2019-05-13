import os

leps = [
"Electron",
"Muon",
]

for lep in leps:

  out = open('histmap_'+lep+'.txt','w')

  os.system('ls -1 '+lep+'*.root &> tmp.txt')
  lines = open('tmp.txt').readlines()
  os.system('rm tmp.txt')
  for line in lines:
    #Electron_HNPairVeryLoose_FakeRates.root
    line = line.strip('\n')
    thisfilename = line

    line = line.replace(lep+'_','')
    line = line.replace('_FakeRates.root','')
    ID = line
    out.write(ID+'\t'+thisfilename+'\n')
    #print line

  out.close()
