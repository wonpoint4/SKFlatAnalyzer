import os

def GetXSECTable(InputSamples,XsecForEachSample):

  out = '------------------------------------------------\n'
  for i in range(0,len(InputSamples)):

    InputSample = InputSamples[i]
    XSEC = XsecForEachSample[i]
    out += InputSample+'\t'+str(XSEC)+' pb\n'

  out += '------------------------------------------------'
  return out
