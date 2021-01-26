import os,json

#f = open('RunBCDEF_SF_ID.json').read()
f = open('RunBCDEF_SF_ISO.json').read()

data = json.loads(f)

for ID in data:

  if "High" in ID:
    continue

  binnings = data[ID]

  print '  if(ID=="'+ID+'"){'

  for binning in binnings:

    etas = data[ID][binning]

    eta_ifst = "if"
    for eta in etas:
      ## eta = "abseta:[1.20,2.10]"
      temp_eta = eta.replace("abseta:[","")
      temp_eta = temp_eta.replace("]","")
      eta_lr = temp_eta.split(',')
      eta_l = float(eta_lr[0])
      eta_r = float(eta_lr[1])
      #print str(eta_l)+'~'+str(eta_r)

      print "    "+eta_ifst+'''(eta>={0} && eta<{1}){{'''.format(str(eta_l),str(eta_r))
      eta_ifst = "else if"

      pts = data[ID][binning][eta]

      pt_ifst = "if"
      for pt in pts:
        temp_pt = pt.replace("pt:[","")
        temp_pt = temp_pt.replace("]","")
        pt_lr = temp_pt.split(',')
        pt_l = float(pt_lr[0])
        pt_r = float(pt_lr[1])

        value = data[ID][binning][eta][pt]['value']
        error = data[ID][binning][eta][pt]['error']


        print "      "+pt_ifst+'''(pt>={0} && pt<{1}){{'''.format(str(pt_l),str(pt_r))
        print "        value = "+str(value)+";"
        print "        error = "+str(error)+";"
        print "      }"
        pt_ifst = "else if"

        #print eta+'\t'+pt+' : '+str(value)

      print '    }' ## end eta

    print '  }' ## END ID


    break ## hmm just print first binning
