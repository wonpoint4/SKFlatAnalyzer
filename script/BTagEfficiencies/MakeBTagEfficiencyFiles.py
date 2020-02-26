import os,ROOT

sample = 'TTLL_powheg'

SKFlat_WD = os.environ['SKFlat_WD']
SKFlatV = os.environ['SKFlatV']
DATA_DIR = os.environ['DATA_DIR']
SKFlatOutputDir = os.environ['SKFlatOutputDir']

Years = [
"2016",
"2017",
"2018",
]

Taggers = [
'DeepCSV',
'DeepJet',
'CSVv2',
]

WPs = [
'Loose',
'Medium',
'Tight',
]

Flavours = [
'B',
'C',
'Light',
]
cond_Flavours = [
'JetFlavor==5',
'JetFlavor==4',
'JetFlavor<=3',
]

out = open('JetTagEfficiencies.C','w')
print>>out,'''double MCCorrection::GetMCJetTagEff(JetTagging::Tagger tagger, JetTagging::WP wp, int JetFlavor, double JetPt, double JetEta){

  if(JetPt<20) JetPt = 20.;
  if(JetPt>=1000.) JetPt = 999.;

'''

for Year in Years:

  f = ROOT.TFile(SKFlatOutputDir+SKFlatV+'/MeasureJetTaggingEfficiency/'+Year+'/MeasureJetTaggingEfficiency_'+sample+'.root')

  out.write('  if(DataYear==%s){\n'%Year)

  for it_Flavour in range(0,len(Flavours)):

    Flavour = Flavours[it_Flavour]
    cond_Flavour = cond_Flavours[it_Flavour]

    h_Den = f.Get('Jet_'+Year+'_eff_'+Flavour+'_denom')

    out.write('    if('+cond_Flavour+'){\n')

    for Tagger in Taggers:

      print>>out,'''      if(tagger==JetTagging::{0}){{'''.format(Tagger)

      for WP in WPs:

        print>>out,'''        if(wp==JetTagging::{0}){{'''.format(WP)

        h_Eff = f.Get('Jet_'+Year+'_'+Tagger+'_'+WP+'_eff_'+Flavour+'_num')
        if not h_Eff:
          out.write('        }\n')
          continue
        h_Eff.Divide(h_Den)

        for iy in range(0,h_Eff.GetYaxis().GetNbins()):

          y_l = h_Eff.GetYaxis().GetBinLowEdge(iy+1)
          y_r = h_Eff.GetYaxis().GetBinUpEdge(iy+1)

          print>>out,'''          if(JetPt >= {0} && JetPt < {1}){{'''.format(y_l,y_r)

          for ix in range(0,h_Eff.GetXaxis().GetNbins()):

            x_l = h_Eff.GetXaxis().GetBinLowEdge(ix+1)
            x_r = h_Eff.GetXaxis().GetBinUpEdge(ix+1)

            this_eff_str = '%.5f'%(h_Eff.GetBinContent(ix+1,iy+1))

            print>>out,'''            if(fabs(JetEta) >= {0} && fabs(JetEta) < {1}) return {2};'''.format(x_l,x_r,this_eff_str)

            #print '%s\t%s\t%s\t\t[%f,%f]\t[%f,%f]'%(Tagger,Flavour,WP,x_l,x_r,y_l,y_r)
          out.write('          }\n') ## End of this pt range

        out.write('        }\n') ## end of this working point

      out.write('      }\n') ## end of this tagger

    out.write('    }\n') ## end of this flavour

  out.write('  }\n') ## end of this Year

print>>out,'''  cout << "[MCCorrection::GetMCJetTagEff] No eff found for tagger = " << tagger << ", wp = " << wp << endl;
  cout << "[MCCorrection::GetMCJetTagEff] Or, wrong pt and eta range : pt = " << JetPt << ", eta = " << JetEta << endl;'''.format(Year,Flavour)
out.write('  return 1.;\n')
out.write('}\n') ## end of this flavour

out.close()
