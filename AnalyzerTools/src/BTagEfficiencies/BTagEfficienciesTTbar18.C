/* These Efficiencies have been derived for Summer15ttbar events and should 
   be used only for the same MC samples or for events with similar topology */  



//################################### FIXME THESE NUMBERS they are old 2016

float BTagSFUtil::TagEfficiencyB_2018(float JetPt, float JetEta) {
  
  
  if (TaggerOP=="CSVv2M") { 
    if (JetPt > 20.00000 && JetPt <= 40.00000){  
      if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.59234 ;  
      else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.57961;  
      else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.52772;  
      else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.43535;  
    } 
    else if (JetPt > 40.00000 && JetPt <= 60.00000){  
      if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.66022 ;  
      else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.64971;  
      else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.60206;  
      else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.52218;  
    } 
    else if (JetPt > 60.00000 && JetPt <= 80.00000){  
      if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.68184 ;  
      else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.67571;  
      else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.62020;  
      else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.53328;  
    } 
    else if (JetPt > 80.00000 && JetPt <= 100.00000){  
      if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.68998 ;  
      else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.68509;  
      else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.63066;  
      else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.53861;  
    } 
    else if (JetPt > 100.00000 && JetPt <= 120.00000){  
      if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.68730 ;  
      else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.68463;  
      else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.61622;  
      else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.52498;  
    } 
    else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
      if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.66520 ;  
      else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.66545;  
      else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.59499;  
      else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.51355;  
    } 
  } 
  if (TaggerOP=="CSVv2T") { 
    if (JetPt > 20.00000 && JetPt <= 40.00000){  
      if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.59234 ;  
      else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.57961;  
      else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.52772;  
      else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.43535;  
    } 
    else if (JetPt > 40.00000 && JetPt <= 60.00000){  
      if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.66022 ;  
      else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.64971;  
      else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.60206;  
      else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.52218;  
    } 
    else if (JetPt > 60.00000 && JetPt <= 80.00000){  
      if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.68184 ;  
      else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.67571;  
      else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.62020;  
      else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.53328;  
    } 
    else if (JetPt > 80.00000 && JetPt <= 100.00000){  
      if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.68998 ;  
      else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.68509;  
      else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.63066;  
      else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.53861;  
    } 
    else if (JetPt > 100.00000 && JetPt <= 120.00000){  
      if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.68730 ;  
      else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.68463;  
      else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.61622;  
      else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.52498;  
    } 
    else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
      if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.66520 ;  
      else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.66545;  
      else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.59499;  
      else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.51355;  
    } 
  } 
  
  return 1.;
}




float BTagSFUtil::TagEfficiencyC_2018(float JetPt, float JetEta) {
  
  if (JetPt > 0. && JetEta < 5.) return 1.;

  return 1.;
}
