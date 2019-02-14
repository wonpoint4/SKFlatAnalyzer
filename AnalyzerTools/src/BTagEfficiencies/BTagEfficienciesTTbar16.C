/* These Efficiencies have been derived for 2016 ttbar events and should 
be used only for the same MC samples or for events with similar topology */  


float BTagSFUtil::TagEfficiencyB_2016(float JetPt, float JetEta) {


if (TaggerOP=="DeepCSV_Loose") { 
   if (JetPt > 20.00000 && JetPt <= 40.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.83433 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.83204;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.80885;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.77735;  
  } 
   else if (JetPt > 40.00000 && JetPt <= 60.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.85844 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.85639;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.83447;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.80207;  
  } 
   else if (JetPt > 60.00000 && JetPt <= 80.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.87174 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.86972;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.84844;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.80978;  
  } 
   else if (JetPt > 80.00000 && JetPt <= 100.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.87931 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.87682;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.85382;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.81212;  
  } 
   else if (JetPt > 100.00000 && JetPt <= 120.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.88245 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.88024;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.85618;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.81299;  
  } 
   else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.88486 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.88051;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.85705;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.81617;  
  } 
 } 
if (TaggerOP=="DeepCSV_Medium") { 
   if (JetPt > 20.00000 && JetPt <= 40.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.68042 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.67330;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.62360;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.55299;  
  } 
   else if (JetPt > 40.00000 && JetPt <= 60.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.71796 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.71099;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.66505;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.59497;  
  } 
   else if (JetPt > 60.00000 && JetPt <= 80.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.73582 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.72891;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.68198;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.60412;  
  } 
   else if (JetPt > 80.00000 && JetPt <= 100.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.74356 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.73623;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.68466;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.59570;  
  } 
   else if (JetPt > 100.00000 && JetPt <= 120.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.74563 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.73732;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.68166;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.58905;  
  } 
   else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.73766 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.72384;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.66359;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.56082;  
  } 
 } 
if (TaggerOP=="DeepCSV_Tight") { 
   if (JetPt > 20.00000 && JetPt <= 40.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.50359 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.49275;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.42951;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.35782;  
  } 
   else if (JetPt > 40.00000 && JetPt <= 60.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.55217 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.54136;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.47978;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.40704;  
  } 
   else if (JetPt > 60.00000 && JetPt <= 80.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.57344 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.56172;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.50024;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.42079;  
  } 
   else if (JetPt > 80.00000 && JetPt <= 100.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.57973 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.56775;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.50292;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.41202;  
  } 
   else if (JetPt > 100.00000 && JetPt <= 120.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.57784 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.56519;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.49637;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.40274;  
  } 
   else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.55623 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.53897;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.46423;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.36610;  
  } 
 } 
 return 1.;
 
 } 
float BTagSFUtil::TagEfficiencyC_2016(float JetPt, float JetEta) {


if (TaggerOP=="DeepCSV_Loose") { 
   if (JetPt > 20.00000 && JetPt <= 40.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.39653 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.40356;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.40823;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.41842;  
  } 
   else if (JetPt > 40.00000 && JetPt <= 60.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.38426 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.39607;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.40289;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.40673;  
  } 
   else if (JetPt > 60.00000 && JetPt <= 80.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.37659 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.38728;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.39743;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.39841;  
  } 
   else if (JetPt > 80.00000 && JetPt <= 100.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.37808 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.39002;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.40032;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.39434;  
  } 
   else if (JetPt > 100.00000 && JetPt <= 120.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.38748 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.39377;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.40124;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.39782;  
  } 
   else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.40511 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.40876;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.42282;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.41951;  
  } 
 } 
if (TaggerOP=="DeepCSV_Medium") { 
   if (JetPt > 20.00000 && JetPt <= 40.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.12826 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.13366;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.13181;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.11913;  
  } 
   else if (JetPt > 40.00000 && JetPt <= 60.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.11673 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.12393;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.12633;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.11763;  
  } 
   else if (JetPt > 60.00000 && JetPt <= 80.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.11010 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.11675;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.12018;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.10958;  
  } 
   else if (JetPt > 80.00000 && JetPt <= 100.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.11052 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.11754;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.11826;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.10721;  
  } 
   else if (JetPt > 100.00000 && JetPt <= 120.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.11433 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.11834;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.12059;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.10769;  
  } 
   else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.12654 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.12752;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.12661;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.10845;  
  } 
 } 
if (TaggerOP=="DeepCSV_Tight") { 
   if (JetPt > 20.00000 && JetPt <= 40.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.02395 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.02583;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.02604;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.02269;  
  } 
   else if (JetPt > 40.00000 && JetPt <= 60.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.02176 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.02406;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.02472;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.02308;  
  } 
   else if (JetPt > 60.00000 && JetPt <= 80.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.02083 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.02314;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.02394;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.02148;  
  } 
   else if (JetPt > 80.00000 && JetPt <= 100.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.02160 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.02363;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.02371;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.02020;  
  } 
   else if (JetPt > 100.00000 && JetPt <= 120.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.02334 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.02409;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.02530;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.02056;  
  } 
   else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.02652 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.02677;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.02605;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.02135;  
  } 
 } 
 return 1.;
 
 } 
float BTagSFUtil::TagEfficiencyLight_2016(float JetPt, float JetEta) {


if (TaggerOP=="DeepCSV_Loose") { 
   if (JetPt > 20.00000 && JetPt <= 40.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.10141 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.12004;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.15840;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.19580;  
  } 
   else if (JetPt > 40.00000 && JetPt <= 60.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.06982 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.08729;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.12566;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.15805;  
  } 
   else if (JetPt > 60.00000 && JetPt <= 80.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.05297 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.06685;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.09833;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.12789;  
  } 
   else if (JetPt > 80.00000 && JetPt <= 100.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.04927 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.06347;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.09164;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.12302;  
  } 
   else if (JetPt > 100.00000 && JetPt <= 120.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.04985 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.06329;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.09141;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.12327;  
  } 
   else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.06612 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.08012;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.11768;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.15626;  
  } 
 } 
if (TaggerOP=="DeepCSV_Medium") { 
   if (JetPt > 20.00000 && JetPt <= 40.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.01167 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.01462;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.01962;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.02147;  
  } 
   else if (JetPt > 40.00000 && JetPt <= 60.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00758 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00975;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.01439;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01733;  
  } 
   else if (JetPt > 60.00000 && JetPt <= 80.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00568 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00719;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.01013;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01286;  
  } 
   else if (JetPt > 80.00000 && JetPt <= 100.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00542 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00648;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00896;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01211;  
  } 
   else if (JetPt > 100.00000 && JetPt <= 120.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00548 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00670;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00909;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01175;  
  } 
   else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00726 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00861;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.01311;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01509;  
  } 
 } 
if (TaggerOP=="DeepCSV_Tight") { 
   if (JetPt > 20.00000 && JetPt <= 40.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00100 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00122;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00162;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.00181;  
  } 
   else if (JetPt > 40.00000 && JetPt <= 60.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00069 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00090;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00130;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.00153;  
  } 
   else if (JetPt > 60.00000 && JetPt <= 80.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00060 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00075;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00093;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.00121;  
  } 
   else if (JetPt > 80.00000 && JetPt <= 100.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00061 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00075;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00085;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.00103;  
  } 
   else if (JetPt > 100.00000 && JetPt <= 120.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00060 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00079;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00087;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.00110;  
  } 
   else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00089 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00105;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00142;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.00153;  
  } 
 } 
 return 1.;
 
 } 
