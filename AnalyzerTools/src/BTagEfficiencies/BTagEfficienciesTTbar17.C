/* These Efficiencies have been derived for TTLJ_powheg events and should 
be used only for the same MC samples or for events with similar topology */  


float BTagSFUtil::TagEfficiencyB_2017(float JetPt, float JetEta) {


if (TaggerOP=="CSVv2_Loose") { 
   if (JetPt > 20.00000 && JetPt <= 40.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.86183 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.86738;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.83541;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.81806;  
  } 
   else if (JetPt > 40.00000 && JetPt <= 60.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.87383 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.87660;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.85353;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.84212;  
  } 
   else if (JetPt > 60.00000 && JetPt <= 80.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.88351 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.88668;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.86381;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.85072;  
  } 
   else if (JetPt > 80.00000 && JetPt <= 100.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.88677 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.88965;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.86689;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.85098;  
  } 
   else if (JetPt > 100.00000 && JetPt <= 120.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.88425 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.88463;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.85895;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.84030;  
  } 
   else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.87608 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.87440;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.84867;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.82376;  
  } 
 } 
if (TaggerOP=="CSVv2_Medium") { 
   if (JetPt > 20.00000 && JetPt <= 40.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.68326 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.66623;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.63343;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.59709;  
  } 
   else if (JetPt > 40.00000 && JetPt <= 60.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.71239 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.69737;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.66263;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.63348;  
  } 
   else if (JetPt > 60.00000 && JetPt <= 80.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.72528 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.71529;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.67098;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.64119;  
  } 
   else if (JetPt > 80.00000 && JetPt <= 100.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.72695 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.71828;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.67259;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.63525;  
  } 
   else if (JetPt > 100.00000 && JetPt <= 120.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.72129 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.71476;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.66100;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.62264;  
  } 
   else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.69700 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.68720;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.62061;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.57762;  
  } 
 } 
if (TaggerOP=="CSVv2_Tight") { 
   if (JetPt > 20.00000 && JetPt <= 40.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.45999 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.41640;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.40955;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.37676;  
  } 
   else if (JetPt > 40.00000 && JetPt <= 60.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.51009 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.47578;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.44107;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.41670;  
  } 
   else if (JetPt > 60.00000 && JetPt <= 80.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.50865 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.48503;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.43879;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.41821;  
  } 
   else if (JetPt > 80.00000 && JetPt <= 100.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.50871 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.48657;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.43412;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.40226;  
  } 
   else if (JetPt > 100.00000 && JetPt <= 120.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.51139 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.49915;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.42178;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.38847;  
  } 
   else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.47269 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.45939;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.36784;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.33376;  
  } 
 } 
if (TaggerOP=="DeepCSV_Loose") { 
   if (JetPt > 20.00000 && JetPt <= 40.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.86914 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.87686;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.85382;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.83420;  
  } 
   else if (JetPt > 40.00000 && JetPt <= 60.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.89085 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.89627;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.87535;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.85429;  
  } 
   else if (JetPt > 60.00000 && JetPt <= 80.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.90224 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.90740;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.88629;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.86678;  
  } 
   else if (JetPt > 80.00000 && JetPt <= 100.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.90727 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.91225;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.89327;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.87449;  
  } 
   else if (JetPt > 100.00000 && JetPt <= 120.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.90979 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.91401;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.89613;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.87858;  
  } 
   else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.90605 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.90812;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.89105;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.87057;  
  } 
 } 
if (TaggerOP=="DeepCSV_Medium") { 
   if (JetPt > 20.00000 && JetPt <= 40.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.70924 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.71813;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.67234;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.63319;  
  } 
   else if (JetPt > 40.00000 && JetPt <= 60.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.75066 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.75571;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.71197;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.67192;  
  } 
   else if (JetPt > 60.00000 && JetPt <= 80.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.77254 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.77679;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.73074;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.68698;  
  } 
   else if (JetPt > 80.00000 && JetPt <= 100.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.78135 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.78415;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.73533;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.68346;  
  } 
   else if (JetPt > 100.00000 && JetPt <= 120.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.78415 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.78379;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.73099;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.67508;  
  } 
   else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.76593 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.75739;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.69782;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.63049;  
  } 
 } 
if (TaggerOP=="DeepCSV_Tight") { 
   if (JetPt > 20.00000 && JetPt <= 40.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.53696 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.53676;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.46895;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.42905;  
  } 
   else if (JetPt > 40.00000 && JetPt <= 60.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.59332 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.59049;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.52333;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.48077;  
  } 
   else if (JetPt > 60.00000 && JetPt <= 80.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.61888 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.61658;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.54717;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.49991;  
  } 
   else if (JetPt > 80.00000 && JetPt <= 100.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.62471 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.62059;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.54901;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.49538;  
  } 
   else if (JetPt > 100.00000 && JetPt <= 120.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.62257 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.61596;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.54117;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.48271;  
  } 
   else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.58590 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.57047;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.49185;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.42590;  
  } 
 } 
 return 1.;
 
 } 
float BTagSFUtil::TagEfficiencyC_2017(float JetPt, float JetEta) {


if (TaggerOP=="CSVv2_Loose") { 
   if (JetPt > 20.00000 && JetPt <= 40.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.43268 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.45897;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.44058;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.43400;  
  } 
   else if (JetPt > 40.00000 && JetPt <= 60.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.41051 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.43372;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.43222;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.43938;  
  } 
   else if (JetPt > 60.00000 && JetPt <= 80.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.41277 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.43438;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.43097;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.43950;  
  } 
   else if (JetPt > 80.00000 && JetPt <= 100.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.41227 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.43598;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.43533;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.44202;  
  } 
   else if (JetPt > 100.00000 && JetPt <= 120.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.40727 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.42990;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.42441;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.43106;  
  } 
   else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.39958 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.42756;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.42579;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.43935;  
  } 
 } 
if (TaggerOP=="CSVv2_Medium") { 
   if (JetPt > 20.00000 && JetPt <= 40.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.13979 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.13061;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.14294;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.13069;  
  } 
   else if (JetPt > 40.00000 && JetPt <= 60.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.12666 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.12040;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.13687;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.13352;  
  } 
   else if (JetPt > 60.00000 && JetPt <= 80.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.12458 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.12186;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.13126;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.13044;  
  } 
   else if (JetPt > 80.00000 && JetPt <= 100.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.12158 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.12173;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.13223;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.13027;  
  } 
   else if (JetPt > 100.00000 && JetPt <= 120.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.12389 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.12774;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.13539;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.13290;  
  } 
   else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.11731 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.12460;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.12799;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.13311;  
  } 
 } 
if (TaggerOP=="CSVv2_Tight") { 
   if (JetPt > 20.00000 && JetPt <= 40.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.02310 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.01822;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.02735;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.02570;  
  } 
   else if (JetPt > 40.00000 && JetPt <= 60.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.02381 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.02049;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.02545;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.02649;  
  } 
   else if (JetPt > 60.00000 && JetPt <= 80.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.02050 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.01927;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.02293;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.02418;  
  } 
   else if (JetPt > 80.00000 && JetPt <= 100.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.02055 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.01996;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.02195;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.02137;  
  } 
   else if (JetPt > 100.00000 && JetPt <= 120.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.02392 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.02475;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.02336;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.02309;  
  } 
   else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.02321 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.02489;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.02052;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.02332;  
  } 
 } 
if (TaggerOP=="DeepCSV_Loose") { 
   if (JetPt > 20.00000 && JetPt <= 40.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.41463 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.44872;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.45259;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.43635;  
  } 
   else if (JetPt > 40.00000 && JetPt <= 60.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.40910 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.44513;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.45274;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.43506;  
  } 
   else if (JetPt > 60.00000 && JetPt <= 80.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.41237 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.44783;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.45775;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.44171;  
  } 
   else if (JetPt > 80.00000 && JetPt <= 100.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.41326 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.45127;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.46696;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.45368;  
  } 
   else if (JetPt > 100.00000 && JetPt <= 120.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.42313 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.45963;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.47635;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.46303;  
  } 
   else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.43083 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.46293;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.47988;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.46804;  
  } 
 } 
if (TaggerOP=="DeepCSV_Medium") { 
   if (JetPt > 20.00000 && JetPt <= 40.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.10297 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.12258;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.13184;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.12306;  
  } 
   else if (JetPt > 40.00000 && JetPt <= 60.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.09935 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.11917;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.13017;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.12504;  
  } 
   else if (JetPt > 60.00000 && JetPt <= 80.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.10361 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.12345;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.13365;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.12799;  
  } 
   else if (JetPt > 80.00000 && JetPt <= 100.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.10949 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.12861;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.13880;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.13239;  
  } 
   else if (JetPt > 100.00000 && JetPt <= 120.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.11702 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.13562;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.14252;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.13483;  
  } 
   else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.12402 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.13674;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.14120;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.12993;  
  } 
 } 
if (TaggerOP=="DeepCSV_Tight") { 
   if (JetPt > 20.00000 && JetPt <= 40.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.01783 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.02229;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.02352;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.02263;  
  } 
   else if (JetPt > 40.00000 && JetPt <= 60.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.01811 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.02257;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.02488;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.02419;  
  } 
   else if (JetPt > 60.00000 && JetPt <= 80.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.02030 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.02481;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.02674;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.02630;  
  } 
   else if (JetPt > 80.00000 && JetPt <= 100.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.02159 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.02629;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.02837;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.02658;  
  } 
   else if (JetPt > 100.00000 && JetPt <= 120.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.02382 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.02802;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.02944;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.02842;  
  } 
   else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.02530 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.02807;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.02810;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.02679;  
  } 
 } 
 return 1.;
 
 } 
float BTagSFUtil::TagEfficiencyLight_2017(float JetPt, float JetEta) {


  if (TaggerOP=="CSVv2_Loose") { 
    if (JetPt > 20.00000 && JetPt <= 40.00000){  
      if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.06246 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.08930;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.12301;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.15619;  
  } 
   else if (JetPt > 40.00000 && JetPt <= 60.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.03838 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.05578;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.08927;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.12962;  
  } 
   else if (JetPt > 60.00000 && JetPt <= 80.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.03229 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.04729;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.08047;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.12114;  
  } 
   else if (JetPt > 80.00000 && JetPt <= 100.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.03177 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.04638;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.08232;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.12559;  
  } 
   else if (JetPt > 100.00000 && JetPt <= 120.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.02989 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.04364;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.07680;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.12485;  
  } 
   else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.03822 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.05551;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.09945;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.15429;  
  } 
 } 
if (TaggerOP=="CSVv2_Medium") { 
   if (JetPt > 20.00000 && JetPt <= 40.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00624 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00760;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.01342;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01798;  
  } 
   else if (JetPt > 40.00000 && JetPt <= 60.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00416 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00525;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00875;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01454;  
  } 
   else if (JetPt > 60.00000 && JetPt <= 80.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00348 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00437;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00707;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01222;  
  } 
   else if (JetPt > 80.00000 && JetPt <= 100.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00345 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00437;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00706;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01223;  
  } 
   else if (JetPt > 100.00000 && JetPt <= 120.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00319 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00419;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00689;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01251;  
  } 
   else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00433 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00595;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00906;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01730;  
  } 
 } 
if (TaggerOP=="CSVv2_Tight") { 
   if (JetPt > 20.00000 && JetPt <= 40.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00075 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00075;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00122;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.00206;  
  } 
   else if (JetPt > 40.00000 && JetPt <= 60.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00044 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00050;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00071;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.00160;  
  } 
   else if (JetPt > 60.00000 && JetPt <= 80.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00044 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00047;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00058;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.00108;  
  } 
   else if (JetPt > 80.00000 && JetPt <= 100.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00050 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00057;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00062;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.00111;  
  } 
   else if (JetPt > 100.00000 && JetPt <= 120.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00059 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00073;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00062;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.00116;  
  } 
   else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00085 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00107;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00094;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.00183;  
  } 
 } 
if (TaggerOP=="DeepCSV_Loose") { 
   if (JetPt > 20.00000 && JetPt <= 40.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.05482 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.08186;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.12651;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.15130;  
  } 
   else if (JetPt > 40.00000 && JetPt <= 60.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.04590 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.07012;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.10830;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.12851;  
  } 
   else if (JetPt > 60.00000 && JetPt <= 80.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.04018 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.06237;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.09850;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.11956;  
  } 
   else if (JetPt > 80.00000 && JetPt <= 100.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.04018 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.06146;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.09851;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.12229;  
  } 
   else if (JetPt > 100.00000 && JetPt <= 120.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.04275 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.06500;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.10209;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.12742;  
  } 
   else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.05520 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.07708;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.12036;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.15160;  
  } 
 } 
if (TaggerOP=="DeepCSV_Medium") { 
   if (JetPt > 20.00000 && JetPt <= 40.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00463 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00733;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.01220;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01438;  
  } 
   else if (JetPt > 40.00000 && JetPt <= 60.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00405 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00647;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.01076;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01318;  
  } 
   else if (JetPt > 60.00000 && JetPt <= 80.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00359 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00561;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00942;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01216;  
  } 
   else if (JetPt > 80.00000 && JetPt <= 100.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00378 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00568;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00926;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01187;  
  } 
   else if (JetPt > 100.00000 && JetPt <= 120.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00406 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00599;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00953;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01217;  
  } 
   else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00541 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00683;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.01111;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01425;  
  } 
 } 
if (TaggerOP=="DeepCSV_Tight") { 
   if (JetPt > 20.00000 && JetPt <= 40.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00065 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00090;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00136;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.00156;  
  } 
   else if (JetPt > 40.00000 && JetPt <= 60.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00044 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00065;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00106;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.00119;  
  } 
   else if (JetPt > 60.00000 && JetPt <= 80.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00045 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00061;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00093;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.00115;  
  } 
   else if (JetPt > 80.00000 && JetPt <= 100.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00050 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00068;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00097;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.00113;  
  } 
   else if (JetPt > 100.00000 && JetPt <= 120.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00056 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00073;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00098;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.00111;  
  } 
   else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
     if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00079 ;  
     else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00088;  
     else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00113;  
     else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.00127;  
  } 
 } 
 return 1.;
 
 } 
