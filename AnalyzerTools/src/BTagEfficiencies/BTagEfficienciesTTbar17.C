/* These Efficiencies have been derived for 2017 JetTagging_TTLL_powheg.root events and should 
be used only for the same MC samples or for events with similar topology */  


float BTagSFUtil::TagEfficiencyB_2017(float JetPt, float JetEta) {


	 if (TaggerOP=="CSVv2_Loose") { 
		 if (JetPt > 20.00000 && JetPt <= 40.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.83171 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.83576;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.79668;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.77934;  
			 } 
		 else if (JetPt > 40.00000 && JetPt <= 60.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.87148 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.87422;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.84933;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.83858;  
			 } 
		 else if (JetPt > 60.00000 && JetPt <= 80.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.88388 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.88619;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.86424;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.85102;  
			 } 
		 else if (JetPt > 80.00000 && JetPt <= 100.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.88775 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.89052;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.86815;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.85160;  
			 } 
		 else if (JetPt > 100.00000 && JetPt <= 120.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.88612 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.88669;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.86097;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.84144;  
			 } 
		 else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.87895 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.87670;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.85112;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.82657;  
			 } 
		 } 
	 if (TaggerOP=="CSVv2_Medium") { 
		 if (JetPt > 20.00000 && JetPt <= 40.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.63876 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.61792;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.57643;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.53836;  
			 } 
		 else if (JetPt > 40.00000 && JetPt <= 60.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.70910 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.69380;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.65622;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.62465;  
			 } 
		 else if (JetPt > 60.00000 && JetPt <= 80.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.72449 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.71302;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.66783;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.63435;  
			 } 
		 else if (JetPt > 80.00000 && JetPt <= 100.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.72711 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.71749;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.67048;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.63058;  
			 } 
		 else if (JetPt > 100.00000 && JetPt <= 120.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.72216 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.71479;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.66111;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.61767;  
			 } 
		 else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.69840 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.68925;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.62157;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.57769;  
			 } 
		 } 
	 if (TaggerOP=="CSVv2_Tight") { 
		 if (JetPt > 20.00000 && JetPt <= 40.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.41772 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.37257;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.35870;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.32352;  
			 } 
		 else if (JetPt > 40.00000 && JetPt <= 60.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.50616 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.47165;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.43395;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.40770;  
			 } 
		 else if (JetPt > 60.00000 && JetPt <= 80.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.50623 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.48115;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.43438;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.41065;  
			 } 
		 else if (JetPt > 80.00000 && JetPt <= 100.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.50732 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.48478;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.43061;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.39712;  
			 } 
		 else if (JetPt > 100.00000 && JetPt <= 120.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.51038 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.49692;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.42013;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.38191;  
			 } 
		 else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.47171 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.45932;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.36676;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.33157;  
			 } 
		 } 
	 if (TaggerOP=="DeepCSV_Loose") { 
		 if (JetPt > 20.00000 && JetPt <= 40.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.84366 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.85016;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.82409;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.80497;  
			 } 
		 else if (JetPt > 40.00000 && JetPt <= 60.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.88838 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.89379;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.87035;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.84781;  
			 } 
		 else if (JetPt > 60.00000 && JetPt <= 80.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.90268 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.90705;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.88511;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.86281;  
			 } 
		 else if (JetPt > 80.00000 && JetPt <= 100.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.90869 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.91339;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.89391;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.87237;  
			 } 
		 else if (JetPt > 100.00000 && JetPt <= 120.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.91242 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.91639;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.89796;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.87712;  
			 } 
		 else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.91174 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.91336;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.89593;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.87537;  
			 } 
		 } 
	 if (TaggerOP=="DeepCSV_Medium") { 
		 if (JetPt > 20.00000 && JetPt <= 40.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.66809 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.67499;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.62431;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.58066;  
			 } 
		 else if (JetPt > 40.00000 && JetPt <= 60.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.74640 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.75139;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.70453;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.66210;  
			 } 
		 else if (JetPt > 60.00000 && JetPt <= 80.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.77173 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.77504;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.72724;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.67878;  
			 } 
		 else if (JetPt > 80.00000 && JetPt <= 100.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.78250 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.78440;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.73397;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.67849;  
			 } 
		 else if (JetPt > 100.00000 && JetPt <= 120.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.78670 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.78599;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.73197;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.67100;  
			 } 
		 else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.77353 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.76540;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.70406;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.63515;  
			 } 
		 } 
	 if (TaggerOP=="DeepCSV_Tight") { 
		 if (JetPt > 20.00000 && JetPt <= 40.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.48662 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.48405;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.41418;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.37179;  
			 } 
		 else if (JetPt > 40.00000 && JetPt <= 60.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.58765 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.58440;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.51423;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.46951;  
			 } 
		 else if (JetPt > 60.00000 && JetPt <= 80.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.61631 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.61186;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.54137;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.49049;  
			 } 
		 else if (JetPt > 80.00000 && JetPt <= 100.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.62352 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.61821;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.54560;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.48791;  
			 } 
		 else if (JetPt > 100.00000 && JetPt <= 120.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.62296 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.61460;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.53862;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.47629;  
			 } 
		 else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.59132 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.57624;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.49524;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.42693;  
			 } 
		 } 
	 return 1.; 
} 
float BTagSFUtil::TagEfficiencyC_2017(float JetPt, float JetEta) {


	 if (TaggerOP=="CSVv2_Loose") { 
		 if (JetPt > 20.00000 && JetPt <= 40.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.31723 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.34742;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.33357;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.33893;  
			 } 
		 else if (JetPt > 40.00000 && JetPt <= 60.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.27863 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.31327;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.34142;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.37417;  
			 } 
		 else if (JetPt > 60.00000 && JetPt <= 80.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.26731 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.30164;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.33196;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.36891;  
			 } 
		 else if (JetPt > 80.00000 && JetPt <= 100.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.27100 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.30781;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.34316;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.37883;  
			 } 
		 else if (JetPt > 100.00000 && JetPt <= 120.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.27830 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.31204;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.34175;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.37916;  
			 } 
		 else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.32315 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.35996;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.37733;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.41261;  
			 } 
		 } 
	 if (TaggerOP=="CSVv2_Medium") { 
		 if (JetPt > 20.00000 && JetPt <= 40.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.10732 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.10329;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.10083;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.09397;  
			 } 
		 else if (JetPt > 40.00000 && JetPt <= 60.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.08898 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.09152;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.10576;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.11064;  
			 } 
		 else if (JetPt > 60.00000 && JetPt <= 80.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.08309 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.08794;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.09660;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.10557;  
			 } 
		 else if (JetPt > 80.00000 && JetPt <= 100.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.08384 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.09219;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.09887;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.10872;  
			 } 
		 else if (JetPt > 100.00000 && JetPt <= 120.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.08672 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.09422;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.10256;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.11389;  
			 } 
		 else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.10032 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.11205;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.11257;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.12662;  
			 } 
		 } 
	 if (TaggerOP=="CSVv2_Tight") { 
		 if (JetPt > 20.00000 && JetPt <= 40.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.02926 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.02360;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.02578;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.02272;  
			 } 
		 else if (JetPt > 40.00000 && JetPt <= 60.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.02051 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.01918;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.02238;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.02469;  
			 } 
		 else if (JetPt > 60.00000 && JetPt <= 80.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.01653 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.01714;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.01882;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.02149;  
			 } 
		 else if (JetPt > 80.00000 && JetPt <= 100.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.01731 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.01848;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.01924;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.02178;  
			 } 
		 else if (JetPt > 100.00000 && JetPt <= 120.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.01930 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.02165;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.02031;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.02254;  
			 } 
		 else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.02366 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.02617;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.02036;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.02523;  
			 } 
		 } 
	 if (TaggerOP=="DeepCSV_Loose") { 
		 if (JetPt > 20.00000 && JetPt <= 40.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.33049 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.36373;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.37278;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.37025;  
			 } 
		 else if (JetPt > 40.00000 && JetPt <= 60.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.30328 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.34303;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.37117;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.37411;  
			 } 
		 else if (JetPt > 60.00000 && JetPt <= 80.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.29777 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.33721;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.36751;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.37307;  
			 } 
		 else if (JetPt > 80.00000 && JetPt <= 100.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.30628 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.34651;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.37795;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.38262;  
			 } 
		 else if (JetPt > 100.00000 && JetPt <= 120.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.32504 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.36096;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.39025;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.39794;  
			 } 
		 else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.36829 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.39636;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.41975;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.41965;  
			 } 
		 } 
	 if (TaggerOP=="DeepCSV_Medium") { 
		 if (JetPt > 20.00000 && JetPt <= 40.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.09688 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.10977;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.10773;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.09776;  
			 } 
		 else if (JetPt > 40.00000 && JetPt <= 60.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.08358 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.09934;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.11077;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.10526;  
			 } 
		 else if (JetPt > 60.00000 && JetPt <= 80.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.08527 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.10079;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.10931;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.10427;  
			 } 
		 else if (JetPt > 80.00000 && JetPt <= 100.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.09155 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.10784;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.11414;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.10876;  
			 } 
		 else if (JetPt > 100.00000 && JetPt <= 120.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.10124 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.11465;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.12003;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.11338;  
			 } 
		 else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.11469 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.12296;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.12393;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.11422;  
			 } 
		 } 
	 if (TaggerOP=="DeepCSV_Tight") { 
		 if (JetPt > 20.00000 && JetPt <= 40.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.02579 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.02785;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.02502;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.02156;  
			 } 
		 else if (JetPt > 40.00000 && JetPt <= 60.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.01802 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.02210;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.02340;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.02247;  
			 } 
		 else if (JetPt > 60.00000 && JetPt <= 80.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.01919 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.02248;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.02382;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.02253;  
			 } 
		 else if (JetPt > 80.00000 && JetPt <= 100.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.02046 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.02399;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.02546;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.02398;  
			 } 
		 else if (JetPt > 100.00000 && JetPt <= 120.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.02272 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.02553;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.02567;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.02512;  
			 } 
		 else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.02362 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.02619;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.02506;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.02308;  
			 } 
		 } 
	 return 1.; 
} 
float BTagSFUtil::TagEfficiencyLight_2017(float JetPt, float JetEta) {


	 if (TaggerOP=="CSVv2_Loose") { 
		 if (JetPt > 20.00000 && JetPt <= 40.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.10332 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.12538;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.13746;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.15632;  
			 } 
		 else if (JetPt > 40.00000 && JetPt <= 60.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.07312 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.08918;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.11685;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.15874;  
			 } 
		 else if (JetPt > 60.00000 && JetPt <= 80.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.05999 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.07360;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.11864;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.15895;  
			 } 
		 else if (JetPt > 80.00000 && JetPt <= 100.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.05317 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.06711;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.12539;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.16556;  
			 } 
		 else if (JetPt > 100.00000 && JetPt <= 120.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.04608 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.05805;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.10537;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.15632;  
			 } 
		 else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.04810 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.06451;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.12250;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.17714;  
			 } 
		 } 
	 if (TaggerOP=="CSVv2_Medium") { 
		 if (JetPt > 20.00000 && JetPt <= 40.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.01066 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.01151;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.01541;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01678;  
			 } 
		 else if (JetPt > 40.00000 && JetPt <= 60.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00843 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00921;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.01490;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01979;  
			 } 
		 else if (JetPt > 60.00000 && JetPt <= 80.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00746 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00806;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.01242;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01637;  
			 } 
		 else if (JetPt > 80.00000 && JetPt <= 100.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00690 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00747;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.01157;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01519;  
			 } 
		 else if (JetPt > 100.00000 && JetPt <= 120.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00618 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00663;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.01099;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01471;  
			 } 
		 else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00617 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00750;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.01132;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01840;  
			 } 
		 } 
	 if (TaggerOP=="CSVv2_Tight") { 
		 if (JetPt > 20.00000 && JetPt <= 40.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00131 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00123;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00151;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.00230;  
			 } 
		 else if (JetPt > 40.00000 && JetPt <= 60.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00089 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00086;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00123;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.00245;  
			 } 
		 else if (JetPt > 60.00000 && JetPt <= 80.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00078 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00083;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00105;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.00150;  
			 } 
		 else if (JetPt > 80.00000 && JetPt <= 100.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00080 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00090;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00101;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.00128;  
			 } 
		 else if (JetPt > 100.00000 && JetPt <= 120.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00090 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00100;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00103;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.00131;  
			 } 
		 else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00106 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00123;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00116;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.00195;  
			 } 
		 } 
	 if (TaggerOP=="DeepCSV_Loose") { 
		 if (JetPt > 20.00000 && JetPt <= 40.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.11391 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.13949;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.16831;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.18273;  
			 } 
		 else if (JetPt > 40.00000 && JetPt <= 60.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.09194 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.11569;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.15109;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.16804;  
			 } 
		 else if (JetPt > 60.00000 && JetPt <= 80.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.07613 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.09860;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.13497;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.15401;  
			 } 
		 else if (JetPt > 80.00000 && JetPt <= 100.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.06932 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.09236;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.12948;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.14946;  
			 } 
		 else if (JetPt > 100.00000 && JetPt <= 120.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.06911 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.09073;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.12805;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.14927;  
			 } 
		 else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.07438 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.09421;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.13860;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.16752;  
			 } 
		 } 
	 if (TaggerOP=="DeepCSV_Medium") { 
		 if (JetPt > 20.00000 && JetPt <= 40.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.01447 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.01891;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.02060;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01850;  
			 } 
		 else if (JetPt > 40.00000 && JetPt <= 60.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.01039 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.01424;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.01932;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.02024;  
			 } 
		 else if (JetPt > 60.00000 && JetPt <= 80.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00848 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.01169;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.01616;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01703;  
			 } 
		 else if (JetPt > 80.00000 && JetPt <= 100.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00800 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.01112;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.01509;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01591;  
			 } 
		 else if (JetPt > 100.00000 && JetPt <= 120.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00811 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.01071;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.01461;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01586;  
			 } 
		 else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00834 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.01022;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.01441;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01652;  
			 } 
		 } 
	 if (TaggerOP=="DeepCSV_Tight") { 
		 if (JetPt > 20.00000 && JetPt <= 40.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00158 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00193;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00200;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.00173;  
			 } 
		 else if (JetPt > 40.00000 && JetPt <= 60.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00121 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00160;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00200;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.00188;  
			 } 
		 else if (JetPt > 60.00000 && JetPt <= 80.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00107 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00140;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00174;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.00161;  
			 } 
		 else if (JetPt > 80.00000 && JetPt <= 100.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00101 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00130;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00156;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.00147;  
			 } 
		 else if (JetPt > 100.00000 && JetPt <= 120.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00098 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00122;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00144;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.00145;  
			 } 
		 else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00104 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00118;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00136;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.00139;  
			 } 
		 } 
	 return 1.; 
} 
