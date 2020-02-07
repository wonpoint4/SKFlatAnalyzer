/* These Efficiencies have been derived for 2018 JetTagging_TTLL_powheg.root events and should 
be used only for the same MC samples or for events with similar topology */  


float BTagSFUtil::TagEfficiencyB_2018(float JetPt, float JetEta) {


	 if (TaggerOP=="DeepCSV_Loose") { 
		 if (JetPt > 20.00000 && JetPt <= 40.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.84026 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.84087;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.81531;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.80086;  
			 } 
		 else if (JetPt > 40.00000 && JetPt <= 60.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.88573 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.88498;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.86131;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.85051;  
			 } 
		 else if (JetPt > 60.00000 && JetPt <= 80.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.90132 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.89971;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.87678;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.86512;  
			 } 
		 else if (JetPt > 80.00000 && JetPt <= 100.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.90882 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.90749;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.88591;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.87473;  
			 } 
		 else if (JetPt > 100.00000 && JetPt <= 120.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.91355 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.91209;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.89053;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.87936;  
			 } 
		 else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.91495 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.91177;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.88930;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.87852;  
			 } 
		 } 
	 if (TaggerOP=="DeepCSV_Medium") { 
		 if (JetPt > 20.00000 && JetPt <= 40.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.68201 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.68024;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.62691;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.57589;  
			 } 
		 else if (JetPt > 40.00000 && JetPt <= 60.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.75782 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.75391;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.70642;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.67134;  
			 } 
		 else if (JetPt > 60.00000 && JetPt <= 80.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.78322 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.77774;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.72840;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.68842;  
			 } 
		 else if (JetPt > 80.00000 && JetPt <= 100.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.79483 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.78774;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.73423;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.68831;  
			 } 
		 else if (JetPt > 100.00000 && JetPt <= 120.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.80048 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.79095;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.73243;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.68206;  
			 } 
		 else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.79025 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.77447;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.70741;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.65078;  
			 } 
		 } 
	 if (TaggerOP=="DeepCSV_Tight") { 
		 if (JetPt > 20.00000 && JetPt <= 40.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.49642 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.48788;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.41910;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.36429;  
			 } 
		 else if (JetPt > 40.00000 && JetPt <= 60.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.59758 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.58630;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.51837;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.47669;  
			 } 
		 else if (JetPt > 60.00000 && JetPt <= 80.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.62744 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.61430;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.54399;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.49937;  
			 } 
		 else if (JetPt > 80.00000 && JetPt <= 100.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.63676 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.62245;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.54685;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.49827;  
			 } 
		 else if (JetPt > 100.00000 && JetPt <= 120.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.63842 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.62178;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.54109;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.48725;  
			 } 
		 else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.61121 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.58933;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.50276;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.44157;  
			 } 
		 } 
	 return 1.; 
} 
float BTagSFUtil::TagEfficiencyC_2018(float JetPt, float JetEta) {


	 if (TaggerOP=="DeepCSV_Loose") { 
		 if (JetPt > 20.00000 && JetPt <= 40.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.33749 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.35796;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.36906;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.40047;  
			 } 
		 else if (JetPt > 40.00000 && JetPt <= 60.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.30994 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.33160;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.34565;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.38641;  
			 } 
		 else if (JetPt > 60.00000 && JetPt <= 80.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.30505 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.32655;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.33955;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.38161;  
			 } 
		 else if (JetPt > 80.00000 && JetPt <= 100.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.31615 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.33703;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.34665;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.38760;  
			 } 
		 else if (JetPt > 100.00000 && JetPt <= 120.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.33290 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.35051;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.36114;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.40364;  
			 } 
		 else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.38392 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.39398;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.39433;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.42191;  
			 } 
		 } 
	 if (TaggerOP=="DeepCSV_Medium") { 
		 if (JetPt > 20.00000 && JetPt <= 40.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.11452 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.12190;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.11313;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.10327;  
			 } 
		 else if (JetPt > 40.00000 && JetPt <= 60.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.09856 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.10703;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.11078;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.11737;  
			 } 
		 else if (JetPt > 60.00000 && JetPt <= 80.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.09930 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.10700;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.10894;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.11526;  
			 } 
		 else if (JetPt > 80.00000 && JetPt <= 100.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.10617 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.11391;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.11068;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.11510;  
			 } 
		 else if (JetPt > 100.00000 && JetPt <= 120.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.11505 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.12055;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.11490;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.11846;  
			 } 
		 else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.13243 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.13116;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.12047;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.11829;  
			 } 
		 } 
	 if (TaggerOP=="DeepCSV_Tight") { 
		 if (JetPt > 20.00000 && JetPt <= 40.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.02894 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.02980;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.02611;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.02200;  
			 } 
		 else if (JetPt > 40.00000 && JetPt <= 60.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.02079 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.02245;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.02333;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.02458;  
			 } 
		 else if (JetPt > 60.00000 && JetPt <= 80.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.02131 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.02307;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.02286;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.02427;  
			 } 
		 else if (JetPt > 80.00000 && JetPt <= 100.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.02328 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.02464;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.02384;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.02394;  
			 } 
		 else if (JetPt > 100.00000 && JetPt <= 120.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.02449 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.02615;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.02422;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.02445;  
			 } 
		 else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.02771 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.02814;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.02447;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.02336;  
			 } 
		 } 
	 return 1.; 
} 
float BTagSFUtil::TagEfficiencyLight_2018(float JetPt, float JetEta) {


	 if (TaggerOP=="DeepCSV_Loose") { 
		 if (JetPt > 20.00000 && JetPt <= 40.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.12510 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.14469;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.17764;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.23014;  
			 } 
		 else if (JetPt > 40.00000 && JetPt <= 60.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.09787 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.11202;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.13511;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.17894;  
			 } 
		 else if (JetPt > 60.00000 && JetPt <= 80.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.08298 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.09683;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.11883;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.15751;  
			 } 
		 else if (JetPt > 80.00000 && JetPt <= 100.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.07849 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.09136;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.11409;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.15093;  
			 } 
		 else if (JetPt > 100.00000 && JetPt <= 120.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.07925 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.09117;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.11431;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.15196;  
			 } 
		 else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.08771 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.09867;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.12628;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.17330;  
			 } 
		 } 
	 if (TaggerOP=="DeepCSV_Medium") { 
		 if (JetPt > 20.00000 && JetPt <= 40.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.02050 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.02314;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.02156;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01937;  
			 } 
		 else if (JetPt > 40.00000 && JetPt <= 60.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.01402 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.01690;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.01967;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.02216;  
			 } 
		 else if (JetPt > 60.00000 && JetPt <= 80.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.01145 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.01368;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.01610;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01852;  
			 } 
		 else if (JetPt > 80.00000 && JetPt <= 100.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.01075 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.01289;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.01515;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01729;  
			 } 
		 else if (JetPt > 100.00000 && JetPt <= 120.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.01076 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.01258;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.01465;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01748;  
			 } 
		 else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.01140 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.01214;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.01456;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01860;  
			 } 
		 } 
	 if (TaggerOP=="DeepCSV_Tight") { 
		 if (JetPt > 20.00000 && JetPt <= 40.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00179 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00199;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00189;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.00156;  
			 } 
		 else if (JetPt > 40.00000 && JetPt <= 60.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00130 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00155;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00172;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.00188;  
			 } 
		 else if (JetPt > 60.00000 && JetPt <= 80.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00117 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00136;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00147;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.00160;  
			 } 
		 else if (JetPt > 80.00000 && JetPt <= 100.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00118 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00137;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00141;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.00142;  
			 } 
		 else if (JetPt > 100.00000 && JetPt <= 120.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00118 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00132;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00132;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.00147;  
			 } 
		 else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00139 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00140;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00140;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.00162;  
			 } 
		 } 
	 return 1.; 
} 
