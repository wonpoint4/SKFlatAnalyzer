/* These Efficiencies have been derived for 2016 JetTagging_TTLL_powheg.root events and should 
be used only for the same MC samples or for events with similar topology */  


float BTagSFUtil::TagEfficiencyB_2016(float JetPt, float JetEta) {


	 if (TaggerOP=="DeepCSV_Loose") { 
		 if (JetPt > 20.00000 && JetPt <= 40.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.80661 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.80283;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.77734;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.74496;  
			 } 
		 else if (JetPt > 40.00000 && JetPt <= 60.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.85604 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.85368;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.83021;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.79445;  
			 } 
		 else if (JetPt > 60.00000 && JetPt <= 80.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.87184 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.86963;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.84629;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.80511;  
			 } 
		 else if (JetPt > 80.00000 && JetPt <= 100.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.88018 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.87781;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.85350;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.80893;  
			 } 
		 else if (JetPt > 100.00000 && JetPt <= 120.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.88483 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.88138;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.85658;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.81362;  
			 } 
		 else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.88858 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.88357;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.86002;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.81864;  
			 } 
		 } 
	 if (TaggerOP=="DeepCSV_Medium") { 
		 if (JetPt > 20.00000 && JetPt <= 40.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.63931 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.62995;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.57548;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.50201;  
			 } 
		 else if (JetPt > 40.00000 && JetPt <= 60.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.71456 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.70730;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.65897;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.58485;  
			 } 
		 else if (JetPt > 60.00000 && JetPt <= 80.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.73479 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.72787;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.67776;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.59509;  
			 } 
		 else if (JetPt > 80.00000 && JetPt <= 100.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.74413 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.73575;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.68231;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.58881;  
			 } 
		 else if (JetPt > 100.00000 && JetPt <= 120.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.74769 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.73773;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.68074;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.58481;  
			 } 
		 else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.74154 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.72711;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.66576;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.56218;  
			 } 
		 } 
	 if (TaggerOP=="DeepCSV_Tight") { 
		 if (JetPt > 20.00000 && JetPt <= 40.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.45517 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.44329;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.37812;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.30672;  
			 } 
		 else if (JetPt > 40.00000 && JetPt <= 60.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.54841 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.53675;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.47336;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.39783;  
			 } 
		 else if (JetPt > 60.00000 && JetPt <= 80.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.57140 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.55978;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.49510;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.41166;  
			 } 
		 else if (JetPt > 80.00000 && JetPt <= 100.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.57859 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.56634;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.49836;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.40512;  
			 } 
		 else if (JetPt > 100.00000 && JetPt <= 120.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.57839 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.56420;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.49360;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.39848;  
			 } 
		 else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.55881 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.54005;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.46457;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.36524;  
			 } 
		 } 
	 return 1.; 
} 
float BTagSFUtil::TagEfficiencyC_2016(float JetPt, float JetEta) {


	 if (TaggerOP=="DeepCSV_Loose") { 
		 if (JetPt > 20.00000 && JetPt <= 40.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.34206 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.35372;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.36309;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.37797;  
			 } 
		 else if (JetPt > 40.00000 && JetPt <= 60.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.30510 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.32183;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.34538;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.36324;  
			 } 
		 else if (JetPt > 60.00000 && JetPt <= 80.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.28532 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.30424;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.32911;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.34586;  
			 } 
		 else if (JetPt > 80.00000 && JetPt <= 100.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.28879 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.30890;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.33023;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.34650;  
			 } 
		 else if (JetPt > 100.00000 && JetPt <= 120.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.30180 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.31728;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.33681;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.34960;  
			 } 
		 else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.35528 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.36509;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.38476;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.38854;  
			 } 
		 } 
	 if (TaggerOP=="DeepCSV_Medium") { 
		 if (JetPt > 20.00000 && JetPt <= 40.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.11123 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.11399;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.10697;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.09249;  
			 } 
		 else if (JetPt > 40.00000 && JetPt <= 60.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.09239 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.09839;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.10106;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.09587;  
			 } 
		 else if (JetPt > 60.00000 && JetPt <= 80.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.08452 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.09055;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.09462;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.08762;  
			 } 
		 else if (JetPt > 80.00000 && JetPt <= 100.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.08594 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.09263;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.09501;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.08640;  
			 } 
		 else if (JetPt > 100.00000 && JetPt <= 120.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.09162 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.09584;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.09687;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.08904;  
			 } 
		 else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.11288 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.11378;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.11314;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.09532;  
			 } 
		 } 
	 if (TaggerOP=="DeepCSV_Tight") { 
		 if (JetPt > 20.00000 && JetPt <= 40.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.02587 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.02697;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.02409;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01863;  
			 } 
		 else if (JetPt > 40.00000 && JetPt <= 60.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.01943 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.02088;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.02144;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01971;  
			 } 
		 else if (JetPt > 60.00000 && JetPt <= 80.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.01769 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.01914;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.01936;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01800;  
			 } 
		 else if (JetPt > 80.00000 && JetPt <= 100.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.01866 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.01908;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.01963;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01762;  
			 } 
		 else if (JetPt > 100.00000 && JetPt <= 120.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.01985 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.02048;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.02064;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01786;  
			 } 
		 else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.02434 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.02463;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.02373;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01934;  
			 } 
		 } 
	 return 1.; 
} 
float BTagSFUtil::TagEfficiencyLight_2016(float JetPt, float JetEta) {


	 if (TaggerOP=="DeepCSV_Loose") { 
		 if (JetPt > 20.00000 && JetPt <= 40.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.16119 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.18085;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.20846;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.23488;  
			 } 
		 else if (JetPt > 40.00000 && JetPt <= 60.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.11066 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.13071;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.17097;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.19327;  
			 } 
		 else if (JetPt > 60.00000 && JetPt <= 80.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.08204 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.09958;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.13518;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.15249;  
			 } 
		 else if (JetPt > 80.00000 && JetPt <= 100.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.07337 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.09115;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.12316;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.14294;  
			 } 
		 else if (JetPt > 100.00000 && JetPt <= 120.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.07078 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.08663;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.11782;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.13980;  
			 } 
		 else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.08189 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.09771;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.13743;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.16930;  
			 } 
		 } 
	 if (TaggerOP=="DeepCSV_Medium") { 
		 if (JetPt > 20.00000 && JetPt <= 40.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.02631 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.02970;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.03149;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.02810;  
			 } 
		 else if (JetPt > 40.00000 && JetPt <= 60.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.01540 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.01892;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.02397;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.02376;  
			 } 
		 else if (JetPt > 60.00000 && JetPt <= 80.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.01114 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.01360;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.01716;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01738;  
			 } 
		 else if (JetPt > 80.00000 && JetPt <= 100.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.01003 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.01200;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.01507;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01528;  
			 } 
		 else if (JetPt > 100.00000 && JetPt <= 120.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00947 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.01132;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.01433;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01460;  
			 } 
		 else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.01078 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.01231;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.01734;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01696;  
			 } 
		 } 
	 if (TaggerOP=="DeepCSV_Tight") { 
		 if (JetPt > 20.00000 && JetPt <= 40.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00210 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00235;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00247;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.00213;  
			 } 
		 else if (JetPt > 40.00000 && JetPt <= 60.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00148 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00176;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00214;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.00218;  
			 } 
		 else if (JetPt > 60.00000 && JetPt <= 80.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00113 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00135;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00154;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.00157;  
			 } 
		 else if (JetPt > 80.00000 && JetPt <= 100.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00107 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00123;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00141;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.00141;  
			 } 
		 else if (JetPt > 100.00000 && JetPt <= 120.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00104 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00120;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00135;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.00130;  
			 } 
		 else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
			 if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00128 ;  
			 else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00141;  
			 else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00180;  
			 else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.00152;  
			 } 
		 } 
	 return 1.; 
} 
