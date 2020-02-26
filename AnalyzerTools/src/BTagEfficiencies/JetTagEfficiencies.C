double MCCorrection::GetMCJetTagEff(JetTagging::Tagger tagger, JetTagging::WP wp, int JetFlavor, double JetPt, double JetEta){

  if(JetPt<20) JetPt = 20.;
  if(JetPt>=1000.) JetPt = 999.;


  if(DataYear==2016){
    if(JetFlavor==5){
      if(tagger==JetTagging::DeepCSV){
        if(wp==JetTagging::Loose){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.77998;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.77645;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.75108;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.71958;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.84008;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.83726;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.81330;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.78016;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.86523;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.86291;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.83947;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.80113;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.87787;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.87549;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.85150;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.80776;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.88578;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.88195;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.85682;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.81460;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.88975;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.88465;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.86137;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.81942;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.88722;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.88174;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.86057;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.82150;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.89026;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.88187;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.86936;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.83630;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.88439;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.88443;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.89013;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.87146;
          }
        }
        if(wp==JetTagging::Medium){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.59791;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.58880;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.53363;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.45864;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.69096;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.68311;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.63313;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.56166;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.72689;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.71977;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.67063;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.59332;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.74154;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.73373;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.68144;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.59110;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.74785;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.73726;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.67877;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.58204;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.74518;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.73090;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.66911;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.56059;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.72638;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.70570;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.63877;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.52330;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.70430;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.67647;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.61817;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.51008;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.64484;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.63054;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.60664;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.55556;
          }
        }
        if(wp==JetTagging::Tight){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.40510;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.39471;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.33214;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.26151;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.51866;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.50696;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.44270;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.37062;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.56306;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.55132;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.48720;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.40843;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.57683;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.56484;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.49803;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.40789;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.57718;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.56229;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.49059;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.39442;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.56440;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.54614;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.46832;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.36255;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.52738;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.50036;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.41589;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.30931;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.47890;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.44286;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.37305;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.28048;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.38867;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.36696;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.34455;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.31808;
          }
        }
      }
      if(tagger==JetTagging::DeepJet){
        if(wp==JetTagging::Loose){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.86215;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.86993;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.85923;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.86572;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.87413;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.87844;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.86092;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.85709;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.87485;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.87852;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.85843;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.84957;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.87692;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.88029;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.85902;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.84967;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.88819;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.89122;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.86985;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.86365;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.90691;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.90968;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.89190;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.88688;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.92522;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.92775;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.91384;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.90892;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.93548;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.93630;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.92941;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.92444;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.94103;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.94733;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.94584;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.95643;
          }
        }
        if(wp==JetTagging::Medium){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.68658;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.69108;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.64258;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.59257;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.72567;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.72846;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.68519;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.64505;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.73482;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.73671;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.69433;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.65947;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.73984;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.74156;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.69843;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.66487;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.75646;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.75765;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.71290;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.68479;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.78456;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.78546;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.74269;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.71440;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.80980;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.80887;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.76591;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.73380;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.80650;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.80079;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.75638;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.71652;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.76666;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.76448;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.71897;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.71895;
          }
        }
        if(wp==JetTagging::Tight){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.49579;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.49548;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.42721;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.35862;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.56683;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.56497;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.49994;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.44449;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.59244;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.59044;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.52744;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.47640;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.60261;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.60120;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.53794;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.48622;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.61920;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.61665;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.55101;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.50407;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.64538;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.64220;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.57475;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.52650;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.66298;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.65618;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.58512;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.53627;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.63458;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.61875;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.54018;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.47754;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.55418;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.53895;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.43749;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.41394;
          }
        }
      }
      if(tagger==JetTagging::CSVv2){
        if(wp==JetTagging::Loose){
        }
        if(wp==JetTagging::Medium){
        }
        if(wp==JetTagging::Tight){
        }
      }
    }
    if(JetFlavor==4){
      if(tagger==JetTagging::DeepCSV){
        if(wp==JetTagging::Loose){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.34547;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.35605;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.36165;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.37606;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.32676;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.34214;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.36182;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.37769;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.29172;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.30982;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.33358;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.35056;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.28676;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.30591;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.32903;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.34596;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.30644;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.32155;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.33955;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.35344;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.33972;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.35260;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.37070;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.38272;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.37799;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.38355;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.40966;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.41026;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.42821;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.43244;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.46733;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.46043;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.49946;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.51099;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.57358;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.50142;
          }
        }
        if(wp==JetTagging::Medium){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.11427;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.11559;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.10544;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.08883;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.10177;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.10750;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.10896;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.09944;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.08702;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.09372;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.09627;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.09075;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.08526;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.09118;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.09378;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.08721;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.09388;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.09711;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.09842;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.08966;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.10704;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.11134;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.10983;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.09701;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.12151;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.12129;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.12112;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.09645;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.14109;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.13419;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.13441;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.09975;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.15762;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.14278;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.16890;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.10826;
          }
        }
        if(wp==JetTagging::Tight){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.02748;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.02811;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.02413;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.01783;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.02199;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.02380;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.02372;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02024;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.01812;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.01961;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01993;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.01862;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.01824;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.01897;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01934;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.01801;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.02018;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.02088;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.02092;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.01831;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.02270;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.02394;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.02331;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02030;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.02692;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.02647;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.02623;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.01869;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.03054;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.02966;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.02557;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.01723;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.03496;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.02588;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.02899;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.01994;
          }
        }
      }
      if(tagger==JetTagging::DeepJet){
        if(wp==JetTagging::Loose){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.49820;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.53148;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.55322;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.62199;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.41267;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.44081;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.45948;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.50768;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.32929;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.35723;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.38286;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.42647;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.29339;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.31889;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.34406;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.39350;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.30038;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.32490;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.34128;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.39168;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.34541;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.36915;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.38489;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.43853;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.42393;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.44777;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.46177;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.50386;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.50654;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.52649;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.55899;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.58415;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.58932;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.62900;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.70569;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.67806;
          }
        }
        if(wp==JetTagging::Medium){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.19609;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.20972;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.19611;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.19084;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.13853;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.15461;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.15925;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.16354;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.09810;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.10971;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.11522;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.12700;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.08592;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.09715;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.10071;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.11209;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.08984;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.09883;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.09997;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.11259;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.10589;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.11659;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.12139;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.13328;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.14109;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.15196;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.15675;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.16724;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.17685;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.18505;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.18660;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.18535;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.21252;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.23017;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.23300;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.19658;
          }
        }
        if(wp==JetTagging::Tight){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.06275;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.06689;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.05741;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.04545;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.03900;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.04437;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.04228;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.03794;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.02741;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.03081;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.02970;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02846;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.02455;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.02744;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.02608;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02552;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.02456;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.02718;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.02592;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02584;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.02849;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.03127;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.03161;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.03164;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.03706;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.03975;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.04100;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.04095;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.04766;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.04632;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.04615;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.04294;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.05190;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.05232;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.05072;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.03704;
          }
        }
      }
      if(tagger==JetTagging::CSVv2){
        if(wp==JetTagging::Loose){
        }
        if(wp==JetTagging::Medium){
        }
        if(wp==JetTagging::Tight){
        }
      }
    }
    if(JetFlavor<=3){
      if(tagger==JetTagging::DeepCSV){
        if(wp==JetTagging::Loose){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.16503;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.18289;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.20612;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.23285;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.14089;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.16320;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.20283;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.22715;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.09229;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.11049;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.14792;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.16650;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.07551;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.09332;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.12657;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.14498;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.07083;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.08657;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.11805;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.14196;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.07571;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.09059;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.12779;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.16265;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.09141;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.10772;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.15560;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.19612;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.13442;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.15410;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.22253;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.25451;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.23121;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.25395;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.35659;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.32475;
          }
        }
        if(wp==JetTagging::Medium){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.02772;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.03041;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.03076;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02691;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.02117;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.02551;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.03087;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02890;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.01254;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.01545;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01951;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.01965;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.01032;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.01244;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01561;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.01579;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00943;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.01124;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01440;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.01491;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00996;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.01139;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01606;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.01676;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.01172;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.01325;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01972;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.01841;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.01815;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.01968;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.03008;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02194;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.03200;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.02928;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.04806;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.03142;
          }
        }
        if(wp==JetTagging::Tight){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00217;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00238;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00236;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00195;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00182;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00216;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00260;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00246;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00127;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00149;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00176;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00184;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00108;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00126;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00144;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00143;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00104;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00121;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00135;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00134;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00116;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00131;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00164;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00154;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00149;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00151;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00220;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00158;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00237;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00238;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00337;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00190;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00373;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00351;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00708;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00303;
          }
        }
      }
      if(tagger==JetTagging::DeepJet){
        if(wp==JetTagging::Loose){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.39123;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.42728;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.48013;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.60114;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.27506;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.29925;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.34908;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.42945;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.15159;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.17425;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.23662;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.31044;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.10188;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.11894;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.17332;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.24947;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.08359;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.09752;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.14507;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.22070;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.07993;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.09461;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.14494;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.23055;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.10313;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.12426;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.18484;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.27818;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.18168;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.21344;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.29802;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.37336;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.33473;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.38117;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.52047;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.54182;
          }
        }
        if(wp==JetTagging::Medium){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.10588;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.11312;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.11331;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.13119;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.04855;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.05833;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.07789;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.10192;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.01590;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.02042;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.03153;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.05309;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.01053;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.01330;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01942;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.03536;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00904;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.01152;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01694;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.03110;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00962;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.01233;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.02012;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.03715;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.01356;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.01739;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.02957;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.05179;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.02619;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.03070;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.05167;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.07475;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.06680;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.07323;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.10775;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.11639;
          }
        }
        if(wp==JetTagging::Tight){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.01579;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.01751;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01614;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.01557;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00545;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00675;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00798;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00956;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00229;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00283;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00332;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00434;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00157;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00195;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00233;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00325;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00131;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00168;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00210;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00328;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00146;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00177;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00275;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00509;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00199;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00255;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00447;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00846;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00436;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00513;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00868;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.01277;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.01320;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.01309;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01317;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.01249;
          }
        }
      }
      if(tagger==JetTagging::CSVv2){
        if(wp==JetTagging::Loose){
        }
        if(wp==JetTagging::Medium){
        }
        if(wp==JetTagging::Tight){
        }
      }
    }
  }
  if(DataYear==2017){
    if(JetFlavor==5){
      if(tagger==JetTagging::DeepCSV){
        if(wp==JetTagging::Loose){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.81892;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.82658;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.79971;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.78155;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.87375;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.87960;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.85566;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.83479;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.89712;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.90195;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.87893;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.85608;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.90730;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.91177;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.89164;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.86990;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.91330;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.91686;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.89867;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.87788;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.91383;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.91551;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.89784;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.87557;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.90460;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.90308;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.88393;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.86311;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.88828;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.88417;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.86757;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.84988;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.87256;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.86076;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.85090;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.84161;
          }
        }
        if(wp==JetTagging::Medium){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.62596;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.63451;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.58261;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.53685;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.72026;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.72630;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.67900;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.63761;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.76186;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.76580;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.71846;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.67358;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.78000;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.78228;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.73285;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.67932;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.78715;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.78558;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.73004;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.66710;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.78089;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.77166;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.70815;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.63178;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.74526;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.72609;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.65144;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.56854;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.67340;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.64797;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.57424;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.50626;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.54038;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.52269;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.45794;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.44720;
          }
        }
        if(wp==JetTagging::Tight){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.43325;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.43209;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.36478;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.32257;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.55396;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.55109;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.48107;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.43821;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.60601;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.60213;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.53169;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.48437;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.62220;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.61715;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.54540;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.48960;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.62182;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.61234;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.53508;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.47046;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.60213;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.58565;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.50014;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.42273;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.53934;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.51266;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.41976;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.33880;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.41415;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.38535;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.30501;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.24675;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.22300;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.21149;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.18220;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.15217;
          }
        }
      }
      if(tagger==JetTagging::DeepJet){
        if(wp==JetTagging::Loose){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.91796;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.92407;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.91814;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.92902;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.92691;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.93156;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.92307;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.92995;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.92597;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.93018;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.92146;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.92817;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.92577;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.92950;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.92017;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.92900;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.93242;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.93545;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.92501;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.93631;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.94469;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.94686;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.93698;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.94618;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.95560;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.95700;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.94702;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.95249;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.95915;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.95866;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.95157;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.95471;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.95896;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.95904;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.95506;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.98447;
          }
        }
        if(wp==JetTagging::Medium){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.77568;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.78422;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.75544;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.75859;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.79739;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.80339;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.77881;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.78658;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.80033;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.80382;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.77722;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.79057;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.80317;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.80532;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.77548;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.79467;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.81755;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.81845;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.78578;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.80755;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.84304;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.84234;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.80939;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.82635;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.86321;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.86098;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.82363;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.83373;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.85899;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.85110;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.80104;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.80430;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.82419;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.80797;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.72912;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.79193;
          }
        }
        if(wp==JetTagging::Tight){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.56953;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.57310;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.51990;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.51489;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.63359;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.63455;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.58260;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.58629;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.65673;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.65610;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.60252;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.61233;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.66644;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.66494;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.60841;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.62340;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.68294;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.67943;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.61907;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.63754;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.71082;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.70473;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.64103;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.65442;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.72683;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.71544;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.64280;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.65175;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.69521;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.67657;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.57979;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.58701;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.61703;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.58362;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.46128;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.53727;
          }
        }
      }
      if(tagger==JetTagging::CSVv2){
        if(wp==JetTagging::Loose){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.80005;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.80493;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.76116;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.74379;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.86134;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.86497;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.83648;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.82341;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.87990;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.88238;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.85846;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.84631;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.88733;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.88973;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.86736;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.85153;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.88719;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.88727;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.86130;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.83987;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.88140;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.87851;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.85314;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.82512;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.86223;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.85503;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.82575;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.79804;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.82488;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.82210;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.79113;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.77782;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.78885;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.80066;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.73064;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.72050;
          }
        }
        if(wp==JetTagging::Medium){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.59054;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.56919;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.52255;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.48491;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.68929;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.67219;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.63667;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.60472;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.72030;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.70665;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.66635;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.63165;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.72737;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.71698;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.66924;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.63097;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.72217;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.71413;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.66019;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.61466;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.70365;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.69347;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.62336;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.57269;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.65976;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.64574;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.53906;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.49107;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.57324;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.56835;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.47029;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.45358;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.46452;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.47134;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.35378;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.40994;
          }
        }
        if(wp==JetTagging::Tight){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.37208;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.32852;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.31076;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.27576;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.47769;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.43749;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.41421;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.38548;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.51001;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.48020;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.44025;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.41438;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.50712;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.48354;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.43027;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.39943;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.50885;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.49451;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.41774;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.37714;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.47571;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.46295;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.36529;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.32380;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.41817;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.40475;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.26774;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.23484;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.31608;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.31188;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.20695;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.20602;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.22307;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.22785;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.11995;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.17391;
          }
        }
      }
    }
    if(JetFlavor==4){
      if(tagger==JetTagging::DeepCSV){
        if(wp==JetTagging::Loose){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.33403;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.36516;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.36809;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.36421;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.31714;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.35456;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.37841;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.38014;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.29758;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.33869;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.36890;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.37010;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.30313;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.34290;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.37314;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.37957;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.33063;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.36674;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.39576;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.40227;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.36827;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.40034;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.42482;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.42174;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.38924;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.41359;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.43506;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.42875;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.39397;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.39919;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.41097;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.42137;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.34370;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.36166;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.37481;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.46975;
          }
        }
        if(wp==JetTagging::Medium){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.10068;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.11138;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.10466;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.09268;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.08735;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.10395;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.11246;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.10692;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.08409;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.09913;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.11025;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.10378;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.08924;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.10551;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.11253;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.10737;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.10433;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.11753;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.12209;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.11556;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.11958;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.13076;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.13287;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.11971;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.12325;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.12818;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.12612;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.11164;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.10160;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.09693;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.08731;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.07719;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.04887;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.04081;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.03379;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.03559;
          }
        }
        if(wp==JetTagging::Tight){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.02828;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.02963;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.02517;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02089;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.01988;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.02358;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.02425;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02267;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.01852;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.02199;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.02313;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02246;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.02004;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.02366;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.02519;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02345;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.02306;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.02641;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.02618;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02536;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.02565;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.02912;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.02808;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02549;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.02493;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.02616;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.02435;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02015;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.01657;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.01524;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01143;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00911;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00413;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00187;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00154;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00000;
          }
        }
      }
      if(tagger==JetTagging::DeepJet){
        if(wp==JetTagging::Loose){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.57724;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.60951;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.64674;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.69873;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.48654;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.52005;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.55769;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.61360;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.38626;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.42168;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.47169;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.53487;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.34534;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.37823;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.42305;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.49216;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.35280;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.38186;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.41970;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.49268;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.40560;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.43407;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.46511;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.52533;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.48953;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.51504;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.53545;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.58762;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.56178;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.58399;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.59768;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.63942;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.61409;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.64395;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.69508;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.74377;
          }
        }
        if(wp==JetTagging::Medium){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.23965;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.25869;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.26076;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.28934;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.16524;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.18697;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.21349;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.24974;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.11642;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.13224;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.15527;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.18798;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.10263;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.11528;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.13055;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.16184;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.10445;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.11545;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.12780;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.16241;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.12511;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.13892;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.14881;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.18444;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.16616;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.18028;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.18980;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.21859;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.20337;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.21103;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.20620;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.22728;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.22691;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.25047;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.22427;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.27046;
          }
        }
        if(wp==JetTagging::Tight){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.07158;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.07663;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.07255;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.07735;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.04333;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.05027;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.05509;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.06494;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.03085;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.03564;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.03922;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.04772;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.02689;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.03126;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.03290;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.03960;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.02705;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.03023;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.03211;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.03999;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.03251;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.03562;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.03528;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.04706;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.04129;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.04463;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.04733;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.05192;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.04929;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.04834;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.04324;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.05012;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.05363;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.05391;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.03379;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.06406;
          }
        }
      }
      if(tagger==JetTagging::CSVv2){
        if(wp==JetTagging::Loose){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.31238;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.33912;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.31795;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.31929;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.30752;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.34292;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.35748;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.37970;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.27133;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.30673;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.33484;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.36432;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.26999;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.30594;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.33745;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.37359;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.28580;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.32048;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.34795;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.38583;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.31607;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.35227;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.37598;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.41089;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.34334;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.38024;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.39215;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.42865;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.36163;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.40291;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.40345;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.43739;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.34783;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.39611;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.37250;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.46619;
          }
        }
        if(wp==JetTagging::Medium){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.10747;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.10134;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.09310;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.08426;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.10013;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.10069;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.11295;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.11348;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.08554;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.08923;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.10138;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.10544;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.08340;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.09109;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.09674;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.10684;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.08885;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.09761;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.10627;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.11869;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.09815;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.10961;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.11418;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.12560;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.10773;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.11897;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.11157;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.12858;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.11132;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.12555;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.11261;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.12848;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.09965;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.13366;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.08218;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.13879;
          }
        }
        if(wp==JetTagging::Tight){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.03191;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.02518;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.02516;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02103;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.02297;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.02019;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.02537;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02583;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.01844;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.01806;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.02081;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02296;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.01712;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.01798;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01865;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02113;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.02039;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.02233;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.02119;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02359;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.02326;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.02561;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.02031;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02532;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.02535;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.02835;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01878;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02531;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.02390;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.02854;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01913;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02525;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.03650;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.04081;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01075;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00712;
          }
        }
      }
    }
    if(JetFlavor<=3){
      if(tagger==JetTagging::DeepCSV){
        if(wp==JetTagging::Loose){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.11572;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.14148;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.16771;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.18285;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.10475;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.12945;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.16391;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.17829;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.08268;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.10568;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.14142;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.15972;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.07096;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.09371;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.13070;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.15022;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.06966;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.09069;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.12889;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.15103;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.07205;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.09197;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.13592;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.16429;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.07809;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.09839;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.14653;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.18372;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.09361;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.11065;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.16911;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.21721;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.12618;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.13445;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.21238;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.26976;
          }
        }
        if(wp==JetTagging::Medium){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.01483;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.01908;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01941;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.01683;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.01265;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.01711;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.02197;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02145;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00913;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.01262;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01734;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.01847;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00815;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.01122;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01532;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.01610;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00810;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.01064;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01453;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.01594;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00823;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.01015;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01465;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.01695;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00880;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.01016;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01438;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.01677;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00926;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00971;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01337;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.01553;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00619;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00613;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00879;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00810;
          }
        }
        if(wp==JetTagging::Tight){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00169;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00198;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00191;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00159;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00134;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00176;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00214;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00198;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00112;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00148;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00187;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00175;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00103;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00132;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00160;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00149;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00099;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00122;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00142;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00142;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00104;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00117;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00139;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00144;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00112;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00127;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00140;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00139;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00109;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00100;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00109;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00126;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00032;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00021;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00039;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00000;
          }
        }
      }
      if(tagger==JetTagging::DeepJet){
        if(wp==JetTagging::Loose){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.44486;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.47981;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.56296;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.65045;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.32008;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.34512;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.41368;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.50948;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.17082;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.19518;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.27538;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.38937;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.11260;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.13153;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.19914;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.30357;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.09120;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.10663;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.16407;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.25928;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.08428;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.10185;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.15929;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.25509;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.10369;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.12835;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.19548;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.29043;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.17563;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.20631;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.29140;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.37942;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.27997;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.32627;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.45748;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.51587;
          }
        }
        if(wp==JetTagging::Medium){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.11313;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.12436;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.13251;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.18653;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.05629;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.06743;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.09756;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.14951;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.01903;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.02468;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.04301;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.08443;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.01186;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.01529;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.02515;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.05157;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00951;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.01238;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01993;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.03965;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00953;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.01244;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.02105;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.03965;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.01264;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.01631;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.02767;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.04764;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.01941;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.02394;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.03980;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.06019;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.03838;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.03979;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.06601;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.08364;
          }
        }
        if(wp==JetTagging::Tight){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.01803;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.02048;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01931;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02745;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00672;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00847;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01184;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02118;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00275;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00345;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00500;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00939;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00176;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00220;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00306;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00553;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00134;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00165;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00236;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00434;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00128;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00164;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00247;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00467;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00152;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00195;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00327;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00628;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00239;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00279;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00450;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00788;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00499;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00448;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00479;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00876;
          }
        }
      }
      if(tagger==JetTagging::CSVv2){
        if(wp==JetTagging::Loose){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.10294;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.12442;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.13411;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.15104;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.09185;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.11232;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.13558;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.16508;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.06614;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.08071;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.11086;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.15491;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.05529;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.06895;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.12535;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.16366;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.04606;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.05849;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.10700;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.15819;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.04479;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.05982;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.12101;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.17730;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.05046;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.07093;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.13723;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.19668;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.06832;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.09592;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.15170;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.20833;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.10331;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.13367;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.16299;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.17648;
          }
        }
        if(wp==JetTagging::Medium){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.01042;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.01119;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01390;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.01474;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.01004;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.01097;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01742;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02107;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00786;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00857;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01334;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.01737;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00709;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00765;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01167;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.01519;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00609;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00665;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01119;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.01521;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00583;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00682;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01078;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.01791;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00620;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00842;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01042;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02052;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00882;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.01230;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01572;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02931;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.01369;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.01847;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01690;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.01511;
          }
        }
        if(wp==JetTagging::Tight){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00147;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00135;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00147;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00203;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00099;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00096;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00147;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00287;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00083;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00085;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00112;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00172;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00079;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00086;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00101;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00129;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00092;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00101;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00106;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00142;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00097;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00108;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00106;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00185;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00117;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00147;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00118;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00223;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00178;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00223;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00191;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00381;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00286;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00384;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00231;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00109;
          }
        }
      }
    }
  }
  if(DataYear==2018){
    if(JetFlavor==5){
      if(tagger==JetTagging::DeepCSV){
        if(wp==JetTagging::Loose){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.81731;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.81863;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.79552;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.79256;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.87172;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.87088;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.84788;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.84135;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.89547;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.89416;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.87079;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.86038;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.90725;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.90566;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.88332;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.87228;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.91481;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.91295;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.89129;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.87963;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.91697;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.91382;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.89073;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.87871;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.90948;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.90434;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.87989;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.86597;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.89942;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.89184;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.87180;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.85849;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.88976;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.88176;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.85958;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.85932;
          }
        }
        if(wp==JetTagging::Medium){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.64435;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.64309;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.59196;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.55903;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.73412;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.73033;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.68408;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.65461;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.77390;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.76902;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.72054;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.68557;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.79242;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.78562;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.73324;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.68873;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.80141;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.79080;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.73104;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.67705;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.79652;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.78059;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.71171;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.64617;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.76451;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.74039;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.66093;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.58451;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.70597;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.67299;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.59465;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.53480;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.58540;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.56283;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.50755;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.45627;
          }
        }
        if(wp==JetTagging::Tight){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.44707;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.43956;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.37763;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.34294;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.56567;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.55440;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.48884;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.45441;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.61752;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.60486;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.53563;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.49554;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.63506;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.62114;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.54706;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.49939;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.63774;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.62012;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.53880;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.48036;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.62106;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.59908;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.50910;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.43463;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.56297;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.53198;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.43462;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.35127;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.44853;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.41239;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.32560;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.26858;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.25862;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.23737;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.19849;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.18251;
          }
        }
      }
      if(tagger==JetTagging::DeepJet){
        if(wp==JetTagging::Loose){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.90350;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.89771;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.88185;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.88890;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.92639;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.91967;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.90268;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.90322;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.93698;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.93072;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.91409;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.90974;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.94226;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.93624;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.91988;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.91407;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.94530;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.93966;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.92307;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.91687;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.94768;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.94219;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.92600;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.92099;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.94867;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.94295;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.92720;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.92193;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.94922;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.94275;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.92973;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.92230;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.95019;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.94326;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.93963;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.95437;
          }
        }
        if(wp==JetTagging::Medium){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.75455;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.73487;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.67707;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.65673;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.81519;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.79801;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.74967;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.73139;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.84287;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.82739;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.78230;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.76125;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.85616;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.84210;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.79798;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.77521;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.86380;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.85060;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.80624;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.78379;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.86701;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.85440;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.81078;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.79006;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.86344;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.84938;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.80072;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.78314;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.84862;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.82986;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.77567;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.75056;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.80714;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.78402;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.71818;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.76046;
          }
        }
        if(wp==JetTagging::Tight){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.51269;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.48877;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.41524;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.39275;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.62874;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.60240;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.52760;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.50520;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.68400;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.65949;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.58546;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.55872;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.70940;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.68748;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.61437;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.58469;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.72322;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.70249;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.62912;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.60009;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.72747;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.70898;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.63592;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.60692;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.71412;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.69143;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.61088;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.58413;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.66348;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.63253;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.54069;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.51472;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.54865;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.51718;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.42257;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.50951;
          }
        }
      }
      if(tagger==JetTagging::CSVv2){
        if(wp==JetTagging::Loose){
        }
        if(wp==JetTagging::Medium){
        }
        if(wp==JetTagging::Tight){
        }
      }
    }
    if(JetFlavor==4){
      if(tagger==JetTagging::DeepCSV){
        if(wp==JetTagging::Loose){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.34276;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.36197;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.37065;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.40591;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.32172;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.34342;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.35767;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.39445;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.30528;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.32682;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.33967;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.37890;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.31244;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.33344;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.34370;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.38576;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.34223;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.35831;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.36743;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.40426;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.38091;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.39125;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.39605;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.42428;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.40770;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.41477;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.40268;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.43533;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.41373;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.41834;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.40360;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.43191;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.41026;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.39487;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.42683;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.45136;
          }
        }
        if(wp==JetTagging::Medium){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.11888;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.12403;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.11306;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.10987;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.10403;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.11392;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.11633;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.12072;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.09783;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.10616;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.10843;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.11525;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.10417;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.11128;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.10987;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.11419;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.11986;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.12357;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.11748;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.12018;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.13715;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.13806;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.12825;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.12359;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.14108;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.13687;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.11886;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.11394;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.12144;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.11176;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.09416;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.08750;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.07486;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.04402;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.04228;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.03502;
          }
        }
        if(wp==JetTagging::Tight){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.03171;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.03173;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.02706;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02408;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.02276;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.02482;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.02506;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02531;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.02075;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.02233;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.02248;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02410;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.02279;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.02418;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.02349;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02380;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.02598;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.02728;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.02537;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02545;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.03057;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.03064;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.02700;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02455;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.02781;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.02897;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.02251;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02134;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.02020;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.01812;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01424;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.01130;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00686;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00513;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00000;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00389;
          }
        }
      }
      if(tagger==JetTagging::DeepJet){
        if(wp==JetTagging::Loose){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.50905;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.50457;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.50338;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.55451;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.45301;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.44987;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.44693;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.48213;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.39407;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.39312;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.39222;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.41866;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.37494;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.37452;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.36822;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.38745;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.38372;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.38139;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.37177;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.38296;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.41544;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.40692;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.39896;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.40606;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.46343;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.44958;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.43480;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.45227;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.52115;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.50880;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.49839;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.51596;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.56456;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.58419;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.62114;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.64981;
          }
        }
        if(wp==JetTagging::Medium){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.19986;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.18723;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.15892;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.15608;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.16601;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.16135;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.15017;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.15309;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.14254;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.13721;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.12882;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.12935;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.13957;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.13531;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.12422;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.12092;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.14585;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.14034;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.12916;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.12310;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.16261;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.15671;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.14600;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.13572;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.18726;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.17851;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.16041;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.15322;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.21809;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.20806;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.17867;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.16995;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.23962;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.21282;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.21545;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.15564;
          }
        }
        if(wp==JetTagging::Tight){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.04096;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.03830;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.03120;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02930;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.03130;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.03001;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.02643;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02671;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.02804;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.02631;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.02359;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02353;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.02833;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.02771;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.02360;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02349;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.03063;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.02959;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.02704;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02636;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.03483;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.03399;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.03139;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02844;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.04096;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.03971;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.03357;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.03308;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.04993;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.04625;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.03895;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.03910;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.05427;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.04530;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.03902;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.01167;
          }
        }
      }
      if(tagger==JetTagging::CSVv2){
        if(wp==JetTagging::Loose){
        }
        if(wp==JetTagging::Medium){
        }
        if(wp==JetTagging::Tight){
        }
      }
    }
    if(JetFlavor<=3){
      if(tagger==JetTagging::DeepCSV){
        if(wp==JetTagging::Loose){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.12958;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.15092;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.18250;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.23031;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.11057;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.12591;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.14901;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.18949;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.08838;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.10256;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.12403;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.16018;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.07925;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.09265;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.11466;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.15089;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.07990;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.09173;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.11447;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.15344;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.08443;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.09491;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.12083;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.16863;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.09421;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.10456;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.13490;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.19391;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.11941;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.12881;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.17347;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.24416;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.17106;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.17961;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.24572;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.31142;
          }
        }
        if(wp==JetTagging::Medium){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.02135;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.02350;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.02126;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02121;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.01748;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.02058;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.02254;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02365;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.01218;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.01481;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01726;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.01932;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.01091;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.01305;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01529;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.01731;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.01087;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.01249;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01452;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.01750;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.01116;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.01183;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01411;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.01834;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.01186;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.01218;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01462;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.01966;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.01380;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.01329;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01759;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02162;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.01171;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00967;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01448;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.01510;
          }
        }
        if(wp==JetTagging::Tight){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00191;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00209;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00193;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00177;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00148;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00172;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00189;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00198;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00122;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00143;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00156;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00167;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00117;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00136;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00141;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00145;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00121;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00133;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00134;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00151;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00138;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00142;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00143;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00161;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00157;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00152;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00138;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00168;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00175;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00138;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00140;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00154;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00061;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00049;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00048;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00067;
          }
        }
      }
      if(tagger==JetTagging::DeepJet){
        if(wp==JetTagging::Loose){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.38928;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.39652;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.42683;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.50429;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.26790;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.28002;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.29937;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.34601;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.14340;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.15187;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.17547;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.21146;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.10197;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.10708;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.12427;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.14883;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.08479;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.08801;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.10486;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.12405;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.07670;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.07930;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.10050;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.12512;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.08170;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.08460;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.11724;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.15674;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.11430;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.12325;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.18902;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.25016;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.22868;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.26852;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.40100;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.44625;
          }
        }
        if(wp==JetTagging::Medium){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.06301;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.06032;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.04914;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.05329;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.03676;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.03787;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.03857;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.04690;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.01635;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.01727;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01827;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02080;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.01089;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.01148;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01244;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.01400;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00864;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00913;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01027;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.01263;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00799;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00837;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01054;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.01414;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00900;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00926;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.01343;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.02026;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.01312;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.01424;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.02428;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.03442;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.03233;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.03513;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.06221;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.06286;
          }
        }
        if(wp==JetTagging::Tight){
          if(JetPt >= 20.0 && JetPt < 30.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00453;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00415;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00290;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00284;
          }
          if(JetPt >= 30.0 && JetPt < 50.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00273;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00273;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00240;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00272;
          }
          if(JetPt >= 50.0 && JetPt < 70.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00138;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00142;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00138;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00164;
          }
          if(JetPt >= 70.0 && JetPt < 100.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00100;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00105;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00105;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00122;
          }
          if(JetPt >= 100.0 && JetPt < 140.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00082;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00085;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00090;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00119;
          }
          if(JetPt >= 140.0 && JetPt < 200.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00075;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00085;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00098;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00151;
          }
          if(JetPt >= 200.0 && JetPt < 300.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00089;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00096;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00131;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00224;
          }
          if(JetPt >= 300.0 && JetPt < 600.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00139;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00163;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00265;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00403;
          }
          if(JetPt >= 600.0 && JetPt < 1000.0){
            if(fabs(JetEta) >= 0.0 && fabs(JetEta) < 0.6) return 0.00296;
            if(fabs(JetEta) >= 0.6 && fabs(JetEta) < 1.2) return 0.00350;
            if(fabs(JetEta) >= 1.2 && fabs(JetEta) < 1.8) return 0.00652;
            if(fabs(JetEta) >= 1.8 && fabs(JetEta) < 2.4) return 0.00755;
          }
        }
      }
      if(tagger==JetTagging::CSVv2){
        if(wp==JetTagging::Loose){
        }
        if(wp==JetTagging::Medium){
        }
        if(wp==JetTagging::Tight){
        }
      }
    }
  }
  cout << "[MCCorrection::GetMCJetTagEff] No eff found for tagger = " << tagger << ", wp = " << wp << endl;
  cout << "[MCCorrection::GetMCJetTagEff] Or, wrong pt and eta range : pt = " << JetPt << ", eta = " << JetEta << endl;
  return 1.;
}
