#ifndef __GLOBAL_H__
#define __GLOBAL_H__

#include<iostream>
#include<vector>
#include"TH2D.h"
#include"TFile.h"
#include"TKey.h"
#include"TString.h"
#include"THStack.h"
#include"TLegend.h"
#include"TCanvas.h"
#include"TLine.h"
#include"TRegexp.h"
#include"TPaveText.h"
#include"TLegendEntry.h"
using namespace std;

int DEBUG=3; //0:quiet 1:error 2:warning 3:info 4:debug 5:all

/////////////////////////////////////////////////////////////////////////////
///////////////////////////// struct and enum ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////
TString GetStringEColor(EColor color){
  switch(color){
  case kBlack: return "kBlack";
  case kRed: return "kRed";
  case kGreen: return "kGreen";
  case kBlue: return "kBlue";
  case kYellow: return "kYellow";
  case kMagenta: return "kMagenta";
  case kOrange: return "kOrange";
  case kGray: return "kGray";
  default: return "UNKNOWN";
  }
}
enum Channel{MUON,ELECTRON};
TString GetStringChannel(Channel channel){
  switch(channel){
  case MUON: return "muon";
  case ELECTRON: return "electron";
  default: return "###ERROR### Bad Channel";
  }
}

#endif
