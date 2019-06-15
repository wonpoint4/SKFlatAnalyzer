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

int DEBUG=4; //0:quiet 1:error 2:warning 3:info 4:debug

/////////////////////////////////////////////////////////////////////////////                                                                               ///////////////////////////// struct and enum ///////////////////////////////
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
  default: return "UNKNOWN";
  }
}
enum SystematicType{ENVELOPE,GAUSSIAN,HESSIAN,MULTI};
TString GetStringSystematicType(SystematicType type){
  switch(type){
  case ENVELOPE: return "ENVELOPE";
  case GAUSSIAN: return "GAUSSIAN";
  case HESSIAN: return "HESSIAN";
  case MULTI: return "MULTI";
  default: return "###ERROR### Bad SystematicType";
  }
}
struct Systematic{
  TString name;
  SystematicType type;
  vector<TString> suffixes;
  int sysbit;
  int varibit;
};
enum Channel{MUON,ELECTRON};
TString GetStringChannel(Channel channel){
  switch(channel){
  case MUON: return "muon";
  case ELECTRON: return "electron";
  default: return "###ERROR### Bad Channel";
  }
}


#endif
