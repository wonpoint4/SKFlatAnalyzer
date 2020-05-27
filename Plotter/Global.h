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


/////////////////////////////////////////////////////////////////////////////
///////////////////////////// Message ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
enum VERBOSITY{QUIET,ERROR,WARNING,INFO,DEBUG,ALL};
int Verbosity=VERBOSITY::INFO; //0:quiet 1:error 2:warning 3:info 4:debug 5:all
int _depth=0;
void PMessage(VERBOSITY level,const char* msg){
  if(level>Verbosity) return;
  if(level==VERBOSITY::QUIET) return;
  if(_depth<0) _depth=0;
  TString indent=std::string(_depth,' ');
  if(level==VERBOSITY::ERROR)        std::cout<<indent+"## ERROR ## "+msg<<endl;
  else if(level==VERBOSITY::WARNING) std::cout<<indent+"# WARNING # "+msg<<endl;
  else if(level==VERBOSITY::INFO)    std::cout<<indent+"## INFO ### "+msg<<endl;
  else if(level==VERBOSITY::DEBUG)   std::cout<<indent+"## DEBUG ## "+msg<<endl;
  else if(level==VERBOSITY::ALL)     std::cout<<indent+"### ALL ### "+msg<<endl;
}
void PError(const char* msg){
  PMessage(VERBOSITY::ERROR,msg);
}
void PWarning(const char* msg){
  PMessage(VERBOSITY::WARNING,msg);
}
void PInfo(const char* msg){
  PMessage(VERBOSITY::INFO,msg);
}
void PDebug(const char* msg){
  PMessage(VERBOSITY::DEBUG,msg);
}
void PAll(const char* msg){
  PMessage(VERBOSITY::ALL,msg);
}
  
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
