#ifndef __STYLE_CC__
#define __STYLE_CC__
class Style{
public:
  int fillcolor=-1;
  int fillstyle=-1;
  int linecolor=-1;
  int linestyle=-1;
  int linewidth=-1;
  int markercolor=-1;
  int markerstyle=-1;
  double markersize=-1;
  TString drawoption;
  Style(int color=-1,int marker=-1,int fillstyle_=-1,TString drawoption_="");
  void Apply(TH1* hist) const;
  void Print() const;
};
Style::Style(int color,int marker,int fillstyle_,TString drawoption_){
  if(color>0){
    fillcolor=0;
    linecolor=color;
    markercolor=color;
    drawoption="e hist";
  }
  if(marker>0){
    markerstyle=marker;
    markersize=0.7;
    drawoption="e p";
  }
  if(fillstyle_>0){
    fillstyle=fillstyle_;
    fillcolor=linecolor;
    drawoption.ReplaceAll("e ","e2 ");
  }
  if(drawoption_!="") drawoption=drawoption_;
}
    
void Style::Apply(TH1* hist) const {
  PAll("[Style::Apply(TH1* hist)]");
  if(hist){
    if(strstr(hist->ClassName(),"THStack")==NULL){
      if(linecolor>=0) hist->SetLineColor(linecolor);
      if(linestyle>=0) hist->SetLineStyle(linestyle);
      if(linewidth>=0) hist->SetLineWidth(linewidth);
      if(fillcolor>=0) hist->SetFillColor(fillcolor);
      if(fillstyle>=0) hist->SetFillStyle(fillstyle);
      if(markercolor>=0) hist->SetMarkerColor(markercolor);
      if(markerstyle>=0) hist->SetMarkerStyle(markerstyle);
      if(markersize>=0) hist->SetMarkerSize(markersize);
      hist->SetOption(drawoption);
    }
  }
}
void Style::Print() const{
  if(fillcolor!=-1) cout<<"FC:"<<fillcolor<<" ";
  if(fillstyle!=-1) cout<<"FS:"<<fillstyle<<" ";
  if(linecolor!=-1) cout<<"LC:"<<linecolor<<" ";
  if(linestyle!=-1) cout<<"LS:"<<linestyle<<" ";
  if(linewidth!=-1) cout<<"LW:"<<linewidth<<" ";
  if(markercolor!=-1) cout<<"MC:"<<markercolor<<" ";
  if(markerstyle!=-1) cout<<"MS:"<<markerstyle<<" ";
  if(markersize!=-1) cout<<"MW:"<<markersize<<" ";
  if(drawoption!="") cout<<"DO:'"<<drawoption<<"' ";
  cout<<endl;
}

#endif
