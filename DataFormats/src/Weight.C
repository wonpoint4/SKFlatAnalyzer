#include "Weight.h"
Weight::Weight(){}
Weight::~Weight(){}
Weight::Weight(const double& w){
  operator=(w);
}
Weight::operator double() const{
  if(nzero>0) return 0.0;
  else if(nzero<0){
    if(value>=0)
      return std::numeric_limits<double>::infinity();
    else
      return -std::numeric_limits<double>::infinity();
  }
  return value;
}
Weight& Weight::operator=(const double& w){
  if(w==0.0){
    nzero=1;
    value=1.0;
  }else if(w==std::numeric_limits<double>::infinity()){
    nzero=-1;
    value=1.0;
  }else if(w==-std::numeric_limits<double>::infinity()){
    nzero=-1;
    value=-1.0;
  }else{
    nzero=0;
    value=w;
  }
  return *this;
}

Weight& Weight::operator*=(const double& w){
  if(w==0.0){
    nzero+=1;
  }else if(w==std::numeric_limits<double>::infinity()){
    nzero-=1;
  }else if(w==-std::numeric_limits<double>::infinity()){
    nzero-=1;
    value*=-1.0;
  }else{
    value*=w;
  }
  return *this;
}
Weight& Weight::operator*=(const Weight& w){
  value*=w.value;
  nzero+=w.nzero;
  return *this;
}
Weight& Weight::operator/=(const double& w){
  if(w==0.0){
    nzero-=1;
  }else if(w==std::numeric_limits<double>::infinity()){
    nzero+=1;
  }else if(w==-std::numeric_limits<double>::infinity()){
    nzero+=1;
    value*=-1.0;
  }else{
    value/=w;
  }
  return *this;
}
Weight& Weight::operator/=(const Weight& w){
  value/=w.value;
  nzero-=w.nzero;
  return *this;
}
