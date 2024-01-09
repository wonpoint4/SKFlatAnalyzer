#ifndef Weight_h
#define Weight_h

#include <limits>
class Weight{
public:
  double value=1.0;
  int nzero=0;

  Weight();
  ~Weight();
  Weight(const double& w);
  operator double() const;

  Weight& operator=(const double& w);
  
  Weight& operator*=(const double& w);
  Weight& operator*=(const Weight& w);
  Weight& operator/=(const double& w);
  Weight& operator/=(const Weight& w);

  Weight& operator*=(const float& w){ return operator*=((const double&)w);}
  Weight& operator*=(const int& w){ return operator*=((const double&)w);}
  Weight& operator/=(const float& w){ return operator/=((const double&)w);}
  Weight& operator/=(const int& w){ return operator/=((const double&)w);}

  const Weight operator*(const double& w) const {return Weight(*this) *= w;}
  const Weight operator*(const float& w) const {return Weight(*this) *= w;}
  const Weight operator*(const int& w) const {return Weight(*this) *= w;}
  const Weight operator*(const Weight& w) const {return Weight(*this) *= w;}
  const Weight operator/(const double& w) const {return Weight(*this) /= w;}
  const Weight operator/(const float& w) const {return Weight(*this) /= w;}
  const Weight operator/(const int& w) const {return Weight(*this) /= w;}
  const Weight operator/(const Weight& w) const {return Weight(*this) /= w;}
};
#endif
