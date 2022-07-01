#ifndef __Track__h
#define __Track__h
#include<stdint.h>
#include<iostream>
#include<iomanip>
class Track : public TObject {
  public:
  Track(): xPos_(0.), yPos_(0.), dxdz_(0.), dydz_(0.), chi2_(0.), ndof_(0.), iden_(0){}
  Track(const double xp, const double yp, const double dxdz, const double dydz, 
        const double chi2, const double ndof, const int identity) 
  : xPos_(xp), yPos_(yp), dxdz_(dxdz), dydz_(dydz), chi2_(chi2), ndof_(ndof), iden_(identity) 
  {}

  double x() { return xPos_; }
  double y() { return yPos_; }
  double dxdz() { return dxdz_; }
  double dydz() { return dydz_; }
  double chi2() { return chi2_; }
  double ndof() { return ndof_; }
  int iden() { return iden_; }
  virtual ~Track() {}

    friend std::ostream &operator<<(std::ostream& os, const Track& s) {
      os << "X= "  << std::setw(4) << s.xPos_
         << "\tY= "  << std::setw(4) << s.yPos_
         << std::endl;
      return os;
    }
  private:
    double  xPos_;
    double  yPos_;
    double  dxdz_;
    double  dydz_;
    double  chi2_;
    double  ndof_;
    int     iden_;
};
#endif
