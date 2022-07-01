#ifndef __Stub__h
#define __Stub__h
#include<stdint.h>
#include<iostream>
#include<iomanip>
class stub : public TObject {
  public:
  stub():center__(0),column__(0),bend_(0.){}
  stub(const float center_, const uint8_t column_, const float bend):center__(center_), column__(column_), bend_(bend) {}
    virtual ~stub() {}
    float center() const { return center__; }//Position in half strip precision
    unsigned int column() const  { return column__; }//column 
    float bend() const         { return bend_; }//
    float realPos() const {return (center__-508.5)*0.09;}

    friend std::ostream &operator<<(std::ostream& os, const stub& s) {
      os << "center_= "  << std::setw(4) << s.center()
         << " column_= " << std::setw(4) << (int)s.column()
         << " Bend= " << std::setw(4) << s.bend()
         << std::endl;
      return os;
    }
  private:
    float center__;
    uint8_t column__;
    float bend_;
};
#endif
