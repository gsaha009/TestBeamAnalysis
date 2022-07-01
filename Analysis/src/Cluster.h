#ifndef __Cluster__h
#define __Cluster__h
#include<stdint.h>
#include<iostream>
#include<iomanip>
class cluster {
 public:
  cluster() : firstStrip_(9999), column_(999), size_(999){}
    
  cluster(unsigned int fStrip, unsigned int col, uint16_t s) :
    firstStrip_(fStrip), column_(col), size_(s) {}
    virtual ~cluster(){}
    float center() const { return float(firstStrip()) + 0.5*(float(size()) - 1.); }
    std::pair< float, float > barycenter() const { return std::make_pair(column(), center()); }
    unsigned int firstStrip() const { return firstStrip_; }
    unsigned int column() const     { return column_; }
    uint16_t size() const           { return size_; }
    
    friend std::ostream &operator<<(std::ostream& os, const cluster& cls) {    
      os << "FStrip="  << std::setw(4) << cls.firstStrip() 
	 << " Column=" << std::setw(4) << cls.column() 
	 << " Size="   << std::setw(6) <<  cls.size() 
	 << " Center="    << std::setw(4) << cls.center() 
	 << " Barycenter(" << std::setw(4) << cls.barycenter().first << "," << cls.barycenter().second << ")" 
	 << std::endl; 
      return os;
    }
    
  private:
    unsigned int firstStrip_;
    unsigned int column_;
    uint16_t size_;

};
#endif
