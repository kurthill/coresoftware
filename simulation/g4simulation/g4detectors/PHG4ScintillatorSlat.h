#ifndef PHG4SCINTILLATORSLAT_H
#define PHG4SCINTILLATORSLAT_H

#include "PHG4ScintillatorSlatDefs.h"

#include <g4main/PHG4HitDefs.h>
#include <phool/PHObject.h>
#include <cmath>
#include <map>

class PHG4ScintillatorSlat : public PHObject
{
 public:
  
  virtual ~PHG4ScintillatorSlat(){}

  virtual void identify(std::ostream& os = std::cout) const {
    os << "PHG4ScintillatorSlat base class" << std::endl;
  }
  
  virtual void add_edep(const double edep, const double e, const double light_yield) {return;}

  virtual void set_key(const PHG4ScintillatorSlatDefs::keytype) {return;}
  virtual void add_hit_key(PHG4HitDefs::keytype) {return;}

  virtual short get_column() const {return 0xFFFF;}
  virtual short get_row() const {return 0xFFFF;}
  virtual PHG4ScintillatorSlatDefs::keytype get_key() const {return 0xFFFFFFFF;}

  virtual double get_edep() const {return NAN;}
  virtual double get_eion() const {return NAN;}
  virtual double get_light_yield() const  {return NAN;}


  
 protected:

  PHG4ScintillatorSlat() {}
  ClassDef(PHG4ScintillatorSlat,1)
};

#endif
