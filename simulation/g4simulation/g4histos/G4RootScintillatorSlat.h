#ifndef G4RootScintillatorSlat_H_
#define G4RootScintillatorSlat_H_

#include "phool/PHObject.h"

class PHG4ScintillatorSlat;

class G4RootScintillatorSlat : public PHObject {

 public:
  G4RootScintillatorSlat();
  G4RootScintillatorSlat(const PHG4ScintillatorSlat &slat);
  virtual ~G4RootScintillatorSlat() {}

  void Reset();
  int isValid() const;
  void identify(std::ostream& os=std::cout) const;

  int get_row() const {return row;}
  int get_column() const {return column;}

  double get_edep() const {return edep;}
  double get_eion() const {return eion;}
  double get_light_yield() const {return light_yield;}

 protected:
  short row;
  short column;
  double edep;
  double eion;
  double light_yield;

  ClassDef(G4RootScintillatorSlat,1)
};
 
#endif /* G4RootScintillatorSlat_H_ */
