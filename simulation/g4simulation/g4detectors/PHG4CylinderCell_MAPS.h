#ifndef PHG4CYLINDERCELL_MAPS_H
#define PHG4CYLINDERCELL_MAPS_H

#include "PHG4CylinderCellv2.h"
#include <cmath>
#include <map>
#include <iostream>

class PHG4CylinderCell_MAPS : public PHG4CylinderCellv2
{
 public:

  PHG4CylinderCell_MAPS();
  virtual ~PHG4CylinderCell_MAPS(){}

  void identify(std::ostream& os = std::cout) const;

  void set_stave_index(const  int si) {stave_number = si;}
  int get_stave_index() const  {return stave_number;}

  void set_half_stave_index(const int i) {half_stave_number = i;}
  int get_half_stave_index() const {return half_stave_number;}

  void set_module_index(const int i) {module_number = i;}
  int get_module_index() const {return module_number;}

  void set_chip_index(const int i) {chip_number = i;}
  int get_chip_index() const {return chip_number;}

  void set_pixel_index(const int i) {pixel_number = i;}
  int get_pixel_index() const {return pixel_number;}

  
 protected:

  int stave_number;
  int half_stave_number;
  int module_number;
  int chip_number;
  int pixel_number;

  
  ClassDef(PHG4CylinderCell_MAPS,1)
};

#endif
