#ifndef PHG4CylinderGeomContainer_H__
#define PHG4CylinderGeomContainer_H__

#include <phool/PHObject.h>

#include <map>

class PHG4CylinderGeom;

class PHG4CylinderGeomContainer: public PHObject
{
 public:
  typedef std::map<int,PHG4CylinderGeom *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;

  PHG4CylinderGeomContainer();
  virtual ~PHG4CylinderGeomContainer();

  void identify(std::ostream& os = std::cout) const;

  int AddLayerGeom(const int i, PHG4CylinderGeom *mygeom);
  int AddLayerGeom(PHG4CylinderGeom *mygeom);
  PHG4CylinderGeom *GetLayerGeom(const int i);
  PHG4CylinderGeom *GetFirstLayerGeom();
  int get_NLayers() const {return layergeoms.size();}
  std::pair<std::map<int,PHG4CylinderGeom *>::const_iterator, std::map<int,PHG4CylinderGeom *>::const_iterator> get_begin_end() const {return std::make_pair(layergeoms.begin(), layergeoms.end());}

 protected:
  std::map<int,PHG4CylinderGeom *> layergeoms ;
  float magfield;
  ClassDef(PHG4CylinderGeomContainer,1)
};

#endif
