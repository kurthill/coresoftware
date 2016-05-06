#ifndef PHG4SCINTILLATORSLATCONTAINER_H__
#define PHG4SCINTILLATORSLATCONTAINER_H__

#include "PHG4ScintillatorSlat.h"

#include <phool/PHObject.h>

#include <map>
#include <set>

class PHG4ScintillatorSlatContainer: public PHObject
{

  public:
  typedef std::map<PHG4ScintillatorSlatDefs::keytype,PHG4ScintillatorSlat *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;
  typedef std::set<int>::const_iterator LayerIter;
  typedef std::pair<LayerIter, LayerIter> LayerRange;

  PHG4ScintillatorSlatContainer();

  virtual ~PHG4ScintillatorSlatContainer() {}

  void Reset();

  void identify(std::ostream& os = std::cout) const;

  ConstIterator AddScintillatorSlat(const PHG4ScintillatorSlatDefs::keytype key, PHG4ScintillatorSlat *newscintillatorSlat);
  
  //! preferred removal method, key is currently the slat id
  void RemoveScintillatorSlat(PHG4ScintillatorSlatDefs::keytype key) {
    slatmap.erase(key);
  }

  //! inefficent, use key where possible instead
  void RemoveScintillatorSlat(PHG4ScintillatorSlat *slat)
  {
    Iterator its = slatmap.begin();
    Iterator last = slatmap.end();
    for (; its != last;)
      {
	if (its->second == slat)
	  {
	    slatmap.erase(its++);
	  }
	else
	  {
	    ++its;
	  }
      }
  }


  //! return all scintillatorSlats matching a given detid
  ConstRange getScintillatorSlats(const short icolumn) const;

  //! return all hist
  ConstRange getScintillatorSlats( void ) const;

  PHG4ScintillatorSlat* findScintillatorSlat(PHG4ScintillatorSlatDefs::keytype key);

  unsigned int size( void ) const
  { return slatmap.size(); }

  double getTotalEdep() const;

 protected:
  Map slatmap;

  ClassDef(PHG4ScintillatorSlatContainer,1)
};

#endif
