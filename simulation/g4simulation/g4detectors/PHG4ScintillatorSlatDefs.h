#ifndef PHG4SCINTILLATORSLATGEFS_H
#define PHG4SCINTILLATORSLATGEFS_H

#include <map>

namespace PHG4ScintillatorSlatDefs
{
  typedef unsigned int keytype;
  static int columnbits = 16; // upper bits used for searching, we want columns
                              // which allows us to create towers quickly by
                              // bunching them up in groups of 5
  static int rowbits = 32-columnbits;
  keytype genkey(const short irow, const short icolumn);
  std::pair<short,short> getrowcol(const keytype key);
}


#endif
