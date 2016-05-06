#include "PHG4ScintillatorSlatv1.h"

using namespace std;

PHG4ScintillatorSlatv1::PHG4ScintillatorSlatv1():
  key(0xFFFFFFFF),
  edep(0),
  eion(0),
  light_yield(0)
{}

short
PHG4ScintillatorSlatv1::get_row() const
{
  return (key&0xFFFF);
}

short
PHG4ScintillatorSlatv1::get_column() const
{
  return (key>>16);
}

void
PHG4ScintillatorSlatv1::identify(std::ostream& os) const
{
  os << "row " << get_row() << " ";
  os << " column " << get_column() << " ";
  os << " energy deposition " << get_edep();
  os << " ionization energy " << get_eion();
  os << " light yield " << get_light_yield();
  os << endl;
}
