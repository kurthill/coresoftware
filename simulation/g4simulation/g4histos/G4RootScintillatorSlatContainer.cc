#include "G4RootScintillatorSlatContainer.h"

#include "G4RootScintillatorSlat.h"

#include <g4detectors/PHG4ScintillatorSlat.h>
#include <TClonesArray.h>

using namespace std;

ClassImp(G4RootScintillatorSlatContainer)

static const int NMAX = 1000;

G4RootScintillatorSlatContainer::G4RootScintillatorSlatContainer():
  etotal(NAN),
  eion(NAN),
  leakage(NAN),
  event(0)
{
 SnglSlats = new TClonesArray("G4RootScintillatorSlat",NMAX);
}

G4RootScintillatorSlatContainer::~G4RootScintillatorSlatContainer()
{
  SnglSlats->Clear();
  delete SnglSlats;
}

void
G4RootScintillatorSlatContainer::Reset()
{
  etotal = NAN;
  leakage = NAN;
  event = 0;
  SnglSlats->Clear();
  if (SnglSlats->GetSize() > NMAX)
    {
      SnglSlats->Expand(NMAX);
    }
  return;
}



G4RootScintillatorSlat *
G4RootScintillatorSlatContainer::AddSlat(const PHG4ScintillatorSlat &slat)
{
  TClonesArray &cl = *SnglSlats;
  int nextindex = SnglSlats->GetLast() + 1;
  if (nextindex == SnglSlats->GetSize())
    {
      SnglSlats->Expand(SnglSlats->GetSize() + 10000);
    }
  new(cl[nextindex]) G4RootScintillatorSlat(slat);
  return (static_cast<G4RootScintillatorSlat *> (cl[nextindex]));
}

void
G4RootScintillatorSlatContainer::identify(ostream& os ) const
{
  os << "Number of Hits: " << SnglSlats->GetLast() << endl;
  return;
}
