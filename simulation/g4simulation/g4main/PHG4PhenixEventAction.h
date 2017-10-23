#ifndef PHG4VUserEventAction_h
#define PHG4VUserEventAction_h

#include <phool/PHTimeServer.h>

#include <Geant4/G4UserEventAction.hh>

#include <list>


class G4Event;
class PHG4EventAction;
class PHCompositeNode;

class PHG4PhenixEventAction : public G4UserEventAction
{

  public:
  PHG4PhenixEventAction( void );

  virtual ~PHG4PhenixEventAction();

  //! register an action. This is called in PHG4Reco::Init based on which actions are found on the tree
  void AddAction( PHG4EventAction* action )
  { actions_.push_back( action ); }

  void BeginOfEventAction(const G4Event*);

  void EndOfEventAction(const G4Event*);

  private:

  //! list of subsystem specific Event actions
  typedef std::list<PHG4EventAction*> ActionList;
  ActionList actions_;

  //! module timer.
  PHTimeServer::timer _timer;
};


#endif
