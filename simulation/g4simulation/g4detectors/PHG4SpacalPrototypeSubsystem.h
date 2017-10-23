#ifndef PHG4SpacalPrototypeSubsystem_h
#define PHG4SpacalPrototypeSubsystem_h

#include "PHG4DetectorSubsystem.h"

class PHG4SpacalPrototypeDetector;
class PHG4SteppingAction;

class PHG4SpacalPrototypeSubsystem : public PHG4DetectorSubsystem
{

public:

  //! constructor
  explicit PHG4SpacalPrototypeSubsystem(const std::string &name = "CEMC");

  //! destructor
  virtual ~PHG4SpacalPrototypeSubsystem(void) {}

  //! init
  /*!
   creates the detector_ object and place it on the node tree, under "DETECTORS" node (or whatever)
   reates the stepping action and place it on the node tree, under "ACTIONS" node
   creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
   */
  int InitRunSubsystem(PHCompositeNode *);

  //! event processing
  /*!
   get all relevant nodes from top nodes (namely hit list)
   and pass that to the stepping action
   */
  int process_event(PHCompositeNode *);

  //! accessors (reimplemented)
  virtual PHG4Detector* GetDetector(void) const;
  PHG4SteppingAction* GetSteppingAction(void) const {return steppingAction_;}


  void
  Print(const std::string &what = "ALL") const;

private:
  //! load the default parameter to param
  void SetDefaultParameters();


  //! detector geometry
  /*! derives from PHG4Detector */
  PHG4SpacalPrototypeDetector* detector_;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SteppingAction* steppingAction_;
};

#endif
