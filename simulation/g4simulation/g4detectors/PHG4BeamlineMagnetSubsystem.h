#ifndef PHG4BeamlineMagnetSubsystem_h
#define PHG4BeamlineMagnetSubsystem_h

#include "PHG4DetectorSubsystem.h"

class PHG4BeamlineMagnetDetector;

class PHG4BeamlineMagnetSubsystem: public PHG4DetectorSubsystem
{

  public:

  //! constructor
  PHG4BeamlineMagnetSubsystem( const std::string &name = "CYLINDER", const int layer = 0 );

  //! destructor
  virtual ~PHG4BeamlineMagnetSubsystem( void )
  {}

  //! init runwise stuff
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

  //! Print info (from SubsysReco)
  void Print(const std::string &what = "ALL") const;

  //! accessors (reimplemented)
  PHG4Detector* GetDetector( void ) const;

 private:
  void SetDefaultParameters();

  //! detector geometry
  /*! defives from PHG4Detector */
  PHG4BeamlineMagnetDetector* detector_;

};

#endif
