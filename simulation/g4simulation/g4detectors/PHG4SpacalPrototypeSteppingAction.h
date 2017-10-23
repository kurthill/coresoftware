// $$Id: PHG4SpacalPrototypeSteppingAction.h,v 1.1 2014/03/24 01:36:44 jinhuang Exp $$

/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $$Revision: 1.1 $$
 * \date $$Date: 2014/03/24 01:36:44 $$
 */

#ifndef PHG4SpacalPrototypeSteppingAction_h
#define PHG4SpacalPrototypeSteppingAction_h

#include <g4main/PHG4SteppingAction.h>
#include <string>

class PHG4Hit;
class PHG4HitContainer;
class PHG4Shower;
class PHG4SpacalPrototypeDetector;

class PHG4SpacalPrototypeSteppingAction : public PHG4SteppingAction
{

public:

  //! constructor
  explicit PHG4SpacalPrototypeSteppingAction(PHG4SpacalPrototypeDetector*);

  //! destroctor
  virtual
  ~PHG4SpacalPrototypeSteppingAction()
  {
  }

  //! stepping action
  virtual bool
  UserSteppingAction(const G4Step*, bool);

  //! reimplemented from base class
  virtual void
  SetInterfacePointers(PHCompositeNode*);

  double
  get_zmin();

  double
  get_zmax();
private:

  //! pointer to the detector
  PHG4SpacalPrototypeDetector* detector_;

//! pointer to hit container
  PHG4HitContainer * hits_;
  PHG4HitContainer * absorberhits_;
  PHG4Hit *hit;
  PHG4HitContainer *savehitcontainer;
  PHG4Shower *saveshower;
};

#endif // PHG4SpacalPrototypeSteppingAction_h
