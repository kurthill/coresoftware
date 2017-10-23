#include "PHG4FPbScSteppingAction.h"
#include "PHG4FPbScDetector.h"
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4TrackUserInfoV1.h>

#include <phool/getClass.h>

#include <Geant4/G4Step.hh>

using namespace std;


PHG4FPbScSteppingAction::PHG4FPbScSteppingAction( PHG4FPbScDetector* detector) :
    detector_( detector ), hits_(NULL), hit(NULL)
{
}


bool PHG4FPbScSteppingAction::UserSteppingAction( const G4Step* aStep, bool)
{
  G4VPhysicalVolume* volume = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  
  int layer_id = 0;
  
  if((detector_->isInScintillator(volume)) == false){return false;}
  
  if((aStep->GetTotalEnergyDeposit()/GeV) == 0.){return false;}
  
  layer_id = detector_->getScintillatorLayer(volume);

  G4StepPoint* prePoint = aStep->GetPreStepPoint();
  G4StepPoint* postPoint = aStep->GetPostStepPoint();
  G4Track* aTrack = aStep->GetTrack();
  
  G4double x_center, y_center, z_center;

  int hit_index = detector_->computeIndex(layer_id, 
				      prePoint->GetPosition().x() / cm, 
				      prePoint->GetPosition().y() / cm, 
				      prePoint->GetPosition().z() / cm, 
				      x_center, y_center, z_center);

  PHG4HitContainer::Iterator it = hits_->findOrAddHit(hit_index);
  
  PHG4Hit* mhit = it->second;
  mhit->set_x(0, x_center);
  mhit->set_y(0, y_center);
  mhit->set_z(0, z_center);
  mhit->set_edep(0);
  
  hit_index = detector_->computeIndex(layer_id, 
				      postPoint->GetPosition().x() / cm, 
				      postPoint->GetPosition().y() / cm, 
				      postPoint->GetPosition().z() / cm, 
				      x_center, y_center, z_center);
  
  it = hits_->findOrAddHit(hit_index);

  mhit = it->second;
  mhit->set_x(1, x_center);
  mhit->set_y(1, y_center);
  mhit->set_z(1, z_center);
  mhit->set_edep((aStep->GetTotalEnergyDeposit()/GeV) + it->second->get_edep());
  
  // Also store the flag that we want to keep this track on output
  //
  if ( G4VUserTrackInformation* p = aTrack->GetUserInformation() )
  {
    // Do nothing; was assume this flag has already been set
    if ( PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p) )
    {
      pp->SetWanted(true);
    }
    else
    {
      std::cout << "WARNING: Unknown UserTrackInformation stored in track" << std::endl;
    }
  }
  else
  {
    PHG4TrackUserInfoV1* pv = new PHG4TrackUserInfoV1();
    pv->SetWanted(true);
    aTrack->SetUserInformation(pv);
  }

  //set the track ID
  {
    mhit->set_trkid(aTrack->GetTrackID());
    if ( G4VUserTrackInformation* p = aTrack->GetUserInformation() )
      {
	if ( PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p) )
	  {
	    mhit->set_trkid(pp->GetUserTrackId());
	    mhit->set_shower_id(pp->GetShower()->get_id());
	  }
      }
  }
  
  {
    if ( G4VUserTrackInformation* p = aTrack->GetUserInformation() )
      {
	if ( PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p) )
	  {
	    pp->GetShower()->add_g4hit_id(hits_->GetID(),mhit->get_hit_id());
	  }
      }
  }
  
  return true;
}


void PHG4FPbScSteppingAction::SetInterfacePointers( PHCompositeNode* topNode )
{
  string hitnodename = "G4HIT_" + detector_->GetName();
  //now look for the map and grab a pointer to it.
  hits_ =  findNode::getClass<PHG4HitContainer>( topNode , hitnodename.c_str() );
  
  // if we do not find the node we need to make it.
  if ( ! hits_ )
    { std::cout << "PHG4FPbScSteppingAction::SetTopNode - unable to find " << hitnodename.c_str() << std::endl; }
  
}



