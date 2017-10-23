#include "PHG4EnvelopeSteppingAction.h"
#include "PHG4EnvelopeDetector.h"

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4Shower.h>

#include <g4main/PHG4TrackUserInfoV1.h>

#include <phool/getClass.h>

#include <Geant4/G4Step.hh>

#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
// this is an ugly hack, the gcc optimizer has a bug which
// triggers the uninitialized variable warning which
// stops compilation because of our -Werror
#include <boost/version.hpp> // to get BOOST_VERSION
#if (__GNUC__ == 4 && __GNUC_MINOR__ == 4 && BOOST_VERSION == 105700 )
#pragma GCC diagnostic ignored "-Wuninitialized"
#pragma message "ignoring bogus gcc warning in boost header lexical_cast.hpp"
#include <boost/lexical_cast.hpp>
#pragma GCC diagnostic warning "-Wuninitialized"
#else
#include <boost/lexical_cast.hpp>
#endif

#include <iostream>

using namespace std;

//______________________________________________________________
PHG4EnvelopeSteppingAction::PHG4EnvelopeSteppingAction( PHG4EnvelopeDetector* detector ):
  detector_( detector ),
  hits_(NULL),
  hit(NULL)
{
}

//____________________________________________________________________________..
bool PHG4EnvelopeSteppingAction::UserSteppingAction( const G4Step* aStep, bool )
{
  G4TouchableHandle touch = aStep->GetPreStepPoint()->GetTouchableHandle();
  G4VPhysicalVolume* volume = touch->GetVolume();
	
  int whichactive = detector_->IsInEnvelope(volume);

  if ( !whichactive  )
    {
      return false;
    }
	
  int layer_id = detector_->get_Layer();
  int tower_id = -1;
  tower_id = touch->GetCopyNumber();

  /* Get energy deposited by this step */
  G4double edep = aStep->GetTotalEnergyDeposit() / GeV;
  G4double eion = (aStep->GetTotalEnergyDeposit() - aStep->GetNonIonizingEnergyDeposit()) / GeV;

  /* Get pointer to associated Geant4 track */
  const G4Track* aTrack = aStep->GetTrack();
	
  //This detector is a black hole! Just put all kinetic energy into edep
  edep = aTrack->GetKineticEnergy() / GeV;
  G4Track* killtrack = const_cast<G4Track *> (aTrack);
  killtrack->SetTrackStatus(fStopAndKill);
	
  /* Make sure we are in a volume */
  if ( detector_->IsActive() )
    {
      /* Check if particle is 'geantino' */
      bool geantino = false;
      if (aTrack->GetParticleDefinition()->GetPDGEncoding() == 0 && aTrack->GetParticleDefinition()->GetParticleName().find("geantino") != string::npos)
	{
	  geantino = true;
	}

      /* Get Geant4 pre- and post-step points */
      G4StepPoint * prePoint = aStep->GetPreStepPoint();
      G4StepPoint * postPoint = aStep->GetPostStepPoint();

      switch (prePoint->GetStepStatus())
	{
	case fGeomBoundary:
	case fUndefined:
	  hit = new PHG4Hitv1();
	  //	  hit->set_layer(0);
	  hit->set_scint_id(tower_id);

	  /* Set hit location (tower index) */
	  //hit->set_index_j(idx_j);
	  //hit->set_index_k(idx_k);
	  //hit->set_index_l(idx_l);
	  hit->set_index_j(0);
	  hit->set_index_k(0);
	  hit->set_index_l(0);

	  /* Set hit location (space point) */
	  hit->set_x( 0, prePoint->GetPosition().x() / cm);
	  hit->set_y( 0, prePoint->GetPosition().y() / cm );
	  hit->set_z( 0, prePoint->GetPosition().z() / cm );

	  /* Set momentum */
	  hit->set_x( 0, prePoint->GetMomentum().x() / GeV );
	  hit->set_y( 0, prePoint->GetMomentum().y() / GeV );
	  hit->set_z( 0, prePoint->GetMomentum().z() / GeV );

	  /* Set hit time */
	  hit->set_t( 0, prePoint->GetGlobalTime() / nanosecond );

	  /* set the track ID */
	  {
	    hit->set_trkid(aTrack->GetTrackID());
	    if ( G4VUserTrackInformation* p = aTrack->GetUserInformation() )
	      {
		if ( PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p) )
		  {
		    hit->set_trkid(pp->GetUserTrackId());
		    hit->set_shower_id(pp->GetShower()->get_id());
		  }
	      }
	  }

	  /* set intial energy deposit */
	  hit->set_edep( 0 );
	  hit->set_eion( 0 );

	  hits_->AddHit(layer_id, hit);
	  
	  {
	    if ( G4VUserTrackInformation* p = aTrack->GetUserInformation() )
	      {
		if ( PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p) )
		  {
		    pp->GetShower()->add_g4hit_id(hits_->GetID(),hit->get_hit_id());
		  }
	      }
	  }
	  
	  break;
	default:
	  break;
	}

      /* Update exit values- will be overwritten with every step until
       * we leave the volume or the particle ceases to exist */
      hit->set_x( 1, postPoint->GetPosition().x() / cm );
      hit->set_y( 1, postPoint->GetPosition().y() / cm );
      hit->set_z( 1, postPoint->GetPosition().z() / cm );

      hit->set_t( 1, postPoint->GetGlobalTime() / nanosecond );

      /* sum up the energy to get total deposited */
      hit->set_edep(hit->get_edep() + edep);
      hit->set_eion(hit->get_eion() + eion);

      if (geantino)
	{
	  hit->set_edep(-1); // only energy=0 g4hits get dropped, this way geantinos survive the g4hit compression
	  hit->set_eion(-1);
	}
      if (edep > 0)
	{
	  if ( G4VUserTrackInformation* p = aTrack->GetUserInformation() )
	    {
	      if ( PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p) )
		{
		  pp->SetKeep(1); // we want to keep the track
		}
	    }
	}
      return true;
    }
  else
    {
      return false;
    }
}
	
void PHG4EnvelopeSteppingAction::SetInterfacePointers( PHCompositeNode* topNode )
{
  string hitnodename;

  if (detector_->SuperDetector() != "NONE")
    {
      hitnodename = "G4HIT_ENVELOPE_" + detector_->SuperDetector();
    }
  else
    {
      hitnodename = "G4HIT_ENVELOPE_" + detector_->GetName();
    }
	
  hits_ =  findNode::getClass<PHG4HitContainer>( topNode , hitnodename.c_str() );

  // if we do not find the node it's messed up.
  if ( ! hits_ )
    {
      std::cout << "PHG4CrystalCalorimeterSteppingAction::SetTopNode - unable to find " << hitnodename << std::endl;
    }
	
	
}
