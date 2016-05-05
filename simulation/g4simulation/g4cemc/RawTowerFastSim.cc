#include "RawTowerFastSim.h"
#include "RawTowerContainer.h"
#include "RawTowerGeomContainer.h"
#include "RawTowerv1.h"

#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellContainer.h>
#include <g4detectors/PHG4CylinderCell.h>
#include <g4detectors/PHG4CylinderCellDefs.h>

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHIODataNode.h>
#include <phool/PHRandomSeed.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>

#include <TRandom3.h>

#include <iostream>
#include <stdexcept>
#include <map>
#include <cassert>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;

RawTowerFastSim::RawTowerFastSim(const std::string& name) :
    SubsysReco(name),
    _hits(NULL), _cells(NULL), _calib_towers(NULL), _raw_towers(NULL), //
    rawtowergeom(NULL), //
		truthinfo(NULL), //
    detector("NONE"), //
    _calib_tower_node_prefix("CALIB"), _raw_tower_node_prefix("RAW"), //
    _em_sigma_p0(0), _em_sigma_p1(0),
    _h_sigma_p0(0), _h_sigma_p1(0),
    _zero_suppression_GeV(0), //
		_spread_e_thresh(0),
    _timer(PHTimeServer::get()->insert_new(name))
{
	RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
	seed = PHRandomSeed(); // fixed seed handled in PHRandomSeed()
	gsl_rng_set(RandomGenerator, seed);
}

int
RawTowerFastSim::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode;
  dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode",
      "DST"));
  if (!dstNode)
    {
      std::cout << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
          << "DST Node missing, doing nothing." << std::endl;
      exit(1);
    }

  try
    {
      CreateNodes(topNode);
    }
  catch (std::exception& e)
    {
      std::cout << e.what() << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

int
RawTowerFastSim::process_event(PHCompositeNode *topNode)
{
  if (verbosity)
    {
      std::cout << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__<< "Process event entered" << std::endl;
    }

	//cout << "*************NEW EVENT***************" << endl;

	// Loop over all raw towers and get energy, then set calib towers to new energy 
  RawTowerContainer::ConstRange begin_end = _raw_towers->getTowers();
  RawTowerContainer::ConstIterator rtiter;
  for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
  {
 		const RawTowerDefs::keytype tkey = rtiter->first;
  	RawTower* raw_tower = rtiter->second;
  	assert(raw_tower);

		// Loop over each cell in a tower
		RawTower::CellConstRange cell_begin_end = raw_tower->get_g4cells();
 		RawTower::CellConstIterator citer;
 		for (citer = cell_begin_end.first; citer != cell_begin_end.second; ++citer)
 	  {
	
			const PHG4CylinderCellDefs::keytype ckey = citer->first;
			PHG4CylinderCell* cell = _cells->findCylinderCell(ckey);

			// Loop over each hit in cell and get particle info
			PHG4CylinderCell::EdepConstRange hit_begin_end = cell->get_g4hits();
			PHG4CylinderCell::EdepConstIterator hiter;
 			for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; ++hiter)
  		{
				const PHG4HitDefs::keytype hkey = hiter->first;
				PHG4Hit* hit = _hits->findHit(hkey);
	
				PHG4Particle* particle = truthinfo->GetParticle( hit->get_trkid() );
				const int pid = particle->get_pid();
	
				// This is the total energy the particle had at impact:
				const double edep = hit->get_edep();

				// Right now, no em particle in the hcal
				if ( detector != "CEMC" && is_em_particle(pid) )
					continue;

				// Smear the deposited energy by particle type
				const double edep_smrd = smear_e( edep, pid ); 

				// If particle energy is above spread threshhold
				if ( edep_smrd > _spread_e_thresh )
				{
					if (detector == "CEMC")
					 spread_e_emc( raw_tower, tkey, edep_smrd, pid );
					else
					 spread_e_hcal( raw_tower, tkey, edep_smrd, pid );
				}
				else
					add_calib_tower( raw_tower, tkey, edep_smrd );

			}// for (hiter....
		}// for (citer....

	} //  for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)

  if (verbosity)
  {
    std::cout << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
			<< "input sum energy = " << _raw_towers->getTotalEdep()
			<< ", output sum digitalized value = " << _calib_towers->getTotalEdep()
			<< std::endl;
	}
  return Fun4AllReturnCodes::EVENT_OK;
}

int
RawTowerFastSim::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}


double RawTowerFastSim::smear_e(double e, int pid)
{
	double e_smrd = 0;

	if ( is_em_particle(pid) )
	{
		double em_gaus = gsl_ran_gaussian(RandomGenerator, 
				//edep * (0.11*pow( edep, -0.5) + 0.005) );
				e * (_em_sigma_p0*pow( e, -0.5) + _em_sigma_p1) );
		e_smrd = e + em_gaus;
	}
	else if ( is_hadron(pid) )
	{
		double h_gaus = gsl_ran_gaussian(RandomGenerator, 
				//edep * (0.175*pow( edep, -0.5) + 0.09) );
				e * (_h_sigma_p0*pow( e, -0.5) + _h_sigma_p1) );
		e_smrd = e + h_gaus;
	}
	else
		e_smrd = 0;

	// Force energy to be non-negative
	if (e_smrd < 0)
		e_smrd = 0;

	return e_smrd;
}

// FIXME:: This needs some reworking to make the hadrons realistic
//
//
void RawTowerFastSim::spread_e_emc(
		const RawTower* raw_tower, unsigned int key, double e, int pid)
{
	const int twr1_bineta = raw_tower->get_bineta();
	const int twr1_binphi = raw_tower->get_binphi();
	double twr1_e = e;

	//cout << "TOTAL E = " << e << endl;

	const int nphibins = rawtowergeom->get_phibins();
	const int netabins = rawtowergeom->get_etabins();
	//cout << "N ETABINS = " << netabins << ", N PHIBINS = " << nphibins << endl;
	
	// Iterate through surrounding towers (5x5) and deposit energy
	for ( int i = -2; i < 3; ++i )
	{
		for (int j = -2; j < 3; ++j )
		{
			// Set the absolute bin positions
			int bineta = twr1_bineta + i;
			int binphi = twr1_binphi + j;
	
			// Make sure we don't fall off the edge in eta!!
			if (bineta < 0 || bineta >= netabins)
			{
				//cout << "OVER EDGE!!!!" << endl;
				continue;
			}
	
			// handle the rap-around in phi
			if (binphi >= nphibins)  binphi = binphi - nphibins;
			else if (binphi < 0)     binphi = nphibins + binphi;
	
			// throw random for the tower's energy then subtract for e cons.
			unsigned int abs_eta = TMath::Abs( i );
			unsigned int abs_phi = TMath::Abs( j );
			double twr_e = 0;
			if ( is_em_particle(pid) )
			{
				// prob = .06
				if ( abs_eta + abs_phi == 1)
					twr_e = gsl_ran_flat(RandomGenerator, 0., 0.12*e );
				// prob = .012  
				else if ( abs_eta == 1 && abs_phi == 1)
					twr_e = gsl_ran_flat(RandomGenerator, 0., 0.024*e );
				// prob = .004
				else if ( (abs_eta == 2 && abs_phi == 0) || 
						      (abs_eta == 0 && abs_phi == 2) )
					twr_e = gsl_ran_flat(RandomGenerator, 0., 0.008*e );
				else
					twr_e = 0;
			}
			else if ( is_hadron(pid) )
			{
				// prob = .071
				if ( abs_eta + abs_phi == 1)
					twr_e = gsl_ran_flat(RandomGenerator, 0., 0.142*e );
				// prob = .027 
				else if ( abs_eta == 1 && abs_phi == 1)
					twr_e = gsl_ran_flat(RandomGenerator, 0., 0.054*e );
				// prob = .01
				else if ( (abs_eta == 2 && abs_phi == 0) || 
						      (abs_eta == 0 && abs_phi == 2) )
					twr_e = gsl_ran_flat(RandomGenerator, 0., 0.02*e );
				// prob = .007 
				else if ( (abs_eta == 2 && abs_phi == 1) || 
						      (abs_eta == 1 && abs_phi == 2) )
					twr_e = gsl_ran_flat(RandomGenerator, 0., 0.014*e );
				else
					twr_e = 0;
			}
	
			twr1_e = twr1_e - twr_e;
			//cout << "TWR BIN PHI = " << binphi << ", BIN ETA = " << bineta << ", E = " << twr_e << endl;
		
			// Add Tower to the container
			add_calib_tower( bineta, binphi, twr_e );
		}// for ( j....
	}// for ( i....
	
	//cout << "CENTER TWR BIN PHI = " << twr1_binphi << ", BIN ETA = " << twr1_bineta << ", E = " << twr1_e << endl;
		
	// Add the center bin
	add_calib_tower(raw_tower, key, twr1_e);

	return;
}


void RawTowerFastSim::spread_e_hcal(
		const RawTower* raw_tower, unsigned int key, double e, int pid)
{
	// right now, assuming em don't make it to the hcal
	if ( !is_hadron(pid) )
		return;

	// parameterizations of percent of total energy in surrounding bins
	double prob[5];
	prob[0] = -0.03337*pow(e,-1./2.) + 0.0741;
	prob[1] = -0.01853*pow(e,-1./2.) + 0.0292;
	prob[2] = -0.00865*pow(e,-1./2.) + 0.0118;
	prob[3] = -0.00675*pow(e,-1./2.) + 0.0088;
	prob[4] = -0.00355*pow(e,-1./2.) + 0.0045;

	const int twr1_bineta = raw_tower->get_bineta();
	const int twr1_binphi = raw_tower->get_binphi();
	double twr1_e = e;

	//cout << "TOTAL E = " << e << endl;

	const int nphibins = rawtowergeom->get_phibins();
	const int netabins = rawtowergeom->get_etabins();
	//cout << "N ETABINS = " << netabins << ", N PHIBINS = " << nphibins << endl;
	
	// Iterate through surrounding towers (5x5) and deposit energy
	for ( int i = -2; i < 3; ++i )
	{
		for (int j = -2; j < 3; ++j )
		{
			// Set the absolute bin positions
			int bineta = twr1_bineta + i;
			int binphi = twr1_binphi + j;
	
			// Make sure we don't fall off the edge in eta!!
			if (bineta < 0 || bineta >= netabins)
			{
				//cout << "OVER EDGE!!!!" << endl;
				continue;
			}
	
			// handle the rap-around in phi
			if (binphi >= nphibins)  binphi = binphi - nphibins;
			else if (binphi < 0)     binphi = nphibins + binphi;
	
			// throw random for the tower's energy then subtract for e cons.
			unsigned int abs_eta = TMath::Abs( i );
			unsigned int abs_phi = TMath::Abs( j );
			double twr_e = 0;
			
			if ( abs_eta + abs_phi == 1)
				twr_e = gsl_ran_flat(RandomGenerator, 0., 2*prob[0]*e );
			else if ( abs_eta == 1 && abs_phi == 1 )
				twr_e = gsl_ran_flat(RandomGenerator, 0., 2*prob[1]*e );
			else if ( (abs_eta == 2 && abs_phi == 0) || 
					      (abs_eta == 0 && abs_phi == 2) )
				twr_e = gsl_ran_flat(RandomGenerator, 0., 2*prob[2]*e );
			else if ( (abs_eta == 2 && abs_phi == 1) || 
					      (abs_eta == 1 && abs_phi == 2) )
				twr_e = gsl_ran_flat(RandomGenerator, 0., 2*prob[3]*e );
			else if ( abs_eta == 2 && abs_phi == 2 ) 
				twr_e = gsl_ran_flat(RandomGenerator, 0., 2*prob[4]*e );
			else
				twr_e = 0;
	
			// subtract this amount from center bin
			twr1_e = twr1_e - twr_e;
			//cout << "TWR BIN PHI = " << binphi << ", BIN ETA = " << bineta << ", E = " << twr_e << endl;
		
			// Add Tower to the container
			add_calib_tower( bineta, binphi, twr_e );
		}// for ( j....
	}// for ( i....
	
	//cout << "CENTER TWR BIN PHI = " << twr1_binphi << ", BIN ETA = " << twr1_bineta << ", E = " << twr1_e << endl;
		
	// Add the center bin
	add_calib_tower(raw_tower, key, twr1_e);

	return;
}

void 
RawTowerFastSim::add_calib_tower(
		const RawTower* raw_tower, RawTowerDefs::keytype key, double e)
{
	// Check if the tower exists
	//   if yes:  add the smeared energy to existing
	//   if no:  make new tower and add to container
	if(e < 0)
		return;

	RawTower* got_twr = _calib_towers->getTower( key );
	if ( got_twr )
		got_twr->set_energy( e + got_twr->get_energy() );
	else
	{
		RawTower *calib_tower = new RawTowerv1( *raw_tower );
   	calib_tower->set_energy( e );
   	_calib_towers->AddTower( key, calib_tower );
	}
	return;
}

// Same as above, but adds by bin
void 
RawTowerFastSim::add_calib_tower( int etabin, int phibin, double e )
{
	if(e < 0)
		return;

	RawTower* got_twr = _calib_towers->getTower( etabin, phibin );
	if ( got_twr )
		got_twr->set_energy( e + got_twr->get_energy() );
	else
	{
		RawTower* twr = new RawTowerv1( etabin , phibin );
		twr->set_energy( e );
		_calib_towers->AddTower( etabin, phibin, twr );
	}

	return;
}

bool RawTowerFastSim::is_em_particle(int pid)
{
	if (pid ==  11 ||
			pid == -11 ||
			pid ==  22)
		return true;
	else return false;
}

bool RawTowerFastSim::is_hadron(int pid)
{
	if (pid < -100 || pid > 100)
		return true;
	else return false;
}

void
RawTowerFastSim::CreateNodes(PHCompositeNode *topNode)
{

  PHNodeIterator iter(topNode);
  PHCompositeNode *runNode = static_cast<PHCompositeNode*>(iter.findFirst(
      "PHCompositeNode", "RUN"));
  if (!runNode)
    {
      std::cerr << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
          << "Run Node missing, doing nothing." << std::endl;
      throw std::runtime_error(
          "Failed to find Run node in RawTowerFastSim::CreateNodes");
    }

  TowerGeomNodeName = "TOWERGEOM_" + detector;
  rawtowergeom = findNode::getClass<RawTowerGeomContainer>(topNode,
      TowerGeomNodeName.c_str());
  if (!rawtowergeom)
    {
      std::cerr << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
          << " " << TowerGeomNodeName << " Node missing, doing bail out!"
          << std::endl;
      throw std::runtime_error(
          "Failed to find " + TowerGeomNodeName
              + " node in RawTowerFastSim::CreateNodes");
    }

  if (verbosity >= 1)
    {
      rawtowergeom->identify();
    }

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst(
      "PHCompositeNode", "DST"));
  if (!dstNode)
    {
      std::cerr << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
          << "DST Node missing, doing nothing." << std::endl;
      throw std::runtime_error(
          "Failed to find DST node in RawTowerFastSim::CreateNodes");
    }

  HitNodeName = "G4HIT_" + detector;
	_hits = findNode::getClass<PHG4HitContainer>(topNode , 
			HitNodeName.c_str());
  if (!_hits)
    {
      std::cerr << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
          << " " << HitNodeName << " Node missing, doing bail out!"
          << std::endl;
      throw std::runtime_error(
          "Failed to find " + HitNodeName
              + " node in RawTowerFastSim::CreateNodes");
    }

  CellNodeName = "G4CELL_" + detector;
	_cells = findNode::getClass<PHG4CylinderCellContainer>(topNode , 
			CellNodeName.c_str());
  if (!_cells)
    {
      std::cerr << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
          << " " << CellNodeName << " Node missing, doing bail out!"
          << std::endl;
      throw std::runtime_error(
          "Failed to find " + CellNodeName
              + " node in RawTowerFastSim::CreateNodes");
    }

  RawTowerNodeName = "TOWER_" + _raw_tower_node_prefix + "_" + detector;
  _raw_towers = findNode::getClass<RawTowerContainer>(dstNode,
      RawTowerNodeName.c_str());
  if (!_raw_towers)
    {
      std::cerr << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
          << " " << RawTowerNodeName << " Node missing, doing bail out!"
          << std::endl;
      throw std::runtime_error(
          "Failed to find " + RawTowerNodeName
              + " node in RawTowerFastSim::CreateNodes");
    }

	TruthInfoNodeName = "G4TruthInfo";
  truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode,
			TruthInfoNodeName.c_str());
  if (!truthinfo) 
		{
      std::cerr << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
          << " " << TruthInfoNodeName << " Node missing, doing bail out!"
          << std::endl;
      throw std::runtime_error(
          "Failed to find " + TruthInfoNodeName
              + " node in RawTowerFastSim::CreateNodes");
  }

  // Create the tower nodes on the tree
  PHNodeIterator dstiter(dstNode);
  PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode*>(dstiter.findFirst(
      "PHCompositeNode", detector));
  if (!DetNode)
    {
      DetNode = new PHCompositeNode(detector);
      dstNode->addNode(DetNode);
    }

  _calib_towers = new RawTowerContainer( RawTowerDefs::convert_name_to_caloid( detector ) );
  CaliTowerNodeName = "TOWER_" + _calib_tower_node_prefix + "_" + detector;
  PHIODataNode<PHObject> *towerNode = new PHIODataNode<PHObject>(_calib_towers,
      CaliTowerNodeName.c_str(), "PHObject");
  DetNode->addNode(towerNode);

  return;
}

