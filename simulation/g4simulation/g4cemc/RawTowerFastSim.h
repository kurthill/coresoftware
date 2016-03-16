#ifndef RawTowerFastSim_H__
#define RawTowerFastSim_H__

#include <fun4all/SubsysReco.h>
#include <string>
#include "RawTowerDefs.h"

#include <phool/PHTimeServer.h>

class PHCompositeNode;
class PHG4HitContainer;
class PHG4CylinderCellContainer;
class RawTowerContainer;
class RawTowerGeomContainer;
class PHG4TruthInfoContainer;
class RawTower;

// rootcint barfs with this header so we need to hide it
#ifndef __CINT__
#include <gsl/gsl_rng.h>
#endif

//! default input DST node is TOWER_RAW_DETECTOR
//! default output DST node is TOWER_CALIB_DETECTOR
class RawTowerFastSim : public SubsysReco
{

public:
  RawTowerFastSim(const std::string& name = "RawTowerFastSim");
  virtual
  ~RawTowerFastSim()
  {
  }

  int
  InitRun(PHCompositeNode *topNode);
  int
  process_event(PHCompositeNode *topNode);
  int
  End(PHCompositeNode *topNode);
	
  void
  Detector(const std::string &d)
	{
    detector = d;
	}

  double get_em_sigma_p0() const 
	{
		return _em_sigma_p0;
	}
  void set_em_sigma_p0(double p0) 
	{
		_em_sigma_p0 = p0;
	}

  double get_em_sigma_p1() const 
	{
		return _em_sigma_p1;
	}
	void set_em_sigma_p1(double p1) 
	{
		_em_sigma_p1 = p1;
	}

  double get_h_sigma_p0() const 
	{
		return _h_sigma_p0;
	}
  
	void set_h_sigma_p0(double p0)
	{
		_h_sigma_p0 = p0;
	}

  double get_h_sigma_p1() const 
	{
		return _h_sigma_p1;
	}
  
	void set_h_sigma_p1(double p1)
	{
		_h_sigma_p1 = p1;
	}

  std::string
  get_calib_tower_node_prefix() const
	{
    return _calib_tower_node_prefix;
	}

  void
  set_calib_tower_node_prefix(std::string calibTowerNodePrefix)
	{
    _calib_tower_node_prefix = calibTowerNodePrefix;
	}

  std::string
  get_raw_tower_node_prefix() const
	{
    return _raw_tower_node_prefix;
	}

  void
  set_raw_tower_node_prefix(std::string rawTowerNodePrefix)
	{
    _raw_tower_node_prefix = rawTowerNodePrefix;
	}

  double
  get_zero_suppression_GeV() const
	{
    return _zero_suppression_GeV;
	}

  void
  set_zero_suppression_GeV(double zeroSuppressionGeV)
	{
    _zero_suppression_GeV = zeroSuppressionGeV;
	}

  double
  get_spread_e_thresh() const
	{
    return _spread_e_thresh;
	}

  void
  set_spread_e_thresh(double spread_e_thresh)
	{  
		_spread_e_thresh = spread_e_thresh;
	}

protected:

	PHG4HitContainer* _hits;
	PHG4CylinderCellContainer* _cells;
  RawTowerContainer* _calib_towers;
  RawTowerContainer* _raw_towers;
  RawTowerGeomContainer *rawtowergeom;
	PHG4TruthInfoContainer* truthinfo;

  std::string detector;
  std::string HitNodeName;
  std::string CellNodeName;
  std::string RawTowerNodeName;
  std::string CaliTowerNodeName;
  std::string TowerGeomNodeName;
	std::string TruthInfoNodeName;

  std::string _calib_tower_node_prefix;
  std::string _raw_tower_node_prefix;

  //! em smear sigma function params p0/sqrt(E) + p1
  double _em_sigma_p0;
  double _em_sigma_p1;

  //! hadron smear sigma function params p0/sqrt(E) + p1
  double _h_sigma_p0;
  double _h_sigma_p1;
	
  //! zero suppression in unit of GeV
  double _zero_suppression_GeV;
	
	//! only spread energy into surrounding towers if above this
	double _spread_e_thresh;

  PHTimeServer::timer _timer;

	unsigned int seed;


  /////////////////////// 
	//  private methods ///
	//
	
  void
  CreateNodes(PHCompositeNode *topNode);
	
	void 
	spread_e_emc( 
			const RawTower* raw_tower, 
			unsigned int key, 
			double e, 
			int pid );

	void 
	spread_e_hcal( 
			const RawTower* raw_tower, 
			unsigned int key, 
			double e, 
			int pid );

	void
	add_calib_tower( 
			const RawTower* raw_tower, 
			unsigned int tkey,
			double e );

	void 
	add_calib_tower( 
			int etabin, 
			int phibin, 
			double e );

	double 
	smear_e( 
			double e, 
			int pid );
	
	bool 
	is_em_particle(int pid);
	
	bool 
	is_hadron(int pid);

	#ifndef __CINT__
  gsl_rng *RandomGenerator;
	#endif
};

#endif /* RawTowerFastSim_H__ */
