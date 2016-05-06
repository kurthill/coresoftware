#ifndef __BBCVERTEXFASTSIMRECO_H__
#define __BBCVERTEXFASTSIMRECO_H__

//===========================================================
/// \file BbcVertexFastSimReco.h
/// \brief simple truth vertex smearing algorithm
/// \author Mike McCumber
//===========================================================

#include <fun4all/SubsysReco.h>

// rootcint barfs with this header so we need to hide it
#ifndef __CINT__
#include <gsl/gsl_rng.h>
#endif


class PHCompositeNode;

/// \class BbcVertexFastSimReco
///
/// \brief simple truth vertex smearing algorithm
///
class BbcVertexFastSimReco : public SubsysReco {

 public:
 
  BbcVertexFastSimReco(const std::string &name = "BbcVertexFastSimReco");
  virtual ~BbcVertexFastSimReco();
		
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  void set_t_smearing(const float t_smear) {_t_smear = t_smear;}
  void set_z_smearing(const float z_smear) {_z_smear = z_smear;}

 private:

  int CreateNodes(PHCompositeNode *topNode);

  float _t_smear;
  float _z_smear;

#ifndef __CINT__
  gsl_rng *RandomGenerator;
#endif

};

#endif // __BBCVERTEXFASTSIMRECO_H__
