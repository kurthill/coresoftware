#ifndef PHPy6ParticleTrigger_h
#define PHPy6ParticleTrigger_h

#include "PHPy6GenTrigger.h"
#include <HepMC/GenEvent.h>

namespace HepMC
{
  class GenEvent;
};

/**
 * Fun4All module based in PHPythia8/PHPy8ParticleTrigger
 *
 * Provides trigger for events gerenated by PHPythia6 event generator (Pythia6 wrapper for Fun4All).
 *
 * This trigger makes use of the conversion of Pythia6 events to HepMC format done in PHPythia6.
 *
 */
class PHPy6ParticleTrigger: public PHPy6GenTrigger
{
public:

  //! constructor
  PHPy6ParticleTrigger(const std::string &name = "PHPy6ParticleTriggerger");

  //! destructor
  ~PHPy6ParticleTrigger( void ){}

#ifndef __CINT__
  bool Apply(const HepMC::GenEvent* evt);
#endif

  void AddParticles(std::string particles);
  void AddParticles(int particle);
  void AddParticles(std::vector<int> particles);

  void AddParents(std::string parents);
  void AddParents(int parent);
  void AddParents(std::vector<int> parents);

  void SetPtHigh(double pt);
  void SetPtLow(double pt);
  void SetPtHighLow(double ptHigh, double ptLow);

  void SetPHigh(double p);
  void SetPLow(double p);
  void SetPHighLow(double pHigh, double pLow);

  void SetEtaHigh(double eta);
  void SetEtaLow(double eta);
  void SetEtaHighLow(double etaHigh, double etaLow);

  void SetAbsEtaHigh(double eta);
  void SetAbsEtaLow(double eta);
  void SetAbsEtaHighLow(double etaHigh, double etaLow);

  void SetPzHigh(double pz);
  void SetPzLow(double pz);
  void SetPzHighLow(double pzHigh, double pzLow);

  void PrintConfig();

protected:

  // trigger variables
  std::vector<int> _theParents;
  std::vector<int> _theParticles;

  double _theEtaHigh, _theEtaLow;
  double _thePtHigh, _thePtLow;
  double _thePHigh, _thePLow;
  double _thePzHigh, _thePzLow;

  bool _doEtaHighCut, _doEtaLowCut, _doBothEtaCut;
  bool _doAbsEtaHighCut, _doAbsEtaLowCut, _doBothAbsEtaCut;
  bool _doPtHighCut, _doPtLowCut, _doBothPtCut;
  bool _doPHighCut, _doPLowCut, _doBothPCut;
  bool _doPzHighCut, _doPzLowCut, _doBothPzCut;

};

#endif
