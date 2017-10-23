#ifndef __JETEVALUATOR_H__
#define __JETEVALUATOR_H__

//===============================================
/// \file JetEvaluator.h
/// \brief Compares reconstructed jets to truth jets
/// \author Michael P. McCumber
//===============================================

#include "JetEvalStack.h"

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class TFile;
class TNtuple;

/// \class JetEvaluator
///
/// \brief Compares reconstructed jets to truth jets
///
/// Plan: This module will trace the jet constituents to
/// the greatest contributor Monte Carlo jet and then
/// test one against the other.
///
class JetEvaluator : public SubsysReco {

 public:
 
  JetEvaluator(const std::string &name = "JETEVALUATOR",
	       const std::string &recojetname = "AntiKt_Tower_r0.3",
	       const std::string &truthjetname = "AntiKt_Truth_r0.3",
	       const std::string &filename = "g4eval_jets.root");
  virtual ~JetEvaluator() {};
		
  int Init(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  void set_strict(bool b) {_strict = b;}
  
 private:

  std::string _recojetname;
  std::string _truthjetname;
  
  unsigned long _ievent;

  JetEvalStack* _jetevalstack;
  
  //----------------------------------
  // evaluator output ntuples

  bool _strict;
  unsigned int _errors;
  
  bool _do_recojet_eval;
  bool _do_truthjet_eval;
  
  TNtuple *_ntp_recojet;
  TNtuple *_ntp_truthjet;

  // evaluator output file
  std::string _filename;
  TFile *_tfile;

  // subroutines
  void printInputInfo(PHCompositeNode *topNode);    ///< print out the input object information (debugging upstream components)
  void fillOutputNtuples(PHCompositeNode *topNode); ///< dump the evaluator information into ntuple for external analysis
  void printOutputInfo(PHCompositeNode *topNode);   ///< print out the ancestry information for detailed diagnosis
};

#endif // __JETEVALUATOR_H__
