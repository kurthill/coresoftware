#ifndef G4ScintillatorSlatTTree_H
#define G4ScintillatorSlatTTree_H

#include <fun4all/SubsysReco.h>
#include <string>

class Fun4AllHistoManager;
class PHCompositeNode;
class TH1;
class TH2;

class G4ScintillatorSlatTTree: public SubsysReco
{
 public:
  G4ScintillatorSlatTTree(const std::string &name = "SCINTILLATORSLATTTREE");
  virtual ~G4ScintillatorSlatTTree(){}

  //! full initialization
  int Init(PHCompositeNode *);

  //! event processing method
  int process_event(PHCompositeNode *);

  int End(PHCompositeNode *);

  void Detector(const std::string &det);

  void SaveScintillatorSlats(const int i=1) {saveslats = i;}

  void HistoFileName(const std::string &name) {_histofilename = name;}

 protected:
 std::string _detector;
 std::string _outnodename;
 std::string _slatnodename;
 std::string _histofilename;
 int saveslats;
 int evtno;
 Fun4AllHistoManager *hm;
 TH1 *etot_hist;
};


#endif
