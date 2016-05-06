#ifndef PHG4EICForwardEcalDetector_h
#define PHG4EICForwardEcalDetector_h

#include <g4main/PHG4Detector.h>
#include <PHG4ForwardEcalDetector.h>

#include <Geant4/globals.hh>
#include <Geant4/G4Types.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4RotationMatrix.hh>
#include <Geant4/G4Material.hh>

#include <string>
#include <map>
#include <vector>
#include <set>

class G4AssemblyVolume;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4VSolid;

/**
 * \file ${file_name}
 * \brief Module to build forward sampling Hadron calorimeterr (endcap) in Geant4
 * \author Nils Feege <nils.feege@stonybrook.edu>
 */

class PHG4EICForwardEcalDetector: public PHG4ForwardEcalDetector
{

public:

  //! constructor
  PHG4EICForwardEcalDetector( PHCompositeNode *Node, const std::string &dnam="BLOCK" );

  //! destructor
  virtual ~PHG4EICForwardEcalDetector();

  //! construct
  virtual void Construct( G4LogicalVolume* world );

  virtual void SetTowerDimensions(G4double dx, G4double dy, G4double dz) {
  _tower_dx = dx;
  _tower_dy = dy;
  _tower_dz = dz;
  }

  void SetMaterialScintillator( G4String material ) { _materialScintillator = material; }
  void SetMaterialAbsorber( G4String material ) { _materialAbsorber = material; }

private:

  G4LogicalVolume* ConstructTower();
  int PlaceTower(G4LogicalVolume* envelope , G4LogicalVolume* tower);
  int ParseParametersFromTable();

  struct towerposition {
    G4double x;
    G4double y;
    G4double z;
  } ;

  std::map< std::string, towerposition > _map_tower;

  /* ECAL tower geometry */
  G4double _tower_dx;
  G4double _tower_dy;
  G4double _tower_dz;

  G4String _materialScintillator;
  G4String _materialAbsorber;

};

#endif
