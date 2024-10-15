#ifndef ONEDetectorConstruction_h
#define ONEDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include <fstream>

class G4VPhysicalVolume;
class G4LogicalVolume;

/// Detector construction class to define materials and geometry.

class ONEDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    ONEDetectorConstruction();
    virtual ~ONEDetectorConstruction();

    virtual G4VPhysicalVolume* Construct();

    G4LogicalVolume* GetLV00() const { return fLV00; }

  protected:
    G4LogicalVolume* fLV00      ;
    std::ifstream    ConfigFile ;
    double           PRESSURE   ;

};

#endif

