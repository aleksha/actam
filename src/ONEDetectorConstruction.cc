//------------------------------------------------------------------------------
#include "ONEDetectorConstruction.hh"
//------------------------------------------------------------------------------
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
//------------------------------------------------------------------------------
ONEDetectorConstruction::ONEDetectorConstruction()
: G4VUserDetectorConstruction(),
  fLV00(0)
{

  ConfigFile.open("CONFIG.txt", std::ios::in);
  ConfigFile >> PRESSURE ;
  ConfigFile.close();

}
//------------------------------------------------------------------------------
ONEDetectorConstruction::~ONEDetectorConstruction(){ }
//------------------------------------------------------------------------------
G4VPhysicalVolume* ONEDetectorConstruction::Construct()
{
  G4bool checkOverlaps = true;
  G4NistManager* nist = G4NistManager::Instance();


  G4Material* w_mat = nist->FindOrBuildMaterial("G4_Galactic");
  G4Material *ArGas   = new G4Material("ArGas"  , 18, 39.948*g/mole, PRESSURE*1.784*kg/m3 );

  // World
  G4double w_xy = 200.0*mm;
  G4double w_z  = 200.0*mm;
  G4Box* solidWorld = new G4Box("World", 0.5*w_xy, 0.5*w_xy, 0.5*w_z);

  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, w_mat, "World");

  // LVs

  // Layer (LV) geometrical parameters
  G4double R_xy  =  50.000*mm;
  G4double l00_z = 100.000*mm;
  G4Tubs* solidLV00 =  new G4Tubs("LV00", 0.*mm, R_xy, 0.5*l00_z, 0.*deg, 360.*deg);

  G4LogicalVolume* logicLV00 = new G4LogicalVolume(solidLV00, ArGas , "LV00");


  G4ThreeVector l00_pos; l00_pos.set(0,0,0);

  G4VPhysicalVolume* physWorld =
    new G4PVPlacement(0, G4ThreeVector(), logicWorld,
                         "World", 0, false, checkOverlaps);

  new G4PVPlacement(0, l00_pos, logicLV00, "LV00", logicWorld, false, 0, checkOverlaps);

  fLV00 = logicLV00;

  return physWorld;
}
//------------------------------------------------------------------------------
