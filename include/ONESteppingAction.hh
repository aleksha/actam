#ifndef ONESteppingAction_h
#define ONESteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

class ONEEventAction;

class G4LogicalVolume;

/// Stepping action class
///

class ONESteppingAction : public G4UserSteppingAction
{
  public:
    ONESteppingAction(ONEEventAction* eventAction);
    virtual ~ONESteppingAction();

    // method from the base class
    virtual void UserSteppingAction(const G4Step*);

  private:
    ONEEventAction*  fEventAction;
    std::ofstream myOUT;

    G4LogicalVolume* fLV00;
};
//------------------------------------------------------------------------------
#endif
