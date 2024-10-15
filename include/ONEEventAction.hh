#ifndef ONEEventAction_h
#define ONEEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class ONERunAction;

/// Event action class
///

class ONEEventAction : public G4UserEventAction
{
  public:
    ONEEventAction(ONERunAction* runAction);
    virtual ~ONEEventAction();
    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);
  private:
    ONERunAction* fRunAction;
};
//------------------------------------------------------------------------------
#endif

