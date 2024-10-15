#ifndef ONERunAction_h
#define ONERunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

/// Run action class

class ONERunAction : public G4UserRunAction
{
  public:
    ONERunAction();
    virtual ~ONERunAction();

    // virtual G4Run* GenerateRun();
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);
};

#endif

