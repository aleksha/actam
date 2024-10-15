#ifndef ONEActionInitialization_h
#define ONEActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

/// Action initialization class.

class ONEActionInitialization : public G4VUserActionInitialization
{
  public:
    ONEActionInitialization();
    virtual ~ONEActionInitialization();

    virtual void BuildForMaster() const;
    virtual void Build() const;
};

#endif


