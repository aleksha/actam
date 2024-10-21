#ifndef ONEPrimaryGeneratorAction_h
#define ONEPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"
#include <fstream>

class G4ParticleGun;
class G4Event;
//class G4Box;

/// The primary generator action class with particle gun.

class ONEPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    ONEPrimaryGeneratorAction();
    virtual ~ONEPrimaryGeneratorAction();

    // method from the base class
    virtual void GeneratePrimaries(G4Event*);

    // method to access particle gun
    const G4ParticleGun* GetParticleGun() const { return fParticleGun; }

  private:
    G4ParticleGun*  fParticleGun; // pointer a to G4 gun class

    std::ifstream    ConfigFile ;
    std::ofstream    GenFile    ;
    int              EV_ID      ;
    int              PID        ;
    double           T_R        ;
    double           T_R_MIN    ;
    double           T_R_MAX    ;
    double           THETA_R    ;
    double           THETA_R_MIN;
    double           THETA_R_MAX;
    double           VX,VY,VZ   ;

};
//------------------------------------------------------------------------------
#endif
