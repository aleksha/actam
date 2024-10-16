#include "ONEPrimaryGeneratorAction.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include <iostream>

//------------------------------------------------------------------------------
ONEPrimaryGeneratorAction::ONEPrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0)
{

  ConfigFile.open("CONFIG.txt", std::ios::in);
  double PRESSURE;
  ConfigFile >> PRESSURE >> PID >> T_R_MIN >> T_R_MAX >> THETA_R_MIN >> THETA_R_MAX >> VX >> VY>> VZ;
  ConfigFile.close();

  double D2R = 3.14159265/180.;
  THETA_R_MIN = THETA_R_MIN*D2R;
  THETA_R_MAX = THETA_R_MAX*D2R;


  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;

  G4ParticleDefinition* particle;
  if(PID<1 || PID>4) particle = particleTable->FindParticle( particleName="proton"  );
  if(PID==1)         particle = particleTable->FindParticle( particleName="deuteron");
  if(PID==2)         particle = particleTable->FindParticle( particleName="triton"  );
  if(PID==3)         particle = particleTable->FindParticle( particleName="He3"     );
  if(PID==4)         particle = particleTable->FindParticle( particleName="alpha"   );
  if(PID==11)        particle = particleTable->FindParticle( particleName="e-"      );
  if(PID==13)        particle = particleTable->FindParticle( particleName="mu-"     );
  if(PID==22)        particle = particleTable->FindParticle( particleName="gamma"   );
  if(PID==-11)       particle = particleTable->FindParticle( particleName="e+"      );
  if(PID==-13)       particle = particleTable->FindParticle( particleName="mu+"     );

  fParticleGun->SetParticleDefinition(particle);

}
//------------------------------------------------------------------------------
ONEPrimaryGeneratorAction::~ONEPrimaryGeneratorAction(){ delete fParticleGun; }
//------------------------------------------------------------------------------
void ONEPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of ecah event
  THETA_R = THETA_R_MIN + G4UniformRand()*( THETA_R_MAX - THETA_R_MIN  );
  T_R     = T_R_MIN     + G4UniformRand()*( T_R_MAX     - T_R_MIN      );
  fParticleGun->SetParticleEnergy(T_R*MeV);
  fParticleGun->SetParticleMomentumDirection( G4ThreeVector( 0.,std::sin(THETA_R),std::cos(THETA_R) ) );
  fParticleGun->SetParticlePosition( G4ThreeVector(VX*mm,VY*mm,VZ*mm) );
  fParticleGun->GeneratePrimaryVertex(anEvent);
}
//------------------------------------------------------------------------------

