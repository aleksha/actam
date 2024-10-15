#include "ONEEventAction.hh"
#include "ONERunAction.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
//------------------------------------------------------------------------------
ONEEventAction::ONEEventAction(ONERunAction* runAction)
: G4UserEventAction(),
  fRunAction(runAction)
{}
//------------------------------------------------------------------------------
ONEEventAction::~ONEEventAction(){}
//------------------------------------------------------------------------------
void ONEEventAction::BeginOfEventAction(const G4Event*){}
//------------------------------------------------------------------------------
void ONEEventAction::EndOfEventAction(const G4Event*){}
//------------------------------------------------------------------------------
