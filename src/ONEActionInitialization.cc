//------------------------------------------------------------------------------
#include "ONEActionInitialization.hh"
#include "ONEPrimaryGeneratorAction.hh"
#include "ONERunAction.hh"
#include "ONEEventAction.hh"
#include "ONESteppingAction.hh"
//------------------------------------------------------------------------------
ONEActionInitialization::ONEActionInitialization()
 : G4VUserActionInitialization()
{}
//------------------------------------------------------------------------------
ONEActionInitialization::~ONEActionInitialization(){}
//------------------------------------------------------------------------------
void ONEActionInitialization::BuildForMaster() const
{
  ONERunAction* runAction = new ONERunAction;
  SetUserAction(runAction);
}
//------------------------------------------------------------------------------
void ONEActionInitialization::Build() const
{
  SetUserAction(new ONEPrimaryGeneratorAction);

  ONERunAction* runAction = new ONERunAction;
  SetUserAction(runAction);

  ONEEventAction* eventAction = new ONEEventAction(runAction);
  SetUserAction(eventAction);
  SetUserAction(new ONESteppingAction(eventAction));
}
//------------------------------------------------------------------------------
