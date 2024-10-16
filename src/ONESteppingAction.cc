
#include "ONESteppingAction.hh"
#include "ONEEventAction.hh"
#include "ONEDetectorConstruction.hh"
#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
//------------------------------------------------------------------------------
ONESteppingAction::ONESteppingAction(ONEEventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fLV00(0)
{
 myOUT .open( "out.data" , std::ios::trunc);
}
//------------------------------------------------------------------------------
ONESteppingAction::~ONESteppingAction(){
  myOUT.close();}
//------------------------------------------------------------------------------
void ONESteppingAction::UserSteppingAction(const G4Step* step)
{
  if ( !fLV00 ){
    const ONEDetectorConstruction* detectorConstruction
      = static_cast<const ONEDetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

    fLV00 = detectorConstruction->GetLV00();
  }

  // get volume of the current step
  G4LogicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()
                                           ->GetVolume()->GetLogicalVolume();


  int vol=-1;
  // check if we are in scoring volume
  if (volume == fLV00) vol= 0 ;

//  if (vol!=0 || vol!=5 || vol!=6 || vol!=7 || vol!=8 || vol!=15) return;

  int    ev_id = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID() ;


  G4Track* trk = step->GetTrack();
//  int    tr_c  = trk->GetDefinition()->GetPDGCharge();
  int    tr_id = trk->GetTrackID() ;
//  double tr_m  = trk->GetDefinition()->GetPDGMass()  ;
  G4int p_code = trk->GetDefinition()->GetPDGEncoding();

  int st_id = 0;
  if( step->IsFirstStepInVolume() ){ st_id = 1; }
  if( step->IsLastStepInVolume()  ){ st_id = 2; }
  if( trk->GetKineticEnergy ()==0 ){ st_id = 3; }


  if(st_id>-1){
    G4StepPoint* pre_step  ; pre_step  = step->GetPreStepPoint()  ;
    G4StepPoint* post_step ; post_step = step->GetPostStepPoint() ;
    double tr_ed = step->GetTotalEnergyDeposit() - step->GetNonIonizingEnergyDeposit() ;
    //double tr_rd = step->GetNonIonizingEnergyDeposit() ;

    double tr_post_x   = post_step->GetPosition().x() / mm ;
    double tr_post_y   = post_step->GetPosition().y() / mm ;
    double tr_post_z   = post_step->GetPosition().z() / mm ;
    double tr_post_t   = post_step->GetGlobalTime ()  / ns ;

    double tr_pre_x   = pre_step->GetPosition().x() / mm ;
    double tr_pre_y   = pre_step->GetPosition().y() / mm ;
    double tr_pre_z   = pre_step->GetPosition().z() / mm ;
    double tr_pre_t   = pre_step->GetGlobalTime ()  / ns ;

//    double tr_x   =  0.5 * (tr_pre_x + tr_post_x);
//    double tr_y   =  0.5 * (tr_pre_y + tr_post_y);
//    double tr_z   =  0.5 * (tr_pre_z + tr_post_z);
//    double tr_t   =  0.5 * (tr_pre_t + tr_post_t);

    if( myOUT.is_open() && vol==0 && tr_ed>0 )
       myOUT << ev_id  << " " << tr_id << " " << p_code << " " << tr_ed
             << " " << tr_pre_x  << " " << tr_pre_y  << " " << tr_pre_z  << " " << tr_pre_t
             << " " << tr_post_x << " " << tr_post_y << " " << tr_post_z << " " << tr_post_t
             << G4endl;

 }
//  // collect energy deposited in this step
//  G4double edepStep = step->GetTotalEnergyDeposit();
//  fEventAction->AddEdep(edepStep);
}
//------------------------------------------------------------------------------

