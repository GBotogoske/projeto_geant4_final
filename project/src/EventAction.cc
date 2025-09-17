#include "EventAction.hh"
#include "RunAction.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4AnalysisManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction()
: G4UserEventAction(),
  fHCID(-1)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
HitCollection*
EventAction::GetHitsCollection(G4int HCID, const G4Event* event) const
{
  auto hitsCollection = static_cast<HitCollection*>(event->GetHCofThisEvent()->GetHC(HCID));

  if ( ! hitsCollection ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << HCID;
    G4Exception("B4cEventAction::GetHitsCollection()",
      "MyCode0003", FatalException, msg);
  }
  return hitsCollection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::BeginOfEventAction(const G4Event* event)
{
    G4int nEvt = event->GetEventID();
    auto analysisManager = G4AnalysisManager::Instance();
    if ( fHCID == -1 ) 
    {
        fHCID = G4SDManager::GetSDMpointer()->GetCollectionID("HitCollection");
    }
    // ***Primary Hit Collection*** //
    auto HC = GetHitsCollection(fHCID, event);
    if(HC->entries()>0)
    {
            G4int PhotonDetected_UV = 0;
            G4int PhotonDetected_VIS = 0;

            PhotonDetected_UV = (*HC)[0]->GetPhotonCounter_UV();
            PhotonDetected_VIS = (*HC)[0]->GetPhotonCounter_VIS();

            analysisManager->FillNtupleIColumn(0,0,  nEvt);
            analysisManager->FillNtupleDColumn(0,1,  PhotonDetected_VIS);
            analysisManager->FillNtupleDColumn(0,2,  PhotonDetected_UV);
            analysisManager->AddNtupleRow(0);
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{
  G4cout<<"\t !!! EndOfEvent "<<event->GetEventID()<<"!!!\n"<<G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
