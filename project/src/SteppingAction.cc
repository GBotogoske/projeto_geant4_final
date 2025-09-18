#include "SteppingAction.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4LogicalVolume.hh"
#include "G4EventManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

SteppingAction::SteppingAction() {}
SteppingAction::~SteppingAction() {}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
    auto track = step->GetTrack();
    auto pos = track->GetPosition();
    auto preVolume = step->GetPreStepPoint()->GetPhysicalVolume(); // <-- seguro
    auto eventID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
    auto Energy = step->GetPreStepPoint()->GetTotalEnergy();

    if(preVolume) {
        G4cout << "EventID: " << eventID
               << " TrackID: " << track->GetTrackID()
               << " Volume: " << preVolume->GetName()
               << " Position: " << pos 
               << " Energy: " << Energy << G4endl;
    } else {
        G4cout << "EventID: " << eventID
               << " TrackID: " << track->GetTrackID()
               << " Volume: [NULL - saiu do mundo]"
               << " Position: " << pos << G4endl;
    }
}
