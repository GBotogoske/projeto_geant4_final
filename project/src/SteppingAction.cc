#include "SteppingAction.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4LogicalVolume.hh"
#include "G4EventManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "G4RunManager.hh"
#include "RunAction.hh"

SteppingAction::SteppingAction() {}
SteppingAction::~SteppingAction() {}

G4ThreadLocal int cont_PEN=0;

void SteppingAction::UserSteppingAction(const G4Step* step)
{

    auto track = step->GetTrack();
    G4String particle = track->GetDefinition()->GetParticleName();

    // PreStep volume info
    G4StepPoint *aPrePoint = step->GetPreStepPoint();
    G4VPhysicalVolume *aPrePV = aPrePoint->GetPhysicalVolume();
    G4String PreVolName = "";
    if (aPrePV)
        PreVolName = aPrePV->GetName();

    // PostStep volume info
    G4StepPoint *aPostPoint = step->GetPostStepPoint();
    G4VPhysicalVolume *aPostPV = aPostPoint->GetPhysicalVolume();
    G4String PostVolName = "";
    if (aPostPV)
        PostVolName = aPostPV->GetName();

    if (particle == "opticalphoton")
    {
        if(G4StrUtil::contains(PreVolName,"argon") && G4StrUtil::contains(PostVolName,"PEN"))
        {
            auto Energy = aPostPoint->GetTotalEnergy();
            if(Energy>9*eV && Energy<12*eV)
            {
                auto runAction = const_cast<RunAction*>(static_cast<const RunAction*>(G4RunManager::GetRunManager()->GetUserRunAction()));
                runAction->AddPenCount(1);
            }
        }
    }

    /* auto track = step->GetTrack();
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
    } */
}
