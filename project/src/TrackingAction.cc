#include "TrackingAction.hh"
#include "G4EventManager.hh"

G4ThreadLocal std::map<G4int, TrackInfo> trackMap;

void MyTrackingAction::PreUserTrackingAction(const G4Track* track) 
{
    /* auto vtx = track->GetVertexPosition();
    G4cout << "[Tracking] Event=" 
           << G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID()
           << " TrackID=" << track->GetTrackID()
           << " ParentID=" << track->GetParentID()
           << " Vertex=" << vtx
           << G4endl; */

    trackMap[track->GetTrackID()] = { track->GetVertexPosition(), track->GetParentID() };
}
