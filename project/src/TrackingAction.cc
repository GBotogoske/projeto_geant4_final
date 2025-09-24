#include "TrackingAction.hh"

std::map<G4int, TrackInfo> trackMap;

void MyTrackingAction::PreUserTrackingAction(const G4Track* track) 
{
    trackMap[track->GetTrackID()] = { track->GetVertexPosition(), track->GetParentID() };
}
