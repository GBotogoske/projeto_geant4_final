// TrackingAction.hh
#include "G4UserTrackingAction.hh"
#include "G4Track.hh"
#include <map>

extern std::map<G4int, G4ThreeVector> trackBirthMap;

class MyTrackingAction : public G4UserTrackingAction {
public:
   
    virtual void PreUserTrackingAction(const G4Track* track) override 
    {
        trackBirthMap[track->GetTrackID()] = track->GetVertexPosition();
    }
};
