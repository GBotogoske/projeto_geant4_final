#include "G4UserTrackingAction.hh"
#include "G4Track.hh"
#include <map>
#include "G4ThreeVector.hh"


struct TrackInfo
{
    G4ThreeVector vertex;
    G4int parentID;
};

extern std::map<G4int, TrackInfo> trackMap;

class MyTrackingAction : public G4UserTrackingAction 
{
public:
    virtual void PreUserTrackingAction(const G4Track* track) override;
};
