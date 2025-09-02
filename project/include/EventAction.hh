#ifndef EventAction_h
#define EventAction_h 1

#include "OneHit.hh"
#include "G4UserEventAction.hh"
#include "globals.hh"

class RunAction;

/// Event action class
///

class EventAction : public G4UserEventAction
{
  public:
    EventAction();
    virtual ~EventAction();

    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);

  private:
    HitCollection *GetHitsCollection(G4int HCID,const G4Event* event) const;
    G4int fHCID;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
