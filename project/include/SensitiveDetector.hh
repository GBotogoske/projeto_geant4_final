#ifndef SensitiveDetector_h
#define SensitiveDetector_h 1

#include "DetectorConstruction.hh"
#include "OneHit.hh"
#include "G4VSensitiveDetector.hh"
#include "G4TouchableHistory.hh"
#include "G4OpBoundaryProcess.hh"
#include <vector>

#include "SiPMSpectrum.hh"

class G4Step;
class G4HCofThisEvent;

class SensitiveDetector : public G4VSensitiveDetector
{
public:
  
    SensitiveDetector(const G4String &SDname, const G4String &HitCollectionName);

    virtual ~SensitiveDetector();

    void Initialize(G4HCofThisEvent *HCE);
    G4bool ProcessHits(G4Step *step, G4TouchableHistory *ROhist);
    void EndOfEvent(G4HCofThisEvent *HCE);

    // // Photon Counter
    void SetCounterStatus_UV(G4int p) {fP_uv +=p;};
    void ResetCounterStatus_UV() {fP_uv = 0;};
    G4int GetCounterStatus_UV() const {return fP_uv;};

    void SetCounterStatus_VIS(G4int p) {fP_vis +=p;};
    void ResetCounterStatus_VIS() {fP_vis = 0;};
    G4int GetCounterStatus_VIS() const {return fP_vis;};


    // // // // // //
    G4bool IsAnOpticalPhoton(G4Step* aStep);
    void CleanSDMemory();
    void PrintSDMemoryStatus();

    private:
    HitCollection *fHitCollection;
    G4int         fP_vis;
    G4int         fP_uv;

};

#endif
