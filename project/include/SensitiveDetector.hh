#ifndef SensitiveDetector_h
#define SensitiveDetector_h 1

#include "DetectorConstruction.hh"
#include "OneHit.hh"
#include "G4VSensitiveDetector.hh"
#include "G4TouchableHistory.hh"
#include "G4OpBoundaryProcess.hh"
#include <vector>

class G4Step;
class G4HCofThisEvent;

class SensitiveDetector : public G4VSensitiveDetector
{
public:
  //Constructor
  SensitiveDetector(const G4String &SDname, const G4String &HitCollectionName);
  //Distructor
  virtual ~SensitiveDetector();

  void Initialize(G4HCofThisEvent *HCE);
  G4bool ProcessHits(G4Step *step, G4TouchableHistory *ROhist);
  void EndOfEvent(G4HCofThisEvent *HCE);

  // Add methods
  // // Total deposit Energy
  void AddEdep(G4double edep) {fE += edep;};
  void DeleteTotalE()  { fE = 0.;};
  G4double GetTotalE()  const {return fE;};

  // // Photon Counter
  void SetCounterStatus(G4int p) {fP +=p;};
  void ResetCounterStatus() {fP = 0;};
  G4int GetCounterStatus() const {return fP;};
  
  // // // // // //
  G4bool IsAnOpticalPhoton(G4Step* aStep);
  void CleanSDMemory();
  void PrintSDMemoryStatus();

private:
  HitCollection *fHitCollection;
  G4double      fE;
  G4int         fP;
};

#endif
