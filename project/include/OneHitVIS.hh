#ifndef OneHitVIS_h
#define OneHitVIS_h 1

#include "G4VHit.hh"
#include "G4Allocator.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4UnitsTable.hh"
#include "G4Track.hh"

class OneHitVIS : public G4VHit
{
  public:

    OneHitVIS();
    virtual ~OneHitVIS();

    // set methods
    void SetEDep(G4double E)                              {fEDep = E;};
    void SetPhotonCounter (G4int p)                       {fPhotonCounter = p;};
    //void SetPhotonCounter (G4int d)                       {fDetectorNumber;};

    // get methods
    G4double GetEDep()                         const { return fEDep;};
    G4int GetPhotonCounter()                   const { return fPhotonCounter;};
    G4int GetDetectorNumber()                   const { return fDetectorNumber;};

    private:
    G4double      fEDep;
    G4int         fPhotonCounter;
    G4int         fDetectorNumber;
  };

typedef G4THitsCollection<OneHitVIS> HitVISCollection;


#endif
