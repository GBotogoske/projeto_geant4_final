#ifndef OneHitUV_h
#define OneHitUV_h 1

#include "G4VHit.hh"
#include "G4Allocator.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4UnitsTable.hh"
#include "G4Track.hh"

class OneHitUV : public G4VHit
{
  public:

    OneHitUV();
    virtual ~OneHitUV();

    // set methods
    void SetEDep(G4double E)                              {fEDep = E;};
    void SetTrackLength(G4double L)                       {fTrackLength = L;};
    void SetPhotonCounter (G4int p)                       {fPhotonCounter = p;};

    // get methods
    G4double GetEDep()                         const { return fEDep;};
    G4double GetTrackLength()                  const { return fTrackLength;};
    G4int GetPhotonCounter()                   const { return fPhotonCounter;};

    private:
    G4double      fEDep;
    G4double      fTrackLength;
    G4int         fPhotonCounter;
  };

typedef G4THitsCollection<OneHitUV> HitUVCollection;


#endif
