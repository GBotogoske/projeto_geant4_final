//GammaBari G4Week - By Corrado Altomare & Davide Serini :D
//Ex2-Physics By Corrado
//

#ifndef OneHit_h
#define OneHit_h 1

//Include Native G4 Classes
#include "G4VHit.hh"
#include "G4Allocator.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4UnitsTable.hh"
#include "G4Track.hh"

class OneHit : public G4VHit{
  public:

    OneHit();
    virtual ~OneHit();

    // methods from base class



    void SetPhotonCounter_VIS (G4int p)                       {fPhotonCounter_VIS = p;};
    G4int GetPhotonCounter_VIS()                   const { return fPhotonCounter_VIS;};

    void SetPhotonCounter_UV (G4int p)                       {fPhotonCounter_UV = p;};
    G4int GetPhotonCounter_UV()                   const { return fPhotonCounter_UV;};

  private:
 
    G4int         fPhotonCounter_VIS;
    G4int         fPhotonCounter_UV;
  };

typedef G4THitsCollection<OneHit> HitCollection;


#endif
