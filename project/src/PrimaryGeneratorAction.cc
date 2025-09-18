#include "PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fGPS(0),fParticleGun(0)
{
    //fGPS  = new G4GeneralParticleSource();
    fParticleGun = new G4ParticleGun(1);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    //delete fGPS;
    delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    /*
    fGPS->SetParticlePosition(G4ThreeVector(0.,0.,0.));
    fGPS->GeneratePrimaryVertex(anEvent);
    */

    G4ParticleDefinition* photon = G4ParticleTable::GetParticleTable()->FindParticle("opticalphoton");

    for (int i = 0; i < 1; i++) 
    {
        
        fParticleGun->SetParticleDefinition(photon);
        fParticleGun->SetParticleEnergy(9.7*eV);       // 128 nm
        fParticleGun->SetParticlePosition(G4ThreeVector(0.0*m,2.0*m,2.0*m));

        G4double theta = std::acos(1 - 2*G4UniformRand()); // de 0 a pi
        G4double phi   = 2*M_PI*G4UniformRand();           // de 0 a 2pi
        G4double x = std::sin(theta)*std::cos(phi);
        G4double y = std::sin(theta)*std::sin(phi);
        G4double z = std::cos(theta);
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(x,y,z));
        fParticleGun->SetParticlePolarization(G4ThreeVector(0,0,1));

        fParticleGun->GeneratePrimaryVertex(anEvent);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
