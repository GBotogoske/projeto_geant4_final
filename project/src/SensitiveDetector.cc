#include "SensitiveDetector.hh"
#include "OneHitVIS.hh"
#include "OneHitUV.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4TouchableHistoryHandle.hh"
#include "G4TouchableHistory.hh"
#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4Tubs.hh"
#include "G4String.hh"
#include "G4OpticalPhoton.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4GeometryTolerance.hh"
#include <math.h>

using namespace std;

static const G4double GeometricalTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
static const G4double h = CLHEP::h_Planck/((CLHEP::eV)*(CLHEP::ns));
static const G4double c = CLHEP::c_light/((CLHEP::nm)/(CLHEP::ns));

SensitiveDetector::SensitiveDetector(const G4String &SDname,const G4String &HitCollectionName)
  : G4VSensitiveDetector(SDname),
    fHitVISCollection(NULL),
    fE(0.),
    fP(0)
{
  G4cout<<"Creating SD with name: "<<SDname<<G4endl;
  collectionName.insert(HitCollectionName);
}

SensitiveDetector::~SensitiveDetector(){}

void SensitiveDetector::Initialize(G4HCofThisEvent *HCE)
{
    fHitVISCollection = new HitVISCollection(GetName(),collectionName[0]);
    static G4int HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); //<<-- this is to get an ID for the colletionName[0]
    //G4cout<<"*** "<<fHitCollection->GetName()<<" initialized [ID = "<<HCID<<"]"<<G4endl;
    HCE->AddHitsCollection(HCID, fHitVISCollection);
}


G4bool SensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *)
{
    G4double eDep   = aStep->GetTotalEnergyDeposit();
    G4int counter = 0;

    SensitiveDetector::AddEdep(eDep);

    G4String thisParticle = aStep->GetTrack()->GetParticleDefinition()->GetParticleName();
    G4String thisVolume   = aStep->GetTrack()->GetVolume()->GetName();

    G4Track* track = aStep->GetTrack();
    G4TrackStatus status = track->GetTrackStatus();

    G4int pDetected = 0;
    if(thisParticle=="opticalphoton")
    {
        G4String procName = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
        G4bool checkAbsorption = G4StrUtil::contains(procName,"Absorption");

        G4cout << thisVolume << std::endl;
        if(checkAbsorption==true and G4StrUtil::contains(thisVolume,"sipmVIS"))
        {
            const auto& sipm_spectrum = SiPMSpectrum::get();

            auto eff_vis = sipm_spectrum.get_effVIS();
            auto E_vis = sipm_spectrum.get_EVIS();
            auto n_vis = sipm_spectrum.getNVIS();

            auto Ephoton = aStep->GetTrack()->GetTotalEnergy();
            G4double p;
            if(Ephoton<=E_vis[0])
            {
                p=eff_vis[0];
            }
            else if(Ephoton>=E_vis[n_vis-1])
            {
               p=eff_vis[n_vis-1]; 
            }
            else
            {
                for (int i = 0; i < n_vis - 1; ++i) 
                {
                    if (Ephoton >= E_vis[i] && Ephoton <= E_vis[i + 1])
                    {
                        G4double x0 = E_vis[i];
                        G4double x1 = E_vis[i + 1];
                        G4double y0 = eff_vis[i];
                        G4double y1 = eff_vis[i + 1];
                        p = y0 + (y1 - y0) * (Ephoton - x0) / (x1 - x0);
                    }
                }
            }
            G4double r = G4RandFlat::shoot();  
            if(r<p)
            {
                pDetected +=1;
            }
        }
 
    }
    SensitiveDetector::SetCounterStatus(pDetected);
    return true;
}

void SensitiveDetector::EndOfEvent(G4HCofThisEvent*)
{
    //Fill the hits
    OneHitVIS *aHit = new OneHitVIS();
    G4double TotE = SensitiveDetector::GetTotalE();
    G4int TotP = SensitiveDetector::GetCounterStatus();

    aHit->SetEDep(TotE);
    aHit->SetPhotonCounter(TotP);

    fHitVISCollection->insert(aHit);

    SensitiveDetector::PrintSDMemoryStatus();
    SensitiveDetector::CleanSDMemory();
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
G4bool SensitiveDetector::IsAnOpticalPhoton(G4Step* aStep)
{
    G4bool flag = false;
    G4Track *aTrack = aStep->GetTrack();
    const G4ParticleDefinition* aDef = aTrack->GetParticleDefinition();
    G4String aName = aDef->GetParticleName();
    if(aName=="opticalphoton"){
    flag=true;
    }
    return flag;
}

void SensitiveDetector::CleanSDMemory()
{
    SensitiveDetector::DeleteTotalE();
    SensitiveDetector::ResetCounterStatus();
}

void SensitiveDetector::PrintSDMemoryStatus()
{
    G4double TotE = SensitiveDetector::GetTotalE();
    G4int    TotP = SensitiveDetector::GetCounterStatus();
    G4cout<<"O.o.O.o.O.o.O.o.O.o.O.o.O.o.O.o.O.o.O.o.O"<<G4endl;
    G4cout<<" Total Energy deposited = "<<G4BestUnit(TotE,"Energy")<<G4endl;
    G4cout<<" Total Photon detected  = "<<TotP<<G4endl;
    G4cout<<"O.o.O.o.O.o.O.o.O.o.O.o.O.o.O.o.O.o.O.o.O"<<G4endl;
}
