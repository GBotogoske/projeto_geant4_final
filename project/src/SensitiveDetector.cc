#include "SensitiveDetector.hh"
#include "OneHit.hh"
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
    fHitCollection(NULL),
    fE(0.),
    fP(0)
{
  G4cout<<"Creating SD with name: "<<SDname<<G4endl;
  collectionName.insert(HitCollectionName);
}

SensitiveDetector::~SensitiveDetector(){}

void SensitiveDetector::Initialize(G4HCofThisEvent *HCE)
{
  fHitCollection = new HitCollection(GetName(),collectionName[0]);
  static G4int HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); //<<-- this is to get an ID for the colletionName[0]
  //G4cout<<"*** "<<fHitCollection->GetName()<<" initialized [ID = "<<HCID<<"]"<<G4endl;
  HCE->AddHitsCollection(HCID, fHitCollection);
}


G4bool SensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *)
{
  G4double eDep   = aStep->GetTotalEnergyDeposit();
  G4int counter = 0;

  SensitiveDetector::AddEdep(eDep);
  
  G4String myParticle = "e-";
  G4String thisParticle = aStep->GetTrack()->GetParticleDefinition()->GetParticleName();
  G4String thisVolume   = aStep->GetTrack()->GetVolume()->GetName();
  
  G4int pDetected = 0;
  if(thisParticle=="opticalphoton"){
    G4String procName = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
    G4bool checkAbsorption = G4StrUtil::contains(procName,"Absorption");
    if(checkAbsorption==true and G4StrUtil::contains(thisVolume,"SiPM"))
      pDetected +=1;
    }
  SensitiveDetector::SetCounterStatus(pDetected);
  return true;
}

void SensitiveDetector::EndOfEvent(G4HCofThisEvent*)
{
  //Fill the hits
  OneHit *aHit = new OneHit();
  G4double TotE = SensitiveDetector::GetTotalE();
  G4int TotP = SensitiveDetector::GetCounterStatus();

  aHit->SetEDep(TotE);
  aHit->SetPhotonCounter(TotP);

  fHitCollection->insert(aHit);

  SensitiveDetector::PrintSDMemoryStatus();
  SensitiveDetector::CleanSDMemory();
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
G4bool SensitiveDetector::IsAnOpticalPhoton(G4Step* aStep){
  G4bool flag = false;
  G4Track *aTrack = aStep->GetTrack();
  const G4ParticleDefinition* aDef = aTrack->GetParticleDefinition();
  G4String aName = aDef->GetParticleName();
  if(aName=="opticalphoton"){
    flag=true;
  }
  return flag;
}

void SensitiveDetector::CleanSDMemory(){
  SensitiveDetector::DeleteTotalE();
  SensitiveDetector::ResetCounterStatus();
}

void SensitiveDetector::PrintSDMemoryStatus(){
  G4double TotE = SensitiveDetector::GetTotalE();
  G4int    TotP = SensitiveDetector::GetCounterStatus();
  G4cout<<"O.o.O.o.O.o.O.o.O.o.O.o.O.o.O.o.O.o.O.o.O"<<G4endl;
  G4cout<<" Total Energy deposited = "<<G4BestUnit(TotE,"Energy")<<G4endl;
  G4cout<<" Total Photon detected  = "<<TotP<<G4endl;
  G4cout<<"O.o.O.o.O.o.O.o.O.o.O.o.O.o.O.o.O.o.O.o.O"<<G4endl;
}
