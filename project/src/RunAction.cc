#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "RunActionMessenger.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4GeneralParticleSource.hh"
#include "G4AnalysisManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(): G4UserRunAction(), fOutputFileName("./Data")
{
    fRunMessenger = new RunActionMessenger(this);
    G4RunManager::GetRunManager()->SetPrintProgress(1);
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    G4cout << "Using " << analysisManager->GetType() << G4endl;

    //Create ntuple
    analysisManager->CreateNtuple("Primary_Hit", "Primary_Hit");
    analysisManager->CreateNtupleIColumn(0,"eventID");
    analysisManager->CreateNtupleDColumn(0,"PhotonDetectedVIS");
    analysisManager->CreateNtupleDColumn(0,"PhotonDetectedUV");
    analysisManager->FinishNtuple(0);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
    delete G4AnalysisManager::Instance();
    delete fRunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{
    G4String fileName   = fOutputFileName;
    G4String fileOutput = fileName;
    fileOutput += ".root";
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->OpenFile(fileOutput);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run*)
{
    auto analysisManager = G4AnalysisManager::Instance();
    // save histograms & ntuple
    G4cout<<"[INFO]: Write the output file "<<fOutputFileName<<".root"<<G4endl;
    analysisManager->Write();
    analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
