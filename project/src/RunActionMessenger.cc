#include "globals.hh"
#include "Randomize.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
 
#include "RunAction.hh"
#include "RunActionMessenger.hh"
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RunActionMessenger::RunActionMessenger(RunAction* runAction)
  : fRunAction (runAction)
{
  fSetNameOfOutputFile=
    new G4UIcmdWithAString("/RunManager/NameOfOutputFile",this);
  fSetNameOfOutputFile->SetGuidance("Name of Output File");
  fSetNameOfOutputFile->SetParameterName("Name of Output File",false);
  fSetNameOfOutputFile->AvailableForStates(G4State_PreInit,G4State_Idle);
  fSetNameOfOutputFile->SetToBeBroadcasted(false);
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RunActionMessenger::~RunActionMessenger() {
  delete fSetNameOfOutputFile;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RunActionMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if( command == fSetNameOfOutputFile)
    fRunAction->SetNameOfOutputFile(newValue);
}
 
