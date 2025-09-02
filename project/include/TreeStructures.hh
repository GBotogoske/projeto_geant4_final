#ifndef TreeStructures_h
#define TreeStructures_h 1

#include <vector>

const G4int VectorSizeMax    = 10000;

//////////////////////////
//TREE OUTPUT STRUCTURE //
//////////////////////////

//Event tree Structure
typedef struct{

  //Scalar
  int RunID;
  G4double PosX;
  G4double PosY;
  G4double PosZ;
  G4double DirX;
  G4double DirY;
  G4double DirZ;

  //Vector
  int VectorSize;
  int VectorSizeInStep;
  //Enumerator
  std::vector<std::string> VolumeName;
  G4int    ModuleEnum[VectorSizeMax];
  G4int    VolumeEnum[VectorSizeMax];
  //Energy
  G4double TotalEnergyDep[VectorSizeMax];
  //Secondary
  G4double TotalPhotonGenerated[VectorSizeMax];
  G4double TotalPhotonGeneratedWLS[VectorSizeMax];
  G4double TotalPhotonGeneratedScint[VectorSizeMax];

  G4double NumbOfPhotonsAbsorbed[VectorSizeMax];
  G4double NumbOfPhotonsAbsorbedWLS[VectorSizeMax];
  G4double NumbOfPhotonsAbsorbedScint[VectorSizeMax];
  //In Step Energy
  G4double InStepEnergy[VectorSizeMax];
  G4double InStepPosX[VectorSizeMax];
  G4double InStepPosY[VectorSizeMax];
  G4double InStepPosZ[VectorSizeMax];


}

  TRACKEREVENTN;

//Optical Photons Struct
typedef struct{

  //Scalar
  int RunID;
  //Enumerator
  char VolumeName[100];
  G4int    ModuleEnum;
  G4int    VolumeEnum;
  //Energy
  G4double KEne;
  //Position
  G4double PosX;
  G4double PosY;
  G4double PosZ;
  //Time
  G4double GlobalTime;
  G4double LocalTime;
  //CreatorProc
  char ProcessName[100];
  G4int ProcessNumber;
}

  TRACKEREVENTNP;

//Variables
static TRACKEREVENTN trackereventn;
static TRACKEREVENTNP trackereventnep;

#endif
