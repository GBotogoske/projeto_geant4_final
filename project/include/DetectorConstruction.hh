#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "G4RunManager.hh"
#include "globals.hh"
#include <vector>


#include <fstream>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

class G4VPhysicalVolume;
class G4LogicalVolume;
class OptMaterials;


/// Detector construction class to define materials and geometry.

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    virtual ~DetectorConstruction();
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField() override;  
    
  protected:
    json config;
    G4LogicalVolume* sipmvis_logical;
    G4LogicalVolume* sipmuv_logical;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
