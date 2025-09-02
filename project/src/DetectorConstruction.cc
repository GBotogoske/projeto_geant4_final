#include "DetectorConstruction.hh"
#include "OneHit.hh"
#include "OptFileManager.hh"

#include "SensitiveDetector.hh"
#include "G4LogicalBorderSurface.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4NistManager.hh"
#include <math.h>



// Option to switch on/off checking of volumes overlaps
G4bool checkOverlaps = true;

DetectorConstruction::DetectorConstruction():G4VUserDetectorConstruction()
{
    std::ifstream f("../configuration/detector.json");
    f >> this->config;
}

DetectorConstruction::~DetectorConstruction(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* DetectorConstruction::Construct()
{

    G4NistManager *fNistManager = G4NistManager::Instance();
    fNistManager->SetVerbose(0);
  
    // ----  World -----

    auto config_world = config["world"];
    G4double world_sizeX = config_world["size"][0].get<double>()*cm;
    G4double world_sizeY = config_world["size"][1].get<double>()*cm;
    G4double world_sizeZ = config_world["size"][2].get<double>()*cm;

    G4Material* world_mat = fNistManager->FindOrBuildMaterial(config_world["material"].get<string>()); 
    G4Box* solidWorld = new G4Box("World",0.5*world_sizeX, 0.5*world_sizeY, 0.5*world_sizeZ);
    G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld,world_mat,"World");
    G4VPhysicalVolume* physicalWorld = new G4PVPlacement(0,G4ThreeVector(),logicWorld,"World",0,false,0,checkOverlaps);

    auto mpt_World = new G4MaterialPropertiesTable();
    mpt_World->AddProperty("RINDEX", {1,1}, {config_world["refraction_index"].get<double>(),config_world["refraction_index"].get<double>()}, 2);
    mpt_World->AddProperty("ABSLENGTH", {1,1}, {config_world["abs_length"].get<double>(),config_world["abs_length"].get<double>()}, 2);

    world_mat->SetMaterialPropertiesTable(mpt_World);

    // ----- Cryostat ----

    auto config_cryostat = config["cryostat"];
    G4double cryostat_sizeX = config_cryostat["size"][0].get<double>()*cm;
    G4double cryostat_sizeY = config_cryostat["size"][1].get<double>()*cm;
    G4double cryostat_sizeZ = config_cryostat["size"][2].get<double>()*cm;
    G4double cryostatThickness = config_cryostat["thickness"].get<double>()*cm;

    G4Material* cryostat_mat = fNistManager->FindOrBuildMaterial(config_cryostat["material"].get<string>());
    G4Box* solidExternalCryostat = new G4Box("Cryostat_ext",0.5*cryostat_sizeX, 0.5*cryostat_sizeY, 0.5*cryostat_sizeZ);
    G4Box* solidInternalCryostat = new G4Box("Cryostat_int",
        0.5*cryostat_sizeX-cryostatThickness, 0.5*cryostat_sizeY-cryostatThickness, 0.5*cryostat_sizeZ-cryostatThickness);
    auto solidCryostat = new G4SubtractionSolid("Cryostat", solidExternalCryostat,solidInternalCryostat,nullptr,G4ThreeVector());
    
    G4LogicalVolume* logicCryostat = new G4LogicalVolume(solidCryostat,cryostat_mat,"Cryostat");
    G4VPhysicalVolume* physicalCryostat = new G4PVPlacement(0,G4ThreeVector(),logicCryostat,"Cryostat",logicWorld,false,0,checkOverlaps);
    

    auto mpt_Cryo = new G4MaterialPropertiesTable();
    mpt_Cryo->AddProperty("RINDEX", {1,1}, {config_cryostat["refraction_index"].get<double>(),config_cryostat["refraction_index"].get<double>()}, 2);
    mpt_Cryo->AddProperty("ABSLENGTH", {1,1}, {config_cryostat["abs_length"].get<double>(),config_cryostat["abs_length"].get<double>()}, 2);

    cryostat_mat->SetMaterialPropertiesTable(mpt_Cryo);

    // -----  Liquid Argon Filling --------
    auto config_LAr = config["liquid_argon"];

    G4Material* lar_mat = fNistManager->FindOrBuildMaterial(config_cryostat["filling"].get<string>());
    G4LogicalVolume* logicCryostatFilling = new G4LogicalVolume(solidInternalCryostat,lar_mat,"LAr");
    G4VPhysicalVolume* physicalCryostatFilling = new G4PVPlacement(0,G4ThreeVector(),logicCryostatFilling,"LAr",logicWorld,false,0,checkOverlaps);   

    std::string file_lar_index = config_LAr["refraction_index"].get<string>();
    std::ifstream f_lar_index(std::string("../configuration/")+file_lar_index);
    json data = json::parse(f_lar_index);
    size_t n_rindex = data.size();
    std::vector<G4double> energies_r(n_rindex);
    std::vector<G4double> rindex(n_rindex);
    for (size_t i = 0; i < n_rindex; i++) 
    {
        energies_r[i] = data[i]["E"].get<double>()*eV;  
        rindex[i]   = data[i]["r"].get<double>();
    }
       
    auto abs_lar = config_LAr["abs_length"];
    size_t n_abs = abs_lar["e"].size();
    std::vector<G4double> absE_LAr(n_abs);
    std::vector<G4double> absLen_LAr(n_abs);
    for (size_t i=0; i<n_abs; i++) 
    {
        absE_LAr[i] = abs_lar["e"][i].get<double>()*eV;    
        absLen_LAr[i] = abs_lar["l"][i].get<double>()*cm;   
    }

    auto scint_lar = config_LAr["scintilation"];
    size_t n_scint = scint_lar["e_s"].size();
    std::vector<G4double> scintE_LAr(n_scint);
    std::vector<G4double> scint(n_scint);
    for (size_t i=0; i<n_scint; i++) 
    {
        scintE_LAr[i] = scint_lar["e_s"][i].get<double>()*eV;    
        scint[i] = scint_lar["s"][i].get<double>();   
    }

    std::string file_lar_rayleigh = config_LAr["rayleigh_scattering"].get<string>();
    std::ifstream f_lar_rayleigh(std::string("../configuration/")+file_lar_rayleigh);
    json data_ray = json::parse(f_lar_rayleigh);
    size_t n_rayleigh = data_ray.size();
    std::vector<G4double> energies_ray(n_rayleigh);
    std::vector<G4double> ray_lar(n_rayleigh);
    for (size_t i = 0; i < n_rayleigh; i++) 
    {
        energies_ray[i] = data_ray[i]["E"].get<double>()*eV;  
        ray_lar[i]   = data_ray[i]["l"].get<double>()*cm;
    }

    auto mpt_LAr = new G4MaterialPropertiesTable();
    mpt_LAr->AddProperty("RINDEX", energies_r.data(), rindex.data(), n_abs);
    mpt_LAr->AddProperty("ABSLENGTH", absE_LAr.data(), absLen_LAr.data(), n_abs);
    mpt_LAr->AddProperty("SCINTILLATIONCOMPONENT1",scintE_LAr,scint,n_scint);
    mpt_LAr->AddProperty("SCINTILLATIONCOMPONENT2",scintE_LAr,scint,n_scint);
    mpt_LAr->AddConstProperty("SCINTILLATIONYIELD",scint_lar["LY"].get<double>()/MeV);
    mpt_LAr->AddConstProperty("SCINTILLATIONTIMECONSTANT1",scint_lar["Slow_Comp"].get<double>()*ns);
    mpt_LAr->AddConstProperty("SCINTILLATIONTIMECONSTANT2",scint_lar["Fast_Comp"].get<double>()*ns);
    mpt_LAr->AddConstProperty("SCINTILLATIONYIELD1", scint_lar["Slow_Amp"].get<double>());
    mpt_LAr->AddConstProperty("SCINTILLATIONYIELD2", scint_lar["Fast_Amp"].get<double>());
    mpt_LAr->AddConstProperty("RESOLUTIONSCALE",scint_lar["Resolution"].get<double>());
    mpt_LAr->AddProperty("RAYLEIGH", energies_ray.data(), ray_lar.data(), n_rayleigh);
 
    lar_mat->GetIonisation()->SetBirksConstant(scint_lar["Birks"].get<double>()*cm/MeV); // se comentar essa linha, o valor eh o default do argonio liquido
    lar_mat->SetMaterialPropertiesTable(mpt_LAr);

    // -----  Liquid Argon && Cryostat Interface Boundary --------

    G4OpticalSurface* surface_cryo_lar = new G4OpticalSurface("surface_cryo_lar");
    surface_cryo_lar-> SetModel(unified);
    surface_cryo_lar-> SetType(dielectric_metal);
    surface_cryo_lar-> SetFinish(polished);
    G4MaterialPropertiesTable* mpt_CryoLar_Surface = new G4MaterialPropertiesTable();
    mpt_CryoLar_Surface->AddProperty("REFLECTIVITY",{1,1}, {config_LAr["reflectivity_cryo"].get<double>(),config_LAr["reflectivity_cryo"].get<double>()},2);
    surface_cryo_lar->SetMaterialPropertiesTable(mpt_CryoLar_Surface);

    new G4LogicalBorderSurface("LiquidArgon-->Cryo", physicalCryostatFilling, physicalCryostat , surface_cryo_lar );

    // ----------------------------------------------------------------------------
    //Assembling the field cage
    auto config_FC = config["FC"];
    G4Material* FC_mat = fNistManager->FindOrBuildMaterial(config_FC["material"].get<string>());

    auto mpt_FC = new G4MaterialPropertiesTable();
    mpt_FC->AddProperty("RINDEX", {1,1}, {config_FC["refraction_index"].get<double>(),config_FC["refraction_index"].get<double>()}, 2);
    mpt_FC->AddProperty("ABSLENGTH", {1,1}, {config_FC["abs_length"].get<double>(),config_FC["abs_length"].get<double>()}, 2);
    FC_mat->SetMaterialPropertiesTable(mpt_FC);

    G4double d_cryo = config_FC["distance_to_wall"].get<double>(); // distance from the membrane

        //Starting with the vertical bars
    G4double FCv_sizeX = config_FC["vertical_bar"][0].get<double>()*cm;
    G4double FCv_sizeY = config_FC["vertical_bar"][1].get<double>()*cm;
    G4double FCv_sizeZ = config_FC["vertical_bar"][2].get<double>()*cm;
    G4double FCv_distance = config_FC["vertical_bar_distance"].get<double>()*cm;

    G4Box* solidBARv = new G4Box("vertical_bar",0.5*FCv_sizeX, 0.5*FCv_sizeY, 0.5*FCv_sizeZ);
    
    int n_barv_long = (cryostat_sizeX-2*cryostatThickness-2*d_cryo)/FCv_distance;
    int n_barv_long_2 = (cryostat_sizeX)/FCv_distance;

    G4cout << "#######  " << n_barv_long <<  "  #######  \n";
    G4cout << "#######  " << n_barv_long_2 <<  "  #######  \n";

    G4LogicalVolume* LogicalFC = new G4LogicalVolume(solidBARv,FC_mat,"FC");
    G4VPhysicalVolume* physicalFC = new G4PVPlacement(0,G4ThreeVector(),LogicalFC,"FC",logicCryostatFilling,false,0,checkOverlaps);

    // -----  Liquid Argon && FC Interface Boundary --------

    G4OpticalSurface* surface_FC_lar = new G4OpticalSurface("surface_FC_lar");
    surface_FC_lar-> SetModel(unified);
    surface_FC_lar-> SetType(dielectric_metal);
    surface_FC_lar-> SetFinish(polished);
    G4MaterialPropertiesTable* mpt_FCLar_Surface = new G4MaterialPropertiesTable();
    mpt_FCLar_Surface->AddProperty("REFLECTIVITY",{1,1}, {config_FC["reflectivity"].get<double>(),config_FC["reflectivity"].get<double>()},2);
    surface_FC_lar->SetMaterialPropertiesTable(mpt_FCLar_Surface);

    new G4LogicalBorderSurface("LiquidArgon-->FC", physicalCryostatFilling, physicalFC , surface_FC_lar );

    return physicalWorld;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
