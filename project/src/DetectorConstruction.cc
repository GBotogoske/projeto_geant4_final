#include "DetectorConstruction.hh"
#include "OneHit.hh"
#include "OptFileManager.hh"

#include "SensitiveDetector.hh"
#include "G4LogicalBorderSurface.hh"

#include "G4Material.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4EllipticalTube.hh"
#include "G4ExtrudedSolid.hh"
#include "G4AssemblyVolume.hh"
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

//nota de coordadenas (x,y,z) --> (direcao da tela, direcao do campo, direcao transversal da tela)

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

    // -----  Liquid Argon Filling --------
    auto config_LAr = config["liquid_argon"];

    auto config_cryostat = config["cryostat"];
    G4Material* lar_mat = fNistManager->FindOrBuildMaterial(config_cryostat["filling"].get<string>());
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

    G4Box* solidWorld = new G4Box("World",0.5*world_sizeX, 0.5*world_sizeY, 0.5*world_sizeZ);
    G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld,lar_mat,"World");
    G4VPhysicalVolume* physicalWorld = new G4PVPlacement(0,G4ThreeVector(),logicWorld,"World",0,false,0,checkOverlaps);

    // ----- Cryostat ----

    G4double cryostat_sizeX = config_cryostat["size"][0].get<double>()*cm;
    G4double cryostat_sizeY = config_cryostat["size"][1].get<double>()*cm;
    G4double cryostat_sizeZ = config_cryostat["size"][2].get<double>()*cm;
    G4double cryostatThickness = config_cryostat["thickness"].get<double>()*cm;

    G4Material* cryostat_mat = fNistManager->FindOrBuildMaterial(config_cryostat["material"].get<string>());
    G4Box* solidExternalCryostat = new G4Box("Cryostat_ext",0.5*cryostat_sizeX+cryostatThickness, 0.5*cryostat_sizeY+cryostatThickness, 0.5*cryostat_sizeZ+cryostatThickness);
    G4Box* solidInternalCryostat = new G4Box("Cryostat_int",0.5*cryostat_sizeX, 0.5*cryostat_sizeY, 0.5*cryostat_sizeZ);
    auto solidCryostat = new G4SubtractionSolid("Cryostat", solidExternalCryostat,solidInternalCryostat,nullptr,G4ThreeVector());
    
    G4LogicalVolume* logicCryostat = new G4LogicalVolume(solidCryostat,cryostat_mat,"Cryostat");
    G4VPhysicalVolume* physicalCryostat = new G4PVPlacement(0,G4ThreeVector(),logicCryostat,"Cryostat",logicWorld,false,0,checkOverlaps);
    
    auto mpt_Cryo = new G4MaterialPropertiesTable();
    mpt_Cryo->AddProperty("RINDEX", {1,1}, {config_cryostat["refraction_index"].get<double>(),config_cryostat["refraction_index"].get<double>()}, 2);
    mpt_Cryo->AddProperty("ABSLENGTH", {1,1}, {config_cryostat["abs_length"].get<double>(),config_cryostat["abs_length"].get<double>()}, 2);

    cryostat_mat->SetMaterialPropertiesTable(mpt_Cryo);

    // -----  Liquid Argon && Cryostat Interface Boundary --------

    G4OpticalSurface* surface_cryo_lar = new G4OpticalSurface("surface_cryo_lar");
    surface_cryo_lar-> SetModel(unified);
    surface_cryo_lar-> SetType(dielectric_metal);
    surface_cryo_lar-> SetFinish(polished);
    G4MaterialPropertiesTable* mpt_CryoLar_Surface = new G4MaterialPropertiesTable();
    mpt_CryoLar_Surface->AddProperty("REFLECTIVITY",{1,1}, {config_LAr["reflectivity_cryo"].get<double>(),config_LAr["reflectivity_cryo"].get<double>()},2);
    surface_cryo_lar->SetMaterialPropertiesTable(mpt_CryoLar_Surface);

    new G4LogicalBorderSurface("LiquidArgon-->Cryo", physicalWorld, physicalCryostat , surface_cryo_lar );

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
    
    int n_barv_long = (cryostat_sizeX-2*d_cryo)/FCv_distance;
    G4UnionSolid* FCv_united =  new G4UnionSolid("FCv_longwall_template", solidBARv, solidBARv, 0, G4ThreeVector(FCv_distance,0,0));
    for(int i=1;i<n_barv_long/2;i++)
    {
        FCv_united = new G4UnionSolid("FCv_longwall_template", FCv_united, solidBARv, 0, G4ThreeVector((i+1)*FCv_distance,0,0));
        FCv_united = new G4UnionSolid("FCv_longwall_template", FCv_united, solidBARv, 0, G4ThreeVector((-i)*FCv_distance,0,0));     
    } 

    int n_barv_end = (cryostat_sizeZ-2*d_cryo)/FCv_distance;
    G4UnionSolid* FCv_united_end =  new G4UnionSolid("FCv_endwall_template", solidBARv, solidBARv, 0, G4ThreeVector(0,0,FCv_distance));
    for(int i=1;i<n_barv_end/2;i++)
    {
        FCv_united_end = new G4UnionSolid("FCv_endwall_template", FCv_united_end, solidBARv, 0, G4ThreeVector(0,0,(i+1)*FCv_distance));
        FCv_united_end = new G4UnionSolid("FCv_endwall_template", FCv_united_end, solidBARv, 0, G4ThreeVector(0,0,(-i)*FCv_distance));     
    } 

    double barH_eixoX1 = config_FC["horizontal_bar_1_majorAxis"].get<double>()*cm;
    double barH_eixoY1 = config_FC["horizontal_bar_1_minorAxis"].get<double>()*cm;
    double barH_eixoX2 = config_FC["horizontal_bar_2_majorAxis"].get<double>()*cm;
    double barH_eixoY2 = config_FC["horizontal_bar_2_minorAxis"].get<double>()*cm;
 
    double cut_value1 = config_FC["cut_horizontal_bar1"].get<double>()*cm;
    double cut_value2 = config_FC["cut_horizontal_bar2"].get<double>()*cm;

    double barH_length = cryostat_sizeX-2*d_cryo-FCv_sizeX;
    double barH_length_end = cryostat_sizeZ-2*d_cryo-FCv_sizeZ-2*cut_value1;

    std::vector<G4TwoVector> ellipse1;
    int npoints = 200;
    for (int i=0; i<npoints; i++) 
    {
        double phi = 2*CLHEP::pi*i/npoints;
        double x = (barH_eixoX1/2.0)*cos(phi);
        double y = (barH_eixoY1/2.0)*sin(phi);
        y = std::max(y, -(cut_value1-barH_eixoY1/2));
        ellipse1.push_back(G4TwoVector(x,y));
    }
    std::vector<G4TwoVector> ellipse2;
    for (int i=0; i<npoints; i++) 
    {
        double phi = 2*CLHEP::pi*i/npoints;
        double x = (barH_eixoX2/2.0)*cos(phi);
        double y = (barH_eixoY2/2.0)*sin(phi);
        y = std::max(y, -(cut_value2-barH_eixoY2/2));
        ellipse2.push_back(G4TwoVector(x,y));
    }

    auto solidBARh1 = new G4ExtrudedSolid("solidEllipse1", ellipse1,  barH_length/2, G4TwoVector(0,0), 1.0,   G4TwoVector(0,0), 1.0);  
    auto solidBARh2 = new G4ExtrudedSolid("solidEllipse2", ellipse2,  barH_length/2, G4TwoVector(0,0), 1.0,   G4TwoVector(0,0), 1.0);  
    G4double FCh_distance = config_FC["horizontal_bar_distance"].get<double>()*cm;
    int n_barh_long = (FCv_sizeY)/FCh_distance;    

    G4UnionSolid* FCh_united =  new G4UnionSolid("FCh_longwall_template", solidBARh2, solidBARh2, 0, G4ThreeVector(FCh_distance,0,0));
    double sep_bar_type = config_FC["sep_bar_type"].get<double>()*cm;

    auto solidBARh1_end = new G4ExtrudedSolid("solidEllipse1_end", ellipse1,  barH_length_end/2, G4TwoVector(0,0), 1.0,   G4TwoVector(0,0), 1.0);  
    auto solidBARh2_end = new G4ExtrudedSolid("solidEllipse2_end", ellipse2,  barH_length_end/2, G4TwoVector(0,0), 1.0,   G4TwoVector(0,0), 1.0); 
    G4UnionSolid* FCh_united_end =  new G4UnionSolid("FCh_endwall_template", solidBARh2_end, solidBARh2_end, 0, G4ThreeVector(FCh_distance,0,0));

    for(int i=1;i<n_barh_long/2;i++)
    {
        if((i+1)*FCh_distance - FCh_distance/2 > sep_bar_type)
        {
            FCh_united = new G4UnionSolid("FCv_longwall_template", FCh_united, solidBARh1, 0, G4ThreeVector((i+1)*FCh_distance,0,0));
            FCh_united = new G4UnionSolid("FCv_longwall_template", FCh_united, solidBARh1, 0, G4ThreeVector((-i)*FCh_distance,0,0));
            FCh_united_end = new G4UnionSolid("FCv_longwall_template", FCh_united_end, solidBARh1_end, 0, G4ThreeVector((i+1)*FCh_distance,0,0));
            FCh_united_end = new G4UnionSolid("FCv_longwall_template", FCh_united_end, solidBARh1_end, 0, G4ThreeVector((-i)*FCh_distance,0,0));        
        }
        else
        {
            FCh_united = new G4UnionSolid("FCv_longwall_template", FCh_united, solidBARh2, 0, G4ThreeVector((i+1)*FCh_distance,0,0));
            FCh_united = new G4UnionSolid("FCv_longwall_template", FCh_united, solidBARh2, 0, G4ThreeVector((-i)*FCh_distance,0,0)); 
            FCh_united_end = new G4UnionSolid("FCv_longwall_template", FCh_united_end, solidBARh2_end, 0, G4ThreeVector((i+1)*FCh_distance,0,0));
            FCh_united_end = new G4UnionSolid("FCv_longwall_template", FCh_united_end, solidBARh2_end, 0, G4ThreeVector((-i)*FCh_distance,0,0));     
        }
        
    } 
    G4RotationMatrix* rotTopY = new G4RotationMatrix();
    rotTopY->rotateY(90*deg);
    G4RotationMatrix* rotTopX = new G4RotationMatrix();
    rotTopX->rotateX(90*deg);
    G4RotationMatrix* minus_rotTopX = new G4RotationMatrix();
    minus_rotTopX->rotateX(-90*deg);
    G4RotationMatrix* rotTopZ = new G4RotationMatrix();
    rotTopZ->rotateZ(90*deg);
    G4RotationMatrix* minus_rotTopZ = new G4RotationMatrix();
    minus_rotTopZ->rotateZ(-90*deg);
    G4RotationMatrix* rotFinal1  = new G4RotationMatrix(); 
    *rotFinal1= (*rotTopY) * (*rotTopX);
    G4RotationMatrix* rotFinal2  = new G4RotationMatrix(); 
    *rotFinal2= (*rotTopY) * (*minus_rotTopX);
   
    G4LogicalVolume* LogicalFCv_wall = new G4LogicalVolume(FCv_united,FC_mat,"FC");
    G4LogicalVolume* LogicalFCv_end = new G4LogicalVolume(FCv_united_end,FC_mat,"FC");
    G4LogicalVolume* LogicalFCh_wall = new G4LogicalVolume(FCh_united,FC_mat,"FC");
    G4LogicalVolume* LogicalFCh_end = new G4LogicalVolume(FCh_united_end,FC_mat,"FC");
    G4VisAttributes* visFC = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5)); 
    visFC->SetForceSolid(true); 
    LogicalFCv_wall->SetVisAttributes(visFC);
    LogicalFCv_end->SetVisAttributes(visFC);
    LogicalFCh_wall->SetVisAttributes(visFC);
    LogicalFCh_end->SetVisAttributes(visFC);


    G4VPhysicalVolume* physicalFC1_v = new G4PVPlacement(0,G4ThreeVector(-FCv_distance/2,0,cryostat_sizeZ/2-d_cryo),
        LogicalFCv_wall,"FCv_long_wall",logicWorld,true,1,checkOverlaps);
    G4VPhysicalVolume* physicalFC2_v = new G4PVPlacement(0,G4ThreeVector(-FCv_distance/2,0,-cryostat_sizeZ/2+d_cryo),
        LogicalFCv_wall,"FCv_long_wall",logicWorld,true,2,checkOverlaps);
    G4VPhysicalVolume* physicalFC3_v = new G4PVPlacement(0,G4ThreeVector(cryostat_sizeX/2-d_cryo,0,-FCv_distance/2),
        LogicalFCv_end,"FCv_end_wall",logicWorld,true,1,checkOverlaps);
    G4VPhysicalVolume* physicalFC4_v = new G4PVPlacement(0,G4ThreeVector(-cryostat_sizeX/2+d_cryo,0,-FCv_distance/2),
        LogicalFCv_end,"FCv_end_wall",logicWorld,true,2,checkOverlaps);

    G4VPhysicalVolume* physicalFC1_h = new G4PVPlacement( rotFinal1 
        ,G4ThreeVector(0,-FCh_distance/2, cryostat_sizeZ/2-d_cryo-FCv_sizeZ/2-(cut_value1-barH_eixoY1/2)),
        LogicalFCh_wall,"FCh_long_wall",logicWorld,true,1,checkOverlaps);
    G4VPhysicalVolume* physicalFC2_h = new G4PVPlacement( rotFinal2 
        ,G4ThreeVector(0,-FCh_distance/2, -(cryostat_sizeZ/2-d_cryo-FCv_sizeZ/2-(cut_value1-barH_eixoY1/2))),
        LogicalFCh_wall,"FCh_long_wall",logicWorld,true,2,checkOverlaps);
    G4VPhysicalVolume* physicalFC3_h = new G4PVPlacement( rotTopZ 
        ,G4ThreeVector(-(cryostat_sizeX/2-d_cryo-FCv_sizeX/2-(cut_value1-barH_eixoY1/2)),-FCh_distance/2, 0),
        LogicalFCh_end,"FCh_end_wall",logicWorld,true,1,checkOverlaps);
    G4VPhysicalVolume* physicalFC4_h = new G4PVPlacement( minus_rotTopZ 
        ,G4ThreeVector((cryostat_sizeX/2-d_cryo-FCv_sizeX/2-(cut_value1-barH_eixoY1/2)),-FCh_distance/2, 0),
        LogicalFCh_end,"FCh_end_wall",logicWorld,true,2,checkOverlaps);            
    

    // -----  Liquid Argon && FC Interface Boundary --------

    G4OpticalSurface* surface_FC_lar = new G4OpticalSurface("surface_FC_lar");
    surface_FC_lar-> SetModel(unified);
    surface_FC_lar-> SetType(dielectric_metal);
    surface_FC_lar-> SetFinish(polished);
    G4MaterialPropertiesTable* mpt_FCLar_Surface = new G4MaterialPropertiesTable();
    mpt_FCLar_Surface->AddProperty("REFLECTIVITY",{1,1}, {config_FC["reflectivity"].get<double>(),config_FC["reflectivity"].get<double>()},2);
    surface_FC_lar->SetMaterialPropertiesTable(mpt_FCLar_Surface);

    new G4LogicalBorderSurface("LiquidArgon-->FC1v", physicalWorld, physicalFC1_v , surface_FC_lar );
    new G4LogicalBorderSurface("LiquidArgon-->FC2v", physicalWorld, physicalFC2_v , surface_FC_lar );
    new G4LogicalBorderSurface("LiquidArgon-->FC3v", physicalWorld, physicalFC3_v , surface_FC_lar );
    new G4LogicalBorderSurface("LiquidArgon-->FC4v", physicalWorld, physicalFC4_v , surface_FC_lar );
    new G4LogicalBorderSurface("LiquidArgon-->FC1h", physicalWorld, physicalFC1_h , surface_FC_lar );
    new G4LogicalBorderSurface("LiquidArgon-->FC2h", physicalWorld, physicalFC2_h , surface_FC_lar );
    new G4LogicalBorderSurface("LiquidArgon-->FC3h", physicalWorld, physicalFC3_h , surface_FC_lar );
    new G4LogicalBorderSurface("LiquidArgon-->FC4h", physicalWorld, physicalFC4_h , surface_FC_lar );

    //Inserting the Acrylic --------------------------------------
    auto config_Acrylic = config["Acrylic"];

    auto density = config_Acrylic["density"].get<double>()*g/cm3;
    std::vector<G4int>    natoms_acrylic = config_Acrylic["n_elements"].get<std::vector<G4int>>();
    std::vector<G4String> elements_acrylic = config_Acrylic["elements"].get<std::vector<G4String>>();

    auto mpt_Acry = new G4MaterialPropertiesTable();
    std::string file_acrylic_index = config_Acrylic["refraction_index"].get<string>();
    std::ifstream f_acrylic_index(std::string("../configuration/")+file_acrylic_index);
    
    json data_acry = json::parse(f_acrylic_index);
    size_t n_rindex_acry = data_acry.size();
    std::vector<G4double> energies_r_acry(n_rindex_acry);
    std::vector<G4double> rindex_acry(n_rindex_acry);
   
    for (size_t i = 0; i < n_rindex_acry; i++) 
    {
        energies_r_acry[i] = data_acry[i]["E"].get<double>()*eV;  
        rindex_acry[i]   = data_acry[i]["r"].get<double>();
    }
    std::string file_acrylic_length = config_Acrylic["abs_length"].get<string>();
    std::ifstream f_acrylic_length(std::string("../configuration/")+file_acrylic_length);
    json data_acry_length = json::parse(f_acrylic_length);
    size_t n_abs_acry = 3;
    std::vector<G4double> absE_acry(n_abs_acry);
    std::vector<G4double> absLen_acry(n_abs_acry);
    for (size_t i=0; i<n_abs_acry; i++) 
    {
        absE_acry[i] = data_acry_length[i]["E"].get<double>()*eV;    
        absLen_acry[i] = data_acry_length[i]["l"].get<double>()*cm;   
    }
    G4Material* Acry_mat = new G4Material("My_Acrylic",density=density,n_abs_acry);
    for (size_t i=0; i<n_abs_acry; i++) 
    {
        G4Element* elem = fNistManager->FindOrBuildElement(elements_acrylic[i]);
        Acry_mat->AddElement(elem,natoms_acrylic[i]);
    }

    mpt_Acry->AddProperty("RINDEX", energies_r_acry.data(), rindex_acry.data(), n_rindex_acry);
    mpt_Acry->AddProperty("ABSLENGTH", absE_acry.data(), absLen_acry.data(), n_abs_acry); 
    Acry_mat->SetMaterialPropertiesTable(mpt_Acry);

    double acrylic_X = cryostat_sizeX-2*d_cryo-FCv_sizeX-2*cut_value1;
    double acrylic_Y = FCv_sizeY;
    double acrylic_Z = cryostat_sizeZ-2*d_cryo-FCv_sizeZ-2*cut_value1;
    double acry_thickness = config_Acrylic["thickness"].get<double>()*cm;
    auto acrylic_wall_long = new G4Box("Acrylic_long", 0.5*acrylic_X, 0.5*acrylic_Y, 0.5*acry_thickness);
    auto acrylic_wall_end = new G4Box("Acrylic_end",  0.5*acry_thickness, 0.5*acrylic_Y, 0.5*acrylic_Z);
    auto logical_Acrylic_wall_long = new G4LogicalVolume(acrylic_wall_long, Acry_mat , "Acrylic_long");
    auto logical_Acrylic_wall_end = new G4LogicalVolume(acrylic_wall_end, Acry_mat , "Acrylic_end");

    G4VPhysicalVolume* physicalAcrylical1 = new G4PVPlacement(0, G4ThreeVector(0,0,cryostat_sizeZ/2-d_cryo-0.5*FCv_sizeZ-cut_value1-0.5*acry_thickness), 
        logical_Acrylic_wall_long, "Acrylical_long", logicWorld, true, 1, checkOverlaps);
    G4VPhysicalVolume* physicalAcrylical2 = new G4PVPlacement(0, G4ThreeVector(0,0,-(cryostat_sizeZ/2-d_cryo-0.5*FCv_sizeZ-cut_value1-0.5*acry_thickness)), 
        logical_Acrylic_wall_long, "Acrylical_long", logicWorld, true, 2, checkOverlaps);
    G4VPhysicalVolume* physicalAcrylical3 = new G4PVPlacement(0, G4ThreeVector(cryostat_sizeX/2-d_cryo-0.5*FCv_sizeX-cut_value1-0.5*acry_thickness,0,0), 
        logical_Acrylic_wall_end, "Acrylical_end", logicWorld, true, 1, checkOverlaps);
    G4VPhysicalVolume* physicalAcrylical4 = new G4PVPlacement(0, G4ThreeVector(-(cryostat_sizeX/2-d_cryo-0.5*FCv_sizeX-cut_value1-0.5*acry_thickness),0,0), 
        logical_Acrylic_wall_end, "Acrylical_end", logicWorld, true, 2, checkOverlaps);
    
    G4VisAttributes* visAcry = new G4VisAttributes(G4Colour(0.8, 0.1, 0.1)); 
    visAcry->SetForceSolid(true); 
    logical_Acrylic_wall_long->SetVisAttributes(visAcry);
    logical_Acrylic_wall_end->SetVisAttributes(visAcry);
    // interface between acrylic, metal and argon

    G4OpticalSurface* surface_acrylic_lar = new G4OpticalSurface("surface_acry_lar");
    surface_acrylic_lar->SetModel(unified);
    surface_acrylic_lar->SetType(dielectric_dielectric);
    surface_acrylic_lar->SetFinish(polished);
    G4MaterialPropertiesTable* mpt_AcryLar_Surface = new G4MaterialPropertiesTable();
    mpt_AcryLar_Surface->AddProperty("REFLECTIVITY",{0,15}, {1,1},2);
    surface_acrylic_lar->SetMaterialPropertiesTable(mpt_AcryLar_Surface);

    G4OpticalSurface* surface_acrylic_cryo = new G4OpticalSurface("surface_acry_metal");
    surface_acrylic_cryo-> SetModel(unified);
    surface_acrylic_cryo-> SetType(dielectric_metal);
    surface_acrylic_cryo-> SetFinish(polished);
    G4MaterialPropertiesTable* mpt_AcryMetal_Surface = new G4MaterialPropertiesTable();
    mpt_AcryMetal_Surface->AddProperty("REFLECTIVITY",{1,1}, {config_LAr["reflectivity_cryo"].get<double>(),config_LAr["reflectivity_cryo"].get<double>()},2);
    surface_acrylic_cryo->SetMaterialPropertiesTable(mpt_AcryMetal_Surface);

    new G4LogicalBorderSurface("LiquidArgon-->Cryo", physicalWorld, physicalCryostat , surface_cryo_lar );


    new G4LogicalBorderSurface("Acrylic1-->Argon", physicalAcrylical1, physicalWorld , surface_acrylic_lar );
    new G4LogicalBorderSurface("Acrylic2-->Argon", physicalAcrylical2, physicalWorld , surface_acrylic_lar );
    new G4LogicalBorderSurface("Argon-->Acrylic1", physicalWorld ,  physicalAcrylical1,surface_acrylic_lar );
    new G4LogicalBorderSurface("Argon-->Acrylic2", physicalWorld , physicalAcrylical2, surface_acrylic_lar ); 
    new G4LogicalBorderSurface("Acrylic3-->Argon", physicalAcrylical3, physicalWorld , surface_acrylic_lar );
    new G4LogicalBorderSurface("Acrylic4-->Argon", physicalAcrylical4, physicalWorld , surface_acrylic_lar );
    new G4LogicalBorderSurface("Argon-->Acrylic3", physicalWorld ,  physicalAcrylical3,surface_acrylic_lar );
    new G4LogicalBorderSurface("Argon-->Acrylic4", physicalWorld , physicalAcrylical4, surface_acrylic_lar ); 
   
    //Inserting PEN

    auto config_PEN = config["PEN"];
    auto density_PEN = config_PEN["density"].get<double>()*g/cm3;
    std::vector<G4int>  natoms_PEN = config_PEN["n_elements"].get<std::vector<G4int>>();
    std::vector<G4String> elements_PEN = config_PEN["elements"].get<std::vector<G4String>>();
    int n_abs_pen = 3;
    G4Material* PEN_mat = new G4Material("My_PEN",density=density_PEN,n_abs_pen);
    for (size_t i=0; i<n_abs_pen; i++) 
    {
        G4Element* elem = fNistManager->FindOrBuildElement(elements_PEN[i]);
        PEN_mat->AddElement(elem,natoms_PEN[i]);
    }
    
    double pen_X = acrylic_X-2*acry_thickness;
    double pen_Y = acrylic_Y;
    double pen_Z = acrylic_Z-2*acry_thickness;
    double pen_thickness = config_PEN["thickness"].get<double>()*cm;
    auto pen_wall_long = new G4Box("PEN_long", 0.5*pen_X, 0.5*pen_Y, 0.5*pen_thickness);
    auto pen_wall_end = new G4Box("PEN_end",  0.5*pen_thickness, 0.5*pen_Y, 0.5*pen_Z);

    auto logical_pen_wall_long = new G4LogicalVolume(pen_wall_long, PEN_mat , "PEN_long");
    auto logical_pen_wall_end = new G4LogicalVolume(pen_wall_end, PEN_mat , "PEN_end");

    G4VPhysicalVolume* physicalPEN1 = new G4PVPlacement(0, G4ThreeVector(0,0,cryostat_sizeZ/2-d_cryo-0.5*FCv_sizeZ-cut_value1-acry_thickness-0.5*pen_thickness), 
        logical_pen_wall_long, "PEN_long", logicWorld, true, 1, checkOverlaps);
    G4VPhysicalVolume* physicalPEN2 = new G4PVPlacement(0, G4ThreeVector(0,0,-(cryostat_sizeZ/2-d_cryo-0.5*FCv_sizeZ-cut_value1-acry_thickness-0.5*pen_thickness)), 
        logical_pen_wall_long, "PEN_long", logicWorld, true, 2, checkOverlaps);
    G4VPhysicalVolume* physicalPEN3 = new G4PVPlacement(0, G4ThreeVector(cryostat_sizeX/2-d_cryo-0.5*FCv_sizeX-cut_value1-acry_thickness-0.5*pen_thickness,0,0), 
        logical_pen_wall_end, "PEN_end", logicWorld, true, 1, checkOverlaps);
    G4VPhysicalVolume* physicalPEN4 = new G4PVPlacement(0, G4ThreeVector(-(cryostat_sizeX/2-d_cryo-0.5*FCv_sizeX-cut_value1-acry_thickness-0.5*pen_thickness),0,0), 
        logical_pen_wall_end, "PEN_end", logicWorld, true, 2, checkOverlaps);
    
    G4VisAttributes* visPEN = new G4VisAttributes(G4Colour(1, 1, 0.1)); 
    visPEN->SetForceSolid(true); 
    logical_pen_wall_long->SetVisAttributes(visPEN);
    logical_pen_wall_end->SetVisAttributes(visPEN);


    // ---- Inserting Cathode grid -------- 

    /* auto config_Cathode = config["Cathode_Grid"];

    double cathode_X = cryostat_sizeX-2*d_cryo-FCv_sizeX-cut_value1;
    double cathode_Y = cryostat_sizeZ-2*d_cryo-FCv_sizeZ-cut_value1;
    double cathode_hole_X = config_Cathode["hole"].get<double>()*cm;
    double cathode_hole_Y = config_Cathode["hole"].get<double>()*cm;
    double cathode_separation = config_Cathode["separation"].get<double>()*cm;


    auto cathode = new G4Box("Cathode", 0.5*cathode_X, 0.5*cryostatThickness, 0.5*cathode_Y);
    auto logicalCathode = new G4LogicalVolume(cathode, FC_mat, "Cathode");
    auto cathodehole = new G4Box("Cathode_Hole", 0.5*cathode_hole_X, 1.1*0.5*cryostatThickness, 0.5*cathode_hole_Y);
    auto logicalHole = new G4LogicalVolume(cathodehole, lar_mat, "Cathode_Hole");
    int n_cathode_x = (cathode_X+cathode_separation)/(cathode_hole_X+cathode_separation);
    int n_cathode_y = (cathode_Y+cathode_separation)/(cathode_hole_Y+cathode_separation);

    G4VPhysicalVolume* physicalCathode = new G4PVPlacement(0, G4ThreeVector(0,0,0), logicalCathode, "Cathode", logicCryostatFilling, false, 0);

    
    for (int ix = 0; ix < n_cathode_x; ix++)
    {
        for (int iz = 0; iz < n_cathode_y; iz++) 
        {
            double x_pos = (ix - (n_cathode_x-1)/2.0)*(cathode_hole_X + cathode_separation);
            double z_pos = (iz - (n_cathode_y-1)/2.0)*(cathode_hole_Y + cathode_separation);
            G4ThreeVector pos(x_pos, 0, z_pos);

            auto physicalHole = new G4PVPlacement(0, pos, logicalHole,"Cathode_Hole_PV", logicalCathode, false, ix*n_cathode_y + iz);

            new G4LogicalBorderSurface( "LAr_in_Hole-->Cathode", physicalHole, physicalCathode, surface_FC_lar );
        }
    }
    new G4LogicalBorderSurface("LiquidArgon-->Cathode", physicalCryostatFilling, physicalCathode , surface_FC_lar );

    G4VisAttributes* visCathode = new G4VisAttributes(G4Colour(0.9,0.9,0.9));
    visCathode->SetForceSolid(true);
    logicalCathode->SetVisAttributes(visCathode); */


    /* auto conjwindows = new G4Box("Window", 0.5*cryolength,0.5*(cryoheight-4*lduw),0.5*winthick);
    auto conjwindowsfb = new G4Box("Window Front/Bottom", 0.5*winthick,0.5*(cryoheight-4*lduw),0.5*cryowidth-fcwdist-fchbarw-penthick-winthick-2*profheight);
    auto logwindow =    new G4LogicalVolume(conjwindows,  Acrylic,"Logical Window");
    auto logwindowfb =  new G4LogicalVolume(conjwindowsfb,Acrylic,"Logical Window FB");
    logwindow->SetVisAttributes( new G4VisAttributes( G4Color::Blue() ) );
    G4VisAttributes * WinVisAtt = new G4VisAttributes(G4Color::Blue());
    WinVisAtt->SetForceSolid(true);
    //logwindowfb->SetVisAttributes( new G4VisAttributes( G4Color::Blue()) );
    logwindowfb->SetVisAttributes(WinVisAtt);
    G4VPhysicalVolume* phywindow1 = new G4PVPlacement(0,{0,0,cryowidth/2-fcwdist-fcvbarw-2*profheight},logwindow, "World/Window", logicalWorld,false,0);
    G4VPhysicalVolume* phywindow2 = new G4PVPlacement(0,{0,0,-cryowidth/2+fcwdist+fcvbarw+2*profheight},logwindow, "World/Window", logicalWorld,false,0);
    G4VPhysicalVolume* phywindowf = new G4PVPlacement(0,{0.5*(cryolength-4*lduw)-fchbarw-2*profheight,0,0},logwindowfb, "World/Window Front", logicalWorld,false,0);
    G4VPhysicalVolume* phywindowb = new G4PVPlacement(0,{-0.5*(cryolength-4*lduw)+fchbarw+2*profheight,0,0},logwindowfb, "World/Window Bottom", logicalWorld,false,0);
 */

    return physicalWorld;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
