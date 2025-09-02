#ifndef OptFileManager_h
#define OptFileManager_h 1

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

class OptFileManager {
  public:
    OptFileManager();
    virtual ~OptFileManager();

    void GetSpectrumFromData (string fin,vector<G4double>* v1,vector<G4double>* v2);
    void readPropertiesFromData(string fin, G4double* abslength, G4double* rindex, G4double* yield, G4double* resolution, G4double* birksconstant, G4double* scint_decaytime ,string* spectrum_data);
    G4double encodingUnit(double num, std::string unit,std::string type);
  private:
};

#endif
