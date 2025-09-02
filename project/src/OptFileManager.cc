#include "OptFileManager.hh"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;
OptFileManager::OptFileManager(){}
OptFileManager::~OptFileManager(){}

//Others FUnctions
void OptFileManager::GetSpectrumFromData (std::string path,std::vector<G4double>* v1,std::vector<G4double>* v2){
  std::ifstream datafile;
  std::string line,line2;
  datafile.open(path);
  int c=0;
  while(getline (datafile,line,'\n')){
    G4String test = (G4String)line;
    if (!G4StrUtil::contains(test,"#")){
    std::stringstream input_line(line);
    while(getline(input_line,line2,':')){
      if(line2!=""){
	  if(c%2==0){
	    G4double conv=std::stod(line2)*eV;
	    v1->push_back(conv);
	  }
	  else {
	    G4double conv=std::stod(line2);
	    v2->push_back(conv);
	  }
      }
	c++;
    }
    }
  }
  datafile.close();};

void OptFileManager::readPropertiesFromData(string fin, G4double* abslength, G4double* rindex, G4double* yield, G4double* resolution, G4double* birksconstant, G4double* scint_decaytime ,string* spectrum_data){
  std::ifstream datafile;
  std::string line,line1,line2,line3,line4;
  datafile.open(fin);
  if (datafile.fail()){
    G4cout<<"[ERROR] "<<fin<<" Does not exist (default values are setted)! "<<G4endl;
    * abslength       = 1;
    * rindex          = 1;
    * yield           = 0.000000001;
    * resolution      = 1;
    * birksconstant   = 1;
    * scint_decaytime = 1;
    * spectrum_data   = "None";
  }
  else{
    int k0,k1,k2;
    k0 = 0;
    while(getline(datafile,line,'\n')){
      //G4cout<<"\t [DEBUG]: process line: "<<line<<G4endl;
      k0 ++;
      //G4cout<<"[DEBUG0]: k0 = "<<k0<<" "<<line<<G4endl;
      std::stringstream input_line(line);
      k1 = 0;
      if(k0>1){
	std::stringstream input_line(line);
	while(getline(input_line,line1,':')){
	  k1++;
	  //G4cout<<"[DEBUG1]: k1 = "<<k1<<" "<<line1<<G4endl;
	  std::stringstream input_line1(line1);
	  k2 = 0;
	  if(k1>1){
	    while(getline(input_line1,line2,' ')){
	      k2++;
	      if(k2==1){
		if(k0==2)*abslength=std::stod(line2);
		if(k0==3)*rindex=std::stod(line2);
		if(k0==4)*yield=std::stod(line2);
		if(k0==5)*resolution=std::stod(line2);
		if(k0==6)*birksconstant=std::stod(line2);
		if(k0==7)*scint_decaytime=std::stod(line2);
		if(k0==8)*spectrum_data=line2;
	      }
	      if(k2==2){
		//G4cout<<" HERE"<<G4endl;
		if(k0==2)*abslength = encodingUnit(*abslength,line2,"Length");
		if(k0==7)*scint_decaytime = encodingUnit(*scint_decaytime,line2,"Time");
	      }
	      //G4cout<<"[DEBUG2]: k2 = "<<k2<<" "<<line2<<G4endl;
	    } 
	  }
	}
      }
      //G4cout<<"[DEBUG_Final]: "<<line<<" "<<k0<<" "<<k1<<" "<<k2<<G4endl;
    }
    G4cout<<"[INFO]: Reading Optical Properties from "<< fin<<
      "\n\t * Absorption Length = "<<G4BestUnit(*abslength,"Length")<<
      "\n\t * Refraction Index = "<<*rindex<<
      "\n\t * Scintillation yeld = "<<*yield<<
      "\n\t * Resolution Scale = "<<*resolution<<
      "\n\t * Birk constant (mm/MeV) = "<<*birksconstant<<
      "\n\t * Decay Time constant = "<<G4BestUnit(*scint_decaytime,"Time")<<
      "\n\t * Emission spectrum = "<<*spectrum_data<<G4endl;
  }
};

G4double OptFileManager::encodingUnit(double num, std::string unit,std::string type){
  double conv=0;
  if(type=="Length"){
    //* length units 
    if (unit=="m") conv=0;
    else if (unit=="dm") conv=-1;
    else if (unit=="cm") conv=-2;
    else if (unit=="mm") conv=-3;
    else if (unit=="um") conv=-6;
    else if (unit=="nm") conv=-9;
    else if (unit=="pm") conv=-12;
    else if (unit=="fm") conv=-15;
    else if (unit=="am") conv=-18;
    else if (unit=="zm") conv=-21;
    else if (unit=="ym") conv=-24;
    else if (unit=="km") conv=3;
    else if (unit=="Mm") conv=6;
    else if (unit=="Gm") conv=9;
    else if (unit=="Tm") conv=12;
    else if (unit=="Pm") conv=15;
    else if (unit=="Em") conv=18;
    else if (unit=="Zm") conv=21;
    else if (unit=="Ym") conv=24;
    else{
      conv=0;
      G4cout<<"[ERROR]: Unit not recognized: default meter was setted"<<G4endl;
    }
    return num*(pow(10,conv))*m;
  }
  //* Time unit
  else if(type=="Time"){
    if (unit=="s") conv=0;
    else if (unit=="ds") conv=-1;
    else if (unit=="cs") conv=-2;
    else if (unit=="ms") conv=-3;
    else if (unit=="us") conv=-6;
    else if (unit=="ns") conv=-9;
    else if (unit=="ps") conv=-12;
    else if (unit=="fs") conv=-15;
    else if (unit=="as") conv=-18;
    else if (unit=="zs") conv=-21;
    else if (unit=="ys") conv=-24;
    else if (unit=="ks") conv=3;
    else if (unit=="Ms") conv=6;
    else if (unit=="Gs") conv=9;
    else if (unit=="Ts") conv=12;
    else if (unit=="Ps") conv=15;
    else if (unit=="Es") conv=18;
    else if (unit=="Zs") conv=21;
    else if (unit=="Ys") conv=24;
    else{
      conv=0;
      G4cout<<"[ERROR]: Unit not recognized: default second was setted"<<G4endl;
    }
    return num*(pow(10,conv))*s;
  }
  else{
    G4cout<<"[ERROR]: Unit not recognized!"<<G4endl;
    return 1;
  }
};
