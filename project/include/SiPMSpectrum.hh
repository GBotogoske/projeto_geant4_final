// Settings.h
#pragma once
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "G4RunManager.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include <nlohmann/json.hpp>

using json = nlohmann::json;

class SiPMSpectrum 
{
    private:
        SiPMSpectrum()
        { 
            loadFromFile("../configuration/detector.json"); 
        }

        std::vector<G4double> eff_uv;
        std::vector<G4double> E_uv;
        int n_uv=0;

        std::vector<G4double> eff_vis;
        std::vector<G4double> E_vis;
        int n_vis=0;

        void loadFromFile(const std::string& filename) 
        {
            std::ifstream file(filename);
            if (!file) return; 
               
            json config;

            file >> config;
            auto config_SiPM = config["SiPM"];

            std::string file_uv = config_SiPM["efficiency_UV"].get<std::string>();
            std::ifstream file_uv_path(std::string("../configuration/")+file_uv);
            json data_uv = json::parse(file_uv_path);
            this->n_uv = data_uv.size();
            this->E_uv.resize(this->n_uv);
            this->eff_uv.resize(this->n_uv);
            for (size_t i = 0; i < this->n_uv; i++) 
            {
                this->E_uv[i] = data_uv[i]["E"].get<double>()*eV;  
                this->eff_uv[i] = data_uv[i]["p"].get<double>();
            }

            std::string file_vis = config_SiPM["efficiency_VIS"].get<std::string>();
            std::ifstream file_vis_path(std::string("../configuration/")+file_vis);
            json data_vis = json::parse(file_vis_path);
            this->n_vis = data_vis.size();
            this->E_vis.resize(this->n_vis);
            this->eff_vis.resize(this->n_vis);
            for (size_t i = 0; i < this->n_vis; i++) 
            {
                this->E_vis[i] = data_vis[i]["E"].get<double>()*eV;  
                this->eff_vis[i] = data_vis[i]["p"].get<double>();
            }
           
        }


    public:
        static SiPMSpectrum& get() 
        {
            static SiPMSpectrum instance; 
            return instance;
        }

        const std::vector<G4double> get_effUV() const
        {
            return this->eff_uv;
        }   
        const std::vector<G4double> get_EUV() const
        {
            return this->E_uv;
        }  
        int getNUV() const 
        { 
            return this->n_uv; 
        }

        const std::vector<G4double> get_effVIS() const
        {
            return this->eff_vis;
        }   
        const std::vector<G4double> get_EVIS() const
        {
            return this->E_vis;
        }  
        int getNVIS() const 
        { 
            return this->n_vis; 
        }

        SiPMSpectrum(const SiPMSpectrum&) = delete;
        SiPMSpectrum& operator=(const SiPMSpectrum&) = delete;
};
