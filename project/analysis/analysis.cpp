
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>

void analysis()
{
    // Open file and get tree
    TFile *file = TFile::Open("../build/Data.root");
    TTree *tree = (TTree*)file->Get("Primary_Hit");

    // Create 2D histogram
    TH2D *h2 = new TH2D("hvis","VIS(z,y)", 26, -6.5, 6.5, 26, -6.5, 6.5);
    TH2D *h3 = new TH2D("huv","VUV(z,y)", 26, -6.5, 6.5, 26, -6.5, 6.5);
    TH2D *h4 = new TH2D("hall","VUV+VIS(z,y)", 26, -6.5, 6.5, 26, -6.5, 6.5);

    // Fill the histogram
    Double_t z, y, f, f2;
    tree->SetBranchAddress("Z", &z);
    tree->SetBranchAddress("Y", &y);
    tree->SetBranchAddress("PhotonDetectedVIS", &f);
    tree->SetBranchAddress("PhotonDetectedUV", &f2);

    float gain=40000.0/25000.0;
    float ratio=1.0;

    Long64_t nentries = tree->GetEntries();
    for (Long64_t i=0; i<nentries; i++) 
    {
        tree->GetEntry(i);

        if(z<(-6) || z>6)
        {
            ratio=gain;
        }
        else
        {
            ratio=1.0;
        }

        if(f>=1)
        {
            h2->Fill(z, y,f*ratio);
        }
            
        if(f2>=1) 
        {
            h3->Fill(z,y,f2*ratio);
        } // fill with weight = f
        
        h4->Fill(z,y,(f+f2)*ratio);
    }

    // Draw heatmap
    TCanvas *c1 = new TCanvas();
    h2->Draw("COLZ");  // COLZ = color map
    c1->SaveAs("heatmapVIS.png");

    TCanvas *c2 = new TCanvas();
    h3->Draw("COLZ");  // COLZ = color map
    c2->SaveAs("heatmapUV.png");

    TCanvas *c5 = new TCanvas();
    h4->Draw("COLZ");  // COLZ = color map
    c5->SaveAs("heatmapAll.png");

        // --- Converter h2 para matriz ---
    int nx = h2->GetNbinsX();
    int ny = h2->GetNbinsY();
    std::vector<std::vector<double>> mat_uv(nx, std::vector<double>(ny,0));
    std::vector<std::vector<double>> mat_vis(nx, std::vector<double>(ny,0));
    std::vector<std::vector<double>> mat_all(nx, std::vector<double>(ny,0));
    std::vector<std::vector<double>> buffer(nx, std::vector<double>(ny,0));
    std::vector<std::vector<double>> myratio(nx, std::vector<double>(ny,0));

    double lightyield_avg=0;
    double lightyield_min=10000;

    for(int i=1; i<=nx; i++)        // ROOT bins começam em 1
    {
        for(int j=1; j<=ny; j++)
        {
            mat_vis[i-1][j-1] = h2->GetBinContent(i,j);
            mat_uv[i-1][j-1] = h3->GetBinContent(i,j);
            mat_all[i-1][j-1] = h4->GetBinContent(i,j);
            buffer[i-1][j-1] =   mat_uv[i-1][j-1]/(  mat_uv[i-1][j-1] + mat_vis[i-1][j-1] );
            myratio[i-1][j-1] =   mat_uv[i-1][j-1]/(mat_vis[i-1][j-1] );
        }
    }

    int cont=0;
    for(int i=0; i<nx; i++)        
    {
        for(int j=1; j<ny-1; j++)
        {
            lightyield_avg+=mat_all[i][j];
            cont++;
        }
    }
    lightyield_avg=lightyield_avg/cont;
    std::cout << "AVG Light Yield: " << lightyield_avg << std::endl;


    /* std::cout << "----- VIS ------" << std::endl;
    // --- Printar na tela ---
    for(int j=ny-1; j>=0; j--)  // imprimir Y crescente de cima para baixo
    {
        for(int i=0; i<nx; i++)
        {
            std::cout << mat_vis[i][j] << " ";
        }
        std::cout << std::endl;
    }

     std::cout << "----- UV ------" << std::endl;
       // --- Printar na tela ---
    for(int j=ny-1; j>=0; j--)  // imprimir Y crescente de cima para baixo
    {
        for(int i=0; i<nx; i++)
        {
            std::cout << mat_uv[i][j] << " ";
        }
        std::cout << std::endl;
    } */

    // Criar histograma para o buffer
    TH2D *hbuffer = new TH2D("hbuffer","UV/(UV+VIS)", nx, -6.5, 6.5, ny, -6.5, 6.5);
    TH2D *hratio = new TH2D("hratio","UV/VIS", nx, -6.5, 6.5, ny, -6.5, 6.5);


    // Preencher com os valores da matriz buffer
    for (int i = 1; i <= nx; i++) {
        for (int j = 1; j <= ny; j++) {
            hbuffer->SetBinContent(i, j, buffer[i-1][j-1]);
            hratio->SetBinContent(i, j, myratio[i-1][j-1]);
        }
    }

    // Plotar
    TCanvas *c3 = new TCanvas();
    hbuffer->GetZaxis()->SetRangeUser(0.0,1.0); // já que é uma fração
    hbuffer->Draw("COLZ");
    c3->SaveAs("heatmapBuffer.png");

    TCanvas *c4 = new TCanvas();
    hratio->GetZaxis()->SetRangeUser(0.0,1.0); // já que é uma fração
    hratio->Draw("COLZ");
    c4->SaveAs("heatmapRatio.png");


}

