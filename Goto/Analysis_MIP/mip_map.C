#include <stdlib.h>
#include <fstream>
#include <math.h>
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TCut.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TH1.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLine.h"
#include "TLatex.h"

using namespace::std;

void mip_map(){

  Int_t MaxSlab = 5;
  Int_t MaxChip = 16;
  Int_t MaxChannel = 64;
  Int_t MaxMemory = 15;

  TString mapchip_path = "../map_chip/map_chip_fev13.dat";
  std::ifstream reading_file(mapchip_path, std::ios::in);
  if(!reading_file){
    cout << "map_chip.dat is not found" << endl;
  }
  Float_t map_pointX[16][64], map_pointY[16][64];
  Int_t htmp_chip = 0, htmp_channel = 0;
  Float_t tmp_x0 = 0, tmp_y0 = 0, tmp_x = 0, tmp_y = 0;
  while(reading_file){
    reading_file >> htmp_chip >> tmp_x0 >> tmp_y0 >> htmp_channel >> tmp_x >> tmp_y;
    map_pointX[htmp_chip][htmp_channel] = -tmp_x;
    map_pointY[htmp_chip][htmp_channel] = tmp_y;
  }

  // Fill to Tree Data
  Double_t mip_mean[5][16][64] = {0};
  Double_t mip_error[5][16][64] = {0};
  Double_t mip_chi2ndf[5][16][64] = {0};

  // Create New TTree
  TFile *fout = new TFile("MIP_Map.root" , "RECREATE");
  TTree *MIP_Tree = new TTree("MIP_Tree", "MIP_Tree");
  MIP_Tree->Branch("mip_mean", mip_mean,
		        "mip_mean[5][16][64]/D");
  MIP_Tree->Branch("mip_error", mip_error,
		        "mip_error[5][16][64]/D");
  MIP_Tree->Branch("mip_chi2ndf", mip_chi2ndf,
		        "mip_chi2ndf[5][16][64]/D");


  TString mip_path = "./Slab_Result/mip_mean.txt";
  std::ifstream reading_mip(mip_path, std::ios::in);
  if(!reading_mip){
    cout << "mip_mean.txt is not found" << endl;
  }
  Int_t tmp_slab = 0, tmp_chip = 0, tmp_channel = 0;
  Double_t tmp_mip_mean = 0, tmp_mip_error = 0, tmp_mip_chi = 0, tmp_mip_ndf = 0;
  while(reading_mip){
    reading_mip >> tmp_slab >> tmp_chip >> tmp_channel >> tmp_mip_mean >> tmp_mip_error >> tmp_mip_chi >> tmp_mip_ndf;
    if(tmp_mip_ndf==0) continue;
    if(mip_mean[tmp_slab][tmp_chip][tmp_channel]==0){
      mip_mean[tmp_slab][tmp_chip][tmp_channel] = tmp_mip_mean;
      mip_error[tmp_slab][tmp_chip][tmp_channel] = tmp_mip_error;
      mip_chi2ndf[tmp_slab][tmp_chip][tmp_channel] = tmp_mip_chi/tmp_mip_ndf;
      cout << tmp_mip_mean << " " << tmp_mip_chi/tmp_mip_ndf << endl;
    }
    else if(mip_mean[tmp_slab][tmp_chip][tmp_channel]!=0 && (tmp_mip_chi/tmp_mip_ndf)<mip_chi2ndf[tmp_slab][tmp_chip][tmp_channel]){
      mip_mean[tmp_slab][tmp_chip][tmp_channel] = tmp_mip_mean;
      mip_error[tmp_slab][tmp_chip][tmp_channel] = tmp_mip_error;
      mip_chi2ndf[tmp_slab][tmp_chip][tmp_channel] = tmp_mip_chi/tmp_mip_ndf;
    }
  }

  fout->cd();
  MIP_Tree->Fill();
  MIP_Tree->Write();

  // Make Mip Map
  Int_t number_of_nodata[5] = {0};
  Int_t number_of_data[5] = {0};
  Int_t number_of_channels[5] = {0};
  TH2F *mappp[MaxSlab];
  TH2F *chip_map = new TH2F("chip_map", "Chip Map;X[mm];Y[mm]", 32, -90, 90, 32, -90, 90);
  TH2F *asic_map = new TH2F("asic_map", "Asic Map;X[mm];Y[mm]", 4, -90, 90, 4, -90, 90);
  for(Int_t islab=0; islab<MaxSlab; islab++){
    TString mapname = TString::Format("Hit Map Slab%i", islab);
    mappp[islab] = new TH2F(mapname, mapname + ";X[mm];Y[mm]", 32, -90, 90, 32, -90, 90);
    for(Int_t ichip=0; ichip<MaxChip; ichip++){
      if(islab==0) asic_map->Fill(map_pointX[ichip][0], map_pointY[ichip][0], ichip);
      for(Int_t ichn=0; ichn<MaxChannel; ichn++){
	number_of_channels[islab]++;
	if(mip_mean[islab][ichip][ichn]==0){
	  cout << "slab " << islab << " chip " << ichip << " channel " << ichn << " mean " << tmp_mip_mean << endl;
	  number_of_nodata[islab]++;
	}
	else if(mip_chi2ndf[islab][ichip][ichn]<3){
          mappp[islab]->Fill(map_pointX[ichip][ichn], map_pointY[ichip][ichn], mip_mean[islab][ichip][ichn]);
	  number_of_data[islab]++;
	}
        if(islab==0) chip_map->Fill(map_pointX[ichip][ichn], map_pointY[ichip][ichn], ichn);
      }
    }
  }
  for(Int_t islab=0; islab<MaxSlab; islab++){
    cout << "in slab " << islab << " : " << number_of_channels[islab] << " channels exist" << endl;
    cout << "in slab " << islab << " : " << number_of_nodata[islab] << " channels are nothing" << endl;
    cout << "in slab " << islab << " : " << number_of_data[islab] << " channels can be used " << endl;
  }

  // Print Hit Map
  TCanvas *canvas1 = new TCanvas("canvas1", "canvas1");
  canvas1->Divide(3, 2);
  for(Int_t islab=0; islab<MaxSlab; islab++){
    canvas1->cd(islab+1);
    gStyle->SetOptStat(0);
    mappp[islab]->GetZaxis()->SetRangeUser(0, 200);
    mappp[islab]->Draw("colz");
    chip_map->Draw("text same");
    asic_map->Draw("text same");

    TLine *border[2][3];
    for(Int_t i=0; i<2; i++){
      for(Int_t j=0; j<3; j++){
        if(i==0)
	  border[i][j] = new TLine(-45+45*j,-90,-45+45*j,90);
        if(i==1)
	  border[i][j] = new TLine(-90,-45+45*j,90,-45+45*j);
        border[i][j]->SetLineColor(kGray);
        border[i][j]->SetLineWidth(1);
        border[i][j]->SetLineStyle(2);
        border[i][j]->Draw("same");
      }
    }
    if(1){
      TLatex latex;
      latex.SetTextAlign(22);  //centered
      latex.SetTextColor(17);
      for(Int_t i=0; i<MaxChip; i++){
        if(i<8)
	  latex.DrawLatex(67.5-45*((Int_t)i/2),-22.5-45*((i+1)%2),Form("#bf{#scale[2]{%d}}",i));
        else
	  latex.DrawLatex(67.5-45*((Int_t)(i-8)/2),67.5-45*((i+1)%2),Form("#bf{#scale[2]{%d}}",i));
      }
    }
  }
  canvas1->Print("Mip_Map.pdf", "pdf");
}
