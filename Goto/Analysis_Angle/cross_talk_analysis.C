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
#include "TH3F.h"
#include "TSpectrum.h"
#include "TPie.h"
#include "TLatex.h"
#include "TLine.h"

using namespace::std;

void cross_talk_analysis(std::string str){
  // Cut File name
  int path_i = str.find_last_of("/")+1;
  int ext_i = str.find_last_of(".");
  std::string pathname = str.substr(0, path_i);
  std::string filename_ = str.substr(path_i, ext_i-path_i);
  //std::string extroot = filename_.substr(ext_i, str.size()-ext_i);

  TString fullpath = str.data();
  TString filename = pathname + filename_;
  TString resultname = "Cross_Talk/" + filename_;



  // ================================================================================================================================================ //
  // TTree open  
  TFile *file = TFile::Open(fullpath); 
  TTree *Event_Tree = (TTree*)file->Get("Event_Tree");

  // ================================================================================================================================================ //
  // Set Loop Max
  Int_t MaxEvent = Event_Tree->GetEntries();
  //Int_t NOver5Hits = Event_Tree->GetEntries("nhit_channel>5");
  Int_t MaxSlab = 5;
  Int_t MaxChip = 16;
  Int_t MaxChannel = 64;
  Int_t MaxMemory = 15;
  Int_t NMAX = MaxSlab*MaxChip*MaxChannel*MaxSlab;

  // ================================================================================================================================================ //
  // Reading ADC data from Tree
  Int_t acqNumber;
  Event_Tree->SetBranchAddress("acqNumber", &acqNumber);
  Int_t nhit_channel;
  Event_Tree->SetBranchAddress("nhit_channel", &nhit_channel);
  Int_t bcid_hit;
  Event_Tree->SetBranchAddress("bcid_hit", &bcid_hit);
  Int_t bcid_original[NMAX];
  Event_Tree->SetBranchAddress("bcid_original", bcid_original);
  Int_t slab_number[NMAX];
  Event_Tree->SetBranchAddress("slab_number", slab_number);
  Int_t chip_number[NMAX];
  Event_Tree->SetBranchAddress("chip_number", chip_number);
  Int_t sca_number[NMAX];
  Event_Tree->SetBranchAddress("sca_number", sca_number);
  Int_t channel_number[NMAX];
  Event_Tree->SetBranchAddress("channel_number", channel_number);
  Float_t hit_x[NMAX];
  Event_Tree->SetBranchAddress("hit_x", hit_x);
  Float_t hit_y[NMAX];
  Event_Tree->SetBranchAddress("hit_y", hit_y);
  Float_t hit_z[NMAX];
  Event_Tree->SetBranchAddress("hit_z", hit_z);
  Int_t charge_hiGain_tdc[NMAX];
  Event_Tree->SetBranchAddress("charge_hiGain_tdc", charge_hiGain_tdc);
  Int_t charge_lowGain_adc[NMAX];
  Event_Tree->SetBranchAddress("charge_lowGain_adc", charge_lowGain_adc);

  // Read Chip Map
  TString mapchip_path = "./map_chip/map_chip_fev13.dat";
  std::ifstream reading_file(mapchip_path, std::ios::in);
  if(!reading_file) cout << "map_chip.dat is not found" << endl;
  Float_t map_pointX[MaxChip][MaxChannel], map_pointY[MaxChip][MaxChannel];
  Int_t tmp_chip = 0, tmp_channel = 0;
  Float_t tmp_x0 = 0, tmp_y0 = 0, tmp_x = 0, tmp_y = 0;
  while(reading_file){
    reading_file >> tmp_chip >> tmp_x0 >> tmp_y0 >> tmp_channel >> tmp_x >> tmp_y;
    map_pointX[tmp_chip][tmp_channel] = -tmp_x;
    map_pointY[tmp_chip][tmp_channel] = tmp_y;
  }

  // ================================================================================================================================================ //
  // Make Cross Talk Energy Histogram
  TH1D *total_E_hist[5][16];
  Int_t hit_number[5][16][64] = {0};
  for(Int_t islab=0; islab<MaxSlab; islab++){
    for(Int_t ichip=0; ichip<MaxChip; ichip++){
      TString histname = TString::Format("Total_Energy_Histogram_slab_%i_chip_%i", islab, ichip);
      total_E_hist[islab][ichip] = new TH1D(histname, histname, 2048, 0, 2048);
    }
  }
  Int_t Ntotal[5][16] = {0};
  Int_t Ncrosstalk[5][16] = {0};
  Int_t cutNcrosstalk[5][16] = {0};
  for(Int_t ievent=0; ievent<MaxEvent;ievent++){
    Event_Tree->GetEntry(ievent);
    if(float(ievent%2000)==0) cout << "Event" << ievent << "End" << endl;

    for(Int_t ihit=0; ihit<nhit_channel; ihit++){
      Ntotal[slab_number[ihit]][chip_number[ihit]]++;
      hit_number[slab_number[ihit]][chip_number[ihit]][channel_number[ihit]]++;
      for(Int_t jhit=0; jhit<ihit; jhit++){
        if(slab_number[ihit]!=slab_number[jhit] || bcid_original[ihit]!=bcid_original[jhit]) continue;
        if((std::pow(hit_x[ihit]-hit_x[jhit], 2.0)+std::pow(hit_y[ihit]-hit_y[jhit], 2.0))<64){ // adjoin channel
	  Int_t total_energy = charge_lowGain_adc[ihit] + charge_lowGain_adc[jhit];
	  total_E_hist[slab_number[ihit]][chip_number[ihit]]->Fill(total_energy);
	  Ncrosstalk[slab_number[ihit]][chip_number[ihit]]++;
	  if(total_energy>100 && total_energy<250) cutNcrosstalk[slab_number[ihit]][chip_number[ihit]]++;
	} 
      }
    }
  }
   
  // Make Hit Map and NHit
  TH2F *mappp[MaxSlab];
  TH2F *chip_map = new TH2F("chip_map", "Chip Map;X[mm];Y[mm]", 32, -90, 90, 32, -90, 90);
  TH2F *asic_map = new TH2F("asic_map", "Asic Map;X[mm];Y[mm]", 4, -90, 90, 4, -90, 90);
  for(Int_t islab=0; islab<MaxSlab; islab++){
    TString mapname = TString::Format("Hit Map Slab%i", islab);
    mappp[islab] = new TH2F(mapname, mapname + ";X[mm];Y[mm]", 32, -90, 90, 32, -90, 90);
    for(Int_t ichip=0; ichip<MaxChip; ichip++){
      if(islab==0) asic_map->Fill(map_pointX[ichip][0], map_pointY[ichip][0], ichip);
      for(Int_t ichn=0; ichn<MaxChannel; ichn++){
        mappp[islab]->Fill(map_pointX[ichip][ichn], map_pointY[ichip][ichn], hit_number[islab][ichip][ichn]);
        if(islab==0) chip_map->Fill(map_pointX[ichip][ichn], map_pointY[ichip][ichn], ichn);
      }
    }
  }


  // ================================================================================================================================================ //
  // Print pdf
  TCanvas *canvas1 = new TCanvas("canvas1", "canvas1");
  canvas1->Print(resultname + "_analysis.pdf[", "pdf");
  canvas1->Divide(4, 4);
  Int_t chip_position[16] = {15, 13, 11, 9, 14, 12, 10, 8, 7, 5, 3, 1, 6, 4, 2, 0};
  Int_t colors[2] = {1, 2};
  Double_t npie[5][16][2];

  // Print Hit Map
  TCanvas *canvas2 = new TCanvas("canvas2", "canvas2");
  canvas2->Divide(3, 2);
  for(Int_t islab=0; islab<MaxSlab; islab++){
    canvas2->cd(islab+1);
    gStyle->SetOptStat(0);
    mappp[islab]->Draw("colz");
    chip_map->Draw("text same");
    asic_map->Draw("text same");
    TLine *border[2][3];
    for(Int_t i=0; i<2; i++){ 
      for(Int_t j=0; j<3; j++){ 
        if(i==0) border[i][j] = new TLine(-45+45*j,-90,-45+45*j,90);
        if(i==1) border[i][j] = new TLine(-90,-45+45*j,90,-45+45*j);
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
        if(i<8) latex.DrawLatex(67.5-45*((Int_t)i/2),-22.5-45*((i+1)%2),Form("#bf{#scale[2]{%d}}",i));
        else latex.DrawLatex(67.5-45*((Int_t)(i-8)/2),67.5-45*((i+1)%2),Form("#bf{#scale[2]{%d}}",i));
      }
    }
  }
  canvas2->Print(resultname + "_analysis.pdf", "pdf");

  for(Int_t islab=0; islab<MaxSlab; islab++){
    for(Int_t ichip=0; ichip<MaxChip; ichip++){
      npie[islab][ichip][0] = double(Ntotal[islab][ichip]-Ncrosstalk[islab][ichip]);
      npie[islab][ichip][1] = double(Ncrosstalk[islab][ichip]);
    }
  }
  for(Int_t islab=0; islab<MaxSlab; islab++){
    for(Int_t ichip=0; ichip<MaxChip; ichip++){
      canvas1->cd(ichip+1);
      gStyle->SetOptStat(1111);
      total_E_hist[islab][chip_position[ichip]]->Draw();
    }
    canvas1->Print(resultname + "_analysis.pdf", "pdf");
    TPie *pie[16];
    for(Int_t ichip=0; ichip<MaxChip; ichip++){
      canvas1->cd(ichip+1);
      TString piename = TString::Format("Cross Talk Pie Chart Slab %i Chip %i NHits", islab, chip_position[ichip]);
      pie[chip_position[ichip]] = new TPie(piename, piename, 2, npie[islab][chip_position[ichip]], colors);
      pie[chip_position[ichip]]->SetEntryRadiusOffset(1,.05);
      pie[chip_position[ichip]]->SetEntryFillStyle(1,3030);
      pie[chip_position[ichip]]->SetCircle(.5,.45,.3);
      pie[chip_position[ichip]]->Draw("rsc");
    }
    canvas1->Print(resultname + "_analysis.pdf", "pdf");
  }
  canvas1->Print(resultname + "_analysis.pdf]", "pdf");
}
