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
#include "TLegend.h"
#include "TPie.h"

using namespace::std;

void ped_shift(){

  // Define
  Int_t MaxSlab = 4;
  Int_t MaxChip = 16;
  Int_t MaxChannel = 64;
  Int_t MaxMemory = 15;
  Int_t chip_mapper[64] = {16,  7, 19, 25, 31, 35, 38, 47,
                           10, 22, 13, 28, 32, 41, 44, 53,
  		            9,  3, 15, 27, 33, 45, 50, 60,
  		           21, 12, 18, 30, 36, 48, 54, 39,
  		            6,  5,  4, 20, 40, 42, 58, 56,
  		            2,  8, 11, 23, 37, 49, 63, 57,
		            1, 24, 17, 29, 43, 52, 55, 51,
                            0, 14, 26, 34, 46, 59, 62, 61};

  // Pedestal Tree open 
  TFile *f_ADC[4];
  TTree *Pedestal_Tree_ADC[4];
  TFile *f_TDC[4];
  TTree *Pedestal_Tree_TDC[4];
  Double_t pedestal_mean_shift[4][16][15][64];

  for(Int_t islab=0; islab<MaxSlab; islab++){
    TString slabname = "";
    if(islab==0) slabname = "1_1_1";
    if(islab==1) slabname = "1_1_2";
    if(islab==2) slabname = "1_1_3";
    if(islab==3) slabname = "1_1_5";

    TString pedestal_filename_adc = "./ADC/adc_dif_" + slabname + "_PedestalMap.root";
    f_ADC[islab] = TFile::Open(pedestal_filename_adc); 
    Pedestal_Tree_ADC[islab] = (TTree*)f_ADC[islab]->Get("Pedestal_Tree");

    TString pedestal_filename_tdc = "./TDC/tdc_dif_" + slabname + "_PedestalMap.root";
    f_TDC[islab] = TFile::Open(pedestal_filename_tdc); 
    Pedestal_Tree_TDC[islab] = (TTree*)f_TDC[islab]->Get("Pedestal_Tree");

    // Reading Calibration data from Tree
    Double_t pedestal_mean_adc[16][15][64];
    Pedestal_Tree_ADC[islab]->SetBranchAddress("pedestal_mean", pedestal_mean_adc);
    Double_t pedestal_mean_tdc[16][15][64];
    Pedestal_Tree_TDC[islab]->SetBranchAddress("pedestal_mean", pedestal_mean_tdc);

    for(Int_t ievent=0; ievent<1; ievent++){
      Pedestal_Tree_ADC[islab]->GetEntry(ievent);
      Pedestal_Tree_TDC[islab]->GetEntry(ievent);
      for(Int_t ichip=0; ichip<MaxChip; ichip++){
        for(Int_t imemory=0; imemory<MaxMemory; imemory++){
	  for(Int_t ichn=0; ichn<MaxChannel; ichn++){
	    pedestal_mean_shift[islab][ichip][imemory][ichn]
	    = pedestal_mean_adc[ichip][imemory][ichn]-pedestal_mean_tdc[ichip][imemory][ichn];
	  }
        }
      }
    }
  }

  // Create New TTree
  TFile *fout = new TFile("Pedestal_Shift.root" , "RECREATE");
  TTree *Pedestal_Shift_Tree = new TTree("Pedestal_Shift_Tree", "Pedestal_Shift_Tree");
  Pedestal_Shift_Tree->Branch("pedestal_mean_shift", pedestal_mean_shift,
		              "pedestal_mean_shift[4][16][15][64]/D");

  Pedestal_Shift_Tree->Fill();
  Pedestal_Shift_Tree->Write();
  fout->Close();
}
