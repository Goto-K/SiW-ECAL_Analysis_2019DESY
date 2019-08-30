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
#include "TSpectrum.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"

using namespace::std;

void ped_time_memory(){

  Int_t MaxSlab = 5;
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
  TFile *f[5];
  TTree *Pedestal_Tree[5];
  Double_t pedestal_mean_copy[5][16][15][64];

  for(Int_t islab=0; islab<MaxSlab; islab++){
    TString slabname = "";
    if(islab==0) slabname = "P1";
    if(islab==1) slabname = "P2";
    if(islab==2) slabname = "P3";
    if(islab==3) slabname = "K1";
    if(islab==4) slabname = "K2";

    TString pedestal_filename = "../pedestal/" + slabname + "_Pedestal.root";
    f[islab] = TFile::Open(pedestal_filename); 
    Pedestal_Tree[islab] = (TTree*)f[islab]->Get("Pedestal_Tree");

    // Reading Calibration data from Tree
    Double_t pedestal_mean[16][15][64];
    Pedestal_Tree[islab]->SetBranchAddress("pedestal_mean", pedestal_mean);

    for(Int_t ievent=0; ievent<1; ievent++){
      Pedestal_Tree[islab]->GetEntry(ievent);
      for(Int_t ichip=0; ichip<MaxChip; ichip++){
        for(Int_t imemory=0; imemory<MaxMemory; imemory++){
	  for(Int_t ichn=0; ichn<MaxChannel; ichn++){
	    pedestal_mean_copy[islab][ichip][imemory][ichn]
	    = TMath::Nint(pedestal_mean[ichip][imemory][ichn]);
	  }
        }
      }
    }
  }

  TGraph *ped_var[5][16][15][64];

  for(Int_t islab=0; islab<MaxSlab; islab++){
    for(Int_t ichip=0; ichip<MaxChip; ichip++){
      for(Int_t imemory=0; imemory<1; imemory++){
        for(Int_t ichn=0; ichn<MaxChannel; ichn++){
	  ped_var[islab][ichip][imemory][ichn] = new TGraph(0);
	}
      }
    }
  }

  // Read Time
  TString time_path = "../run_data/run_second.txt";
  std::ifstream reading_time(time_path, std::ios::in);
  if(!reading_time){
    cout << "run_second is not found" << endl;
  }
  TString tmp_run;
  Float_t tmp_first = 0, tmp_final = 0, tmp_acqtime = 0;
  Int_t MaxRun[5][15][64] = {0};

  while(reading_time){
    reading_time >> tmp_run >> tmp_first >> tmp_final >> tmp_acqtime;

    // TTree open
    TString filename = "./Memory_Result/" + tmp_run + "_merge_PedestalMap_all.root";
    TFile *file = TFile::Open(filename); 
    if(file){
    
      TTree *Pedestal_Tree = (TTree*)file->Get("Pedestal_Tree");

      Int_t MaxEvent = Pedestal_Tree->GetEntries();
      Float_t runtime = tmp_first + tmp_acqtime/2;

      // Read data from Pedestal_Tree
      Double_t pedestal_mean[5][16][1][64];
      Pedestal_Tree->SetBranchAddress("pedestal_mean", pedestal_mean);
      Double_t pedestal_error[5][16][1][64];
      Pedestal_Tree->SetBranchAddress("pedestal_error", pedestal_error);

      for(Int_t ievent=0; ievent<MaxEvent; ievent++){
        Pedestal_Tree->GetEntry(ievent);
        for(Int_t islab=0; islab<MaxSlab; islab++){
          for(Int_t ichip=0; ichip<MaxChip; ichip++){
            for(Int_t imemory=0; imemory<1; imemory++){
	      for(Int_t ichn=0; ichn<MaxChannel; ichn++){
		Float_t adc = pedestal_mean[islab][ichip][imemory][ichn] - pedestal_mean_copy[islab][ichip][imemory][ichn];
		if(pedestal_error[islab][ichip][imemory][ichn]>1 || pedestal_error[islab][ichip][imemory][ichn]==0 || 
		   pedestal_mean[islab][ichip][imemory][ichn]<10) continue;
	        ped_var[islab][ichip][imemory][ichn]->SetPoint(ped_var[islab][ichip][imemory][ichn]->GetN(), runtime, adc);
	      }
	    }
	  }
        }
      }
    }
  }
 
  Int_t imemory = 0;
  Int_t Color[8] = {1, 2, 3, 4, 6, 7, 8, 9};
  TCanvas *canvas_pedestal_map[5][16];
  TLegend *legend[5][16][8];
  TText *t[5][15][64];
  for(Int_t islab=0; islab<MaxSlab; islab++){
    for(Int_t ichip=0; ichip<MaxChip; ichip++){
      canvas_pedestal_map[islab][ichip]= new TCanvas(TString::Format("pedestal_map_chip%i",ichip),
		                                     TString::Format("pedestal_map_chip%i",ichip), 1200, 800);
      canvas_pedestal_map[islab][ichip]->Print(TString::Format("./Time_Dependance/slab_%i/pedestal_vs_time_chip%i.pdf[", islab, ichip));
      for(Int_t ichn=0; ichn<8; ichn++){
        legend[islab][ichip][ichn] = new TLegend( 0.6, 0.6, 0.90, 0.90);
        canvas_pedestal_map[islab][ichip]->DrawFrame(0, -40, 835000, 40);
        for(Int_t jchn=0; jchn<8; jchn++){
          ped_var[islab][ichip][imemory][ichn*8+jchn]->GetXaxis()->SetTitle("Time[second]");
          ped_var[islab][ichip][imemory][ichn*8+jchn]->GetYaxis()->SetTitle("Pedestal Mean");
          ped_var[islab][ichip][imemory][ichn*8+jchn]->SetMarkerColor(Color[jchn]);
          ped_var[islab][ichip][imemory][ichn*8+jchn]->SetMarkerSize(5);
          ped_var[islab][ichip][imemory][ichn*8+jchn]->SetLineColor(Color[jchn]);
	  ped_var[islab][ichip][imemory][ichn*8+jchn]->Draw("PLsame");
          legend[islab][ichip][ichn]->AddEntry(ped_var[islab][ichip][imemory][ichn*8+jchn],
			                       TString::Format("ped_slab%i_chip%i_chn%i_sca%i", islab, ichip, ichn*8+jchn, imemory), "p");
	  gPad->Modified();
          gPad->Update();
        }
        legend[islab][ichip][ichn]->Draw("same");
        canvas_pedestal_map[islab][ichip]->RedrawAxis();
        canvas_pedestal_map[islab][ichip]->Print(TString::Format("./Time_Dependance/slab_%i/pedestal_vs_time_chip%i.pdf", islab, ichip));
      }
      canvas_pedestal_map[islab][ichip]->Print(TString::Format("./Time_Dependance/slab_%i/pedestal_vs_time_chip%i.pdf]", islab, ichip));
    }
  }

}


