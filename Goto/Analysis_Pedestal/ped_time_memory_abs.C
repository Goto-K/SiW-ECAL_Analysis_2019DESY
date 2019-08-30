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
#include "TDatime.h"
#include "TGaxis.h"

using namespace::std;

void ped_time_memory_abs(){

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
  UInt_t Offset = 788918400;

  // Pedestal Tree open 
  TFile *f[5];
  TTree *Pedestal_Tree[5];
  Double_t pedestal_mean_copy[5][16][15][64];

  TGraph *ped_var[5][16][15][64];
  TGraph *temp_var[5];

  for(Int_t islab=0; islab<MaxSlab; islab++){
    temp_var[islab] = new TGraph(0);
    for(Int_t ichip=0; ichip<MaxChip; ichip++){
      for(Int_t imemory=0; imemory<1; imemory++){
        for(Int_t ichn=0; ichn<MaxChannel; ichn++){
	  ped_var[islab][ichip][imemory][ichn] = new TGraph(0);
	}
      }
    }
  }

  // Read Temperature
  TString temp_path = "../run_data/run_temperature.txt";
  std::ifstream reading_temp(temp_path, std::ios::in);
  if(!reading_temp){
    cout << "run_temperature is not found" << endl;
  }
  Float_t tmp_temp;
  Int_t tmp_time[5] = {0}, tmp_slab;

  while(reading_temp){
    reading_temp >> tmp_time[0] >> tmp_time[1] >> tmp_time[2] >> tmp_time[3] >> tmp_time[4] >> tmp_slab >> tmp_temp;
    TDatime *Temptime = new TDatime(2019, tmp_time[0], tmp_time[1], tmp_time[2], tmp_time[3], tmp_time[4]);
    if(tmp_temp>100.0 || tmp_temp<0) continue;
    temp_var[tmp_slab-1]->SetPoint(temp_var[tmp_slab-1]->GetN(), Temptime->Convert()-Offset, tmp_temp+240);
  }

  temp_var[0]->Draw("P");
  
  // Read Time
  TString time_path = "../run_data/run_time.txt";
  std::ifstream reading_time(time_path, std::ios::in);
  if(!reading_time){
    cout << "run_time is not found" << endl;
  }
  TString tmp_run;
  Int_t tmp_begin[6] = {0}, tmp_end[6] = {0};
  Int_t MaxRun = 0;
  TDatime *testBegin = new TDatime();
  TDatime *testEnd = new TDatime();

  while(reading_time){
    reading_time >> tmp_run >> 
	    	    tmp_begin[0] >> tmp_begin[1] >> tmp_begin[2] >> tmp_begin[3] >> tmp_begin[4] >> tmp_begin[5] >>
		    tmp_end[0] >> tmp_end[1] >> tmp_end[2] >> tmp_end[3] >> tmp_end[4] >> tmp_end[5];

    TDatime *acqBegin = new TDatime(tmp_begin[0], tmp_begin[1], tmp_begin[2], tmp_begin[3], tmp_begin[4], tmp_begin[5]);
    if(MaxRun==0) *testBegin = TDatime(tmp_begin[0], tmp_begin[1], tmp_begin[2], tmp_begin[3], tmp_begin[4], tmp_begin[5]);
    *testEnd = TDatime(tmp_end[0], tmp_end[1], tmp_end[2], tmp_end[3], tmp_end[4], tmp_end[5]);
    MaxRun++;
    // TTree open
    TString filename = "./Memory_Result/" + tmp_run + "_merge_PedestalMap_all.root";
    TFile *file = TFile::Open(filename); 
    if(file){
    
      TTree *Pedestal_Tree = (TTree*)file->Get("Pedestal_Tree");

      Int_t MaxEvent = Pedestal_Tree->GetEntries();

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
		Float_t adc = pedestal_mean[islab][ichip][imemory][ichn];
		if(pedestal_error[islab][ichip][imemory][ichn]>1 || pedestal_error[islab][ichip][imemory][ichn]==0 || 
		   pedestal_mean[islab][ichip][imemory][ichn]<10) continue;
	        ped_var[islab][ichip][imemory][ichn]->SetPoint(ped_var[islab][ichip][imemory][ichn]->GetN(), acqBegin->Convert()-Offset, adc);
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
  gStyle->SetTitleXOffset(1.4);
  gStyle->SetPadGridX(1);

  for(Int_t islab=0; islab<MaxSlab; islab++){
    for(Int_t ichip=0; ichip<MaxChip; ichip++){
      canvas_pedestal_map[islab][ichip]= new TCanvas(TString::Format("pedestal_map_chip%i",ichip),
		                                     TString::Format("pedestal_map_chip%i",ichip), 1200, 800);
      canvas_pedestal_map[islab][ichip]->Print(TString::Format("./Time_Dependance/slab_%i/pedestal_vs_time_chip%i.pdf[", islab, ichip));
      for(Int_t ichn=0; ichn<8; ichn++){
        legend[islab][ichip][ichn] = new TLegend( 0.101, 0.7, 0.899, 0.899);
        legend[islab][ichip][ichn]->SetNColumns(3);
        legend[islab][ichip][ichn]->AddEntry(temp_var[islab],
			                     TString::Format("temperature_slab%i", islab), "p");

	TH1F *frame = (TH1F*)canvas_pedestal_map[islab][ichip]->DrawFrame(testBegin->Convert()-Offset, 200, testEnd->Convert()-Offset+15000, 400);
        frame->GetXaxis()->SetTitle("Time [yy/mm/dd]");
        frame->GetYaxis()->SetTitle("Pedestal Mean [ADC ch]");
        frame->GetXaxis()->SetNdivisions(510);
        frame->GetXaxis()->SetTimeDisplay(1);
        frame->GetXaxis()->SetTimeFormat("#splitline{%y/%m/%d}{%H:%M}");
        frame->GetXaxis()->SetLabelOffset(0.02);
	frame->GetXaxis()->SetLabelSize(0.03);

        TGaxis *tgaxis = new TGaxis(testEnd->Convert()-Offset+15000, 200, testEnd->Convert()-Offset+15000, 400, -40, 160, 510, "+L");
	tgaxis->SetTitle("Temperature [degree]");
	tgaxis->SetLabelSize(0.035);
	tgaxis->SetLabelFont(42);
	tgaxis->SetTitleSize(0.035);
	tgaxis->SetTitleFont(42);
        tgaxis->Draw();

	gPad->Modified();
        gPad->Update();

        for(Int_t jchn=0; jchn<8; jchn++){
          temp_var[islab]->SetMarkerColor(46);
          temp_var[islab]->SetMarkerSize(1);
          temp_var[islab]->SetMarkerStyle(2);
          temp_var[islab]->Draw("Psame");

          ped_var[islab][ichip][imemory][ichn*8+jchn]->SetMarkerColor(Color[jchn]);
          ped_var[islab][ichip][imemory][ichn*8+jchn]->SetMarkerSize(5);
          ped_var[islab][ichip][imemory][ichn*8+jchn]->SetMarkerStyle(7);
          ped_var[islab][ichip][imemory][ichn*8+jchn]->SetLineColor(Color[jchn]);
	  ped_var[islab][ichip][imemory][ichn*8+jchn]->Draw("Psame");
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


