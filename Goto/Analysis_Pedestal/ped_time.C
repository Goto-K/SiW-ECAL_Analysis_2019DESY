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

void ped_time(){

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
    temp_var[tmp_slab-1]->SetPoint(temp_var[tmp_slab-1]->GetN(), Temptime->Convert()-Offset, tmp_temp);
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
		Float_t adc = pedestal_mean[islab][ichip][imemory][ichn] - pedestal_mean_copy[islab][ichip][imemory][ichn];
		//Float_t adc = pedestal_mean[islab][ichip][imemory][ichn];
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
  TPad *pad1[5][16][8];
  TPad *pad2[5][16][8];
  TLegend *legend[5][16][8];
  TText *t[5][15][64];
  gStyle->SetTitleXOffset(1.4);
  gStyle->SetPadGridX(1);

  //---- Adust five following parameters--------
  Double_t gap=0.000;            // gap width between plots
  Double_t height1=0.7, height2=1.0-height1;
  Double_t top_margin=0.050;     // top margin <= bottom margin
  Double_t bottom_margin=0.075;  // bottom margin
  Double_t left_margin=0.05;
  Double_t right_margin=0.10;
  //------ corresponding to pad[1](top), pad[2](second top....--------
  Double_t apron=0.05; // pad aprons at top and bottom for label overlaps
  Double_t plotspace=1.0-gap-top_margin-bottom_margin;
  Double_t x1=0.0, x2=1.0, y1=0.0, y2=1.0;

  for(Int_t islab=0; islab<MaxSlab; islab++){
    for(Int_t ichip=0; ichip<MaxChip; ichip++){
      canvas_pedestal_map[islab][ichip]= new TCanvas(TString::Format("pedestal_map_chip%i",ichip),
		                                     TString::Format("pedestal_map_chip%i",ichip), 1200, 800);
      canvas_pedestal_map[islab][ichip]->Print(TString::Format("./Time_Dependance/slab_%i/pedestal_vs_time_chip%i.pdf[", islab, ichip));
      for(Int_t ichn=0; ichn<8; ichn++){
	
	pad1[islab][ichip][ichn] = new TPad("pad1", "pad1", 0.0, 0.5, 1.0, 1.0);
	pad1[islab][ichip][ichn]->SetTopMargin(0.3);
        pad1[islab][ichip][ichn]->SetBottomMargin(0.05);
        //pad1[islab][ichip][ichn]->SetLeftMargin(left_margin);
        //pad1[islab][ichip][ichn]->SetRightMargin(right_margin);
	pad1[islab][ichip][ichn]->SetFillColorAlpha(0,0.00);
        pad1[islab][ichip][ichn]->Draw();
        pad1[islab][ichip][ichn]->cd();
	
	TH1F *frame1 = (TH1F*)pad1[islab][ichip][ichn]->DrawFrame(testBegin->Convert()-Offset, temp_var[islab]->GetHistogram()->GetMinimum()*0.8,
			                                          testEnd->Convert()-Offset+15000, temp_var[islab]->GetHistogram()->GetMaximum()*1.7);
        frame1->GetYaxis()->SetTitle("Temperature");
	frame1->GetYaxis()->SetTitleSize(0.05);
        frame1->GetXaxis()->SetNdivisions(510);
        frame1->GetXaxis()->SetTimeDisplay(1);
	frame1->GetXaxis()->SetLabelOffset(100);

        canvas_pedestal_map[islab][ichip]->cd();
	pad2[islab][ichip][ichn] = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.5);
	pad2[islab][ichip][ichn]->SetTopMargin(0.0);
        pad2[islab][ichip][ichn]->SetBottomMargin(0.25);
        //pad2[islab][ichip][ichn]->SetLeftMargin(left_margin);
        //pad2[islab][ichip][ichn]->SetRightMargin(right_margin);
	pad2[islab][ichip][ichn]->SetFillColorAlpha(1,0.00);
        pad2[islab][ichip][ichn]->Draw();
        pad2[islab][ichip][ichn]->cd();

	Double_t maximum_ped = ped_var[islab][ichip][imemory][ichn*8]->GetHistogram()->GetMaximum();
	Double_t minimum_ped = ped_var[islab][ichip][imemory][ichn*8]->GetHistogram()->GetMinimum();
        for(Int_t jchn=0; jchn<8; jchn++){
	  if(maximum_ped<ped_var[islab][ichip][imemory][ichn*8+jchn]->GetHistogram()->GetMaximum()){
	    maximum_ped = ped_var[islab][ichip][imemory][ichn*8+jchn]->GetHistogram()->GetMaximum();
	  }
	  if(minimum_ped>ped_var[islab][ichip][imemory][ichn*8+jchn]->GetHistogram()->GetMinimum()){
	    minimum_ped = ped_var[islab][ichip][imemory][ichn*8+jchn]->GetHistogram()->GetMinimum();
	  }
	}
        
	TH1F *frame2 = (TH1F*)pad2[islab][ichip][ichn]->DrawFrame(testBegin->Convert()-Offset, minimum_ped-5, testEnd->Convert()-Offset+15000, maximum_ped+5);
        frame2->GetXaxis()->SetTitle("Time [yy/mm/dd]");
        frame2->GetXaxis()->SetTitleOffset(2.0);
	frame2->GetXaxis()->SetTitleSize(0.05);
        frame2->GetYaxis()->SetTitle("Pedestal [ADC counts]");
	frame2->GetYaxis()->SetTitleSize(0.05);
        frame2->GetXaxis()->SetNdivisions(510);
        frame2->GetXaxis()->SetTimeDisplay(1);
        frame2->GetXaxis()->SetTimeFormat("#splitline{%y/%m/%d}{%H:%M}");
        frame2->GetXaxis()->SetLabelOffset(0.02);
	frame2->GetXaxis()->SetLabelSize(0.04);

        pad1[islab][ichip][ichn]->cd();

        gPad->Update();
	Double_t yleg1 = (40*0.65+30 - gPad->GetY1())/(gPad->GetY2()-gPad->GetY1());
	Double_t yleg2 = (40*0.90+30 - gPad->GetY1())/(gPad->GetY2()-gPad->GetY1());

        temp_var[islab]->SetMarkerColor(46);
        temp_var[islab]->SetMarkerSize(1);
        temp_var[islab]->SetMarkerStyle(2);
        temp_var[islab]->Draw("P");

	legend[islab][ichip][ichn] = new TLegend( 0.101, 0.5, 0.899, 0.7);
        legend[islab][ichip][ichn]->SetNColumns(3);
        legend[islab][ichip][ichn]->AddEntry(temp_var[islab],
			                     TString::Format("temperature_slab%i", islab), "p");

	gPad->Modified();
        gPad->Update();

        canvas_pedestal_map[islab][ichip]->cd();
        pad2[islab][ichip][ichn]->cd();
        for(Int_t jchn=0; jchn<8; jchn++){
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
        canvas_pedestal_map[islab][ichip]->cd();
        pad1[islab][ichip][ichn]->cd();
        legend[islab][ichip][ichn]->Draw("same");
        canvas_pedestal_map[islab][ichip]->cd();
        pad2[islab][ichip][ichn]->RedrawAxis();
        canvas_pedestal_map[islab][ichip]->Print(TString::Format("./Time_Dependance/slab_%i/pedestal_vs_time_chip%i.pdf", islab, ichip));
        pad1[islab][ichip][ichn]->cd();
        pad1[islab][ichip][ichn]->Clear();
        pad2[islab][ichip][ichn]->cd();
        pad2[islab][ichip][ichn]->Clear();
        canvas_pedestal_map[islab][ichip]->cd();
      }
      canvas_pedestal_map[islab][ichip]->Print(TString::Format("./Time_Dependance/slab_%i/pedestal_vs_time_chip%i.pdf]", islab, ichip));
    }
  }

}
