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

void tdc_analysis(std::string str){

  TString filename = str.data();

  // TTree open  
  TFile *file = TFile::Open(filename); 
  TTree *TDC_Tree = (TTree*)file->Get("TDC_Tree");

  Int_t MaxEvent = TDC_Tree->GetEntries();
  Int_t MaxChip = 16;
  Int_t MaxChannel = 64;

  Int_t nhit_chan = 5*16*64*15;

  Int_t nhit_channel;
  TDC_Tree->SetBranchAddress("nhit_channel", &nhit_channel);

  Int_t bcid_hit;
  TDC_Tree->SetBranchAddress("bcid_hit", &bcid_hit);

  Int_t charge_hiGain_tdc[nhit_chan];
  TDC_Tree->SetBranchAddress("charge_hiGain_tdc", charge_hiGain_tdc);

  Int_t charge_lowGain_adc[nhit_chan];
  TDC_Tree->SetBranchAddress("charge_lowGain_adc", charge_lowGain_adc);

  Int_t slab_number[nhit_chan];
  TDC_Tree->SetBranchAddress("slab_number", slab_number);

  Int_t chip_number[nhit_chan];
  TDC_Tree->SetBranchAddress("chip_number", chip_number);

  Int_t channel_number[nhit_chan];
  TDC_Tree->SetBranchAddress("channel_number", channel_number);

  /*TGraph *tdc0_vs_tdc1[MaxChip][MaxChannel];
  TGraph *tdc0_vs_tdc2[MaxChip][MaxChannel];
  TGraph *tdc0_vs_tdc3[MaxChip][MaxChannel];
  TGraph *tdc0_vs_tdc4[MaxChip][MaxChannel];
  TGraph *tdc1_vs_tdc2[MaxChip][MaxChannel];
  TGraph *tdc1_vs_tdc3[MaxChip][MaxChannel];
  TGraph *tdc1_vs_tdc4[MaxChip][MaxChannel];
  TGraph *tdc2_vs_tdc3[MaxChip][MaxChannel];
  TGraph *tdc2_vs_tdc4[MaxChip][MaxChannel];
  TGraph *tdc3_vs_tdc4[MaxChip][MaxChannel];*/

  TGraph *tdc_vs_tdc[9][MaxChip][MaxChannel];

  for(Int_t ichip=0; ichip<MaxChip; ichip++){
    for(Int_t ichn=0; ichn<MaxChannel; ichn++){
      /*tdc0_vs_tdc1[ichip][ichn] = new TGraph();
      tdc0_vs_tdc2[ichip][ichn] = new TGraph();
      tdc0_vs_tdc3[ichip][ichn] = new TGraph();
      tdc0_vs_tdc4[ichip][ichn] = new TGraph();
      tdc1_vs_tdc2[ichip][ichn] = new TGraph();
      tdc1_vs_tdc3[ichip][ichn] = new TGraph();
      tdc1_vs_tdc4[ichip][ichn] = new TGraph();
      tdc2_vs_tdc3[ichip][ichn] = new TGraph();
      tdc2_vs_tdc4[ichip][ichn] = new TGraph();
      tdc3_vs_tdc4[ichip][ichn] = new TGraph();*/
    }
  }

  canvas->cd(); 
  for(Int_t ievent=0; ievent<MaxEvent; ievent++){
    TDC_Tree->GetEntry(ievent);
    if(slab_number[0]==0&&slab_number[1]==1&&slab_number[2]==2&&slab_number[3]==3&&slab_number[4]==4&&nhit_channel<6&&bcid_hit%2==0){
      for(Int_t ichip=0; i<MaxChip; ichip++){
        if(chip_number[0]==ichip&&chip_number[1]==ichip&&chip_number[2]==ichip&&chip_number[3]==ichip&&chip_number[4]==ichip){
          for(Int_t ichn=0; i<MaxChnnel; ichn++){
	    if(channel_number[0]==ichn&&channel_number[1]==ichn&&channel_number[2]==ichn&&channel_number[3]==ichn&&channel_number[4]==ichn){
	      tdc0_vs_tdc1[ichip][ichn]->SetPoint(tdc0_vs_tdc1[i]->GetN(), charge_hiGain_tdc[0], charge_hiGain_tdc[1]);
	      tdc0_vs_tdc2[ichip][ichn]->SetPoint(tdc0_vs_tdc2[i]->GetN(), charge_hiGain_tdc[0], charge_hiGain_tdc[2]);
	      tdc0_vs_tdc3[ichip][ichn]->SetPoint(tdc0_vs_tdc3[i]->GetN(), charge_hiGain_tdc[0], charge_hiGain_tdc[3]);
	      tdc0_vs_tdc4[ichip][ichn]->SetPoint(tdc0_vs_tdc4[i]->GetN(), charge_hiGain_tdc[0], charge_hiGain_tdc[4]);
	      tdc1_vs_tdc2[ichip][ichn]->SetPoint(tdc1_vs_tdc2[i]->GetN(), charge_hiGain_tdc[1], charge_hiGain_tdc[2]);
	      tdc1_vs_tdc3[ichip][ichn]->SetPoint(tdc1_vs_tdc3[i]->GetN(), charge_hiGain_tdc[1], charge_hiGain_tdc[3]);
	      tdc1_vs_tdc4[ichip][ichn]->SetPoint(tdc1_vs_tdc4[i]->GetN(), charge_hiGain_tdc[1], charge_hiGain_tdc[4]);
	      tdc2_vs_tdc3[ichip][ichn]->SetPoint(tdc2_vs_tdc3[i]->GetN(), charge_hiGain_tdc[2], charge_hiGain_tdc[3]);
	      tdc2_vs_tdc4[ichip][ichn]->SetPoint(tdc2_vs_tdc4[i]->GetN(), charge_hiGain_tdc[2], charge_hiGain_tdc[4]);
	      tdc3_vs_tdc4[ichip][ichn]->SetPoint(tdc3_vs_tdc4[i]->GetN(), charge_hiGain_tdc[3], charge_hiGain_tdc[4]);
	    }
	  }
	}
      }
    }
  }

  TCanvas *canvas[9][8];
  TLegend *legend[9][8];

  for(Int_t ichip=0; ichip<MaxChip; ichip++){
    for(Int_t ichns=0; ichns<8; ichns++){
      canvas[ichip][ichns] = new TCanvas("canvas", "canvas", 1200, 1200);
      legend[ichip][ichns] = new TLegend(0.70, 0.01, 0.99, 0.39);
    }
  }

  for(Int_t ichip=0; ichip<MaxChip; ichip++){
    for(Int_t ichns=0; ichns<8; ichns++){
      for(Int_t i=0; i<8; i++){
      TString GraphName = TString::Format("chip_%d_channel_%d", ichip, ichns*8+i);
      tdc0_vs_tdc1[ichip][ichns*8+i]->SetMarkerStyle(4);
      tdc0_vs_tdc1[ichip][ichns*8+i]->SetMarkerColor(i+1);
      if(i==0){
        tdc0_vs_tdc1[ichip][ichns*8+i]->GetXaxis()->SetTitle("TDC_P1");
        tdc0_vs_tdc1[ichip][ichns*8+i]->GetYaxis()->SetTitle("TDC_P2");
	tdc0_vs_tdc1[ichip][ichns*8+i]->Draw("AP");
      else{
        tdc0_vs_tdc1[ichip][ichns*8+i]->Draw("PSAME");
      }
    legend->AddEntry(tdc_vs_tdc[i], GraphName, "p");

  }
  gPad->Modified();
  gPad->Update();
  legend->Draw();

  canvas->Print("Test.pdf");

}


