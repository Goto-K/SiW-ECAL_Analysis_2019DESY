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

void tdc_multi_channel_v2(std::string str){

  TString filename = str.data();

  // TTree open  
  TFile *file = TFile::Open(filename); 
  TTree *TDC_Tree = (TTree*)file->Get("TDC_Tree");

  Int_t MaxEvent = TDC_Tree->GetEntries();
  Int_t Chip = 13;
  Int_t Channel[] = {39, 42, 45, 48, 50, 54, 56, 58, 60};

  Int_t nhit_chan = 5*16*64*15;

  Int_t nhit_channel;
  TDC_Tree->SetBranchAddress("nhit_channel", &nhit_channel);

  Int_t bcid_hit;
  TDC_Tree->SetBranchAddress("bcid_hit", &bcid_hit);

  Int_t charge_hiGain_tdc[nhit_chan];
  TDC_Tree->SetBranchAddress("charge_hiGain_tdc", charge_hiGain_tdc);

  Int_t slab_number[nhit_chan];
  TDC_Tree->SetBranchAddress("slab_number", slab_number);

  Int_t chip_number[nhit_chan];
  TDC_Tree->SetBranchAddress("chip_number", chip_number);

  Int_t channel_number[nhit_chan];
  TDC_Tree->SetBranchAddress("channel_number", channel_number);

  TGraph *tdc_vs_tdc[9];

  for(Int_t i=0; i<9; i++){
    tdc_vs_tdc[i] = new TGraph();
  }

  TCanvas *canvas = new TCanvas("canvas", "canvas", 1200, 1200);
  TLegend *legend = new TLegend(0.70, 0.01, 0.99, 0.39);

  canvas->cd(); 
  Int_t i = 0;
  for(Int_t ievent=0; ievent<MaxEvent; ievent++){
    TDC_Tree->GetEntry(ievent);
    for(Int_t i=0; i<9; i++){
      Int_t ich = Channel[i];
      if(slab_number[0]==0&&slab_number[1]==1&&slab_number[2]==2&&slab_number[3]==3&&slab_number[4]==4&&
	 chip_number[0]==13&&chip_number[1]==13&&chip_number[2]==13&&chip_number[3]==13&&chip_number[4]==13&&
	 channel_number[0]==ich&&channel_number[1]==ich&&channel_number[2]==ich&&channel_number[3]==ich&&channel_number[4]==ich&&nhit_channel<6&&bcid_hit%2==0){
         
	 tdc_vs_tdc[i]->SetPoint(tdc_vs_tdc[i]->GetN(), charge_hiGain_tdc[0], charge_hiGain_tdc[1]);

      }
    }
  }

  for(Int_t i=0; i<9; i++){
    TString GraphName = TString::Format("chip_%d_channel_%d", Chip, Channel[i]);
    tdc_vs_tdc[i]->SetMarkerStyle(4);
    tdc_vs_tdc[i]->SetMarkerColor(i+1);
    tdc_vs_tdc[i]->GetXaxis()->SetTitle("TDC_P1");
    tdc_vs_tdc[i]->GetYaxis()->SetTitle("TDC_P2");
    if(i==0) tdc_vs_tdc[i]->Draw("AP");
    else{
      tdc_vs_tdc[i]->Draw("PSAME");
    }
    legend->AddEntry(tdc_vs_tdc[i], GraphName, "p");

  }
  gPad->Modified();
  gPad->Update();
  legend->Draw();

  canvas->Print("Test.pdf");

}


