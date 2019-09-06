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

void tdc_multi_channel(std::string str){

  TString filename = str.data();

  // TTree open  
  TFile *file = TFile::Open(filename); 
  TTree *TDC_Tree = (TTree*)file->Get("TDC_Tree");

  Int_t MaxEvent = TDC_Tree->GetEntries();
  Int_t Chip = 13;
  Int_t Channel[] = {39, 42, 45, 48, 50, 54, 56, 58, 60};

  Int_t Max_nhit_chan = 5*16*64*15;

  Int_t nhit_channel;
  TDC_Tree->SetBranchAddress("nhit_channel", nhit_channel)

  Double_t charge_hiGain_tdc[nhit_chan];
  TDC_Tree->SetBranchAddress("charge_hiGain_tdc", charge_hiGain_tdc);

  TGraph *tdc_vs_tdc[9];

  TCanvas *canvas = new TCanvas("canvas", "canvas");
  TLegend *legend = new TLegend(0.6, 0.6, 0.99, 0.99);

  canvas->cd(); 
  Int_t i = 0;
  for(Int_t ich : Channel){
    //TString cut_condition = TString::Format("slab_number[0]==0&&slab_number[2]==1&&slab_number[2]==2&&slab_number[3]==3&&slab_number[4]==4&&chip_number[0]==13&&chip_number[1]==13&&chip_number[2]==13&&chip_number[3]==13&&chip_number[4]==13&&channel_number[0]==%d&&channel_number[1]==%d&&channel_number[2]==%d&&channel_number[3]==%d&&channel_number[4]==%d&&nhit_channel<6&&bcid_hit%2==0", ich, ich, ich, ich, ich);
    TString cut_condition = TString::Format("slab_number[0]==0&&slab_number[2]==1&&slab_number[2]==2&&slab_number[3]==3&&slab_number[4]==4&&chip_number[0]==13&&chip_number[1]==13&&chip_number[2]==13&&chip_number[3]==13&&chip_number[4]==13&&channel_number[0]==%d&&channel_number[1]==%d&&channel_number[2]==%d&&channel_number[3]==%d&&channel_number[4]==%d&&nhit_channel<6", ich, ich, ich, ich, ich);
    TString GraphName = TString::Format("chip_%d_channel_%d", Chip, ich);

    tdc_vs_tdc[i] = new TGraph();
    tdc_vs_tdc[i]->SetName(GraphName);
    tdc_vs_tdc[i]->SetTitle(GraphName);

    TDC_Tree->Draw(TString::Format("charge_hiGain_tdc[1]:charge_hiGain_tdc[0]>>%s", GraphName.Data()), cut_condition + "&&bcid_hit%2==0", "PSAME");
    canvas->Print(TString::Format("Test_%d.pdf", i));

    tdc_vs_tdc[i]->SetMarkerStyle(4);
    tdc_vs_tdc[i]->SetMarkerColor(i);

    legend->AddEntry(tdc_vs_tdc[i], GraphName, "p");

    i++;
  }
  legend->Draw();

  canvas->Print("Test.pdf");

}


