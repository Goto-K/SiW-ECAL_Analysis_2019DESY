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

using namespace::std;

void tdc_analysis(std::string str, Int_t most_hit_chip, Int_t tdc_limit = 4096){
  // Cut File name
  int path_i = str.find_last_of("/")+1;
  int ext_i = str.find_last_of(".");
  std::string pathname = str.substr(0, path_i);
  std::string filename_ = str.substr(path_i, ext_i-path_i);
  //std::string extroot = filename_.substr(ext_i, str.size()-ext_i);

  TString fullpath = str.data();
  TString filename = pathname + filename_;

  // ================================================================================================================================================ //
  // TTree open  
  TFile *file = TFile::Open(fullpath); 
  TTree *Event_Tree = (TTree*)file->Get("Event_Tree");

  // ================================================================================================================================================ //
  // Set Loop Max
  Int_t MaxEvent = Event_Tree->GetEntries();
  Int_t NOver5Hits = Event_Tree->GetEntries("nhit_channel>5");
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
  Int_t slab_number[NMAX];
  Event_Tree->SetBranchAddress("slab_number", slab_number);
  Int_t chip_number[NMAX];
  Event_Tree->SetBranchAddress("chip_number", chip_number);
  Int_t sca_number[NMAX];
  Event_Tree->SetBranchAddress("sca_number", sca_number);
  Int_t channel_number[NMAX];
  Event_Tree->SetBranchAddress("channel_number", channel_number);
  Int_t charge_hiGain_tdc[NMAX];
  Event_Tree->SetBranchAddress("charge_hiGain_tdc", charge_hiGain_tdc);
  Int_t charge_lowGain_adc[NMAX];
  Event_Tree->SetBranchAddress("charge_lowGain_adc", charge_lowGain_adc);

  // ================================================================================================================================================ //
  // Make graph
  TGraph* tdc_vs_adc_fo[64];
  TGraph* tdc_vs_adc_fe[64];
  TGraph* tdc_vs_adc_bo[64];
  TGraph* tdc_vs_adc_be[64];
  TGraph* tdc_vs_tdc_o[64];
  TGraph* tdc_vs_tdc_e[64];
  for(Int_t i=0; i<64; i++){
    TString nametvto = TString::Format("chip_%i_channel_%i_tdc_vs_tdc_odd", most_hit_chip, i);
    TString nametvte = TString::Format("chip_%i_channel_%i_tdc_vs_tdc_even", most_hit_chip, i);
    TString nametvafo = TString::Format("chip_%i_channel_%i_cut_front_odd", most_hit_chip, i);
    TString nametvafe = TString::Format("chip_%i_channel_%i_cut_front_even", most_hit_chip, i);
    TString nametvabo = TString::Format("chip_%i_channel_%i_cut_back_odd", most_hit_chip, i);
    TString nametvabe = TString::Format("chip_%i_channel_%i_cut_back_even", most_hit_chip, i);
    tdc_vs_tdc_o[i] = new TGraph(0);
    tdc_vs_tdc_o[i]->SetTitle(nametvto);
    tdc_vs_tdc_o[i]->SetName(nametvto);
    tdc_vs_tdc_e[i] = new TGraph(0);
    tdc_vs_tdc_e[i]->SetTitle(nametvte);
    tdc_vs_tdc_e[i]->SetName(nametvte);
    tdc_vs_adc_fo[i] = new TGraph(0);
    tdc_vs_adc_fo[i]->SetTitle(nametvafo);
    tdc_vs_adc_fo[i]->SetName(nametvafo);
    tdc_vs_adc_fe[i] = new TGraph(0);
    tdc_vs_adc_fe[i]->SetTitle(nametvafe);
    tdc_vs_adc_fe[i]->SetName(nametvafe);
    tdc_vs_adc_bo[i] = new TGraph(0);
    tdc_vs_adc_bo[i]->SetTitle(nametvabo);
    tdc_vs_adc_bo[i]->SetName(nametvabo);
    tdc_vs_adc_be[i] = new TGraph(0);
    tdc_vs_adc_be[i]->SetTitle(nametvabe);
    tdc_vs_adc_be[i]->SetName(nametvabe);
  }
  // ================================================================================================================================================ //
  // Fill Event Display
  bool odd = false;
  for (Int_t ievent=0; ievent<MaxEvent;ievent++){
    Event_Tree->GetEntry(ievent);
    if(float(ievent%2000)==0) cout << "Event" << ievent << "End" << endl;

    if(bcid_hit%2==1) odd = true;
    if(bcid_hit%2==0) odd = false;
    if(nhit_channel>5){

      // =============================================================================================== //
      // Cut condition loop
      Int_t nhit[5][16][64] = {0};
      Int_t tdc[5][16][64] = {0};
      Int_t adc[5][16][64] = {0};
      //Int_t nchips[16] = {0};
      //Int_t nchannels[64] = {0};
      bool bhit[64] = {false};
      for(Int_t ihit=0; ihit<nhit_channel; ihit++){
	nhit[slab_number[ihit]][chip_number[ihit]][channel_number[ihit]]++;
	tdc[slab_number[ihit]][chip_number[ihit]][channel_number[ihit]] = charge_hiGain_tdc[ihit];
	adc[slab_number[ihit]][chip_number[ihit]][channel_number[ihit]] = charge_lowGain_adc[ihit];
      }
      for(Int_t ichn=0; ichn<MaxChannel; ichn++){
        if(nhit[0][most_hit_chip][ichn]>0&&nhit[1][most_hit_chip][ichn]>0&&
	   nhit[2][most_hit_chip][ichn]>0&&nhit[3][most_hit_chip][ichn]>0&&
	   nhit[4][most_hit_chip][ichn]>0) bhit[ichn] = true;
      }

      // =============================================================================================== //
      for(Int_t islab=0; islab<MaxSlab; islab++){
        for(Int_t ichn=0; ichn<MaxChannel; ichn++){
          if(tdc[0][most_hit_chip][ichn]>tdc_limit || !bhit[ichn]) continue;
	  if(odd){
	    tdc_vs_tdc_o[ichn]->SetPoint(tdc_vs_tdc_o[ichn]->GetN(), tdc[0][most_hit_chip][ichn], tdc[1][most_hit_chip][ichn]);
	  }
	  if(!odd){
	    tdc_vs_tdc_e[ichn]->SetPoint(tdc_vs_tdc_e[ichn]->GetN(), tdc[0][most_hit_chip][ichn], tdc[1][most_hit_chip][ichn]);
	  }
	  if(adc[0][most_hit_chip][ichn]>100 && adc[0][most_hit_chip][ichn]<200 && odd){
	    tdc_vs_adc_fo[ichn]->SetPoint(tdc_vs_adc_fo[ichn]->GetN(), adc[1][most_hit_chip][ichn], tdc[0][most_hit_chip][ichn]-tdc[1][most_hit_chip][ichn]);
	  }
	  if(adc[0][most_hit_chip][ichn]>100 && adc[0][most_hit_chip][ichn]<200 && !odd){
	    tdc_vs_adc_fe[ichn]->SetPoint(tdc_vs_adc_fe[ichn]->GetN(), adc[1][most_hit_chip][ichn], tdc[0][most_hit_chip][ichn]-tdc[1][most_hit_chip][ichn]);
	  }
	  if(adc[1][most_hit_chip][ichn]>100 && adc[1][most_hit_chip][ichn]<200 && odd){
	    tdc_vs_adc_bo[ichn]->SetPoint(tdc_vs_adc_bo[ichn]->GetN(), adc[0][most_hit_chip][ichn], tdc[0][most_hit_chip][ichn]-tdc[1][most_hit_chip][ichn]);
	  }
	  if(adc[1][most_hit_chip][ichn]>100 && adc[1][most_hit_chip][ichn]<200 && !odd){
	    tdc_vs_adc_be[ichn]->SetPoint(tdc_vs_adc_be[ichn]->GetN(), adc[0][most_hit_chip][ichn], tdc[0][most_hit_chip][ichn]-tdc[1][most_hit_chip][ichn]);
	  }
	}
      }
    }

  }
  // ================================================================================================================================================ //
  // Write Event Display
  TFile *fout = new TFile(filename+"_Time_walk.root" , "RECREATE");
  TCanvas *c = new TCanvas();
  c->Divide(2, 3);
  c->Print(filename + "_Time_walk.pdf[", "pdf");
  for(Int_t ichn=0; ichn<MaxChannel; ichn++) {
    c->cd(1);
    tdc_vs_tdc_o[ichn]->GetXaxis()->SetTitle("TDC charge_hiGain slab P1");
    tdc_vs_tdc_o[ichn]->GetYaxis()->SetTitle("TDC charge_hiGain slab P2");
    tdc_vs_tdc_o[ichn]->Write();
    tdc_vs_tdc_o[ichn]->Draw("AP");
    c->cd(2);
    tdc_vs_tdc_e[ichn]->GetXaxis()->SetTitle("TDC charge_hiGain slab P1");
    tdc_vs_tdc_e[ichn]->GetYaxis()->SetTitle("TDC charge_hiGain slab P2");
    tdc_vs_tdc_e[ichn]->Write();
    tdc_vs_tdc_e[ichn]->Draw("AP");
    c->cd(3);
    tdc_vs_adc_fo[ichn]->GetXaxis()->SetTitle("ADC charge_lowGain slab P2");
    tdc_vs_adc_fo[ichn]->GetYaxis()->SetTitle("TDC difference charge_hiGain slab P1 - slab P2");
    tdc_vs_adc_fo[ichn]->Write();
    tdc_vs_adc_fo[ichn]->Draw("AP");
    c->cd(4);
    tdc_vs_adc_fe[ichn]->GetXaxis()->SetTitle("ADC charge_lowGain slab P2");
    tdc_vs_adc_fe[ichn]->GetYaxis()->SetTitle("TDC difference charge_hiGain slab P1 - slab P2");
    tdc_vs_adc_fe[ichn]->Write();
    tdc_vs_adc_fe[ichn]->Draw("AP");
    c->cd(5);
    tdc_vs_adc_bo[ichn]->GetXaxis()->SetTitle("ADC charge_lowGain slab P1");
    tdc_vs_adc_bo[ichn]->GetYaxis()->SetTitle("TDC difference charge_hiGain slab P1 - slab P2");
    tdc_vs_adc_bo[ichn]->Write();
    tdc_vs_adc_bo[ichn]->Draw("AP");
    c->cd(6);
    tdc_vs_adc_be[ichn]->GetXaxis()->SetTitle("ADC charge_lowGain slab P1");
    tdc_vs_adc_be[ichn]->GetYaxis()->SetTitle("TDC difference charge_hiGain slab P1 - slab P2");
    tdc_vs_adc_be[ichn]->Write();
    tdc_vs_adc_be[ichn]->Draw("AP");
    c->Print(filename + "_Time_walk.pdf", "pdf");
  }
  c->Print(filename + "_Time_walk.pdf]", "pdf");
  fout->cd();
  fout->Close();
}
