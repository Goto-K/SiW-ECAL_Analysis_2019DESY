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

void mip_sn(){

  // Define
  Int_t MaxSlab = 5;
  Int_t MaxChip = 16;
  Int_t MaxChannel = 64;
  Int_t MaxMemory = 15;

  Double_t mip_mean[5][16][64] = {0};
  Double_t mip_error[5][16][64] = {0};
  Double_t mip_chi2ndf[5][16][64] = {0};

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
      //cout << tmp_mip_mean << " " << tmp_mip_chi/tmp_mip_ndf << endl;
    }
    else if(mip_mean[tmp_slab][tmp_chip][tmp_channel]!=0 && (tmp_mip_chi/tmp_mip_ndf)<mip_chi2ndf[tmp_slab][tmp_chip][tmp_channel]){
      mip_mean[tmp_slab][tmp_chip][tmp_channel] = tmp_mip_mean;
      mip_error[tmp_slab][tmp_chip][tmp_channel] = tmp_mip_error;
      mip_chi2ndf[tmp_slab][tmp_chip][tmp_channel] = tmp_mip_chi/tmp_mip_ndf;
    }
  }


  // Pedestal Tree open 
  TFile *f[5];
  TTree *Pedestal_Tree[5];
  Double_t pedestal_width_copy[5][16][15][64] = {0};
  Double_t pedestal_width_sum[5][16][64] = {0};
  Double_t s_over_n[5][16][64] = {0};
  Double_t s_over_n_slab[5] = {0};
  Double_t s_over_n_slab_sum[5] = {0};

  for(Int_t islab=0; islab<MaxSlab; islab++){
    TString slabname = "";
    if(islab==0) slabname = "P1";
    if(islab==1) slabname = "P2";
    if(islab==2) slabname = "P3";
    if(islab==3) slabname = "K1";
    if(islab==4) slabname = "K2";

    TString pedestal_filename = "/Users/kiichigoto/Desktop/laboratory/ILC/Test_Beam_2019/Analysis/Pedestal_Map/" + slabname + "_Pedestal.root";
    f[islab] = TFile::Open(pedestal_filename); 
    Pedestal_Tree[islab] = (TTree*)f[islab]->Get("Pedestal_Tree");

    // Reading Calibration data from Tree
    Double_t pedestal_width[16][15][64];
    Pedestal_Tree[islab]->SetBranchAddress("pedestal_width", pedestal_width);

    for(Int_t ievent=0; ievent<1; ievent++){
      Pedestal_Tree[islab]->GetEntry(ievent);
      for(Int_t ichip=0; ichip<MaxChip; ichip++){
        for(Int_t imemory=0; imemory<MaxMemory; imemory++){
	  for(Int_t ichn=0; ichn<MaxChannel; ichn++){
	    pedestal_width_copy[islab][ichip][imemory][ichn]
	    = pedestal_width[ichip][imemory][ichn];
	  }
        }
      }
    }
  }

  // Calculating number of nhits and mip
  Int_t count[5] = {0};
  for(Int_t islab=0; islab<MaxSlab; islab++){
    for(Int_t ichip=0; ichip<MaxChip; ichip++){
      for(Int_t ichn=0; ichn<MaxChannel; ichn++){
        if(mip_chi2ndf[islab][ichip][ichn]==0 || mip_chi2ndf[islab][ichip][ichn]>3) continue;
        if(islab!=2 && mip_mean[islab][ichip][ichn]<100) continue;
        if(islab==2 && mip_mean[islab][ichip][ichn]<30) continue;
	Int_t i = 0;
        for(Int_t imemory=0; imemory<MaxMemory; imemory++){
	  if(pedestal_width_copy[islab][ichip][imemory][ichn]!=0) i++;
	  pedestal_width_sum[islab][ichip][ichn] = pedestal_width_sum[islab][ichip][ichn] + pedestal_width_copy[islab][ichip][imemory][ichn];
	}
	if(i==0) continue;
	s_over_n[islab][ichip][ichn] = mip_mean[islab][ichip][ichn]/(pedestal_width_sum[islab][ichip][ichn]/i);
	s_over_n_slab_sum[islab] = s_over_n_slab_sum[islab] + s_over_n[islab][ichip][ichn];
	count[islab]++;
      }
    }
    s_over_n_slab[islab] = s_over_n_slab_sum[islab]/count[islab];
    cout << "slab " << islab << " S/N ratio : " << s_over_n_slab[islab] << " for " << count[islab] << " channels" << endl;
  }

  ofstream outtext(".S_over_N_ratio/s_over_n.txt", ios::app);
  for(Int_t islab=0; islab<MaxSlab; islab++){
    for(Int_t ichip=0; ichip<MaxChip; ichip++){
      for(Int_t ichn=0; ichn<MaxChannel; ichn++){
        outtext<<TString::Format("%i %i %i %f\n", islab, ichip, ichn, s_over_n[islab][ichip][ichn]);
      }
    }
  }
  outtext.close(); 
}
