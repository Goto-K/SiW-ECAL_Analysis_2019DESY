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

using namespace::std;

void build_event_goto(std::string str){

  // Cut File name
  int path_i = str.find_last_of("/")+1;
  int ext_i = str.find_last_of(".");
  std::string pathname = str.substr(0, path_i);
  std::string extroot = str.substr(ext_i, str.size()-ext_i);
  std::string filename_ = str.substr(path_i, ext_i-path_i);

  TString fullpath = str.data();
  TString filename = pathname + filename_;

  // TTree open  
  TFile *file = TFile::Open(fullpath); 
  TTree *fev10 = (TTree*)file->Get("fev10");

  // Set Loop Max
  Int_t MaxEvent = fev10->GetEntries();
  Int_t MaxSlab = 5;
  Int_t MaxChip = 16;
  Int_t MaxChannel = 64;
  Int_t MaxMemory = 15;

  // Reading ADC data from Tree
  Int_t acqNumber;
  fev10->SetBranchAddress("acqNumber", &acqNumber);

  Int_t bcid[5][16][15];
  fev10->SetBranchAddress("bcid", bcid);

  Int_t corrected_bcid[5][16][15];
  fev10->SetBranchAddress("corrected_bcid", corrected_bcid);

  Int_t badbcid[5][16][15];
  fev10->SetBranchAddress("badbcid", badbcid);

  Int_t nhits[5][16][15];
  fev10->SetBranchAddress("nhits", nhits);

  Int_t charge_hiGain[5][16][15][64];
  fev10->SetBranchAddress("charge_hiGain", charge_hiGain);

  Int_t charge_lowGain[5][16][15][64];
  fev10->SetBranchAddress("charge_lowGain", charge_lowGain);

  Int_t gain_hit_high[5][16][15][64];
  fev10->SetBranchAddress("gain_hit_high", gain_hit_high);

  // Read Chip Map
  TString mapchip_path = "./map_chip_fev13.dat";
  std::ifstream reading_file(mapchip_path, std::ios::in);
  if(!reading_file){
    cout << "map_chip.dat is not found" << endl;
  }
  Float_t map_pointX[16][64], map_pointY[16][64];
  Float_t map_pointZ[5] = {0.5, 2.5, 4.5, 6.5, 8.5};
  Int_t tmp_chip = 0, tmp_channel = 0;
  Float_t tmp_x0 = 0, tmp_y0 = 0, tmp_x = 0, tmp_y = 0;
  while(reading_file){
    reading_file >> tmp_chip >> tmp_x0 >> tmp_y0 >> tmp_channel >> tmp_x >> tmp_y;
    map_pointX[tmp_chip][tmp_channel] = -tmp_x;
    map_pointY[tmp_chip][tmp_channel] = tmp_y;
  }



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

    TString pedestal_filename = "./pedestal/" + slabname + "_Pedestal.root";
//  TFile *file_1 = TFile::Open(slab.data() + "_Pedestal.root"); 
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
  //// Define //---------------------------------------------------------------//


  Int_t NMAX = 16*64*15*5;
  Int_t nhit_channel;
  Int_t bcid_hit;
  Int_t bcid_original[NMAX];
  Int_t slab_number[NMAX];
  Int_t chip_number[NMAX];
  Int_t channel_number[NMAX];
  Int_t sca_number[NMAX];
  Float_t hit_x[NMAX];
  Float_t hit_y[NMAX];
  Float_t hit_z[NMAX];
  Int_t charge_hiGain_tdc[NMAX];
  Int_t charge_lowGain_adc[NMAX];

  // Create New TTree
  TFile *fout = new TFile(filename + "_BUILD_EVENT.root" , "RECREATE");
  TTree *Event_Tree = new TTree("Event_Tree", "Event_Tree");
  Event_Tree->Branch("acqNumber", &acqNumber,
		        "acqNumber/I");
  Event_Tree->Branch("nhit_channel", &nhit_channel,
		        "nhit_channel/I");
  Event_Tree->Branch("bcid_hit", &bcid_hit,
		        "bcid_hit/I");
  Event_Tree->Branch("bcid_original", bcid_original,
		        "bcid_original[nhit_channel]/I");
  Event_Tree->Branch("slab_number", slab_number,
		        "slab_number[nhit_channel]/I");
  Event_Tree->Branch("chip_number", chip_number,
		        "chip_number[nhit_channel]/I");
  Event_Tree->Branch("sca_number", sca_number,
		        "sca_number[nhit_channel]/I");
  Event_Tree->Branch("hit_x", hit_x,
		        "hit_x[nhit_channel]/F");
  Event_Tree->Branch("hit_y", hit_y,
		        "hit_y[nhit_channel]/F");
  Event_Tree->Branch("hit_z", hit_z,
		        "hit_z[nhit_channel]/F");
  Event_Tree->Branch("channel_number", channel_number,
		        "channel_number[nhit_channel]/I");
  Event_Tree->Branch("charge_hiGain_tdc", charge_hiGain_tdc,
		        "charge_hiGain_tdc[nhit_channel]/I");
  Event_Tree->Branch("charge_lowGain_adc", charge_lowGain_adc,
		        "charge_lowGain_adc[nhit_channel]/I");

  //// Main Data Analysis Loop //----------------------------------------------//

  fout->cd();

  Int_t nhit_all[MaxEvent];
  for(Int_t ievent=0; ievent<MaxEvent; ievent++){
    if(float(ievent%500)==0) cout << "Event" << ievent << "End" << endl;
    fev10->GetEntry(ievent);

    // Fill bcid
    std::vector<std::vector<Int_t> > bcid_all;
    bcid_all.resize(4096);
    for(Int_t ibcid=0; ibcid<4096; ibcid++){
      for(Int_t islab=0; islab<MaxSlab; islab++){
        for(Int_t ichip=0; ichip<MaxChip; ichip++){
          for(Int_t imemory=0; imemory<MaxMemory; imemory++){
	    if(badbcid[islab][ichip][imemory]!=0) continue;
	    else if(bcid[islab][ichip][imemory]==ibcid){ // bcid loop
	      for(Int_t ichn=0; ichn<MaxChannel; ichn++){
	      if(gain_hit_high[islab][ichip][imemory][ichn]==0||gain_hit_high[islab][ichip][imemory][ichn]<0||charge_hiGain[islab][ichip][imemory][ichn]<0) continue;
	      bcid_all.at(ibcid).push_back(ibcid);
	      bcid_all.at(ibcid).push_back(islab);
	      bcid_all.at(ibcid).push_back(ichip);
	      bcid_all.at(ibcid).push_back(ichn);
	      bcid_all.at(ibcid).push_back(imemory);
	      bcid_all.at(ibcid).push_back(map_pointX[ichip][ichn]);
	      bcid_all.at(ibcid).push_back(map_pointY[ichip][ichn]);
	      bcid_all.at(ibcid).push_back(map_pointZ[islab]);
	      bcid_all.at(ibcid).push_back(charge_hiGain[islab][ichip][imemory][ichn]);
	      bcid_all.at(ibcid).push_back(charge_lowGain[islab][ichip][imemory][ichn] - pedestal_mean_copy[islab][ichip][imemory][ichn]);
	      nhit_all[ievent]++;
	      }
	    }
	  }
	}
      }
    }

    // merge bcid < 2
    std::vector<std::vector<Int_t> > bcid_merge;
    bcid_merge.resize(4096);
    for(Int_t ibcid=1; ibcid<4095; ibcid++){
      if(!bcid_all[ibcid].empty()){
	bcid_merge[ibcid] = bcid_all[ibcid];
	if(!bcid_all[ibcid+1].empty() && bcid_all[ibcid].size()>bcid_all[ibcid+1].size()){
	  bcid_merge[ibcid].insert(bcid_merge[ibcid].end(), bcid_all[ibcid+1].begin(), bcid_all[ibcid+1].end());
	  bcid_all[ibcid+1].clear();
	  if(!bcid_all[ibcid-1].empty() && bcid_all[ibcid].size()>bcid_all[ibcid-1].size()){
	    bcid_merge[ibcid].insert(bcid_merge[ibcid].begin(), bcid_all[ibcid-1].begin(), bcid_all[ibcid-1].end());
	    bcid_all[ibcid-1].clear();
	  }
	}
      }
    }

    // Fill Tree
    for(Int_t ibcid=1; ibcid<4095; ibcid++){
      if(!bcid_merge[ibcid].empty()){
	Int_t Nsize = bcid_merge[ibcid].size();
	nhit_channel = Nsize/10;
	
        for(Int_t i=0; i<nhit_channel; i++){
	  bcid_hit = ibcid;
	  bcid_original[i] = bcid_merge[ibcid][i*10+0];
	  slab_number[i] = bcid_merge[ibcid][i*10+1];
	  chip_number[i] = bcid_merge[ibcid][i*10+2];
	  channel_number[i] = bcid_merge[ibcid][i*10+3];
	  sca_number[i] = bcid_merge[ibcid][i*10+4];
	  hit_x[i] = bcid_merge[ibcid][i*10+5];
	  hit_y[i] = bcid_merge[ibcid][i*10+6];
	  hit_z[i] = bcid_merge[ibcid][i*10+7];
	  charge_hiGain_tdc[i] = bcid_merge[ibcid][i*10+8];
	  charge_lowGain_adc[i] = bcid_merge[ibcid][i*10+9];
	}
      Event_Tree->Fill();
      }
    }
  }
//  cout << MaxEvent << endl;
  Event_Tree->Write();
  fout->Close();
}
