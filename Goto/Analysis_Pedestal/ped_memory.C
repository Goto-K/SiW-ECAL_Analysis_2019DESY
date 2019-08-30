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
#include "TSpectrum.h"

using namespace::std;

void ped_memory(std::string str){

  // Cut File name
  int path_i = str.find_last_of("/")+1;
  int ext_i = str.find_last_of(".");
  std::string pathname = str.substr(0, path_i);
  std::string filename_ = str.substr(path_i, ext_i-path_i);

  TString fullpath = str.data();
  TString filename = pathname + filename_;
  TString resultname = "./Memory_Result/" + filename_;

  // TTree_open  
  TFile *file = TFile::Open(fullpath); 
  TTree *fev10 = (TTree*)file->Get("fev10");

  // Define
  Int_t MaxEvent = fev10->GetEntries();
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

  cout << "event_number" << MaxEvent << endl;
  Int_t corrected_bcid[MaxSlab][MaxChip][MaxMemory];
  fev10->SetBranchAddress("corrected_bcid", corrected_bcid);
  Int_t charge_hiGain[MaxSlab][MaxChip][MaxMemory][MaxChannel];
  fev10->SetBranchAddress("charge_hiGain", charge_hiGain); 
  Int_t charge_lowGain[MaxSlab][MaxChip][MaxMemory][MaxChannel];
  fev10->SetBranchAddress("charge_lowGain", charge_lowGain);   
  Int_t gain_hit_high[MaxSlab][MaxChip][MaxMemory][MaxChannel];
  fev10->SetBranchAddress("gain_hit_high", gain_hit_high);
  Int_t gain_hit_low[MaxSlab][MaxChip][MaxMemory][MaxChannel];
  fev10->SetBranchAddress("gain_hit_low", gain_hit_low);
  Int_t badbcid[MaxSlab][MaxChip][MaxMemory];
  fev10->SetBranchAddress("badbcid", badbcid);
  Int_t nhits[MaxSlab][MaxChip][MaxMemory];
  fev10->SetBranchAddress("nhits", nhits);

  // Read Chip Map
  TString mapchip_path = "../map_chip/map_chip_fev13.dat";
  std::ifstream reading_file(mapchip_path, std::ios::in);
  if(!reading_file){
    cout << "map_chip.dat is not found" << endl;
  }
  
  Float_t map_pointX[16][64], map_pointY[16][64];
  Int_t tmp_chip = 0, tmp_channel = 0;
  Float_t tmp_x0 = 0, tmp_y0 = 0, tmp_x = 0, tmp_y = 0;
  while(reading_file){
    reading_file >> tmp_chip >> tmp_x0 >> tmp_y0 >> tmp_channel >> tmp_x >> tmp_y;
    map_pointX[tmp_chip][tmp_channel] = -tmp_x;
    map_pointY[tmp_chip][tmp_channel] = tmp_y;
  }

  /*// Calculating number of hit
  Double_t hit_number[5][16][64]={0};
  Double_t nhits_number[5][64]={0};
  for(Int_t ievent = 0; ievent<MaxEvent; ievent++){
    fev10->GetEntry(ievent);
    for(Int_t islab=0; islab<MaxSlab; islab++){
      for(Int_t ichip=0; ichip<MaxChip; ichip++){
        for(Int_t imemory=0; imemory<MaxMemory; imemory++){
	  if(badbcid[islab][ichip][imemory]==1)continue;
	  for(Int_t ichn=0; ichn<MaxChannel; ichn++){
	    if(gain_hit_high[islab][ichip][imemory][ichn]<0)continue;
	    hit_number[islab][ichip][ichn] += gain_hit_high[islab][ichip][imemory][ichn];
	  }
	}
      }
    }
  }

  // Calculating most hit chips and channels
  Int_t most_hit_chips[5] = {0};
  Int_t most_hit_channels[5] = {0};
  Int_t max[5] = {0};
  for(Int_t islab=0; islab<MaxSlab; islab++){
    for(Int_t ichip=0; ichip<MaxChip; ichip++){
      for(Int_t ichn=0; ichn<MaxChannel; ichn++){
        if(hit_number[islab][ichip][ichn]>max[islab]){
	  max[islab] = hit_number[islab][ichip][ichn];
	  most_hit_chips[islab] = ichip;
	  most_hit_channels[islab] = ichn;
	}
      }
    }
  }*/

  TH2F* pedestal_map[5];
  TH2F* pedestal_width_map[5];
  TH2F* pedestal_error_map[5];

  for(Int_t islab=0; islab<MaxSlab; islab++){
    pedestal_map[islab]= new TH2F(TString::Format("pedestal_map_slab%i",islab),
		                  TString::Format("pedestal_map_slab%i;X[mm];Y[mm]",islab),32,-90,90,32,-90,90);
    pedestal_width_map[islab]= new TH2F(TString::Format("pedestal_width_map_slab%i",islab),
		                        TString::Format("pedestal_width_map_slab%i;X[mm];Y[mm]",islab),32,-90,90,32,-90,90);
    pedestal_error_map[islab]= new TH2F(TString::Format("pedestal_error_map_slab%i",islab),
		                        TString::Format("pedestal_error_map_slab%i;X[mm];Y[mm]",islab),32,-90,90,32,-90,90);
  }

  // Fill to Tree Data
  Double_t pedestal_mean[5][16][1][64] = {-1};
  Double_t pedestal_width[5][16][1][64] = {-1};
  Double_t pedestal_error[5][16][1][64] = {-1};
  Double_t pedestal_chi2ndf[5][16][1][64] = {-1};

  // Create New TTree
  TFile *fout = new TFile(resultname + "_PedestalMap_all.root" , "RECREATE");
  TTree *Pedestal_Tree = new TTree("Pedestal_Tree", "Pedestal_Tree");
  Pedestal_Tree->Branch("pedestal_mean", pedestal_mean,
		        "pedestal_mean[5][16][1][64]/D");
  Pedestal_Tree->Branch("pedestal_width", pedestal_width,
		        "pedestal_width[5][16][1][64]/D");
  Pedestal_Tree->Branch("pedestal_error", pedestal_error,
		        "pedestal_error[5][16][1][64]/D");
  Pedestal_Tree->Branch("pedestal_chi2ndf", pedestal_chi2ndf,
		        "pedestal_chi2ndf[5][16][1][64]/D");

  TDirectory *hist_sca[5];

  // --------------------
  // all sca
  std::vector<std::vector<std::vector<std::vector<TH1F*> > > > ped_sca;

  for(Int_t islab=0; islab<MaxSlab; islab++){
    std::vector<std::vector<std::vector<TH1F*> > > pedtemp_sca1;
    for(Int_t ichip=0; ichip<MaxChip; ichip++){
      std::vector<std::vector<TH1F*> >pedtemp_sca;
      for(Int_t ichn=0; ichn<MaxChannel; ichn++){
        std::vector<TH1F*> pedtemp_sca2;
        for(Int_t imemory=0; imemory<1; imemory++){
	  TH1F *ped_sca2 = new TH1F(TString::Format("ped_slab%i_chip%i_chn%i_sca%i", islab, ichip, ichn, imemory),
	 		            TString::Format("ped_slab%i_chip%i_chn%i_sca%i", islab, ichip, ichn, imemory), 1000, 0.5, 1000.5);
	  pedtemp_sca2.push_back(ped_sca2);
        }
        pedtemp_sca.push_back(pedtemp_sca2);
      }
      pedtemp_sca1.push_back(pedtemp_sca);
    }
    ped_sca.push_back(pedtemp_sca1);
  }

  // -----------------------------------------------------------------------------------------------------   
  // SCA analysis
  for(Int_t ievent=0; ievent<MaxEvent; ievent++){
    fev10->GetEntry(ievent);
    for(Int_t islab=0; islab<MaxSlab; islab++){
      for(Int_t ichip=0; ichip<MaxChip; ichip++){
        for(Int_t imemory=0; imemory<1; imemory++){
          if(badbcid[islab][ichip][imemory]==0){
	    for(Int_t ichn=0; ichn<MaxChannel; ichn++){
	      if(charge_lowGain[islab][ichip][imemory][ichn]>10 && gain_hit_low[islab][ichip][imemory][ichn]==0){
	        ped_sca.at(islab).at(ichip).at(ichn).at(imemory)->Fill(charge_lowGain[islab][ichip][imemory][ichn]);
	      }
	    }
	  }
	}
      }//imemory
    }//ichip 
  }  // end first loop analysis to fill pedestal historgrams

  std::vector<std::vector<std::vector<std::vector<Double_t> > > > ped_mean;
  std::vector<std::vector<std::vector<std::vector<Double_t> > > > ped_error;
  std::vector<std::vector<std::vector<std::vector<Double_t> > > > ped_width;

  //initialize pedestal vectors
  for(Int_t i=0; i<5; i++){
    std::vector<std::vector<std::vector<Double_t> > > slab_ped_mean;
    std::vector<std::vector<std::vector<Double_t> > > slab_ped_error;
    std::vector<std::vector<std::vector<Double_t> > > slab_ped_width;

    for(Int_t j=0; j<16; j++){
      std::vector<std::vector<Double_t> > chip_ped_mean;
      std::vector<std::vector<Double_t> > chip_ped_error;
      std::vector<std::vector<Double_t> > chip_ped_width;

      for(Int_t k=0; k<64; k++){
        std::vector<Double_t> chn_ped_mean;
        std::vector<Double_t> chn_ped_error;
        std::vector<Double_t> chn_ped_width;
        chip_ped_mean.push_back(chn_ped_mean);
        chip_ped_error.push_back(chn_ped_error);
        chip_ped_width.push_back(chn_ped_width);
      }  
      slab_ped_mean.push_back(chip_ped_mean);
      slab_ped_error.push_back(chip_ped_error);
      slab_ped_width.push_back(chip_ped_width);
    }
    ped_mean.push_back(slab_ped_mean);
    ped_error.push_back(slab_ped_error);
    ped_width.push_back(slab_ped_width);
  }

  // do pedestal (chip/channel/sca based) analysis
  for(Int_t islab=0; islab<MaxSlab; islab++){
    cout << "SLAB " << islab << " START" << endl;
    TString TDirectory_Name = Form("SLAB_%i_Pedestal", islab);
    hist_sca[islab] = fout->mkdir(TDirectory_Name);
    hist_sca[islab]->cd();
    for(Int_t ichip=0; ichip<MaxChip; ichip++){
      for(Int_t ichn=0; ichn<MaxChannel; ichn++){
        for(Int_t imemory=0; imemory<1; imemory++){

	  ped_sca.at(islab).at(ichip).at(ichn).at(imemory)->SetTitle(TString::Format("ped_slab%i_chip%i_chn%i_sca%i", islab, ichip, ichn, imemory));
	  ped_sca.at(islab).at(ichip).at(ichn).at(imemory)->SetName(TString::Format("ped_slab%i_chip%i_chn%i_sca%i", islab, ichip, ichn, imemory));

	  ped_mean.at(islab).at(ichip).at(ichn).push_back(0.);
	  ped_width.at(islab).at(ichip).at(ichn).push_back(0.);
	  ped_error.at(islab).at(ichip).at(ichn).push_back(0.);

	  //====================================================================================================
	  if(ped_sca.at(islab).at(ichip).at(ichn).at(imemory)->GetEntries()>10){ //max_entries/2 ){
	    TSpectrum *s = new TSpectrum();
	    int npeaks = s->Search(ped_sca.at(islab).at(ichip).at(ichn).at(imemory), 2, "", 0.2); 

	    if(npeaks > 0){

              Double_t *mean_peak=s->GetPositionX();
              Double_t *mean_high=s->GetPositionY();
              double mean_peak_higher=0;
              double mean_high_higher=0;
	      int npeak_max=0;

              for(int ipeak=0; ipeak<npeaks; ipeak++){
                if(mean_high[ipeak]>mean_high_higher && mean_high[ipeak]>50){
                  mean_high_higher=mean_high[ipeak];
                  mean_peak_higher=mean_peak[ipeak];
		  npeak_max=ipeak;
                }
              }
	    }

	    TH1 *pedestal_histogram = ped_sca.at(islab).at(ichip).at(ichn).at(imemory);
	    Double_t mean_peak = pedestal_histogram->GetBinCenter(pedestal_histogram->GetMaximumBin());
	    
	    TF1 *f0 = new TF1("f0", "gaus", mean_peak-2*ped_sca.at(islab).at(ichip).at(ichn).at(imemory)->GetRMS(), 
			                    mean_peak+2*ped_sca.at(islab).at(ichip).at(ichn).at(imemory)->GetRMS());
	    ped_sca.at(islab).at(ichip).at(ichn).at(imemory)->Fit("f0", "RQNOC");

	    TF1 *f1 = new TF1("f1", "gaus", f0->GetParameter(1)-2.*f0->GetParameter(2), 
			                    f0->GetParameter(1)+2.*f0->GetParameter(2));
	    ped_sca.at(islab).at(ichip).at(ichn).at(imemory)->Fit("f1","RQME");

  	    gStyle->SetOptFit(1111);
	    ped_sca.at(islab).at(ichip).at(ichn).at(imemory)->Write();
	  
	    ped_mean.at(islab).at(ichip).at(ichn).at(imemory)=f1->GetParameter(1);
	    ped_width.at(islab).at(ichip).at(ichn).at(imemory)=f1->GetParameter(2);
	    ped_error.at(islab).at(ichip).at(ichn).at(imemory)=f1->GetParError(1);

	    pedestal_mean[islab][ichip][imemory][ichn] = f1->GetParameter(1);
	    pedestal_width[islab][ichip][imemory][ichn] = f1->GetParameter(2);
	    pedestal_error[islab][ichip][imemory][ichn] = f1->GetParError(1);

	    if(imemory==0){
	      pedestal_map[islab]->Fill(map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , pedestal_mean[islab][ichip][imemory][ichn]);
	      pedestal_width_map[islab]->Fill(map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , pedestal_width[islab][ichip][imemory][ichn]);
	      pedestal_error_map[islab]->Fill(map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , pedestal_error[islab][ichip][imemory][ichn]);
	    }

	  }  

	}

      }

    }
    cout << "SLAB " << islab << " END" << endl;
  }

  fout->cd();
  Pedestal_Tree->Fill();
  Pedestal_Tree->Write();
 
  TH2F* chip_map = new TH2F("chip_map", "Chip Map;X[mm];Y[mm]", 32, -90, 90, 32, -90, 90);
  TH2F* asic_map = new TH2F("asic_map", "Asic Map;X[mm];Y[mm]", 4, -90, 90, 4, -90, 90);
  for(Int_t chip=0; chip<16; chip++){
    asic_map->Fill(map_pointX[chip][0], map_pointY[chip][0], chip);
    for(Int_t channel=0; channel<64; channel++){
      chip_map->Fill(map_pointX[chip][channel], map_pointY[chip][channel], channel);
    }
  }

  /*TCanvas *canvas_pedestal_map = new TCanvas("pedestal_map", "pedestal_map", 1200, 1200);
  canvas_pedestal_map->Print(resultname + "_pedestal.pdf[", "pdf");
  Int_t imemory = 0;
  canvas_pedestal_map->Divide(2, 3);
  for(Int_t islab=0; islab<MaxSlab; islab++){
    canvas_pedestal_map->cd(islab+1);
    pedestal_map[islab]->SetStats(kFALSE);
    pedestal_map[islab]->SetTitle(TString::Format("pedestal_map_slab%i", islab));
    pedestal_map[islab]->GetXaxis()->SetTitle("x");
    pedestal_map[islab]->GetYaxis()->SetTitle("y");
    pedestal_map[islab]->GetZaxis()->SetRangeUser(200, 400);
    //noise.at(ichip)->SetLineColor(2);
    pedestal_map[islab]->Draw("colz");
    chip_map->Draw("text same");
    asic_map->Draw("text same");
    TLine *border[2][3];
    for(Int_t i=0; i<2; i++){ 
      for(Int_t j=0; j<3; j++){ 
        if(i==0)
 	  border[i][j] = new TLine(-45+45*j,-90,-45+45*j,90);
        if(i==1)
	  border[i][j] = new TLine(-90,-45+45*j,90,-45+45*j);
        border[i][j]->SetLineColor(kGray);
        border[i][j]->SetLineWidth(1);
        border[i][j]->SetLineStyle(2);
        border[i][j]->Draw("same");
      }
    }
    if(1){
      TLatex latex;
      latex.SetTextAlign(22);  //centered
      latex.SetTextColor(17); 
      for(Int_t i=0; i<16; i++){
        if(i<8)
	  latex.DrawLatex(67.5-45*((Int_t)i/2),-22.5-45*((i+1)%2),Form("#bf{#scale[2]{%d}}",i));
        else
	  latex.DrawLatex(67.5-45*((Int_t)(i-8)/2),67.5-45*((i+1)%2),Form("#bf{#scale[2]{%d}}",i));
      }
    }
    pedestal_map[islab]->Write();
  }
  canvas_pedestal_map->Print(resultname + "_pedestal.pdf", "pdf");

  for(Int_t islab=0; islab<MaxSlab; islab++){
    canvas_pedestal_map->cd(islab+1);
    pedestal_width_map[islab]->SetStats(kFALSE);
    pedestal_width_map[islab]->SetTitle(TString::Format("pedestal_width_map_slab%i", islab));
    pedestal_map[islab]->GetXaxis()->SetTitle("x");
    pedestal_width_map[islab]->GetXaxis()->SetTitle("x");
    pedestal_width_map[islab]->GetYaxis()->SetTitle("y");
    pedestal_width_map[islab]->GetZaxis()->SetRangeUser(0,7);
    //noise.at(ichip)->SetLineColor(2);
    pedestal_width_map[islab]->Draw("colz");
    chip_map->Draw("text same");
    asic_map->Draw("text same");
    TLine *border[2][3];
    for(Int_t i=0; i<2; i++){ 
      for(Int_t j=0; j<3; j++){ 
        if(i==0)
 	  border[i][j] = new TLine(-45+45*j,-90,-45+45*j,90);
        if(i==1)
	  border[i][j] = new TLine(-90,-45+45*j,90,-45+45*j);
        border[i][j]->SetLineColor(kGray);
        border[i][j]->SetLineWidth(1);
        border[i][j]->SetLineStyle(2);
        border[i][j]->Draw("same");
      }
    }
    if(1){
      TLatex latex;
      latex.SetTextAlign(22);  //centered
      latex.SetTextColor(17); 
      for(Int_t i=0; i<16; i++){
        if(i<8)
	  latex.DrawLatex(67.5-45*((Int_t)i/2),-22.5-45*((i+1)%2),Form("#bf{#scale[2]{%d}}",i));
        else
	  latex.DrawLatex(67.5-45*((Int_t)(i-8)/2),67.5-45*((i+1)%2),Form("#bf{#scale[2]{%d}}",i));
      }
    } 
    pedestal_width_map[islab]->Write();
  }
  canvas_pedestal_map->Print(resultname + "_pedestal.pdf", "pdf");

  for(Int_t islab=0; islab<MaxSlab; islab++){
    canvas_pedestal_map->cd(islab+1);
    pedestal_error_map[islab]->SetStats(kFALSE);
    pedestal_error_map[islab]->SetTitle(TString::Format("pedestal_error_map_slab%i", islab));
    pedestal_map[islab]->GetXaxis()->SetTitle("x");
    pedestal_error_map[islab]->GetXaxis()->SetTitle("x");
    pedestal_error_map[islab]->GetYaxis()->SetTitle("y");
    pedestal_error_map[islab]->GetZaxis()->SetRangeUser(0,0.1);
    //noise.at(ichip)->SetLineColor(2);
    pedestal_error_map[islab]->Draw("colz");
    chip_map->Draw("text same");
    asic_map->Draw("text same");
    TLine *border[2][3];
    for(Int_t i=0; i<2; i++){ 
      for(Int_t j=0; j<3; j++){ 
        if(i==0)
 	  border[i][j] = new TLine(-45+45*j,-90,-45+45*j,90);
        if(i==1)
	  border[i][j] = new TLine(-90,-45+45*j,90,-45+45*j);
        border[i][j]->SetLineColor(kGray);
        border[i][j]->SetLineWidth(1);
        border[i][j]->SetLineStyle(2);
        border[i][j]->Draw("same");
      }
    }
    if(1){
      TLatex latex;
      latex.SetTextAlign(22);  //centered
      latex.SetTextColor(17); 
      for(Int_t i=0; i<16; i++){
        if(i<8)
	  latex.DrawLatex(67.5-45*((Int_t)i/2),-22.5-45*((i+1)%2),Form("#bf{#scale[2]{%d}}",i));
        else
	  latex.DrawLatex(67.5-45*((Int_t)(i-8)/2),67.5-45*((i+1)%2),Form("#bf{#scale[2]{%d}}",i));
      }
    } 
    pedestal_error_map[islab]->Write();
  }
  canvas_pedestal_map->Print(resultname + "_pedestal.pdf", "pdf");
  canvas_pedestal_map->Print(resultname + "_pedestal.pdf]", "pdf");*/

  fout->Close();
}



