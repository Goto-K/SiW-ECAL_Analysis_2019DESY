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

using namespace::std;

// version 1: caluculating all memory cell pedestal of each ASIC ch
//            and Fitting gaus
//            and Get parameter
//            and Fill to Tree
void ped_calib_fill(std::string str){

  // Cut File name
  int path_i = str.find_last_of("/")+1;
  int ext_i = str.find_last_of(".");
  std::string pathname = str.substr(0, path_i);
  std::string extroot = str.substr(ext_i, str.size()-ext_i);
  std::string filename_i = str.substr(path_i, ext_i-path_i);
  int ext_j = filename_i.find_last_of(".");
  std::string extraw = filename_i.substr(ext_j, filename_i.size()-ext_j);
  std::string true_filename = filename_i.substr(0, ext_j);

  TString fullpath = str.data();
  TString filename = pathname + true_filename;

  // TTree open  
  TFile *file = TFile::Open(fullpath); 
  TTree *fev10 = (TTree*)file->Get("fev10");

  // Read Chip Map
  TString mapchip_path = "map_chip.dat";
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
    map_pointY[tmp_chip][tmp_channel] = -tmp_y;
  }
  
  Int_t MaxEvent = fev10->GetEntries();
  Int_t MaxChip = 16;
  Int_t MaxChannel = 64;
  Int_t MaxSca = 15;

  TH2F* pedestal_map[15];
  TH2F* pedestal_width_map[15];
  TH2F* pedestal_error_map[15];
  TH2F* pedestal_chi2ndf_map[15];
  TH2F* pedestal_npeaks_map[15];
  TH2F* pedestal_entries_map[15];

  for(int isca=0; isca<MaxSca; isca++) {
    pedestal_map[isca]= new TH2F(TString::Format("pedestal_map_sca%i",isca),TString::Format("pedestal_map_sca%i;X[mm];Y[mm]",isca),32,-90,90,32,-90,90);
    pedestal_width_map[isca]= new TH2F(TString::Format("pedestal_width_map_sca%i",isca),TString::Format("pedestal_width_map_sca%i;X[mm];Y[mm]",isca),32,-90,90,32,-90,90);
    pedestal_error_map[isca]= new TH2F(TString::Format("pedestal_error_map_sca%i",isca),TString::Format("pedestal_error_map_sca%i;X[mm];Y[mm]",isca),32,-90,90,32,-90,90);
    pedestal_chi2ndf_map[isca]= new TH2F(TString::Format("pedestal_chi2ndf_map_sca%i",isca),TString::Format("pedestal_chi2ndf_map_sca%i;X[mm];Y[mm]",isca),32,-90,90,32,-90,90);
    pedestal_npeaks_map[isca]= new TH2F(TString::Format("pedestal_npeaks_map_sca%i",isca),TString::Format("pedestal_npeaks_map_sca%i;X[mm];Y[mm]",isca),32,-90,90,32,-90,90);
    pedestal_entries_map[isca]= new TH2F(TString::Format("pedestal_entries_map_sca%i",isca),TString::Format("pedestal_entries_map_sca%i;X[mm];Y[mm]",isca),32,-90,90,32,-90,90);
  }

  // Read data from fev10
  Int_t charge_hiGain[16][15][64];
  //fev10->SetBranchAddress("charge_hiGain", charge_hiGain);
  fev10->SetBranchAddress("charge_lowGain", charge_hiGain);
  Int_t gain_hit_high[16][15][64];
  //fev10->SetBranchAddress("gain_hit_high", gain_hit_high);
  fev10->SetBranchAddress("gain_hit_low", gain_hit_high);

  Int_t badbcid[16][15];
  fev10->SetBranchAddress("badbcid", badbcid);

  // Fill to Tree Data
  Double_t pedestal_mean[16][15][64] = {-1};
  Double_t pedestal_width[16][15][64] = {-1};
  Double_t pedestal_error[16][15][64] = {-1};
  Double_t pedestal_chi2ndf[16][15][64] = {-1};
  Double_t charge_hiGain_ped_calib[16][15][64] = {-1};

  // Create New TTree
  TFile *fout = new TFile(filename + "_PedestalMap.root" , "RECREATE");
  TTree *Pedestal_Tree = new TTree("Pedestal_Tree", "Pedestal_Tree");
//  Pedestal_Tree->Branch("gain_hit_high", gain_hit_high,
//		        "gain_hit_high[16][15][64]/I");
//  Pedestal_Tree->Branch("badbcid", badbcid,
//		        "badbcid[16][15]/I");
  Pedestal_Tree->Branch("pedestal_mean", pedestal_mean,
		        "pedestal_mean[16][15][64]/D");
  Pedestal_Tree->Branch("pedestal_width", pedestal_width,
		        "pedestal_width[16][15][64]/D");
  Pedestal_Tree->Branch("pedestal_error", pedestal_error,
		        "pedestal_error[16][15][64]/D");
  Pedestal_Tree->Branch("pedestal_chi2ndf", pedestal_chi2ndf,
		        "pedestal_chi2ndf[16][15][64]/D");
//  Pedestal_Tree->Branch("charge_hiGain_ped_calib", charge_hiGain_ped_calib,
//             	        "charge_hiGain_ped_calib[16][15][64]/D");

  TDirectory *hist_sca[16];

 
  // --------------------
  // all sca
  //  std::vector<TCanvas*> chip;
  std::vector<std::vector<std::vector<TH1F*> > > ped_sca;

  std::vector<TH1F*> pedestal_chip ;
  std::vector<TH1F*> pedestal_diff_chip ;

  for(int ichip=0; ichip<MaxChip; ichip++) {
    TH1F *ped_chip = new TH1F(TString::Format("ped_chip%i", ichip), 
		    	      TString::Format("ped_chip%i", ichip), 1000, 0.5, 1000.5);
    pedestal_chip.push_back(ped_chip);

    TH1F *ped_diff_chip = new TH1F(TString::Format("ped_diff_chip%i", ichip),
		                   TString::Format("ped_diff_chip%i", ichip), 1002, -500, 500);
    pedestal_diff_chip.push_back(ped_diff_chip);
    
    std::vector<std::vector<TH1F*> >pedtemp_sca;

    for(int ichn=0; ichn<MaxChannel; ichn++) {
      std::vector<TH1F*> pedtemp_sca2;

      for(int isca=0; isca<MaxSca; isca++) {
	TH1F *ped_sca2 = new TH1F(TString::Format("ped_chip%i_chn%i_sca%i", ichip, ichn, isca),
			          TString::Format("ped_chip%i_chn%i_sca%i", ichip, ichn, isca), 1000, 0.5, 1000.5);
	pedtemp_sca2.push_back(ped_sca2);
      }
      pedtemp_sca.push_back(pedtemp_sca2);
    }
    ped_sca.push_back(pedtemp_sca);
  }

  // -----------------------------------------------------------------------------------------------------   
  // SCA analysis
  for (int ievent=0; ievent<MaxEvent; ievent++) {
    fev10->GetEntry(ievent);

    for(int ichip=0; ichip<MaxChip; ichip++) {

      for(int isca=0; isca<MaxSca; isca++) {

	for(int ichn=0; ichn<MaxChannel; ichn++) {

	  //if(charge_hiGain[ichip][isca][ichn]>10 && badbcid[ichip][isca]==0){
	  if(charge_hiGain[ichip][isca][ichn]>10 && gain_hit_high[ichip][isca][ichn]==0){
	    ped_sca.at(ichip).at(ichn).at(isca)->Fill(charge_hiGain[ichip][isca][ichn]);
	  }

	} 

      }//isca

    }//ichip 
   
  }  // end first loop analysis to fill pedestal historgrams


  std::vector<std::vector<std::vector<Double_t> > > ped_mean;
  std::vector<std::vector<std::vector<Double_t> > > ped_error;
  std::vector<std::vector<std::vector<Double_t> > > ped_width;

  //initialize pedestal vectors
  for(int i=0; i<16; i++) {
    std::vector<std::vector<Double_t> > chip_ped_mean;
    std::vector<std::vector<Double_t> > chip_ped_error;
    std::vector<std::vector<Double_t> > chip_ped_width;

    for(int j=0; j<64; j++) {
      std::vector<Double_t> chn_ped_mean;
      std::vector<Double_t> chn_ped_error;
      std::vector<Double_t> chn_ped_width;
      chip_ped_mean.push_back(chn_ped_mean);
      chip_ped_error.push_back(chn_ped_error);
      chip_ped_width.push_back(chn_ped_width);
    }  
    ped_mean.push_back(chip_ped_mean);
    ped_error.push_back(chip_ped_error);
    ped_width.push_back(chip_ped_width);
  }

  // do pedestal (chip/channel/sca based) analysis
  for(int ichip=0; ichip<MaxChip; ichip++) {
    cout << "CHIP " << ichip << " START" << endl;
    TString TDirectory_Name = Form("ALL_PEDESTAL_%02dCHIP", ichip);
    hist_sca[ichip] = fout->mkdir(TDirectory_Name);
    hist_sca[ichip]->cd();
    for(int ichn=0; ichn<MaxChannel; ichn++) {
      for(int isca=0; isca<MaxSca; isca++) {

	ped_sca.at(ichip).at(ichn).at(isca)->SetTitle(TString::Format("ped_chip%i_chn%i_sca%i",ichip,ichn,isca));
	ped_sca.at(ichip).at(ichn).at(isca)->SetName(TString::Format("ped_chip%i_chn%i_sca%i",ichip,ichn,isca));
	ped_sca.at(ichip).at(ichn).at(isca)->Write();

	pedestal_entries_map[isca]->Fill(map_pointX[ichip][ichn], map_pointY[ichip][ichn], 
					   ped_sca.at(ichip).at(ichn).at(isca)->GetEntries());
	
	ped_mean.at(ichip).at(ichn).push_back(0.);
	ped_error.at(ichip).at(ichn).push_back(0.);
	ped_width.at(ichip).at(ichn).push_back(0.);

	if(ped_sca.at(ichip).at(ichn).at(isca)->GetEntries()>10){ //max_entries/2 ) {
	  TSpectrum *s = new TSpectrum();
	  int npeaks = s->Search(ped_sca.at(ichip).at(ichn).at(isca), 2, "", 0.2); 
	  pedestal_npeaks_map[isca]->Fill(map_pointX[ichip][ichn], map_pointY[ichip][ichn], npeaks);

	  if(npeaks > 0) {

            Double_t *mean_peak=s->GetPositionX();
            Double_t *mean_high=s->GetPositionY();
            double mean_peak_higher=0;
            double mean_high_higher=0;
	    int npeak_max=0;

            for(int ipeak=0; ipeak<npeaks; ipeak++) {
              if(mean_high[ipeak]>mean_high_higher && mean_high[ipeak]>50) {
                mean_high_higher=mean_high[ipeak];
                mean_peak_higher=mean_peak[ipeak];
		npeak_max=ipeak;
              }
            }

            for(int ipeak=0; ipeak<npeaks; ipeak++) {
	      if(ipeak != npeak_max) pedestal_diff_chip.at(ichip)->Fill(mean_peak[npeak_max] - mean_peak[ipeak]);
	    }

	  }

	  //	  if(npeaks == 1) {

	  //	    Double_t *mean_peak=s->GetPositionX();
	  TH1 *pedestal_histogram = ped_sca.at(ichip).at(ichn).at(isca);
	  Double_t mean_peak = pedestal_histogram->GetBinCenter(pedestal_histogram->GetMaximumBin());
	    
	    //	    TF1 *f0 = new TF1("f0", "gaus", mean_peak[0]-2*ped_sca.at(ichip).at(ichn).at(isca)->GetRMS(), 
	    //			      mean_peak[0]+2*ped_sca.at(ichip).at(ichn).at(isca)->GetRMS());
	    TF1 *f0 = new TF1("f0", "gaus", mean_peak-2*ped_sca.at(ichip).at(ichn).at(isca)->GetRMS(), 
			      mean_peak+2*ped_sca.at(ichip).at(ichn).at(isca)->GetRMS());
	    ped_sca.at(ichip).at(ichn).at(isca)->Fit("f0", "RQNOC");

	    TF1 *f1 = new TF1("f1", "gaus", f0->GetParameter(1)-2.*f0->GetParameter(2), 
			      f0->GetParameter(1)+2.*f0->GetParameter(2));
	    ped_sca.at(ichip).at(ichn).at(isca)->Fit("f1","RQME");

	    ped_mean.at(ichip).at(ichn).at(isca)=f1->GetParameter(1);
	    ped_error.at(ichip).at(ichn).at(isca)=f1->GetParError(1);
	    ped_width.at(ichip).at(ichn).at(isca)=f1->GetParameter(2);

	    pedestal_chip.at(ichip)->Fill(f1->GetParameter(1));

	    pedestal_map[isca]->Fill(map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , f1->GetParameter(1));
	    pedestal_width_map[isca]->Fill(map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , f1->GetParameter(2));
	    pedestal_error_map[isca]->Fill(map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , f1->GetParError(1));
	    pedestal_chi2ndf_map[isca]->Fill(map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , f1->GetChisquare() / f1->GetNDF());
	      
	    pedestal_mean[ichip][isca][ichn] = f1->GetParameter(1);
	    pedestal_width[ichip][isca][ichn] = f1->GetParameter(2);
	    pedestal_error[ichip][isca][ichn] = f1->GetParError(1);
	    pedestal_chi2ndf[ichip][isca][ichn] = f1->GetChisquare() / f1->GetNDF();
	    //	  }

	}
      }
    }
    cout << "CHIP " << ichip << " END" << endl;
  }

  fout->cd();
  
  Pedestal_Tree->Fill();
  Pedestal_Tree->Write();
  
  // good pedestal events (not tagged events)
  TCanvas *canvas_pedestal_map = new TCanvas("pedestal_map", "pedestal_map", 1200, 1200);
  canvas_pedestal_map->Divide(4, 4);
  for(int isca=0; isca<MaxSca; isca++){
    canvas_pedestal_map->cd(isca+1);
    pedestal_map[isca]->SetStats(kFALSE);
    pedestal_map[isca]->SetTitle(TString::Format("pedestal_map_sca%i", isca));
    pedestal_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_map[isca]->GetYaxis()->SetTitle("y");
    pedestal_map[isca]->GetZaxis()->SetRangeUser(200, 400);
    //noise.at(ichip)->SetLineColor(2);
    pedestal_map[isca]->Draw("colz");
    pedestal_map[isca]->Write();
  }

  TCanvas *canvas_pedestal_width_map = new TCanvas("pedestal_width_map", "pedestal_width_map", 1200, 1200);
  canvas_pedestal_width_map->Divide(4,4);
  for(int isca=0; isca<MaxSca; isca++){
    canvas_pedestal_width_map->cd(isca+1);
    pedestal_width_map[isca]->SetStats(kFALSE);
    pedestal_width_map[isca]->SetTitle(TString::Format("pedestal_width_map_sca%i", isca));
    pedestal_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_width_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_width_map[isca]->GetYaxis()->SetTitle("y");
    pedestal_width_map[isca]->GetZaxis()->SetRangeUser(0,7);
    //noise.at(ichip)->SetLineColor(2);
    pedestal_width_map[isca]->Draw("colz");
    pedestal_width_map[isca]->Write();
  }

  TCanvas *canvas_pedestal_error_map = new TCanvas("pedestal_error_map", "pedestal_error_map", 1200, 1200);
  canvas_pedestal_error_map->Divide(4,4);
  for(int isca=0; isca<MaxSca; isca++){
    canvas_pedestal_error_map->cd(isca+1);
    pedestal_error_map[isca]->SetStats(kFALSE);
    pedestal_error_map[isca]->SetTitle(TString::Format("pedestal_error_map_sca%i", isca));
    pedestal_error_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_error_map[isca]->GetYaxis()->SetTitle("y");
    //noise.at(ichip)->SetLineColor(2);
    pedestal_error_map[isca]->Draw("colz");
    pedestal_error_map[isca]->Write();
  }

  TCanvas *canvas_pedestal_npeaks_map = new TCanvas("pedestal_npeaks_map", "pedestal_npeaks_map", 1200, 1200);
  canvas_pedestal_npeaks_map->Divide(4,4);
  for(int isca=0; isca<MaxSca; isca++){
    canvas_pedestal_npeaks_map->cd(isca+1);
    pedestal_npeaks_map[isca]->SetStats(kFALSE);
    pedestal_npeaks_map[isca]->SetTitle(TString::Format("pedestal_npeaks_map_sca%i", isca));
    pedestal_npeaks_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_npeaks_map[isca]->GetYaxis()->SetTitle("y");
    pedestal_npeaks_map[isca]->GetZaxis()->SetRangeUser(0.5, 4.5);
    //noise.at(ichip)->SetLineColor(2);
    pedestal_npeaks_map[isca]->Draw("colz");
    pedestal_npeaks_map[isca]->Write();
  }

  TCanvas *canvas_pedestal_chi2ndf_map = new TCanvas("pedestal_chi2ndf_map", "pedestal_chi2ndf_map", 1200, 1200);
  canvas_pedestal_chi2ndf_map->Divide(4,4);
  for(int isca=0; isca<MaxSca; isca++){
    canvas_pedestal_chi2ndf_map->cd(isca+1);
    pedestal_chi2ndf_map[isca]->SetStats(kFALSE);
    pedestal_chi2ndf_map[isca]->SetTitle(TString::Format("pedestal_chi2ndf_map_sca%i", isca));
    pedestal_chi2ndf_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_chi2ndf_map[isca]->GetYaxis()->SetTitle("y");
    //noise.at(ichip)->SetLineColor(2);
    pedestal_chi2ndf_map[isca]->Draw("colz");
    pedestal_chi2ndf_map[isca]->Write();
  }

  TCanvas *canvas_pedestal_entries_map = new TCanvas("pedestal_entries_map", "pedestal_entries_map", 1200, 1200);
  canvas_pedestal_entries_map->Divide(4,4);
  for(int isca=0; isca<MaxSca; isca++){
    canvas_pedestal_entries_map->cd(isca+1);
    pedestal_entries_map[isca]->SetStats(kFALSE);
    pedestal_entries_map[isca]->SetTitle(TString::Format("pedestal_entries_map_sca%i", isca));
    pedestal_entries_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_entries_map[isca]->GetYaxis()->SetTitle("y");
    //noise.at(ichip)->SetLineColor(2);
    pedestal_entries_map[isca]->Draw("colz");
    pedestal_entries_map[isca]->Write();
  }
  

  canvas_pedestal_map->Write();
  canvas_pedestal_width_map->Write();
  canvas_pedestal_error_map->Write();
  canvas_pedestal_npeaks_map->Write();
  canvas_pedestal_entries_map->Write();
  canvas_pedestal_chi2ndf_map->Write();

   
  TCanvas *canvas_pedestal_pedestal = new TCanvas("pedestal_average", "pedestal_average", 1200, 1200);
  canvas_pedestal_pedestal->Divide(4,4);
  for(int ichip=0; ichip<MaxChip; ichip++) {
    canvas_pedestal_pedestal->cd(ichip+1);
    //gPad->SetLogy();
    pedestal_chip.at(ichip)->GetXaxis()->SetRangeUser(0, 500);
    pedestal_chip.at(ichip)->SetTitle(TString::Format("Average pedestal, chip-%i", ichip));
    pedestal_chip.at(ichip)->GetXaxis()->SetTitle("ADC");
    pedestal_chip.at(ichip)->GetYaxis()->SetTitle("#");
    //noise.at(ichip)->SetLineColor(2);
    pedestal_chip.at(ichip)->Draw("hs");
    pedestal_chip.at(ichip)->Write();
  }

  canvas_pedestal_pedestal->Write();

  TCanvas *canvas_pedestal_diff = new TCanvas("pedestal_diff", "pedestal_diff", 1200, 1200);
  canvas_pedestal_diff->Divide(4,4);
  for(int ichip=0; ichip<MaxChip; ichip++) {
    canvas_pedestal_diff->cd(ichip+1);
    pedestal_diff_chip.at(ichip)->GetXaxis()->SetRangeUser(-100, 100);
    pedestal_diff_chip.at(ichip)->SetTitle(TString::Format("Pedestal diff, chip-%i", ichip));
    pedestal_diff_chip.at(ichip)->GetXaxis()->SetTitle("ADC");
    pedestal_diff_chip.at(ichip)->GetYaxis()->SetTitle("#");
    pedestal_diff_chip.at(ichip)->Draw("hs");
  }

  canvas_pedestal_diff->Write();

  fout->Close();

}


