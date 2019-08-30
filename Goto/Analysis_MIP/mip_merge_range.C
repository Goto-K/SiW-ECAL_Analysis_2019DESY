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

// ================================================================================================================================================ //
// ================================================================================================================================================ //
 
Double_t langaufun(Double_t *x, Double_t *par) {
     Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
     Double_t mpshift  = -0.22278298;       // Landau maximum location
     Double_t np = 100.0;      // number of convolution steps
     Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
     Double_t xx;
     Double_t mpc;
     Double_t fland;
     Double_t sum = 0.0;
     Double_t xlow,xupp;
     Double_t step;
     Double_t i;
     mpc = par[1] - mpshift * par[0];
     xlow = x[0] - sc * par[3];
     xupp = x[0] + sc * par[3];
     step = (xupp-xlow) / np;
     for(i=1.0; i<=np/2; i++) {
        xx = xlow + (i-.5) * step;
        fland = TMath::Landau(xx,mpc,par[0]) / par[0];
        sum += fland * TMath::Gaus(x[0],xx,par[3]);
        xx = xupp - (i-.5) * step;
        fland = TMath::Landau(xx,mpc,par[0]) / par[0];
        sum += fland * TMath::Gaus(x[0],xx,par[3]);
     }
     return (par[2] * step * sum * invsq2pi / par[3]);
}


TF1 *langaufit(TH1F *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF)
{
  Int_t i;
  Char_t FunName[100];
  sprintf(FunName,"Fitfcn_%s",his->GetName());
  TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
  if (ffitold) delete ffitold;
  TF1 *ffit = new TF1(FunName,langaufun,fitrange[0],fitrange[1],4);
  ffit->SetParameters(startvalues);
  ffit->SetParNames("Width","MP","Area","GSigma");
  for (i=0; i<4; i++) {
     ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
  }
  his->Fit(FunName,"RQOCB0");   // fit within specified range, use ParLimits, do not plot
  ffit->GetParameters(fitparams);    // obtain fit parameters
  for (i=0; i<4; i++) {
     fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
  }
  ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
  NDF[0] = ffit->GetNDF();           // obtain ndf
  return (ffit);              // return fit function
}


Int_t langaupro(Double_t *params, Double_t &maxx, Double_t &FWHM) {
  Double_t p, x = 0.0, fy, fxr, fxl;
  Double_t step;
  Double_t l, lold;
  Int_t i = 0;
  Int_t MAXCALLS = 10000;
  p = params[1] - 0.1 * params[0];
  step = 0.05 * params[0];
  lold = -2.0;
  l    = -1.0;
  while ( (l != lold) && (i < MAXCALLS) ) {
    i++;
    lold = l;
    x = p + step;
    l = langaufun(&x,params);
     if (l < lold)
        step = -step/10;
     p += step;
  }
  if (i == MAXCALLS)
    return (-1);
  maxx = x;
  fy = l/2;
  p = maxx + params[0];
  step = params[0];
  lold = -2.0;
  l    = -1e300;
  i    = 0;
  while ( (l != lold) && (i < MAXCALLS) ) {
    i++;
    lold = l;
    x = p + step;
    l = TMath::Abs(langaufun(&x,params) - fy);
    if (l > lold)
      step = -step/10;
      p += step;
  }
  if (i == MAXCALLS)
    return (-2);
  fxr = x;
  p = maxx - 0.5 * params[0];
  step = -params[0];
  lold = -2.0;
  l    = -1e300;
  i    = 0;
  while ( (l != lold) && (i < MAXCALLS) ) {
    i++;
    lold = l;
    x = p + step;
    l = TMath::Abs(langaufun(&x,params) - fy);
    if (l > lold)
      step = -step/10;
    p += step;
  }
  if (i == MAXCALLS)
    return (-3);
  fxl = x;
  FWHM = fxr - fxl;
  return (0);
}

// ================================================================================================================================================ //
// ================================================================================================================================================ //

void mip_merge_range(std::string str){

  // Cut File name
  int path_i = str.find_last_of("/")+1;
  int ext_i = str.find_last_of(".");
  std::string pathname = str.substr(0, path_i);
  std::string filename_ = str.substr(path_i, ext_i-path_i);

  TString fullpath = str.data();
  TString filename = pathname + filename_;
  TString resultname = "./Result/" + filename_;

  // TTree_open  
  TFile *file = TFile::Open(fullpath); 
  TTree *fev10 = (TTree*)file->Get("fev10");

  // Define
  Int_t MaxEvent = fev10->GetEntries();
  Int_t MaxSlab = 5;
  Int_t MaxChip = 16;
  Int_t MaxChannel = 64;
  Int_t MaxMemory = 15;

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

  TFile *fout = new TFile(filename + "_Mip.root" , "RECREATE");

  // Calculating number of hit
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
  }

  // mip hist initialized
  std::vector<TH1F*> mip_sca;
  for(int islab=0; islab<MaxSlab; islab++){
    Int_t ichip = most_hit_chips[islab];
    Int_t ichn = most_hit_channels[islab];
    TH1F *miptemp_sca2 = new TH1F(TString::Format("mip_slab%i_chip%i_chn%i", islab, ichip, ichn),
			          TString::Format("mip_slab%i_chip%i_chn%i", islab, ichip, ichn), 550, 0, 550);
    mip_sca.push_back(miptemp_sca2);
  }

  // Calculating number of nhits and mip
  for(Int_t ievent = 0; ievent<MaxEvent; ievent++){
    fev10->GetEntry(ievent);
    for(Int_t islab=0; islab<MaxSlab; islab++){
      Int_t ichip = most_hit_chips[islab];
      for(Int_t imemory=0; imemory<MaxMemory; imemory++){
	if(badbcid[islab][ichip][imemory]==1) continue;
	for(Int_t ichn=0; ichn<MaxChannel; ichn++){
	  if(gain_hit_high[islab][ichip][imemory][ichn]<0) continue;
	  nhits_number[islab][nhits[islab][ichip][imemory]]++;
	  if(ichn==most_hit_channels[islab] && charge_lowGain[islab][ichip][imemory][ichn]>10){
	    mip_sca.at(islab)->Fill(charge_lowGain[islab][ichip][imemory][ichn]-pedestal_mean_copy[islab][ichip][imemory][ichn]);
	  }
	}
      }
    }
  }

  fout->cd();
  // do mip analysis
  Double_t mip_mean[5] = {0};
  TCanvas *canvas_mip_sca = new TCanvas("mip_sca", "mip_sca");
  canvas_mip_sca->Divide(2, 3);
  for(Int_t islab=0; islab<MaxSlab; islab++){
    Int_t ichip = most_hit_chips[islab];
    Int_t ichn = most_hit_channels[islab];
    canvas_mip_sca->cd(islab+1);
    Double_t fr[2];
    Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];
    fr[0]=0.6*mip_sca.at(islab)->GetMean();
    fr[1]=2.0*mip_sca.at(islab)->GetMean();
    pllo[0]=0.0; 
    pllo[1]=mip_sca.at(islab)->GetMean() - 70;
    if(islab==2) pllo[1]=mip_sca.at(islab)->GetMean() - 40;
    pllo[2]=0.0; 
    pllo[3]=0.0;;
    plhi[0]=0.0; 
    plhi[1]=mip_sca.at(islab)->GetMean() + 10;
    plhi[2]=0.0; 
    plhi[3]=20;
    sv[0]=5.0;
    sv[1]=mip_sca.at(islab)->GetMean();
    sv[2]=10000000;
    sv[3]=10;
    Double_t chisqr;
    Int_t    ndf;
    Double_t SNRPeak, SNRFWHM;
    langaupro(fp,SNRPeak,SNRFWHM);

    mip_sca.at(islab)->SetTitle(TString::Format("mip_chip%i_chn%i",ichip,ichn));
    mip_sca.at(islab)->SetName(TString::Format("mip_chip%i_chn%i",ichip,ichn));
    mip_sca.at(islab)->Write();

    TF1 *fitsnr = langaufit(mip_sca.at(islab), fr, sv, pllo, plhi, fp, fpe, &chisqr, &ndf);

    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(111);
    mip_sca.at(islab)->Draw();
    fitsnr->Draw("lsame");
    //gPad->SetLogy(1);
    gPad->Modified();
    gPad->Update();
    gPad->Write();
    mip_mean[islab] = fp[1];
  }
  canvas_mip_sca->Print(filename + "_mip.pdf", "pdf");
 
  fout->Close();
}
