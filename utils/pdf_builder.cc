/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : pdf_builder.cc
 * @created     : Monday Dec 16, 2024 18:43:19 CET
 */

#include <iostream>
#include <fstream>  
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "THnBase.h"
#include "THn.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TROOT.h"
#include "TF1.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TCanvas.h"

THnD* build_response_matrix(TTree* t, const TH1* hEnu, const int ecal_nb, const double ecal_x0, const double ecal_x1)
{
  TString name = Form("%s_response", hEnu->GetName());
  TString titl = Form("%s response matrix", hEnu->GetTitle());

  int nbins[2] = {hEnu->GetNbinsX(), ecal_nb}; 
  double xmin[2] = {hEnu->GetXaxis()->GetXmin(), ecal_x0};
  double xmax[2] = {hEnu->GetXaxis()->GetXmax(), ecal_x1};

  THnD* hResponse = new THnD(name, titl, 2, nbins, xmin, xmax);

  TTreeReader reader(t);
  TTreeReaderValue<double> enu(reader, "E");
  TTreeReaderValue<double> erec(reader, "Ecal");

  while (reader.Next()) {
    double point[2] = {*enu, *erec};
    hResponse->Fill(point);
  }

  for (int i=1; i<=nbins[0]; i++) {
    hResponse->GetAxis(0)->SetRange(i, i);
    TH1D* h = hResponse->Projection(1, "h");
    const double integral = h->Integral();
    double norm = 0; 
    if (integral > 0) norm = 1.0 / integral;
    for (int j=1; j<=nbins[1]; j++) {
      int ibin[2] = {i, j};
      hResponse->SetBinContent(ibin, h->GetBinContent(j) * norm);
    }
    delete h;
  }

  hResponse->GetAxis(0)->SetRange(0, 0); 

  return hResponse;
}

int neutrino_pdf_builder(const TString& spectrum_file = "./b8-neutrino.txt")
{
  TH1::AddDirectory(kFALSE);
  std::ifstream spectrum_file_stream(spectrum_file);
  if (!spectrum_file_stream.is_open()) {
    std::cerr << "Error: cannot open file " << spectrum_file << std::endl;
    return 1;
  }

  std::string line;
  TGraph* g_nu = new TGraph();
  while (std::getline(spectrum_file_stream, line)) {
    std::istringstream iss(line);
    double energy, flux, err;
    iss >> energy >> flux;
    g_nu->SetPoint(g_nu->GetN(), energy, flux);
  }

  TF1* fB8 = new TF1("fB8", [g_nu](double* x, double* p) { 
      double y = g_nu->Eval(x[0]); 
      if ( y > 0.0 ) return y; 
      else return 0.0;
      }, 0, 20, 0);

  TH1D* hb8 = new TH1D("hb8", "hb8", 100, 0, 20);
  for (int i = 0; i < 100; ++i) {
    double x0 = hb8->GetBinLowEdge(i+1);
    double x1 = hb8->GetBinLowEdge(i+2);
    hb8->SetBinContent(i+1, fB8->Integral(x0, x1) / (x1 - x0));
  }

  // build response matrix
  TFile* f_flat = new TFile("energy_marley_v1.root"); 
  TTree* t = f_flat->Get<TTree>("energy_tree");
  THnD* h2 = build_response_matrix(t, hb8, 42, -1.0, 20);
  f_flat->Close();

  TH1D* hErec = h2->Projection(1, "hErec"); 
  hErec->SetTitle("Reconstructed energy");
  hErec->Reset(); 
  for (int i=1; i<h2->GetAxis(0)->GetNbins(); i++) {
    h2->GetAxis(0)->SetRange(i, i); 
    TH1D* h = h2->Projection(1, "h");
    const double norm = hb8->GetBinContent(i);
    hErec->Add(h, norm);
    delete h;
  }
  h2->GetAxis(0)->SetRange(0, 0); 

  TCanvas* c_neutrino = new TCanvas("c_neutrino", "c_neutrino", 1200, 600); 
  c_neutrino->Divide(3, 1);
  c_neutrino->cd(1);
  hb8->Draw();
  c_neutrino->cd(2);
  h2->Projection(1, 0)->Draw("colz");
  c_neutrino->cd(3);
  hErec->Draw("hist");
  
  // Experiment with response matrix applied to 2D (E_nu, cos_nadir) distribution
  TFile* file_nadir = new TFile("./nadir.root"); 
  TH1D* h_nadir = file_nadir->Get<TH1D>("cos_nadir"); 
  h_nadir->Rebin(2); 

  int nbins[2]   = {hb8->GetNbinsX(), h_nadir->GetNbinsX()}; 
  double xmin[2] = {hb8->GetXaxis()->GetXmin(), h_nadir->GetXaxis()->GetXmin()};
  double xmax[2] = {hb8->GetXaxis()->GetXmax(), h_nadir->GetXaxis()->GetXmax()};

  THnD* hnTarget = new THnD("hnTarget", "target", 2, nbins, xmin, xmax); 
  Long64_t i = 0; 
  auto* iter = hnTarget->CreateIter(false); 
  int coords[2]; 

  while ( (i = iter->Next(coords)) >=0 )
  {
    double flux = hb8->GetBinContent(coords[0]); 
    double cnad = h_nadir->GetBinContent(coords[1]); 
    hnTarget->SetBinContent(i, flux * cnad); 
  }

  TCanvas* c2D = new TCanvas(); 
  c2D->Divide(2, 1); 
  c2D->cd(1); 
  TH2D* h2Target = hnTarget->Projection(1, 0);
  h2Target->SetEntries(nbins[0] * nbins[1]); 
  h2Target->Draw("colz"); 

  //TFile* boron_file = new TFile("boron_pdf.root", "RECREATE");
  //hb8->Write();
  //h2->Write("nueCCresponseMatrix");
  //boron_file->Close();
  return 0;
}

int background_pdf_builder() 
{
  TFile* file_nadir = new TFile("./nadir.root", "update");
  TH1F* h_cnadir = file_nadir->Get<TH1F>("nadir");

  TF1* f_nadir = new TF1("f_nadir", 
      [h_cnadir](double* x, double* p) { 
      double y = h_cnadir->Interpolate(x[0]); 
      if ( y > 0.0 ) return y; 
      else return 0.0;
      }, -1.0, 1.0, 0);

  TH1D* h_cnadir_rb = new TH1D("cos_nadir", "cos_nadir_surf", 100, -1, 1); 
  for (int ix=1; ix<=h_cnadir_rb->GetNbinsX(); ix++) {
    double x0 = h_cnadir_rb->GetBinLowEdge(ix); 
    double x1 = h_cnadir_rb->GetBinLowEdge(ix+1); 
    h_cnadir_rb->SetBinContent(ix, f_nadir->Integral(x0, x1, 1e-4) ); 
  }
  h_cnadir_rb->Write( h_cnadir_rb->GetName(), TFile::kOverwrite); 


  // background energy spectrum
  TF1* fgg = new TF1("fgg", [](double* x, double* p) { 
      double x0 = x[0];
      double A0 = p[0];
      double mu0 = p[1];
      double sigma0 = p[2];
      double A1 = p[3];
      double mu1 = p[4];
      double sigma1 = p[5];

      if (x0 < 5.0) return 0.0; 
      double y0 = A0 * TMath::Gaus(x0, mu0, sigma0);
      double y1 = A1 * TMath::Gaus(x0, mu1, sigma1);
      return (y0+y1);
      }, 
      -1, 20.0, 6);

  fgg->SetParameters(20.0, 6, 1, 7.0, 8.5, 1.0);

  const int nbins[2] = {42, 20};
  const double xmin[2] = {-1.0, -1.0};
  const double xmax[2] = {20.0,  1.0};

  THnD* hbkg = new THnD("hbkg", "hbkg", 2, nbins, xmin, xmax);

  for (int iy=1;  iy<=nbins[1]; iy++) {
    double nadir_x0 = hbkg->GetAxis(1)->GetBinLowEdge(iy);
    double nadir_x1 = hbkg->GetAxis(1)->GetBinUpEdge(iy);
    double nadir = f_nadir->Integral(nadir_x0, nadir_x1, 1e-3); 
    //printf("nadir range: %f %f nadir=%f\n", nadir_x0, nadir_x1, nadir);

    for (int ix=1; ix<=nbins[0]; ix++) {
      double x0 = hbkg->GetAxis(0)->GetBinLowEdge(ix);
      double x1 = hbkg->GetAxis(0)->GetBinUpEdge(ix);
      double bkg = fgg->Integral(x0, x1, 1e-3) * nadir;
      const int ibin[2] = {ix, iy};
      //printf("\tx0=%f x1=%f nadir=%f bkg=%f\n", x0, x1, nadir, bkg);
      hbkg->SetBinContent(ibin, bkg);
    }
  }
  hbkg->SetEntries( 42*20 );

  TFile* bkg_file = new TFile("background_pdf.root", "RECREATE");
  hbkg->Write("neutron");
  bkg_file->Close();

  file_nadir->Close(); 

  return 0;
}


