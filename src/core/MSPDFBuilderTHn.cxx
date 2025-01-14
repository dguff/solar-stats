// Copyright (C) 2018 Matteo Agostini <matteo.agostini@ph.tum.de>

// This is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 2.1 of the License, or
// (at your option) any later version.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

// c/c++ libs
#include <iostream>

// ROOT libs
#include "TFile.h"
#include "THnBase.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TF1.h"
#include "TTimer.h"
#include "TROOT.h"

// m-stats libs
#include "MSPDFBuilderTHn.h"

// Prob3++
#include "BargerPropagator.h"

namespace mst {

MSPDFBuilderTHn::MSPDFBuilderTHn(const std::string& name): MSObject(name)
{
  fPDFMap = new PDFMap;
  fRespMatrixMap = new RespMatrixMap;
}

MSPDFBuilderTHn::~MSPDFBuilderTHn()
{
  if (fPDFMap != nullptr) {
    for (auto& it : *fPDFMap) {
      if (it.second) delete it.second;
    }
    fPDFMap->clear();
    delete fPDFMap;
  }

  if (fRespMatrixMap) {
    for (auto& it : *fRespMatrixMap) {
      if (it.second) delete it.second; 
    }
    fRespMatrixMap->clear();
  }

  fMarleyGen.clear();

  if (fOscillogram) delete fOscillogram;
  if (fNadirPDF) delete fNadirPDF; 
  if (fNadirFun) delete fNadirFun; 

  delete fTmpPDF;
  delete fRnd;
}

void MSPDFBuilderTHn::RegisterNadirPDF(TH1* pdf) {
  fNadirPDF = dynamic_cast<TH1D*>(pdf); 
  fNadirFun = new TF1("nadir_fun",
      [this](double* x, double* p) {
      return p[0]*this->fNadirPDF->Interpolate(x[0]); 
      }, 
      fNadirPDF->GetXaxis()->GetXmin(), 
      fNadirPDF->GetXaxis()->GetXmax(), 
      1); 
  fNadirFun->SetParameter(0, 1.0); 

  return;
}

void MSPDFBuilderTHn::RegisterPDF(THn* hist) {
   // Check if an hist with the same name was already loaded
  if (fPDFMap->find(hist->GetName()) != fPDFMap->end()) {
    std::cerr << "error: PDF already loaded\n";
    return;
  }

  fPDFMap->insert( PDFPair( hist->GetName(), new MSTHnPDF(hist->GetName(), hist)) );
  return;
}

void MSPDFBuilderTHn::RegisterPDF(MSTHnPDF* pdf) {
   // Check if an hist with the same name was already loaded
  if (fPDFMap->find(pdf->GetName()) != fPDFMap->end()) {
    std::cerr << "error: PDF already loaded\n";
    return;
  }

  fPDFMap->insert( PDFPair( pdf->GetName(), pdf) );
  return;
}

void MSPDFBuilderTHn::RegisterResponseMatrix(THn* hist) {
   // Check if an hist with the same name was already loaded
  if (fRespMatrixMap->find(hist->GetName()) != fRespMatrixMap->end()) {
    std::cerr << "error: Response matrix already loaded\n";
    return;
  }

  fRespMatrixMap->insert( RespMatrixPair( hist->GetName(), hist) );
  return;
}

double MSPDFBuilderTHn::AddHistToPDF(const std::string& histName, double scaling, NeutrinoPropagator* propagator) {
   // define the total rate
   double total_rate = 0;

   // find hist by name
   PDFMap::iterator im = fPDFMap->find(histName);
   if (im == fPDFMap->end()) {
      std::cerr << "error: PDF with name " << histName.c_str() << " not loaded\n";
      exit(EXIT_FAILURE);
   }
   // check if the temporary PDF exist already
   if (!fTmpPDF) {
      fTmpPDF = fHandler.CreateHn();
      fTmpPDF->SetName("privatePDF");
      fTmpPDF->SetTitle("privatePDF");
      ResetPDF();
      size_t iaxis = 0; 
      for (const auto& axis : fHandler.GetAxes()) {
        if (axis.fSetRange) {
          fTmpPDF->GetAxis(iaxis)->SetRangeUser(axis.fRangeMin, axis.fRangeMax);
        }
        iaxis++;
      }
   }

   EPDFType pdfType = im->second->GetPDFType();
   THn* hn = nullptr; 
   bool response_matrix_applied = false; 

   if (pdfType == EPDFType::kNeutrino) {
     MSTHnPDFNeutrino* pdf = static_cast<MSTHnPDFNeutrino*>(im->second);
     if (pdf->ApplyOscillation()) {
       if (propagator == nullptr) {
         fprintf(stderr, "MSPDFBuilderTHn::AddHistToPDF ERROR: No neutrino propagator associated with %s\n", 
             histName.c_str());
         exit(EXIT_FAILURE); 
       } 
       else {
         MSTHnHandler::axis energy_axis_settings;
         energy_axis_settings.fNbins = pdf->GetTHn()->GetAxis(0)->GetNbins();
         energy_axis_settings.fMin = pdf->GetTHn()->GetAxis(0)->GetXmin();
         energy_axis_settings.fMax = pdf->GetTHn()->GetAxis(0)->GetXmax();
         if (fOscillogram == nullptr) fOscillogram = ComputeOscillationProb(propagator, energy_axis_settings); 
       }
     }

     for (auto& channel : pdf->GetChannels()) {
       // Reset normalization in channel attributes
       channel.fNormalization = 0;
       // Apply oscillation to PDF
       THn* hn_osc = ApplyOscillationProb(im->second->GetTHn(), channel.fPDG);

       // get response matrix from the PDF builder
       auto im_resp = fRespMatrixMap->find(channel.fResponeMatrix);
       if (im_resp == fRespMatrixMap->end()) {
         fprintf(stderr, "MSPDFBuilderTHn::AddHistToPDF ERROR: Response matrix %s not found for %s channel %s\n", 
             histName.data(), channel.fName.data(), channel.fResponeMatrix.c_str());
         exit(EXIT_FAILURE);
       }
       THnD* respMatrix = dynamic_cast<THnD*>(im_resp->second);
       hn = ApplyResponseMatrixAndCrossSection(hn_osc, respMatrix, channel );
       response_matrix_applied = true;
       channel.fNormalization *= exposure_conversion;

       fTmpPDF->Add(hn, scaling * exposure_conversion ); 
/*
 *       TTimer* timer = new TTimer("gSystem->ProcessEvents();", 100, kFALSE);
 *       TCanvas* c = new TCanvas("c", "c", 1200, 1000); 
 *       c->Divide(3,1);
 *       c->cd(1); 
 *       TH2D* h2_surv = fOscillogram->Projection(1,0);
 *       h2_surv->SetEntries( h2_surv->GetNbinsX() * h2_surv->GetNbinsY() );
 *       h2_surv->Draw("colz"); gPad->Update();
 *       c->cd(2); 
 *       TH2D* h2_osc = (TH2D*) hn_osc->Projection(1,0);
 *       h2_osc->SetEntries( h2_osc->GetNbinsX() * h2_osc->GetNbinsY() );
 *       h2_osc->Draw("colz"); gPad->Update();
 *       printf("h2_osc integral = %f\n", h2_osc->Integral());
 *       c->cd(3);
 *       TH2D* h2_recosc = (TH2D*) hn->Projection(1,0,"A");
 *       h2_recosc->SetEntries( h2_recosc->GetNbinsX() * h2_recosc->GetNbinsY() );
 *       h2_recosc->Draw("colz"); gPad->Update();
 *       printf("h2_recosc integral = %f\n", h2_recosc->Integral());
 *       printf("rate = %f\n", rate );
 *
 *       timer->TurnOn(); timer->Reset(); 
 *       getchar(); 
 *       timer->TurnOff();
 */

       total_rate += scaling * channel.fNormalization;

       //printf("[%s] %s.%s: normalization: %g, rate: %g\n", 
           //fName.data(), histName.data(), channel.fName.data(), 
           //channel.fNormalization, scaling * channel.fNormalization);
       //getchar();

       delete hn_osc;
       delete hn; 
     }
   }
   else if (pdfType == EPDFType::kComponent) {
     MSTHnPDFComponent* pdf = static_cast<MSTHnPDFComponent*>(im->second);
     if (pdf->GetRespMatrix()) {
       hn = ApplyResponseMatrix(im->second->GetTHn(), pdf->GetRespMatrix());
       response_matrix_applied = true;
     } else {
       hn = im->second->GetTHn();
     }
     total_rate = scaling;
     fTmpPDF->Add(hn, scaling);
     if (response_matrix_applied) delete hn;
   }
   else {
      fprintf(stderr, "MSPDFBuilderTHn::AddHistToPDF ERROR: PDF type %i not allowed for %s\n", 
          pdfType, histName.data()); 
      exit(EXIT_FAILURE);
   }


   return total_rate;
}

double MSPDFBuilderTHn::AddHistToPDF(const std::string& histName, const string& channel, const double scaling, NeutrinoPropagator* propagator) {
   // define the total rate
   double total_rate = 0;

   // find hist by name
   PDFMap::iterator im = fPDFMap->find(histName);
   if (im == fPDFMap->end()) {
      std::cerr << "error: PDF with name " << histName.c_str() << " not loaded\n";
      exit(EXIT_FAILURE);
   }
   // check if the temporary PDF exist already
   if (!fTmpPDF) {
      fTmpPDF = fHandler.CreateHn();
      fTmpPDF->SetName("privatePDF");
      fTmpPDF->SetTitle("privatePDF");
      ResetPDF();
      size_t iaxis = 0; 
      for (const auto& axis : fHandler.GetAxes()) {
        if (axis.fSetRange) {
          fTmpPDF->GetAxis(iaxis)->SetRangeUser(axis.fRangeMin, axis.fRangeMax);
        }
        iaxis++;
      }
   }

   EPDFType pdfType = im->second->GetPDFType();
   THn* hn = nullptr; 
   bool response_matrix_applied = false; 

   if (pdfType == EPDFType::kComponent) {
      fprintf(stderr, "MSPDFBuilderTHn::AddHistToPDF WARNING: Cannot use channel \"%s\" for kComponent PDF %s. Using standard method.\n", 
          channel.c_str(), histName.c_str());
      return AddHistToPDF(histName, scaling, propagator);
   }
   else if (pdfType == EPDFType::kNeutrino) {
     MSTHnPDFNeutrino* pdf = static_cast<MSTHnPDFNeutrino*>(im->second);
     if (pdf->ApplyOscillation()) {
       if (propagator == nullptr) {
         fprintf(stderr, "MSPDFBuilderTHn::AddHistToPDF ERROR: No neutrino propagator associated with %s\n", 
             histName.c_str());
         exit(EXIT_FAILURE); 
       } 
       else {
         MSTHnHandler::axis energy_axis_settings;
         energy_axis_settings.fNbins = pdf->GetTHn()->GetAxis(0)->GetNbins();
         energy_axis_settings.fMin = pdf->GetTHn()->GetAxis(0)->GetXmin();
         energy_axis_settings.fMax = pdf->GetTHn()->GetAxis(0)->GetXmax();
         if (fOscillogram == nullptr) fOscillogram = ComputeOscillationProb(propagator, energy_axis_settings); 
       }
     }

     MSTHnPDFNeutrino::NuIntChannel_t& ch = pdf->GetChannel(channel);
     // Reset normalization in channel attributes
     ch.fNormalization = 0;
     // Apply oscillation to PDF
     THn* hn_osc = ApplyOscillationProb(im->second->GetTHn(), ch.fPDG);

     // get response matrix from the PDF builder
     auto im_resp = fRespMatrixMap->find(ch.fResponeMatrix);
     if (im_resp == fRespMatrixMap->end()) {
       fprintf(stderr, "MSPDFBuilderTHn::AddHistToPDF ERROR: Response matrix %s not found for %s channel %s\n", 
           histName.data(), ch.fName.data(), ch.fResponeMatrix.c_str());
       exit(EXIT_FAILURE);
     }
     THnD* respMatrix = dynamic_cast<THnD*>(im_resp->second);
     hn = ApplyResponseMatrixAndCrossSection(hn_osc, respMatrix, ch);
     response_matrix_applied = true;
     ch.fNormalization *= exposure_conversion;

     fTmpPDF->Add(hn, scaling * exposure_conversion ); 
     /*
      *       TTimer* timer = new TTimer("gSystem->ProcessEvents();", 100, kFALSE);
      *       TCanvas* c = new TCanvas("c", "c", 1200, 1000); 
      *       c->Divide(3,1);
      *       c->cd(1); 
      *       TH2D* h2_surv = fOscillogram->Projection(1,0);
      *       h2_surv->SetEntries( h2_surv->GetNbinsX() * h2_surv->GetNbinsY() );
      *       h2_surv->Draw("colz"); gPad->Update();
      *       c->cd(2); 
      *       TH2D* h2_osc = (TH2D*) hn_osc->Projection(1,0);
      *       h2_osc->SetEntries( h2_osc->GetNbinsX() * h2_osc->GetNbinsY() );
      *       h2_osc->Draw("colz"); gPad->Update();
      *       printf("h2_osc integral = %f\n", h2_osc->Integral());
      *       c->cd(3);
      *       TH2D* h2_recosc = (TH2D*) hn->Projection(1,0,"A");
      *       h2_recosc->SetEntries( h2_recosc->GetNbinsX() * h2_recosc->GetNbinsY() );
      *       h2_recosc->Draw("colz"); gPad->Update();
      *       printf("h2_recosc integral = %f\n", h2_recosc->Integral());
      *       printf("rate = %f\n", rate );
      *
      *       timer->TurnOn(); timer->Reset(); 
      *       getchar(); 
      *       timer->TurnOff();
      */

     total_rate = scaling * ch.fNormalization;

     //printf("[%s] %s.%s: normalization: %g, rate: %g\n", 
     //fName.data(), histName.data(), ch.fName.data(), 
     //ch.fNormalization, scaling * ch.fNormalization);
     //getchar();

     delete hn_osc;
     delete hn; 
   }
   else {
     fprintf(stderr, "MSPDFBuilderTHn::AddHistToPDF ERROR: PDF type %i not allowed for %s\n", 
         pdfType, histName.data()); 
     exit(EXIT_FAILURE);
   }


   return total_rate;
}


THn* MSPDFBuilderTHn::GetPDF (const std::string& objName) { 
   if (!fTmpPDF) return 0;
   THn* clone = (THn*) fTmpPDF->Clone();
   clone->SetName(objName.c_str());
   clone->SetTitle(objName.c_str());
   ResetPDF();
   return clone;
}

THn* MSPDFBuilderTHn::GetMCRealizaton(int ctsNum, bool addPoissonFluctuation, const std::string& procedure_name) {
   // set internal random number generator if set
   TRandom* rndTmpCopy = nullptr;
   if (fRnd != nullptr)  {
      rndTmpCopy = gRandom;
      gRandom = fRnd;
   }

   // optinally add Poission fluctuatoins on the number of cts
   if (addPoissonFluctuation) {
      ctsNum = gRandom->Poisson(ctsNum);
   }

   // Build a THn<int> with the same axis of the PDF's
   const int dim = fTmpPDF->GetNdimensions();
   Int_t bin[dim], first[dim], last[dim];
   Double_t min[dim], max[dim];
   for (int i = 0 ; i < dim; i++) {
      bin[i]   = fTmpPDF->GetAxis(i)->GetNbins();
      first[i] = fTmpPDF->GetAxis(i)->GetFirst();
      last[i]  = fTmpPDF->GetAxis(i)->GetLast();
      min[i]   = fTmpPDF->GetAxis(i)->GetXmin();
      max[i]   = fTmpPDF->GetAxis(i)->GetXmax();
   }

   THn* realization = new THnD (Form("mc_seed_%u",gRandom->GetSeed()),
                                Form("mc_seed_%u",gRandom->GetSeed()),
                                dim, bin, min, max);

   for (int i = 0 ; i < dim; i++) 
      realization->GetAxis(i)->SetRange(first[i],last[i]);

   // Get the method for producing the data realization
   const MCRealizationProcedure procedure = GetMCRealizationProcedure(procedure_name);

   // Fill realizations using n-dimensional method
   if (procedure == MCRealizationProcedure::kSampling) {
     // Note that GetRandom works on the full range of the histograms and cannot
     // be limited to the user range. That's not a problem since the sampling MUST
     // be done on all histogram range in order to preserve the the actual rate.
     // Note that over- and under-flow bins are not considered
     Double_t rndPoint[dim];
     for ( int j = 0; j < ctsNum; j++ ) {
       fTmpPDF->GetRandom(rndPoint, kFALSE);
       realization->Fill(rndPoint);
     }
   }
   else if (procedure == MCRealizationProcedure::kBinSampling) {
     fHandler.NormalizeHn(fTmpPDF, ctsNum);
     auto* iter = realization->CreateIter(true);
     int coords[dim];
     Long64_t i = 0;
     while ( (i = iter->Next(coords)) >= 0 ) {
       const double expctd_bc = fTmpPDF->GetBinContent(coords);
       const double rnd_bc = gRandom->Poisson(expctd_bc);
       realization->SetBinContent(i, rnd_bc);
     }
     delete iter;
   }
   else if (procedure == MCRealizationProcedure::kAsimov) {
     fHandler.NormalizeHn(fTmpPDF, ctsNum);
     auto* iter = realization->CreateIter(true);
     int coords[dim];
     Long64_t i = 0;
     while ( (i = iter->Next(coords)) >= 0 ) {
       const double expctd_bc = fTmpPDF->GetBinContent(coords);
       realization->SetBinContent(i, expctd_bc);
     }
     delete iter;
   }
   else {
     fprintf(stderr, "MSPDFBuilderTHn::GetMCRealizaton ERROR: unknown procedure %s\n", procedure_name.c_str());
     fprintf(stderr, "MSPDFBuilderTHn::GetMCRealizaton ERROR: available procedures are: sampling, bin_sampling, asimov - fall back to sampling\n");
     Double_t rndPoint[dim];
     for ( int j = 0; j < ctsNum; j++ ) {
       fTmpPDF->GetRandom(rndPoint, kFALSE);
       realization->Fill(rndPoint);
     }
   }

   if (rndTmpCopy != nullptr) gRandom = rndTmpCopy;
   return realization;
}

THn* MSPDFBuilderTHn::ComputeOscillationProb(NeutrinoPropagator* propagator, const MSTHnHandler::axis& energy_axis_settings) {
  //printf("MSPDFBuilderTHn::CreateOscillogramHD with oscillation parameters:\n"); 
  //printf("Δm12 = %g\n", propagator->GetDeltaMSq21()); 
  //printf("Δm23 = %g\n", propagator->GetDeltaMSq32()); 
  //printf("sin²θ12 = %g\n", propagator->GetSinSqTheta12()); 
  //printf("sin²θ13 = %g\n", propagator->GetSinSqTheta13()); 
  //printf("sin²θ23 = %g\n", propagator->GetSinSqTheta23()); 
  //getchar(); 
  //
  // check if axis settings include the nadir angle
  bool nadir_found = false;
  int nadir_axis = -1;
  size_t i = 0;
  for (const auto& axis : fHandler.GetAxes()) {
    TString label = axis.fLabel;
    if (label.Contains("nadir", TString::kIgnoreCase)) {
      nadir_found = true;
      nadir_axis = i;
    }
    i++;
  }

  if (nadir_found) {
    const MSTHnHandler::axis& nadir_axis_settings = fHandler.GetAxes().at(nadir_axis);
    int nbins[2]; double xmin[2]; double xmax[2]; 
    nbins[0] = energy_axis_settings.fNbins;
    nbins[1] = nadir_axis_settings.fNbins;
    xmin[0] = energy_axis_settings.fMin;
    xmin[1] = nadir_axis_settings.fMin;
    xmax[0] = energy_axis_settings.fMax;
    xmax[1] = nadir_axis_settings.fMax;

    THnD *oscillogram = new THnD("surv_map", "surv_map", 2, nbins, xmin, xmax);
    TAxis* energy_axis = oscillogram->GetAxis(0);
    TAxis* nadir_axis = oscillogram->GetAxis(1);
    nadir_axis->SetRangeUser( nadir_axis_settings.fRangeMin, nadir_axis_settings.fRangeMax );
    auto* it = oscillogram->CreateIter(true);
    int coords[2];
    Long64_t i = 0;
    while( (i = it->Next(coords)) >= 0 )  {
      const double nu_ene  = energy_axis->GetBinCenter( coords[0] );
      const double cos_nad = nadir_axis->GetBinCenter( coords[1] );
      //const double cos_nad_x0 = nadir_axis->GetBinLowEdge( coords[1] );
      //const double cos_nad_x1 = nadir_axis->GetBinUpEdge( coords[1] );
      //const double dcn = cos_nad_x1 - cos_nad_x0;
      //const double cos_nad_exposure = fNadirFun->Integral( cos_nad_x0, cos_nad_x1 ) * 1000.0 / dcn;
      propagator->SetEnergy( nu_ene*1e-3 ); 
      propagator->DefinePath( cos_nad, 1.47e8 ); 
      propagator->propagate( 1 );
      const double p11 = propagator->GetProb(1, 1); 
      const double p21 = propagator->GetProb(2, 1); 
      const double p31 = propagator->GetProb(3, 1);

      double f_2 = propagator->GetSinSqTheta12Sun();
      double f_3 = propagator->GetSinSqTheta13(); 
      double surv = (1 - f_2 - f_3) * p11 + f_2 * p21 + f_3 * p31;

      oscillogram->SetBinContent( coords, surv ); 
      //printf("E = %g, cos(nad) = %g, p11 = %g, p21 = %g, p31 = %g, f2 = %g, f3 = %g, surv = %g\n", 
          //nu_ene, cos_nad, p11, p21, p31, f_2, f_3, surv);
    }

    delete it;
    return oscillogram;
  }
  else {
    int nbins[1]; double xmin[1]; double xmax[1]; 
    nbins[0] = energy_axis_settings.fNbins;
    xmin[0] = energy_axis_settings.fMin;
    xmax[0] = energy_axis_settings.fMax;

    THnD *oscillogram = new THnD("surv_curve", "surv_curve", 1, nbins, xmin, xmax);
    TAxis* energy_axis = oscillogram->GetAxis(0);
    auto* it = oscillogram->CreateIter(true);
    int coords[1];
    Long64_t i = 0;
    while( (i = it->Next(coords)) >= 0 )  {
      const double nu_ene  = energy_axis->GetBinCenter( coords[0] );
      propagator->SetEnergy( nu_ene*1e-3 ); 
      propagator->DefinePath( 1, 1.47e8 ); 
      propagator->propagate( 1 );
      const double p11 = propagator->GetProb(1, 1); 
      const double p21 = propagator->GetProb(2, 1); 
      const double p31 = propagator->GetProb(3, 1);

      double f_2 = propagator->GetSinSqTheta12Sun();
      double f_3 = propagator->GetSinSqTheta13(); 
      double surv = (1 - f_2 - f_3) * p11 + f_2 * p21 + f_3 * p31;

      oscillogram->SetBinContent( coords, surv ); 
      //if (i%10000 == 0)
      //printf("E = %g, p11 = %g, p21 = %g, p31 = %g, f2 = %g, f3 = %g, surv = %g\n", 
          //nu_ene, p11, p21, p31, f_2, f_3, surv);
    }

    delete it;
    return oscillogram;
  }
}

THn* MSPDFBuilderTHn::BuildNadirPDF() const {
  if (fNadirPDF == nullptr) {
    fprintf(stderr, "MSPDFBuilderTHn::BuildNadirPDF ERROR: Nadir PDF not registered\n");
    exit(EXIT_FAILURE);
  }

  // get the nadir axis settings from the handler
  MSTHnHandler::axis nadir_axis_settings; 
  for (const auto& axis : fHandler.GetAxes()) {
    TString label = axis.fLabel;
    if (label.Contains("nadir", TString::kIgnoreCase)) {
      nadir_axis_settings = axis;
    }
  }

  printf("BuildNadirPDF: nadir_axis_settings\n");
  nadir_axis_settings.print();

  THnD* nadir_pdf = new THnD("nadir_pdf", "nadir_pdf", 1,
      &nadir_axis_settings.fNbins, &nadir_axis_settings.fMin, &nadir_axis_settings.fMax);

  int i_range_min = 1; 
  int i_range_max = nadir_axis_settings.fNbins;

  if (nadir_axis_settings.fSetRange) {
    i_range_min = nadir_pdf->GetAxis(0)->FindBin(nadir_axis_settings.fRangeMin);
    i_range_max = nadir_pdf->GetAxis(0)->FindBin(nadir_axis_settings.fRangeMax);
  }

  for (int i = 1; i <= nadir_axis_settings.fNbins; i++) {
    if (i < i_range_min || i >= i_range_max) continue;
    double x0 = nadir_pdf->GetAxis(0)->GetBinLowEdge(i);
    double x1 = nadir_pdf->GetAxis(0)->GetBinUpEdge(i);
    double prob = fNadirFun->Integral(x0, x1, 1e-3);
    nadir_pdf->SetBinContent(i, prob);
  }

  fHandler.NormalizeHn(nadir_pdf);

  return nadir_pdf;
}

THn* MSPDFBuilderTHn::ApplyOscillationProb(const THn* target, const int pdg) {
  THnD* product = (THnD*)fOscillogram->Clone("hn_osc"); 
  product->SetNameTitle(Form("%s_osc", target->GetName()), Form("%s_osc", target->GetTitle()));
  product->Reset();
  
  if (fOscillogram->GetNdimensions() == 2) {
    TAxis* energy_axis = product->GetAxis(0);
    TAxis* nadir_axis = product->GetAxis(1);
    auto* it = product->CreateIter(true);
    int coords[2];
    Long64_t i = 0;
    while( (i = it->Next(coords)) >= 0 )  {
      //const double nu_ene  = energy_axis->GetBinCenter( coords[0] );
      //const double nu_flux = target->GetBinContent( coords[0] );
      //const double cos_nad = nadir_axis->GetBinCenter( coords[1] );
      //const double cos_nad_x0 = nadir_axis->GetBinLowEdge( coords[1] );
      //const double cos_nad_x1 = nadir_axis->GetBinUpEdge( coords[1] );
      //const double cos_nad_exposure = fNadirFun->Integral( cos_nad_x0, cos_nad_x1 );
      double surv = fOscillogram->GetBinContent( i ); 
      double nu_flux = target->GetBinContent( i ); 

      if (pdg > 0 && pdg != 12) { surv = 1.0 - surv; }

      product->SetBinContent(i, nu_flux * surv );
      //if (i%100 == 0)
      //printf("E = %g, cos(nad) = %g, cos(nad) exposure = %g, nu_flux = %g, surv = %g, product = %g\n", 
          //nu_ene, cos_nad, cos_nad_exposure, nu_flux, surv, nu_flux * surv * cos_nad_exposure);
    }

    delete it;
  }
  else if (fOscillogram->GetNdimensions() == 1) {
    TAxis* energy_axis = product->GetAxis(0);
    auto* it = product->CreateIter(true);
    int coords[1];
    Long64_t i = 0;
    while( (i = it->Next(coords)) >= 0 )  {
      const double nu_ene  = energy_axis->GetBinCenter( coords[0] );
      const double nu_flux = target->GetBinContent( coords[0] );
      double surv = fOscillogram->GetBinContent( i ); 
      if (pdg > 0 && pdg != 12) { surv = 1.0 - surv; }

      product->SetBinContent(i, nu_flux * surv);
      //if (i%100 == 0)
      //printf("E = %g, nu_flux = %g, surv = %g, product = %g\n", 
          //nu_ene, nu_flux, surv, nu_flux * surv);
    }

    delete it;
  }
  else {
    fprintf(stderr, "MSPDFBuilderTHn::ApplyOscillationProb ERROR: Oscillogram has %d dimensions\n", fOscillogram->GetNdimensions());
    exit(EXIT_FAILURE);
  }

  return product;
}

THn* MSPDFBuilderTHn::ApplyResponseMatrix(const THn* target, const  THn* responseMatrix) {
  // TODO: make these defined in configuration file
  const int iaxis_transform_target = 0; 
  const int iaxis_transform_response = 1; 

  THnD* product = dynamic_cast<THnD*>( fHandler.CreateHn() );
  product->SetName(Form("%s_%s", target->GetName(), responseMatrix->GetName())); 

  const int ndim_target = target->GetNdimensions(); 
  const int ndim_product = product->GetNdimensions(); 

  int ibresp[2] = {0, 0}; 
  int ibtarget[2] = {0, 0}; 
  int ibprod[2] = {0, 0}; 

  const int nbins_target = target->GetAxis(iaxis_transform_target)->GetNbins(); 
  const int nbins_product = product->GetAxis(iaxis_transform_target)->GetNbins(); 

  double response = 0.0; 
  double integral = 0.0;
  double tmp = 0.0; 

  for (int inadir=1; inadir<=product->GetAxis(1)->GetNbins(); inadir++) { // nadir bins loop
    for (int ienu = 1; ienu <= nbins_target; ienu++) { // true energy bins loop
      for (int ierec = 1; ierec <= nbins_product; ierec++) { // reco energy bins loop
        ibresp[0] = ienu; ibresp[1] = ierec;
        response = responseMatrix->GetBinContent(ibresp); 
        if (response != 0 ) {
          ibtarget[0] = ienu; ibtarget[1] = inadir;
          double target_val = target->GetBinContent(ibtarget);
          if (target_val != 0) {
            ibprod[0] = ierec; ibprod[1] = inadir;
            tmp = target_val * response;
            integral += tmp;
            product->AddBinContent(ibprod, tmp);
          }
        }
      }
    }
  }

  // Apply axis settings as defined in the THn handler
  size_t idim = 0; 
  for (const auto& axis : fHandler.GetAxes()) {
    if (axis.fSetRange) {
      product->GetAxis(idim)->SetRangeUser(axis.fRangeMin, axis.fRangeMax);
    }
    idim++;
  }

  return product;
}

THn* MSPDFBuilderTHn::ApplyResponseMatrixAndCrossSection(const THn* target, 
    const  THn* responseMatrix, MSTHnPDFNeutrino::NuIntChannel_t& channel) {
  // TODO: make these defined in configuration file
  const int iaxis_transform_target = 0; 
  const int iaxis_transform_response = 1; 

  THnD* product = static_cast<THnD*>(fHandler.CreateHn());
  product->SetName(Form("%s_%s", target->GetName(), responseMatrix->GetName())); 

  const int ndim_target = target->GetNdimensions(); 
  const int ndim_product = product->GetNdimensions(); 

  const std::vector<double>& crossSection = channel.fCrossSection;
  double& normalization = channel.fNormalization;

  if (ndim_product == 2) {
    int ibresp[2] = {0, 0}; 
    int ibtarget[2] = {0, 0}; 
    int ibprod[2] = {0, 0}; 

    const int nbins_target = target->GetAxis(iaxis_transform_target)->GetNbins(); 
    const int nbins_product = product->GetAxis(iaxis_transform_target)->GetNbins(); 

    double response = 0.0; 
    double xsec = 0.0;
    double tmp = 0.0; 

    for (size_t n=1; n<=product->GetAxis(1)->GetNbins(); n++) { // nadir bins loop
      ibtarget[1] = n; 
      ibprod[1] = n; 
      for (size_t i = 1; i <= nbins_target; i++) { // true energy bins loop
        ibresp[0] = i;
        ibtarget[0] = i;
        xsec = crossSection.at(i-1);
        normalization += xsec * target->GetBinContent(ibtarget);
        for (size_t j = 1; j <= nbins_product; j++) { // reco energy bins loop
          ibresp[1] = j; 
          ibprod[0] = j;
          response = responseMatrix->GetBinContent(ibresp); 
          if (response > 0 ) {
            double target_val = target->GetBinContent(ibtarget);
            if (target_val > 0) {
              tmp = target_val * response * xsec;
              product->AddBinContent(ibprod, tmp);
              //if ((i+j)%100 == 0)
                //printf("ibresp = [%d, %d], ibtarget = [%d, %d], ibprod = [%d, %d], response = %g, target_val = %g, xsec = %g\n", 
                    //ibresp[0], ibresp[1], ibtarget[0], ibtarget[1], ibprod[0], ibprod[1], response, target_val, xsec);
            }
          }
        }
      }
    }
  }
  else if (ndim_product == 1) {
    int ibresp[2] = {0, 0}; 
    int ibtarget[1] = {0};
    int ibprod[1] = {0};

    const int nbins_target = target->GetAxis(iaxis_transform_target)->GetNbins(); 
    const int nbins_product = product->GetAxis(iaxis_transform_target)->GetNbins(); 

    double response = 0.0; 
    double xsec = 0.0;
    double tmp = 0.0; 

    for (size_t i = 1; i <= nbins_target; i++) { // true energy bins loop
      ibresp[0] = i;
      ibtarget[0] = i;
      xsec = crossSection.at(i-1);
      normalization += xsec * target->GetBinContent(ibtarget);
      for (size_t j = 1; j <= nbins_product; j++) { // reco energy bins loop
        ibresp[1] = j; 
        ibprod[0] = j;
        response = responseMatrix->GetBinContent(ibresp); 
        if (response > 0 ) {
          double target_val = target->GetBinContent(ibtarget);
          if (target_val > 0) {
            tmp = target_val * response * xsec;
            product->AddBinContent(ibprod, tmp);
            //if ((i+j)%1 == 0)
              //printf("ibresp = [%d, %d], ibtarget = [%d], ibprod = [%d], response = %g, target_val = %g, xsec = %g\n", 
                  //ibresp[0], ibresp[1], ibtarget[0], ibprod[0], response, target_val, xsec);
          }
        }
      }
    }
  }

  // Apply axis settings as defined in the THn handler
  size_t idim = 0; 
  for (const auto& axis : fHandler.GetAxes()) {
    if (axis.fSetRange) {
      product->GetAxis(idim)->SetRangeUser(axis.fRangeMin, axis.fRangeMax);
    }
    idim++;
  }

  return product;
}


} // namespace mst
