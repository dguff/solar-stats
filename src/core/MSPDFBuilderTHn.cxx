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
#include "TH1D.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TF1.h"
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

  if (fNadirPDF) delete fNadirPDF; 
  if (fNadirFun) delete fNadirFun; 

  delete fTmpPDF;
  delete fRnd;
}

void MSPDFBuilderTHn::RegisterNadirPDF(TH1* pdf) {
  fNadirPDF = dynamic_cast<TH1D*>(pdf); 
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

void MSPDFBuilderTHn::AddHistToPDF(const std::string& histName, double scaling, NeutrinoPropagator* propagator) {
   // find hist by name
   PDFMap::iterator im = fPDFMap->find(histName);
   if (im == fPDFMap->end()) {
      std::cerr << "error: PDF with name " << histName.c_str() << " not loaded\n";
      return;
   }
   // check if the temporary PDF exist already
   if (!fTmpPDF) {
      fTmpPDF = fHandler.CreateHn();
      fTmpPDF->SetName("privatePDF");
      fTmpPDF->SetTitle("privatePDF");
      ResetPDF();
   }

   // Add hist to pdf with scaling factor
   THn* hn = nullptr; 
   THn* hn_osc = nullptr;
   if (im->second->ApplyOscillation()) {
      if (propagator == nullptr) {
        fprintf(stderr, "MSPDFBuilderTHn::AddHistToPDF ERROR: No neutrino propagator associated with %s\n", 
            histName.c_str());
        exit(EXIT_FAILURE); 
      } 
      else {
        hn_osc = CreateOscillogramHD(im->second, propagator); 
        new TCanvas(); 
        hn_osc->Projection(1, 0)->Draw("colz"); 
        gPad->Update(); 
        getchar(); 
      }
        

      if ( im->second->GetRespMatrix() ) {
        hn = ApplyResponseMatrix(hn_osc, im->second->GetRespMatrix());
      }

      delete hn_osc;
   }

   if (im->second->GetRespMatrix()) {
      hn = ApplyResponseMatrix(im->second->GetTHn(), im->second->GetRespMatrix());
   } else {
      hn = im->second->GetTHn();
   }

   fTmpPDF->Add(hn, scaling);

   if (im->second->GetRespMatrix()) delete hn;
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

   THn* realization = new THnI (Form("mc_seed_%u",gRandom->GetSeed()),
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

THn* MSPDFBuilderTHn::CreateOscillogramHD(MSTHnPDF* pdf, NeutrinoPropagator* propagator) {
  printf("MSPDFBuilderTHn::CreateOscillogramHD with oscillation parameters:\n"); 
  printf("Δm12 = %g\n", propagator->GetDeltaMSq21()); 
  printf("Δm23 = %g\n", propagator->GetDeltaMSq32()); 
  printf("sin²θ12 = %g\n", propagator->GetSinSqTheta12()); 
  printf("sin²θ13 = %g\n", propagator->GetSinSqTheta13()); 
  printf("sin²θ23 = %g\n", propagator->GetSinSqTheta23()); 
  getchar(); 
  const int ndim = 2; 
  const THnBase* pdfHn = pdf->GetTHn();
  const MSTHnHandler::axis& axis_nadir = fHandler.GetAxes().at(1); 
  const int nbin[ndim] = {pdfHn->GetAxis(0)->GetNbins(), axis_nadir.fNbins}; 
  const double xmin[ndim] = {pdfHn->GetAxis(0)->GetXmin(), axis_nadir.fMin};
  const double xmax[ndim] = {pdfHn->GetAxis(0)->GetXmax(), axis_nadir.fMax}; 

  THnD *oscillogram = new THnD("oscillogram", "oscillogram", ndim, nbin, xmin, xmax);
  
  oscillogram->GetAxis(0)->SetRangeUser(2.5, 20.0); 
  auto* it = oscillogram->CreateIter(true);

  int coords[ndim];
  Long64_t i = 0;
  while( (i = it->Next(coords)) >= 0 ) 
  {
    const double nu_ene  = pdfHn->GetAxis(0)->GetBinCenter( coords[0] );
    const double nu_flux = pdfHn->GetBinContent( &coords[0] );
    const double cos_nad = oscillogram->GetAxis(1)->GetBinCenter( coords[1] );
    const double cos_nad_x0 = oscillogram->GetAxis(1)->GetBinLowEdge( coords[1] );
    const double cos_nad_x1 = oscillogram->GetAxis(1)->GetBinUpEdge( coords[1] );
    const double cos_nad_exposure = fNadirPDF->Integral( cos_nad_x0, cos_nad_x1 ) / (cos_nad_x1 - cos_nad_x0);
    propagator->SetEnergy( nu_ene ); 
    propagator->DefinePath( cos_nad, 1.47e8 ); 
    propagator->propagate( 1 );
    const double p11 = propagator->GetProb(1, 1); 
    const double p21 = propagator->GetProb(2, 1); 
    const double p31 = propagator->GetProb(3, 1);

    double f_2 = propagator->GetSinSqTheta12Sun();
    double f_3 = propagator->GetSinSqTheta13(); 
    double surv = (1 - f_2 - f_3) * p11 + f_2 * p21 + f_3 * p31;

    printf("oscillogram: coords = (%d, %d), nu_ene = %g, nu_flux = %g, cos_nad = %g, cos_nad_exposure = %g, p11 = %g, p21 = %g, p31 = %g, surv = %g\n", 
        coords[0], coords[1], nu_ene, nu_flux, cos_nad, cos_nad_exposure, p11, p21, p31, surv);

    oscillogram->SetBinContent( coords, nu_flux * surv * cos_nad_exposure ); 
  }

  getchar();
  fHandler.NormalizeHn( oscillogram ); 
  delete it;

  return oscillogram;
}

THn* MSPDFBuilderTHn::ApplyResponseMatrix(const THn* target, const  THn* responseMatrix) {
  // TODO: make these defined in configuration file
  const int iaxis_transform_target = 0; 
  const int iaxis_transform_response = 1; 

  THn* product = fHandler.CreateHn();
  product->SetName(Form("%s_%s", target->GetName(), responseMatrix->GetName())); 

  const int ndim_target = target->GetNdimensions(); 
  const int ndim_product = product->GetNdimensions(); 

  int ibresp[2] = {0, 0}; 
  int ibtarget[2] = {0, 0}; 
  int ibprod[2] = {0, 0}; 

  const int nbins_target = target->GetAxis(iaxis_transform_target)->GetNbins(); 
  const int nbins_product = product->GetAxis(iaxis_transform_response)->GetNbins(); 

  double response = 0.0; 

  for (size_t n=1; n<=product->GetAxis(1)->GetNbins(); n++) {
    ibtarget[1] = n; 
    ibprod[1] = n; 
    for (size_t i = 1; i <= nbins_target; i++) {
      ibresp[0] = i; 
      for (size_t j = 1; j <= nbins_product; j++) {
        ibresp[1] = j; 
        response = responseMatrix->GetBinContent(ibresp); 
        if (response > 0) {
          product->AddBinContent(ibprod, target->GetBinContent(ibtarget) * response);
        }
      }
    }
  }

  return product;
}

} // namespace mst
