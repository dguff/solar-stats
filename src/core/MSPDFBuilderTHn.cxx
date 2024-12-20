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
      std::cerr << "error: PDF not loaded\n";
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
   if (im->second->ApplyOscillation()) {
      if (propagator == nullptr) {
        fprintf(stderr, "MSPDFBuilderTHn::AddHistToPDF ERROR: No neutrino propagator given\n");
        exit(EXIT_FAILURE); 
      } 
      

   }

   if (im->second->GetRespMatrix()) {

   }

   fTmpPDF->Add(hn, scaling);
}

THn* MSPDFBuilderTHn::GetPDF (const std::string& objName) { 
   if (!fTmpPDF) return 0;
   THn* clone = (THn*) fTmpPDF->Clone();
   clone->SetName(objName.c_str());
   clone->SetTitle(objName.c_str());
   ResetPDF();
   return clone;
}

THn* MSPDFBuilderTHn::GetMCRealizaton(int ctsNum, bool addPoissonFluctuation) {
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

   // Fill realizations using n-dimensional method
   // Note that GetRandom works on the full range of the histograms and cannot
   // be limited to the user range. That's not a problem since the sampling MUST
   // be done on all histogram range in order to preserve the the actual rate.
   // Note that over- and under-flow bins are not considered
   Double_t rndPoint[dim];
   for ( int j = 0; j < ctsNum; j++ ) {
      fTmpPDF->GetRandom(rndPoint, kFALSE);
      realization->Fill(rndPoint);
   }

   if (rndTmpCopy != nullptr) gRandom = rndTmpCopy;
   return realization;
}

THn* MSPDFBuilderTHn::CreateOscillogramHD(MSTHnPDF* pdf, NeutrinoPropagator* propagator) {
  const int ndim = 2; 
  const THnBase* pdfHn = pdf->GetTHn();
  const MSTHnHandler::axis& axis_nadir = fHandler.GetAxes().at(1); 
  const int nbin[ndim] = {pdfHn->GetAxis(0)->GetNbins(), axis_nadir.fNbins}; 
  const double xmin[ndim] = {pdfHn->GetAxis(0)->GetXmin(), axis_nadir.fMin};
  const double xmax[ndim] = {pdfHn->GetAxis(0)->GetXmax(), axis_nadir.fMax}; 

  THnD *oscillogram = new THnD("oscillogram", "oscillogram", ndim, nbin, xmin, xmax);
  
  oscillogram->GetAxis(0)->SetRangeUser(2.5, 20.0); 
  auto* it = oscillogram->CreateIter(true);

  int* coords = new int[ndim];
  Long64_t i = 0;
  while( (i = it->Next(coords)) >= 0 ) 
  {
    const double nu_ene  = pdfHn->GetAxis(0)->GetBinCenter( coords[0] );
    const double nu_flux = pdfHn->GetBinContent( &coords[0] );
    const double cos_nad = fNadirPDF->GetBinContent( coords[1] );
    propagator->SetEnergy( nu_ene ); 
    propagator->DefinePath( cos_nad, 1.47e8 ); 
    propagator->propagate( 1 );
    const double p11 = propagator->GetProb(1, 1); 
    const double p21 = propagator->GetProb(2, 1); 
    const double p31 = propagator->GetProb(3, 1);

    double f_2 = propagator->GetSinSqTheta12Sun();
    double f_3 = propagator->GetSinSqTheta13(); 
    double surv = (1 - f_2 - f_3) * p11 + f_2 * p21 + f_3 * p31;

    oscillogram->SetBinContent( coords, nu_flux * surv * cos_nad ); 
  }

  fHandler.NormalizeHn( oscillogram ); 
  delete[] coords;

  return oscillogram;
}

THn* MSPDFBuilderTHn::ApplyResponseMatrix(THn* target, THn* responseMatrix) {
  THn* out = fHandler.CreateHn();
  out->SetName(Form("%s_%s", target->GetName(), responseMatrix->GetName())); 

  auto* it = out->CreateIter(true); 
  Long64_t i = 0; 
  int* coords = new int[out->GetNdimensions()];

  while ( (i = it->Next(coords)) >= 0 ) {


  }

  return out; 
}

} // namespace mst
