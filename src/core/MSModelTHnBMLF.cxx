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

// m-stats libs
#include "MSMath.h"
#include "MSModelTHnBMLF.h"

namespace mst {

double MSModelTHnBMLF::NLogLikelihood(double* par, NeutrinoPropagator* propagator)
{
   fPDFBuilder->ResetPDF();

   // retrieve parameters from Minuit and compute the total exposure
   //for (int i =0; i < fParNameList->size(); i++) {
      //const double par_cts = GetMinuitParameter(par, fParNameList->at(i));
      //fPDFBuilder->AddHistToPDF(fParNameList->at(i),  par_cts, propagator);
   //}
   for (const auto& par_itr : *fParameters) {
     if (par_itr.second->IsInput()) {
       const std::string par_name = GetLocalName( par_itr.second->GetName() );
       const double par_cts = GetMinuitParameter(par, par_name);
       fPDFBuilder->AddHistToPDF(par_name, par_cts, propagator);
     }
   }

   const THn* pdf = fPDFBuilder->GetPDF("tmpPDF");
   if (pdf == 0) {
      std::cerr << "NLogLikelihood >> error: PDFBuilder returned unknown object type\n";
      exit(1);
   }
   if (fDataSet == 0) {
      std::cerr << "NLogLikelihood >> error: DataHist not set\n";
      exit(1);
   }

   double logLikelihood = 0.0;
   // loop over dimensions
   auto it = fDataSet->CreateIter(kTRUE);
   Long64_t i = 0;
   int coords[fDataSet->GetNdimensions()];
   while ((i = it->Next(coords)) >= 0) {
     //printf("[%i, %i] -> (%g, %g): bc = %g - pdf = %g\n", 
         //coords[0], coords[1], 
         //fDataSet->GetAxis(0)->GetBinCenter(coords[0]),
         //fDataSet->GetAxis(1)->GetBinCenter(coords[1]),
         //fDataSet->GetBinContent(i), 
         //fExposure*pdf->GetBinContent(i));
     logLikelihood += MSMath::LogPoisson(fDataSet->GetBinContent(i), 
         fExposure*pdf->GetBinContent(i));
   }

   delete pdf;
   delete it;
   return (-logLikelihood);
}

bool MSModelTHnBMLF::AreInputHistsConsistent () 
{
   const THn* pdf = fPDFBuilder->GetPDF("tmpPDF");
   if (pdf == 0) {
      std::cerr << "error: PDFBuilder returned unknown object type\n";
      exit(1);
   }
   if (fDataSet == 0) {
      std::cerr << "error: DataHist not set\n";
      exit(1);
   }

   if (fDataSet->GetNdimensions() != pdf->GetNdimensions()) return false;
   for (int i = 0; i < fDataSet->GetNdimensions(); i++) {
      if (fDataSet->GetAxis(i)->GetNbins() != pdf->GetAxis(i)->GetNbins()) return false;
      if (fDataSet->GetAxis(i)->GetXmin()  != pdf->GetAxis(i)->GetXmin() ) return false;
      if (fDataSet->GetAxis(i)->GetXmax()  != pdf->GetAxis(i)->GetXmax() ) return false;
   }
   return true;

}

} // namespace mst
