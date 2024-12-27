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
//

// c/c++ libs
#include <iostream>

// ROOT libgs
#include "TFile.h"
#include "THnBase.h"
#include "TH1.h"
#include "TH2.h"

// m-stats libs
#include "MSTHnHandler.h"

namespace mst {

THn* MSTHnHandler::LoadHist( const std::string& fileName, 
                             const std::string& histName, 
                             const std::string& newHistName ) 
{
   // Get pointer of the file 
   TFile inputFile(fileName.c_str(), "READ");
   if(inputFile.IsOpen() == kFALSE) {
      std::cerr << "error: input file not found\n";
      exit(1);
   }

   // load new hist checking for the type
   THn* outHist = nullptr;
   TObject* hist = inputFile.Get(histName.c_str());
   if (!hist) {
      std::cerr << "error: hist " << histName << " not found in the file\n";
      exit(1);
   } else {
      // Create local THn
      if (dynamic_cast<THnBase*>(hist)) {
         outHist = THn::CreateHn(newHistName.c_str(), newHistName.c_str(), 
                   dynamic_cast<THnBase*>(hist));
      } else if (dynamic_cast<TH2*>(hist)) {
         outHist = THn::CreateHn(newHistName.c_str(), newHistName.c_str(), 
                   dynamic_cast<TH2*>(hist));
      } else if (dynamic_cast<TH1*>(hist)) {
         outHist = THn::CreateHn(newHistName.c_str(), newHistName.c_str(), 
                   dynamic_cast<TH1*>(hist));
      }
      else {
         std::cerr << "error: hist " << newHistName << " is not of type THnBase or TH1\n";
         exit(1);
      }

      delete hist;
   }
   inputFile.Close();

   outHist->SetTitle(newHistName.c_str());
   outHist->SetName(newHistName.c_str());

   return outHist;
}

THn* MSTHnHandler::BuildHist(const std::string& fileName, 
                             const std::string& histName,
                             const std::string& newHistName,
                             const bool normalize) {

   // Get pointer of the file 
   TFile inputFile(fileName.c_str(), "READ");
   if(inputFile.IsOpen() == kFALSE) {
      std::cerr << "error: input file not found\n";
      exit(1);
   }

   // load new hist checking for the type
   THn* outHist = nullptr;
   TObject* hist = inputFile.Get(histName.c_str());
   if (!hist) {
      std::cerr << "error: hist " << histName << " not found in the file\n";
      exit(1);
   } else {
      // Create local THn
      if (dynamic_cast<THnBase*>(hist)) {
         outHist = THn::CreateHn(newHistName.c_str(), newHistName.c_str(), 
                   dynamic_cast<THnBase*>(hist));
      } else if (dynamic_cast<TH2*>(hist)) {
         outHist = THn::CreateHn(newHistName.c_str(), newHistName.c_str(), 
                   dynamic_cast<TH2*>(hist));
      } else if (dynamic_cast<TH1*>(hist)) {
         outHist = THn::CreateHn(newHistName.c_str(), newHistName.c_str(), 
                   dynamic_cast<TH1*>(hist));
      }
      else {
         std::cerr << "error: hist " << newHistName << " is not of type THnBase or TH1\n";
         exit(1);
      }

      delete hist;
   }
   inputFile.Close();

   // project hist
   if (fProjectID.size()) {
      int IDs[fProjectID.size()];
      for (int i =0; i < fProjectID.size(); i++) IDs[i] = fProjectID[i];

      THn* tmp = outHist->Projection(fProjectID.size(), IDs);
      delete outHist;
      outHist = tmp;
   }

   // consistency check
   if (fAxis.size() > outHist->GetNdimensions()) {
      std::cerr << "error: requested to set parameters for not existing axis.\n";
      exit(1);
   }

   // rebin
   if (fAxis.size()) {

      // check if a rebin is needed
      for (int i = 0; i < fAxis.size(); i++) {
         if (fAxis[i].fNgroup > 1) {

            // a new hist is substitute to the original one because the method THn 
            // doesn't not provide a method that modify the object itself
            int ngroup[outHist->GetNdimensions()]; 
            for (size_t idim = 0; idim < outHist->GetNdimensions(); idim++) ngroup[idim] = 1;
            for (int i = 0; i < fAxis.size(); i++) ngroup[i] = fAxis[i].fNgroup;

            THn* tmp = outHist->Rebin(ngroup);
            delete outHist;
            outHist = tmp;

            // exit from loop
            break; 
         }
      }

      // set range user if changed
      for (int idim = 0; idim < fAxis.size(); idim++) {
        MSTHnHandler::axis axis_temp;
         if (fAxis.at(idim).fSetRange) {
           outHist->GetAxis(idim)->SetRangeUser(fAxis[idim].fMin,fAxis[idim].fMax);
         }
         axis_temp.fMin   = outHist->GetAxis(idim)->GetXmin();
         axis_temp.fMax   = outHist->GetAxis(idim)->GetXmax();
         axis_temp.fNbins = outHist->GetAxis(idim)->GetNbins();

         if (axis_temp == fAxis.at(idim)) continue;
         else {
            std::cerr << "MSTHnHandler::BuildHist() WARNING: axis " << idim << " differs from set\n";
            std::cerr << "initial settings: \n";
            fAxis.at(idim).print();
            std::cerr << "current settings: \n";
            axis_temp.print();

            fAxis.at(idim) = axis_temp;
         }
      }
   }

   // Normalize 
   if (normalize) {
      auto it = outHist->CreateIter(fRespectUserRange);
      Long64_t i = 0;
      double integral = 0;
      while ((i = it->Next()) >= 0) integral += outHist->GetBinContent(i);
      outHist->Scale(1./integral);
   }

  outHist->SetTitle(newHistName.c_str());
  outHist->SetName(newHistName.c_str());
  return outHist;
}

THn* MSTHnHandler::CreateHn() {
  if (fAxis.size() == 0) {
    fprintf(stderr, "MSTHnHandler::CreateHn ERROR: axes template vector is empty.\n"); 
    exit(EXIT_FAILURE);
  }
  std::vector<int>  nBins(fAxis.size());
  std::vector<double> min(fAxis.size());
  std::vector<double> max(fAxis.size());
  int idim = 0; 
  for (const auto& axis : fAxis) {
     nBins.at(idim) = axis.fNbins;
     min.at(idim) = axis.fMin;
     max.at(idim) = axis.fMax;
     //printf("axis %d: %d bins, min %f, max %f\n", idim, axis.fNbins, axis.fMin, axis.fMax);
     idim++; 
  }

  THnD* outHist = new THnD("tmp", "tmp", fAxis.size(), nBins.data(), min.data(), max.data());
  return outHist;
}

void MSTHnHandler::NormalizeHn(THn* hn, const double norm) {
  auto it = hn->CreateIter(fRespectUserRange);
  Long64_t i = 0;
  double integral = 0;
  while ((i = it->Next()) >= 0) integral += hn->GetBinContent(i);
  hn->Scale( norm / integral);

  return;
}

} // namespace mst
