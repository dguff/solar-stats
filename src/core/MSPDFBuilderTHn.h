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

/*!
 * \class mst::MSPDFBuilderTHn
 *
 * \brief 
 * class used for building the pdf
 * 
 * \details 
 * The object stored in an internal map the histograms registered through the
 * function MSPDFBuilderTHn::LoadHist and then add them to a tmp hist (the final
 * PDF) when the method MSPDFBuilderTHn::AddHistToPDF is called. A copy of the
 * internal tmp hist can be retrieved with MSPDFBuilderTHn::GetPDF. The internal
 * hist must be reset through the MSPDFBuilderTHn::Reset method
 *
 * \author Matteo Agostini
 */

#ifndef MST_MSPDFBuilderTHn_H
#define MST_MSPDFBuilderTHn_H

// c/c++ libs
#include <climits>
#include <map>
#include <string>

// ROOT libs
#include "THn.h"
#include "TRandom3.h"

// m-stats libs
#include "MSTHnPDF.h"
#include "MSTHnHandler.h"
#include "MSObject.h"
#include "NeutrinoPropagator.h"

namespace mst {

class MSPDFBuilderTHn : public MSObject
{
 public:
   //! Constructor
   MSPDFBuilderTHn(const std::string& name = "");
   //! Destructor
   virtual ~MSPDFBuilderTHn();

   //! Map of hists
   using PDFPair = std::pair<const std::string, MSTHnPDF*>;
   //! Pair for hist map
   using PDFMap  = std::map <const std::string, MSTHnPDF*>;

   //! Pair of response matrices
   using RespMatrixPair = std::pair<const std::string, THn*>;
   //! Map of response matrices
   using RespMatrixMap  = std::map <const std::string, THn*>;

   //! register PDF 
   void RegisterPDF(THn*);
   //! register PDF
   void RegisterPDF(MSTHnPDF* pdf); 

   //! register nadir pdf
   void RegisterNadirPDF(TH1* pdf); 

   //! register response matrix
   void RegisterResponseMatrix(THn*);

   //! Add scaled histogram to tmp PDF
   void AddHistToPDF(const std::string& histName, double scaling = 1, NeutrinoPropagator* propagator = nullptr);

   //! Set Seed
   void SetSeed(unsigned int seed) { delete fRnd; fRnd = new TRandom3(seed); }

   //! Reset tmp PDF 
   void ResetPDF() { if (fTmpPDF) fTmpPDF->Reset(); }

   //! Get copy of tmp PDF and reset tmpPDF
   THn* GetPDF (const std::string& objName);

   //! Get MC realization extracted by tmpPDF
   THn* GetMCRealizaton(int ctsNum, const bool addPoissonFluctuation = false, const string& procedure_name = "sampling");

   //! Get the response matrix
   inline THn* GetResponseMatrix(const std::string& name) { return fRespMatrixMap->at(name); }

   //! Get the hist handler
   inline MSTHnHandler& GetHistHandler() { return fHandler; }

 protected:
   //! Map of histograms
   PDFMap* fPDFMap {nullptr};
   //! Map of response matrices
   RespMatrixMap* fRespMatrixMap {nullptr};
   //! Nadir Exposure PDF 
   TH1D* fNadirPDF {nullptr};
   //! Nadir Exposure PDF function
   TF1* fNadirFun {nullptr};
   //! Temporary PDF
   THn* fTmpPDF {nullptr};
   //! Pseudo-random number generator
   TRandom* fRnd {nullptr};
   //! Hist handler
   MSTHnHandler fHandler;

   enum class MCRealizationProcedure {
     kUndefined = 0, kSampling = 1, kBinSampling = 2, kAsimov = 3
   };

   inline MCRealizationProcedure GetMCRealizationProcedure(const string& procedure_name) {
     if (procedure_name == "sampling") return MCRealizationProcedure::kSampling;
     else if (procedure_name == "bin_sampling") return MCRealizationProcedure::kBinSampling;
     else if (procedure_name == "asimov") return MCRealizationProcedure::kAsimov;
     else return MCRealizationProcedure::kUndefined;
   }

   // Private methods
   THn* CreateOscillogramHD(MSTHnPDF* pdf, NeutrinoPropagator* propagator);
   THn* ApplyResponseMatrix(const THn* target, const THn* responseMatrix); 
};

} // namespace mst

#endif //MST_MSPDFBuilderTHn_H

