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
#include <THn.h>
#include <TRandom3.h>

// m-stats libs
#include "MSObject.h"

namespace mst {

class MSPDFBuilderTHn : public MSObject
{
 public:
   //! Constructor
   MSPDFBuilderTHn(const std::string& name = "");
   //! Destructor
   virtual ~MSPDFBuilderTHn();

   //! Map of hists
   using HistPair = std::pair<const std::string, THn*>;
   //! Pair for hist map
   using HistMap  = std::map <const std::string, THn*>;

   //! register histogram 
   void RegisterHist(THn*);

   //! Add scaled histogram to tmp PDF
   void AddHistToPDF(const std::string& histName, double scaling = 1);

   //! Set Seed
   void SetSeed(unsigned int seed) { delete fRnd; fRnd = new TRandom3(seed); }

   //! Reset tmp PDF 
   void ResetPDF() { if (fTmpPDF) fTmpPDF->Reset(); }

   //! Get copy of tmp PDF and reset tmpPDF
   THn* GetPDF (const std::string& objName);

   //! Get MC realization extracted by tmpPDF
   THn* GetMCRealizaton(int ctsNum, bool addPoissonFluctuation = false);

 protected:
   // Map of histograms
   HistMap* fHistMap {nullptr};
   THn*     fTmpPDF  {nullptr};
   TRandom* fRnd     {nullptr};
};

} // namespace mst

#endif //MST_MSPDFBuilderTHn_H

