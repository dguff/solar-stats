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
 * \class mst::MSTHnHandler
 *
 * \brief 
 * class used for consistenly handle multiple THn that should undergo the same
 * transforamtion. The class provides method for setting projections, new axis
 * ranges and new bining. The function BuildHist load an histogram from file and
 * return a modified close of it according to the settings requested.
 * 
 * \details 
 *
 * \author Matteo Agostini
 */

#ifndef MST_MSTHnHandler_H
#define MST_MSTHnHandler_H

// c/c++ libs
#include <map>

// ROOT libs
#include <THn.h>

// m-stats libs
#include "MSObject.h"

namespace mst {

class MSTHnHandler : public MSObject
{
   public:
     struct axis {
       bool   fSetRange = {false};
       double fMin      = {0.0};
       double fMax      = {0.0};
       double fRangeMin = {0.0};
       double fRangeMax = {0.0};
       int    fNbins    = {0};
       int    fNgroup   = {1};
       std::string fLabel = {""};
       bool operator==(const axis& a) const {
          return fSetRange == a.fSetRange 
            && fMin == a.fMin       && fMax == a.fMax 
            && fRangeMin == a.fRangeMin && fRangeMax == a.fRangeMax
            && fNbins == a.fNbins   && fNgroup == a.fNgroup;
       }
       void operator=(const axis& a) {
          fSetRange = a.fSetRange;
          fMin = a.fMin;
          fMax = a.fMax;
          fRangeMin = a.fRangeMin;
          fRangeMax = a.fRangeMax;
          fNbins = a.fNbins;
          fNgroup = a.fNgroup;
          fLabel = a.fLabel;
       }
       void print() const {
         printf("axis settings: label: %s, nbins = %d, min = %f, max = %f, range = [%g, %g] ngroup = %d\n", 
                fLabel.data(), fNbins, fMin, fMax, fRangeMin, fRangeMax, fNgroup);
       }
     };

      //! Constructor
      MSTHnHandler(const std::string& name = "") : MSObject (name) {}
      //! Destructor
      virtual ~MSTHnHandler() {}

      //! define sequence of axis to which project the THn
      void ProjectToAxis(const std::vector<int>& axisID) {
         for (const auto& i : axisID) fProjectID.push_back(i);
      }

      //! Get the sequence of axis to which project the THn
      inline const std::vector<int>& GetProjectID() const { return fProjectID; }
      
      //! Set the user range of a specific axis
      inline void SetRange(const int axisID, const double min, const double max) {
         if (axisID >= fAxis.size()) fAxis.resize(axisID+1);
         fAxis.at(axisID).fRangeMin = min;
         fAxis.at(axisID).fRangeMax = max;
         fAxis.at(axisID).fSetRange = true;
      }

      //! Set the label of a specific axis
      void inline SetLabel(const int axisID, const std::string& label) {
         if (axisID >= fAxis.size()) fAxis.resize(axisID+1);
         fAxis.at(axisID).fLabel = label;
      }

      //! Set the axis limits for a specific axis for checking consistency between THns
      inline void SetLimits(const int axisID, const double min, const double max) {
         if (axisID >= fAxis.size()) fAxis.resize(axisID+1);
         fAxis.at(axisID).fMin = min;
         fAxis.at(axisID).fMax = max;
         return;
      }

      //! Set the number of bins for a specific axis
      inline void SetNbins(const int axisID, const int nbins) {
         if (axisID >= fAxis.size()) fAxis.resize(axisID+1);
         fAxis.at(axisID).fNbins = nbins;
      }

      //! Rebin a specific axis
      void Rebin(const int axisID, const int ngroup) { 
         if (axisID >= fAxis.size()) fAxis.resize(axisID+1);
         fAxis.at(axisID).fNgroup = ngroup;
      }

      //! Rebin a THn to match the given axes settings
      static THn* RebinHist(const THn* hn, const std::vector<int>& iaxis_target, const std::vector<TAxis*>& axes_reference, int* bin_group=nullptr);
      
      //! Rebin a THn to match the given axes settings
      static THn* RebinHist(const THn* hn, const std::vector<int>& iaxis_target, const std::vector<axis>& axes_reference, int* bin_group=nullptr);
      
      //! Print THn axis settings
      inline static void PrintAxes(const THn* hn) {
        printf("%s axis settings:\n", hn->GetName());
        for (int i = 0; i < hn->GetNdimensions(); i++) {
          TAxis* axis = hn->GetAxis(i);
          printf("\taxis %d: title: %s,  %d bins, min %f, max %f\n", 
              i, axis->GetTitle(), axis->GetNbins(), axis->GetXmin(), axis->GetXmax());
        }
      }

      //! Integrate a THn
      static double Integral(const THn* hn, const bool respectRange = false) {
        auto it = hn->CreateIter(respectRange);
        Long64_t i = 0;
        double integral = 0;
        while ((i = it->Next()) >= 0) integral += hn->GetBinContent(i);
        delete it;
        return integral;
      }

      //! The normalization of the histogram will be performed in the 
      //! user range if respectAxisUserRange is true. Otherwise by default it
      //! includes all bins, including over- and under-shot bins
      void RespectAxisUserRange(const bool respectAxisUserRange = true) { 
         fRespectUserRange = respectAxisUserRange;
      }

      //! Load histogram from file, withour further manipulation
      THn* LoadHist(const std::string& filename, 
                    const std::string& histName, 
                    const std::string& newHistName, 
                    const bool normalize = false); 

      //! Load histogram from file, manipulate it and return a copy
      THn* BuildHist(const std::string& fileName, 
                     const std::string& histName,
                     const std::string& newHistName,
                     bool  normalize = false);


      //! Create a THn that is the combination of two THns
      THn* FactorizeTHn(const THn* pdf0, const THn* pdf1);

      THn* CreateHn(); 

      inline bool IsConsistent(const THn* hn) const {
        if (hn->GetNdimensions() != fAxis.size()) return false;
        for (int i = 0; i < hn->GetNdimensions(); i++) {
          if (fAxis.at(i).fNbins != hn->GetAxis(i)->GetNbins()) return false;
          if (fAxis.at(i).fMin != hn->GetAxis(i)->GetXmin()) return false;
          if (fAxis.at(i).fMax != hn->GetAxis(i)->GetXmax()) return false;
        }
        return true;
      }

      static void NormalizeHn(THn* hn, const double norm = 1.0, const bool respectUserRange = false);  

      //! Get axes vector
      const std::vector<axis>& GetAxes() const { return fAxis; }
      std::vector<axis>& GetAxes() { return fAxis; }

      //! Reset settings
      void Reset() { fProjectID.clear(); fAxis.clear(); fRespectUserRange = false; }

   protected:
      std::vector<int> fProjectID;
      std::vector<axis> fAxis;

      bool fRespectUserRange = {false};
};

} // namespace mst

#endif //MST_MSTHnHandler_H

