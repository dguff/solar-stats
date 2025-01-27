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
 * \class mst::MSMinimizer
 *
 * \brief 
 * Base class implementing the minimizer
 *
 * \details 
 * Class handling the interface betweeen models and minimizer.
 *
 * \author Matteo Agostini
 */

#ifndef MST_MSMinimizer_H
#define MST_MSMinimizer_H

// ROOT libs
#include <TMinuit.h>

// m-stats libs
#include "MSModel.h"
#include "MSObject.h"

// Prob3++
#include "NeutrinoPropagator.h"

namespace mst {

class MSMinimizer : public MSObject
{
   public:
      //! Type defining a vector of Models used as class member
      using MSModelVector = std::vector<MSModel*>;

      //! Structure to hold the input parameters for the neutrino propagator
      struct PropagatorInputs_t {
        int ix12 = -1;
        int ix13 = -1;
        int ix23 = -1;
        int idcp = -1;
        int idm21 = -1;
        int idm32 = -1;

        double x12 = {};
        double x13 = {};
        double x23 = {};
        double dcp = {};
        double dm21 = {};
        double dm32 = {};

        bool useSinSq = true; 
        int  nubar = 1;

        //! Set oscillation parameters 
        void SetParameters(const double* par) {
          (ix12 > 0) ? x12 = par[ix12] : x12 = 0.0;
          (ix13 > 0) ? x13 = par[ix13] : x13 = 0.0;
          (ix23 > 0) ? x23 = par[ix23] : x23 = 0.0;
          (idcp > 0) ? dcp = par[idcp] : dcp = 0.0;
          (idm21 > 0) ? dm21 = par[idm21] : dm21 = 0.0;
          (idm32 > 0) ? dm32 = par[idm32] : dm32 = 0.0;
        }
      }; 

   public:
      //! Constructor
      MSMinimizer(const std::string& name = "");
      //! Desstructor
      virtual ~MSMinimizer();

      //! Add model (the function does NOT take ownership of the object
      void AddModel(MSModel* model) { fModelVector->push_back(model); }
      //! Get the number of models added to the minimizer
      unsigned int GetNModels() const { return fModelVector->size(); }

      //! Get the internal vector of models
      const MSModelVector* GetModels() const { return fModelVector; }

      //! Get pointer to parameter map (the minimizer keeps ownership of the obj!)
      const MSParameterMap* GetParameterMap() const { return fGlobalParMap; }

      //! Get pointer to parameter map (the minimizer keeps ownership of the obj!)
      MSParameterMap* GetParameterMap() { return fGlobalParMap; }

      //! Get pointer to parameter (the minimizer keeps ownership of the obj!)
      MSParameter* GetParameter(const std::string& parGloablName) {
         MSParameterMap::iterator it = fGlobalParMap->find(parGloablName);
         return it != fGlobalParMap->end() ? it->second : 0;
      }

      //! Get the index of a parameter in the Minuit parameter array
      size_t GetMinuitParameterIndex(const std::string& parGloablName) const {
         MSParameterMap::iterator it = fGlobalParMap->find(parGloablName);
         if (it == fGlobalParMap->end()) {
           fprintf(stderr, "MSMinimizer::GetMinuitParameterIndex: parameter %s not found\n",
                   parGloablName.c_str());
           exit(EXIT_FAILURE);
         }
         size_t index = std::distance(fGlobalParMap->begin(), it);
         return index;
      }


      //! Print summary of the parameters
      void PrintParSummary() const {
         for (const auto& it : *fGlobalParMap) it.second->PrintSummary();
      }

      //! Initialize minuit with custom
      //! -- verbosity can be from -1 (min) to 3 (max)
      //! -- errVal is usually set to 0.5 for NLL or 1 for chi^2
      TMinuit* InitializeMinuit(int verbosity = 0, double errVal = .5);

      //! Set minuit verbosity
      void SetMinuitErrVal (double level);

      //! Set minuit verbosity
      void SetMinuitVerbosity (int level);

      //! Set minuit maxcalls
      void SetMinuitMaxCalls (int maxcalls) { fMinuitMaxCalls = maxcalls; }

      //! Set minuit tolerance
      void SetMinuitTolerance (double tolerance) { fMinuitTolerance = tolerance; }

      //! Set minuit precision
      void SetMinuitPrecision (double precision) { fMinuitPrecision = precision; }

      //! Get the pointer to minuit
      TMinuit* GetMinuit() const { return fMinuit;}

      //! Get the pointer to the neutrino propagator
      NeutrinoPropagator* GetNeutrinoPropagator() const { return fNeutrinoPropagator; }
      //! Update oscillation parameters in the propagator
      void UpdateOscillationParameters();
      //! Get the propagator inputs by reference
      inline PropagatorInputs_t& GetPropagatorInputs() { return fPropagatorInputs; }
      //! Get the propagator inputs by reference
      inline const PropagatorInputs_t& GetPropagatorInputs() const { return fPropagatorInputs; }

      //! Sync parameter info from the model to minuit
      //! Optionally, do not reset the starting values of the fit parameters and 
      //! leave the  values corresponding to the best fit point found in the 
      //! previous interation
      void SyncFitParameters(bool resetFitStartValue  = true);

      //! Call the minimizer
      void Minimize(const std::string& minimizer = "MINIMIZE", 
                    bool resetFitStartValue = true);
      
      //! Call the minimizer
      void Minimize(const std::string& minimizer, bool resetFitStartValue, 
                    int maxcalls, double tolerance) {
         SetMinuitMaxCalls(maxcalls);
         SetMinuitTolerance(tolerance);
         Minimize(minimizer, resetFitStartValue);
      }

      //! Get output status after minuit last call
      int GetMinuitStatus() const { return fMinuitErrorFlag; }

      //! Get output status after minuit last call
      int GetNMinuitFails() const { return fNMinuitFails; }

      //! Get NLL minimum (synced with last migrad call)
      double GetMinNLL() const { return fMinNLL; }

      //! Get estimated vertical distance from minimum
      //!  (synced with last migrad call)
      double GetEDM() const { return fEDM; }

      //! Get minuit quality code of covariance matrix
      //!  (synced with last migrad call)
      double GetCovQual() const { return fCovQual; }


      //! Wrapper function of NLogLikelihood for minuit
      static void  FCNNLLLikelihood(int& npar, double* grad, double& fval,
            double* par, int flag);

   private:

      //! Global pointer for using FCNNLLLikelihood as Minuit FCN
      static MSMinimizer* global_pointer;

      //! Pointer to the model
      MSModelVector* fModelVector {nullptr};
      //! Pointer to the global parameter map
      MSParameterMap* fGlobalParMap {nullptr};
      //! Pointer to a local copy of the parameter map
      MSParameterMap* fLocalParMap {nullptr};

      //! Pointer to minuit
      TMinuit* fMinuit {nullptr};
      //! Argument list used by minuit functions
      double fMinuitArglist[2] {0.0,0.0};
      //! Error flag used by minuit functions
      int fMinuitErrorFlag {0};
      //! Number of times migrad returned an error flag
      int fNMinuitFails {0};

      //! Pointer to the neutrino propagator
      NeutrinoPropagator* fNeutrinoPropagator {nullptr};
      //! Structure to hold the input parameters for the neutrino propagator
      PropagatorInputs_t fPropagatorInputs;

      //! Maximum interations of minuit during minimization
      int fMinuitMaxCalls {2000};

      //! Tolerance on the maximum error during minimization
      double fMinuitTolerance {1e-6};

      //! Numerical precision of the minimizer
      double fMinuitPrecision {1e-16};

      //! Minimum of the negative log likelihood function
      //! Synced with mnstat-fmin
      //!    "the best function value found so far"
      double fMinNLL {0.0};

      //! Distance from minumum
      //! Synced with mnstat-fedm:
      //!    "the estimated vertical distance remaining to minimum"
      double fEDM {0.0};

      //! minuit quality code of covariance matrix
      //! synced with mnstat-istat:
      //!    "a status integer indicating how good is the covariance"
      //!     0 = not calculated at all
      //!     1 = approximation only, not accurate
      //!     2 = full matrix, but forced positive-definite
      //!     3 = full accurate covariance matrix
      int fCovQual {0};
};

} // namespace mst

#endif // MST_MSMinimizer_H
