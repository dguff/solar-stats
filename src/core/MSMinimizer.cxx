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

// root libs
#include <TString.h>

// m-stats libs
#include "MSMinimizer.h"

// Prob3++
#include "BargerPropagator.h"

namespace mst {

MSMinimizer* MSMinimizer::global_pointer = 0;

MSMinimizer::MSMinimizer(const std::string& name) : MSObject(name)
{
   fModelVector        = new MSModelVector();
   fLocalParMap        = new MSParameterMap();
   fNeutrinoPropagator = new BargerPropagator(); 
}

MSMinimizer::~MSMinimizer()
{
  printf("calling MSMinimizer destructor\n");
  getchar(); 
  
   if (fModelVector) {
      for (auto& i : *fModelVector) delete i;
      delete fModelVector;
   }

   if (fLocalParMap) {
      for (auto& it : *fLocalParMap) delete it.second;
      delete fLocalParMap;
   }

   if (fNeutrinoPropagator) delete fNeutrinoPropagator;

   delete fMinuit;
}

TMinuit* MSMinimizer::InitializeMinuit (int verbosity, double errVal)
{
   // Check if module list contains at least one model
   if (fModelVector->size() == 0) {
      std::cerr << "MSMinimizer::InitializeMinuit: module list empty"
                << std::endl;
      exit(1);
   }
   // Set global pointer
   global_pointer = this;

   // Fixed the pointer to the global parameter list
   fGlobalParMap = fModelVector->at(0)->GetParameters();

   // Initialize the propagator inputs container
   size_t i = 0;
   for (const auto& par : *fGlobalParMap) {
     if (par.second->IsOscillation() == false) {i++; continue;}

     TString parName = par.first;
     if ( parName.Contains("Theta12") ) {
       fPropagatorInputs.x12 = par.second->GetFitStartValue();
       fPropagatorInputs.ix12 = i;
     } else if ( parName.Contains("Theta13")) {
       fPropagatorInputs.x13 = par.second->GetFitStartValue();
       fPropagatorInputs.ix13 = i;
     } else if ( parName.Contains("Theta23") ) {
       fPropagatorInputs.x23 = par.second->GetFitStartValue();
       fPropagatorInputs.ix23 = i;
     } else if ( parName.Contains("deltaCP") ) {
       fPropagatorInputs.dcp = par.second->GetFitStartValue();
       fPropagatorInputs.idcp = i;
     } else if ( parName.Contains("dm21")  ) {
       fPropagatorInputs.dm21 = par.second->GetFitStartValue();
       fPropagatorInputs.idm21 = i;
     } else if ( parName.Contains("dm32") ) {
       fPropagatorInputs.dm32 = par.second->GetFitStartValue();
       fPropagatorInputs.idm32 = i;
     }
     i++;
   }

   // Initialize minuit
   delete fMinuit;
   fMinuit = new TMinuit(fGlobalParMap->size());
   fMinuit->SetFCN(&MSMinimizer::FCNNLLLikelihood);

   SetMinuitVerbosity(verbosity);
   SetMinuitErrVal(errVal);

   return fMinuit;
}

void MSMinimizer::SetMinuitErrVal(double errVal)
{
   if (!fMinuit) {
      std::cerr << "MSMinimizer::SetMinuitErrVal: minuit not initialized yet"
                << std::endl;
   }
   fMinuitArglist[0] = errVal;
   fMinuit->mnexcm("SET ERR", fMinuitArglist, 1, fMinuitErrorFlag);
}

void MSMinimizer::SetMinuitVerbosity(int level)
{
   if (!fMinuit) {
      std::cerr << "MSMinimizer::SetMinuitVerbosity: minuit not initialized yet"
                << std::endl;
   }
   fMinuit->SetPrintLevel(level);
}

void MSMinimizer::SyncFitParameters( bool resetParStartVal)
{
   // Check whether minuit has been last sync against this minimizer
   // and define value of global pointer
   bool forceUpdateAll = false;
   if (global_pointer != this) {
      forceUpdateAll = true;
      global_pointer = this;
   }

   // Initialize minuit if not done manually
   if (!fMinuit) InitializeMinuit();

   const MSParameterMap::const_iterator gItB = fGlobalParMap->begin();
   const MSParameterMap::const_iterator gItE = fGlobalParMap->end();

   // Parse parameter to minuit (only if different from previous call)
   for (MSParameterMap::const_iterator gIt = gItB; gIt != gItE; ++gIt) {
      const int d = distance(gItB, gIt);
      MSParameterMap::const_iterator lIt = fLocalParMap->find(gIt->first);

      // Intialize all parameters the first time the function is called or
      // if minuit was synced against a different MSMinimizer
      if (lIt == fLocalParMap->end() || forceUpdateAll) {
         if (fVerbosity) std::cerr << "MSMinimizer::SyncFitParameters: "
                                   << "par[" << d
                                   << "] \"" << gIt->first << "\""
                                   << " -> synced all fields"
                                   << std::endl;

         fMinuit->mnparm ( d, gIt->second->GetName().data(),
               gIt->second->GetFitStartValue(),
               gIt->second->GetFitStartStep(),
               gIt->second->GetRangeMin(),
               gIt->second->GetRangeMax(),
               fMinuitErrorFlag);

         if (gIt->second->IsFixed()) {
            if (fVerbosity) std::cerr << "MSMinimizer::SyncFitParameters: "
                                      << "par[" << d
                                      << "] \"" << gIt->first << "\""
                                      << " -> fixed"
                                      << std::endl;
            fMinuitArglist[0] = d+1;
            fMinuit->mnexcm("FIX", fMinuitArglist , 1, fMinuitErrorFlag);
         }

      // otherwise update only changed field
      } else if (gIt->second->IsFixed()) {

         if (gIt->second->GetFitStartValue() != lIt->second->GetFitStartValue()) {
            if (fVerbosity) std::cerr << "MSMinimizer::SyncFitParameters: "
                                      << "par[" << d
                                      << "] \"" << gIt->first << "\""
                                      << " -> synced starting value"
                                      << std::endl;
            fMinuitArglist[0] = d+1;
            fMinuitArglist[1] = gIt->second->GetFitStartValue();
            fMinuit->mnexcm("SET PAR",fMinuitArglist,2, fMinuitErrorFlag);
         }

         if (!lIt->second->IsFixed()) {
            if (fVerbosity) std::cerr << "MSMinimizer::SyncFitParameters: "
                                      << "par[" << d
                                      << "] \"" << gIt->first << "\""
                                      << " -> fixed"
                                      << std::endl;
            fMinuitArglist[0] = d+1;
            fMinuit->mnexcm("FIX", fMinuitArglist , 1, fMinuitErrorFlag);
         }

      } else {
         if (lIt->second->IsFixed()) {
            if (fVerbosity) std::cerr << "MSMinimizer::SyncFitParameters: "
                                      << "par[" << d
                                      << "] \"" << gIt->first << "\""
                                      << " -> released"
                                      << std::endl;
            fMinuitArglist[0] = d+1;
            fMinuit->mnexcm("RELEASE", fMinuitArglist , 1, fMinuitErrorFlag);
         }

         if ( resetParStartVal ||
              gIt->second->GetFitStartValue() != lIt->second->GetFitStartValue()||
              gIt->second->GetFitStartStep()  != lIt->second->GetFitStartStep() ||
              gIt->second->GetRangeMin()      != lIt->second->GetRangeMin() ||
              gIt->second->GetRangeMax()      != lIt->second->GetRangeMax() ) {
            if (fVerbosity) std::cerr << "MSMinimizer::SyncFitParameters: "
                                      << "par[" << d
                                      << "] \"" << gIt->first << "\""
                                      << " -> synced all fields"
                                      << std::endl;

            fMinuit->mnparm ( d, gIt->second->GetName().data(),
                                 gIt->second->GetFitStartValue(),
                                 gIt->second->GetFitStartStep(),
                                 gIt->second->GetRangeMin(),
                                 gIt->second->GetRangeMax(),
                                 fMinuitErrorFlag);
         }
      }
   }

   // Clear and Make local copy of the gloabal parameter map
   for (MSParameterMap::iterator it = fLocalParMap->begin();
         it != fLocalParMap->end(); ++it)
      delete it->second;
   fLocalParMap->clear();

   for (MSParameterMap::const_iterator gIt = gItB; gIt != gItE; ++gIt) {
      MSParameter* newPar = new MSParameter(*(gIt->second));
      fLocalParMap->insert(MSParameterPair(newPar->GetName(), newPar));
   }

}

void MSMinimizer::Minimize(const std::string& minimizer, bool resetFitStartValue) {
   // Sync parameters
   SyncFitParameters(resetFitStartValue);
   // Set maxcalls
   fMinuitArglist[0] = fMinuitMaxCalls;
   // Set tolerance
   fMinuitArglist[1] = fMinuitTollerance;
   // Run actual minimization
   fMinuit->mnexcm(minimizer.c_str(), fMinuitArglist, 2, fMinuitErrorFlag);

   if (GetMinuitStatus()) fNMinuitFails++;

   // Retrive fit results from minuit and store info
   MSParameterMap::const_iterator gIt0 = fGlobalParMap->begin();
   for (MSParameterMap::const_iterator gIt = fGlobalParMap->begin();
         gIt != fGlobalParMap->end(); ++gIt) {
      int d = distance(gIt0, gIt);
      double fitBestValue, fitBestValueErr;
      double fitRangeMin, fitRangeMax;
      TString name;
      int index;
      fMinuit->mnpout(d, name, fitBestValue, fitBestValueErr,
                      fitRangeMin, fitRangeMax, index);
      gIt->second->SetFitBestValue(fitBestValue);
      gIt->second->SetFitBestValueErr(fitBestValueErr);
   }

   // Retrive info about the status of the minimation
   double errdef;
   int npari, nparx;
   fMinuit->mnstat(fMinNLL,fEDM,errdef,npari,nparx,fCovQual);

   // update oscillation parameters in the propagator
   for (const auto& par : *fGlobalParMap) {
     if (par.second->IsOscillation() == false) continue;

     TString parName = par.first;
     if ( parName.Contains("Theta12") ) {
       fPropagatorInputs.x12 = par.second->GetFitBestValue();
     } else if ( parName.Contains("Theta13")) {
       fPropagatorInputs.x13 = par.second->GetFitBestValue();
     } else if ( parName.Contains("Theta23") ) {
       fPropagatorInputs.x23 = par.second->GetFitBestValue();
     } else if ( parName.Contains("deltaCP") ) {
       fPropagatorInputs.dcp = par.second->GetFitBestValue();
     } else if ( parName.Contains("dm21")  ) {
       fPropagatorInputs.dm21 = par.second->GetFitBestValue();
     } else if ( parName.Contains("dm32") ) {
       fPropagatorInputs.dm32 = par.second->GetFitBestValue();
     }
   }
  
   UpdateOscillationParameters();
}

void MSMinimizer::FCNNLLLikelihood(int & npar, double * /*grad*/,
      double &fval, double * par, int /*flag*/)
{
   fval = 0.0;
   MSModelVector* modelVector = global_pointer->fModelVector;
   NeutrinoPropagator* propagator = global_pointer->fNeutrinoPropagator;
   PropagatorInputs_t& propagator_inputs = global_pointer->fPropagatorInputs;

   propagator_inputs.SetParameters( par );

   propagator->SetMNS(propagator_inputs.x12, propagator_inputs.x13,
                      propagator_inputs.x23, 
                      propagator_inputs.dm21, propagator_inputs.dm32,
                      propagator_inputs.dcp,
                      1.0, // Energy placeholder
                      propagator_inputs.useSinSq, propagator_inputs.nubar);

   for (const auto& i : *modelVector) fval += i->NLogLikelihood(par, propagator);
}

void MSMinimizer::UpdateOscillationParameters() {
  if (!fNeutrinoPropagator) return;

  fNeutrinoPropagator->SetMNS(fPropagatorInputs.x12, fPropagatorInputs.x13,
                              fPropagatorInputs.x23, 
                              fPropagatorInputs.dm21, fPropagatorInputs.dm32,
                              fPropagatorInputs.dcp,
                              1.0, // Energy placeholder
                              fPropagatorInputs.useSinSq, fPropagatorInputs.nubar);
  return;
}

} // namespace mst
