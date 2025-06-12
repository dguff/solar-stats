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
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "Fit/ParameterSettings.h"

// m-stats libs
#include "MSMinimizer.h"

// Prob3++
#include "NeutrinoPropagator.h"
#include "BargerPropagator.h"

namespace mst {

MSMinimizer* MSMinimizer::global_pointer = 0;

MSMinimizer::MSMinimizer(const std::string& name) : MSObject(name)
{
   fModelVector        = new MSModelVector();
   fLocalParMap        = new MSParameterMap();
   fNeutrinoPropagator = new BargerPropagator(); 
   BargerPropagator* p_local = dynamic_cast<BargerPropagator*>(fNeutrinoPropagator);
   p_local->UseMassEigenstates(true); 
   p_local->SetDefaultOctant(23, 2); 
}

MSMinimizer::~MSMinimizer()
{
   if (fModelVector) {
      for (auto& i : *fModelVector) delete i;
      delete fModelVector;
   }

   if (fLocalParMap) {
      for (auto& it : *fLocalParMap) delete it.second;
      delete fLocalParMap;
   }

   if (fNeutrinoPropagator) delete fNeutrinoPropagator;
}

void MSMinimizer::SetupMiminizerOptions(const rapidjson::Value& jopts)
{
  // Check if the json object is an array
  if (!jopts.IsArray()) {
    std::cerr << "MSMinimizer::SetupMiminizerOptions: json object is not an array" << std::endl;
    exit(EXIT_FAILURE);
  }

  for (const auto &mz_itr : jopts.GetArray()) {
    fMinimizers.push_back(MSMinimizerEngine_t());
    auto& engine = fMinimizers.back();
    const auto& mz = mz_itr.GetObject();
    (mz.HasMember("type")) ? engine.fMinimizerOptions.SetMinimizerType( mz["type"].GetString()) 
      : engine.fMinimizerOptions.SetMinimizerType("Minuit");
    (mz.HasMember("algorithm")) ? engine.fMinimizerOptions.SetMinimizerAlgorithm( mz["algorithm"].GetString())
      : engine.fMinimizerOptions.SetMinimizerAlgorithm("Migrad");
    (mz.HasMember("strategy")) ? engine.fMinimizerOptions.SetStrategy( mz["strategy"].GetInt())
      : engine.fMinimizerOptions.SetStrategy(1);
    (mz.HasMember("verbosity")) ? engine.fMinimizerOptions.SetPrintLevel( mz["verbosity"].GetInt())
      : engine.fMinimizerOptions.SetPrintLevel(0);
    (mz.HasMember("errordef")) ? engine.fMinimizerOptions.SetErrorDef( mz["errordef"].GetDouble())
      : engine.fMinimizerOptions.SetErrorDef(0.5);
    (mz.HasMember("precision")) ? engine.fMinimizerOptions.SetPrecision( mz["precision"].GetDouble())
      : engine.fMinimizerOptions.SetPrecision(1e-16);
    (mz.HasMember("maxcalls")) ? engine.fMinimizerOptions.SetMaxFunctionCalls( mz["maxcalls"].GetDouble())
      : engine.fMinimizerOptions.SetMaxFunctionCalls(10000);
    (mz.HasMember("maxiterations")) ? engine.fMinimizerOptions.SetMaxIterations( mz["maxiterations"].GetDouble())
      : engine.fMinimizerOptions.SetMaxIterations(10000);
    (mz.HasMember("tolerance")) ? engine.fMinimizerOptions.SetTolerance( mz["tolerance"].GetDouble())
      : engine.fMinimizerOptions.SetTolerance(1e-4);

    if (engine.fMinimizerOptions.PrintLevel()) {
      printf("MSMinimizer::SetupMiminizerOptions() Minimizer options %ld created with the following options:\n", fMinimizers.size());
      engine.fMinimizerOptions.Print();
    }
   }
 
   return;
}

void MSMinimizer::InitializeMinimizer(const int iengine)
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

   // Initialize minimizer
   auto& minimizer = fMinimizers.at(iengine).fMinimizer;
   auto& options = fMinimizers.at(iengine).fMinimizerOptions;
   minimizer = std::unique_ptr<ROOT::Math::Minimizer>(
      ROOT::Math::Factory::CreateMinimizer(
         options.MinimizerType(),
         options.MinimizerAlgorithm()));
   minimizer->SetOptions( options ); 

   fFcn = new ROOT::Math::Functor(this, &MSMinimizer::FCNNLLLikelihood, fGlobalParMap->size());
   minimizer->SetFunction(*fFcn);

   if (options.PrintLevel() > 0) {
      printf("MSMinimizer::InitializeMinimizer() Minimizer %i created with the following options:\n", iengine);
      options.Print();
   }

   return;
}

void MSMinimizer::SyncFitParameters(const int iengine, bool resetParStartVal, bool forceUpdateAll)
{
  if (fVerbosity) {
    printf("MSMinimizer::SyncFitParameters() Minimizer %i - reset %i - force update: %i\n", 
        iengine, resetParStartVal, forceUpdateAll);
  }
   // Check whether minuit has been last sync against this minimizer
   // and define value of global pointer
   //if (global_pointer != this) {
      //forceUpdateAll = true;
      //global_pointer = this;
   //}
   if (iengine != fCurrentMinimizer) {
     forceUpdateAll = true;
   }

   // Initialize minuit if not done manually
   if (!fMinimizers.at(iengine).fMinimizer) InitializeMinimizer(iengine);
   auto& fMinimizer = fMinimizers.at(iengine).fMinimizer;

   const MSParameterMap::const_iterator gItB = fGlobalParMap->begin();
   const MSParameterMap::const_iterator gItE = fGlobalParMap->end();

   // Parse parameter to minuit (only if different from previous call)
   for (MSParameterMap::const_iterator gIt = gItB; gIt != gItE; ++gIt) {
     const int d = distance(gItB, gIt);
     MSParameterMap::const_iterator lIt = fLocalParMap->find(gIt->first);

     // Intialize all parameters the first time the function is called or
     // if minuit was synced against a different MSMinimizer
     if (lIt == fLocalParMap->end() || forceUpdateAll) {
       if (fVerbosity) std::cerr << "MSMinimizer::SyncFitParameters: forced update "
         << "par[" << d
           << "] \"" << gIt->first << "\""
           << " -> synced all fields"
           << std::endl;

       fMinimizer->SetLimitedVariable(d, 
           gIt->second->GetName().data(),
           gIt->second->GetFitStartValue(),
           gIt->second->GetFitStartStep(),
           gIt->second->GetRangeMin(),
           gIt->second->GetRangeMax());

       if (gIt->second->IsFixed()) {
         if (fVerbosity) std::cerr << "MSMinimizer::SyncFitParameters: forced update "
           << "par[" << d
             << "] \"" << gIt->first << "\""
             << " -> fixed"
             << std::endl;
         fMinimizer->FixVariable(d);
       }
     } // otherwise, if global parameter is fixed...
     else if (gIt->second->IsFixed()) {
       //if (gIt->second->GetFitStartValue() != lIt->second->GetFitStartValue()) {
         //if (fVerbosity) std::cerr << "MSMinimizer::SyncFitParameters: "
           //<< "par[" << d
             //<< "] \"" << gIt->first << "\""
             //<< " -> synced starting value"
             //<< std::endl;
         ////fMinimizerArglist[0] = d+1;
         ////fMinimizerArglist[1] = gIt->second->GetFitStartValue();
         ////fMinimizer->mnexcm("SET PAR",fMinuitArglist,2, fMinuitErrorFlag);
         //fMinimizer->SetVariableValue(d, gIt->second->GetFitStartValue());
       //}
       //if (!lIt->second->IsFixed()) {
         //if (fVerbosity) std::cerr << "MSMinimizer::SyncFitParameters: "
           //<< "par[" << d
             //<< "] \"" << gIt->first << "\""
             //<< " -> fixed"
             //<< std::endl;
         ////fMinimizerArglist[0] = d+1;
         ////fMinimizer->mnexcm("FIX", fMinuitArglist , 1, fMinuitErrorFlag);
         //fMinimizer->FixVariable(d);
       //}
       //else {
         if (fVerbosity) std::cerr << "MSMinimizer::SyncFitParameters: "
           << "par[" << d << "]\"" << gIt->first.data() << "\" -> " 
             << " -> parameter already fixed in local map - update" << std::endl;

         fMinimizer->SetFixedVariable(d, 
             gIt->second->GetName().data(),
             gIt->second->GetFitStartValue());
       //}
     } 
     else { // global parameter is not fixed
       if (lIt->second->IsFixed()) {
         if (fVerbosity) std::cerr << "MSMinimizer::SyncFitParameters: "
           << "par[" << d
             << "] \"" << gIt->first << "\""
             << " -> released"
             << std::endl;
         //fMinimizerArglist[0] = d+1;
         //fMinimizer->mnexcm("RELEASE", fMinuitArglist , 1, fMinuitErrorFlag);
         fMinimizer->ReleaseVariable(d);
       }

       if ( resetParStartVal ) {
         if (fVerbosity) std::cerr << "MSMinimizer::SyncFitParameters: reset "
           << "par[" << d
             << "] \"" << gIt->first << "\""
             << " -> synced all fields"
             << std::endl;

         //fMinimizer->mnparm ( d, gIt->second->GetName().data(),
         //gIt->second->GetFitStartValue(),
         //gIt->second->GetFitStartStep(),
         //gIt->second->GetRangeMin(),
         //gIt->second->GetRangeMax(),
         //fMinimizerErrorFlag);
         fMinimizer->SetLimitedVariable(d,
             gIt->second->GetName().data(),
             gIt->second->GetFitStartValue(),
             gIt->second->GetFitStartStep(),
             gIt->second->GetRangeMin(),
             gIt->second->GetRangeMax());

         if (gIt->second->IsFixed()) {
           fMinimizer->FixVariable(d);
         }
       }
       else /* if (  
           gIt->second->GetFitStartValue() != lIt->second->GetFitStartValue()||
           gIt->second->GetFitStartStep()  != lIt->second->GetFitStartStep() ||
           gIt->second->GetRangeMin()      != lIt->second->GetRangeMin() ||
           gIt->second->GetRangeMax()      != lIt->second->GetRangeMax() ) */
       {
         if (fVerbosity) std::cerr << "MSMinimizer::SyncFitParameters: update "
           << "par[" << d
             << "] \"" << gIt->first << "\""
             << " -> update value from global map"
             << std::endl;

         fMinimizer->SetLimitedVariable(d,
             gIt->second->GetName().data(),
             gIt->second->GetFitStartValue(),
             gIt->second->GetFitStartStep(),
             gIt->second->GetRangeMin(),
             gIt->second->GetRangeMax());

         if (gIt->second->IsFixed()) {
           fMinimizer->FixVariable(d);
         }

       }
     }
   }

   // Print current variables settings
   /*
    *int ipar = 0;
    *printf("Minimizer %i initial state:\n", iengine);
    *for (const auto& param : *fGlobalParMap) {
    * ROOT::Fit::ParameterSettings settings;
    * fMinimizer->GetVariableSettings(ipar, settings);
    * printf("par[%s] index: %i, value: %g, step: %g, min: %g, max: %g, fixed: %i\n",
    *        param.first.c_str(), ipar,
    *        settings.Value(), settings.StepSize(),
    *        settings.LowerLimit(), settings.UpperLimit(),
    *        settings.IsFixed());
    * ipar++;
    *}
    *getchar();
    */

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

void MSMinimizer::Minimize(const int iengine, bool resetFitStartValue, bool forceUpdateAll) {
   // Sync parameters
   SyncFitParameters(iengine, resetFitStartValue);
   // Run actual minimization
   //fMinimizer->mnexcm(minimizer.c_str(), fMinuitArglist, 2, fMinuitErrorFlag);
   auto& fMinimizer = fMinimizers.at(iengine).fMinimizer;
   // Dump current state of the minimizer and the parameters
   if (fVerbosity) {
     printf("MSMinimizer::Minimize() Minimizer %i state before minimization:\n", iengine);
     fMinimizer->PrintResults();
   }

   fMinimizer->Minimize();
   fCurrentMinimizer = iengine;


   if (GetMinuitStatus()) fNMinuitFails++;

   // Retrive fit results from minuit and store info
   MSParameterMap::const_iterator gIt0 = fGlobalParMap->begin();
   for (MSParameterMap::const_iterator gIt = fGlobalParMap->begin();
       gIt != fGlobalParMap->end(); ++gIt) {
     int d = distance(gIt0, gIt);
     gIt->second->SetFitBestValue( fMinimizer->X()[d] );
     gIt->second->SetFitBestValueErr( fMinimizer->Errors()[d] );

     // update oscillation parameters in the propagator
     if (gIt->second->IsOscillation() == false) continue;

     TString parName = gIt->first;
     if ( parName.Contains("Theta12") ) {
       fPropagatorInputs.x12 = gIt->second->GetFitBestValue();
     } else if ( parName.Contains("Theta13")) {
       fPropagatorInputs.x13 = gIt->second->GetFitBestValue();
     } else if ( parName.Contains("Theta23") ) {
       fPropagatorInputs.x23 = gIt->second->GetFitBestValue();
     } else if ( parName.Contains("deltaCP") ) {
       fPropagatorInputs.dcp = gIt->second->GetFitBestValue();
     } else if ( parName.Contains("dm21")  ) {
       fPropagatorInputs.dm21 = gIt->second->GetFitBestValue();
     } else if ( parName.Contains("dm32") ) {
       fPropagatorInputs.dm32 = gIt->second->GetFitBestValue();
     }
   }
   UpdateOscillationParameters();

   // Retrive info about the status of the minimation
   double errdef;
   int npari, nparx;

   fMinNLL = fMinimizer->MinValue();
   fEDM = fMinimizer->Edm();
   fCovQual = fMinimizer->CovMatrixStatus();

  
}

double MSMinimizer::FCNNLLLikelihood(const double* par)
    //(int & npar, double * [>grad<],
      //double &fval, double * par, int [>flag<])
{
   double fval = 0.0;
   MSModelVector* modelVector = global_pointer->fModelVector;
   NeutrinoPropagator* propagator = global_pointer->fNeutrinoPropagator;
   PropagatorInputs_t& propagator_inputs = global_pointer->fPropagatorInputs;

/*
 *   UInt_t ipar = 0;
 *   for (const auto& pp : *fGlobalParMap) {
 *     printf("par[%s] index: %i, value: %f\n", pp.first.c_str(), ipar, par[ipar]);
 *     ipar++;
 *   }
 *
 */
   propagator_inputs.SetParameters( par );


   propagator->SetMNS(propagator_inputs.x12, propagator_inputs.x13,
                      propagator_inputs.x23, 
                      propagator_inputs.dm21, propagator_inputs.dm32,
                      propagator_inputs.dcp,
                      1.0, // Energy placeholder
                      propagator_inputs.useSinSq, propagator_inputs.nubar);

   for (const auto& i : *modelVector) fval += i->NLogLikelihood(par, propagator);
   return fval;
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
