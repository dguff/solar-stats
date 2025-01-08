#ifndef MST_MSHistFit_CXX
#define MST_MSHistFit_CXX

// c/c++ libs
#include <csignal>
#include <cstdlib> 
#include <map>
#include <sstream>
#include <fstream>

// ROOT libs
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THn.h"
#include "TROOT.h"
#include "TFile.h"
#include "TBranch.h"
#include "TTree.h"
#include "TStyle.h"

// rapidjson's DOM-style API
#include "rapidjson/document.h"
#include "rapidjson/istreamwrapper.h"
#include "rapidjson/filewritestream.h"
#include "rapidjson/prettywriter.h"

// MARLEY
#include "marley/RootJSONConfig.hh"

// m-stats libs
#include "MSPDFBuilderTHn.h"
#include "MSModelTHnBMLF.h"
#include "MSModelPulls.h"
#include "MSMinimizer.h"

using namespace std;

namespace mst {

  /* 
   * Load and check integrity of the json config file
   */
  inline rapidjson::Document LoadConfig (const string& configFileName, bool verbose = false) {

    std::ifstream inputStream (configFileName);
    rapidjson::IStreamWrapper inputStreamWrapper (inputStream);
    rapidjson::Document json;
    if (json.ParseStream(inputStreamWrapper).HasParseError()) {
      cerr << "error in json config file: invalid syntax" << endl;
      exit(1);
    }
    inputStream.close();


    // lambda to check if a file exist and has the correct properties
    auto isMemberCorrect = [&verbose] (const rapidjson::Value& json, 
        const string& memberName, 
        const string& memberType, 
        const string& subMemberType = "",  
        const int memberSize = 0) {

      // verbose dump
      if (verbose) cout << "info: checking json config file: member " << memberName << "... ";

      // check existance of the member
      const auto member = json.FindMember(memberName.c_str());
      if (member == json.MemberEnd()) {
        cerr << "error in json config file: member " << memberSize 
          << " not found" << endl;
        exit(1);
      }

      // check member type
      auto isTypeCorrect = [] (const rapidjson::Value& member, const string& memberType){
        bool c = false;
        if      (!strcmp (memberType.c_str(), "Bool"  )) c = member.IsBool();  
        else if (!strcmp (memberType.c_str(), "Number")) c = member.IsNumber();
        else if (!strcmp (memberType.c_str(), "Int"   )) c = member.IsInt();   
        else if (!strcmp (memberType.c_str(), "Double")) c = member.IsDouble();
        else if (!strcmp (memberType.c_str(), "String")) c = member.IsString();
        else if (!strcmp (memberType.c_str(), "Array" )) c = member.IsArray(); 
        else if (!strcmp (memberType.c_str(), "Object")) c = member.IsObject();
        return c;
      };
      if (isTypeCorrect(member->value, memberType) == false) {
        cerr << "error in json config file: member " << memberName
          << " must be of type " << memberType << endl;
        exit(1);
      }

      // check array submembers
      if (strcmp(subMemberType.c_str(), "")) {
        for (const auto& i : member->value.GetArray()) {
          if (isTypeCorrect(i, subMemberType.c_str()) == false) {
            cerr << "error in json config file: members of member " << memberName
              << " must be of type " << subMemberType << endl;
            exit(1);
          }
        }
      }

      // check size
      if (memberSize > 0 && member->value.Size() != memberSize) {
        cerr << "error in json config file: member " << memberName
          << " is not of size " << memberSize << endl;
        exit(1);
      }

      if (verbose) cout << "done." << endl;
      return;
    };

    isMemberCorrect(json,"fittingModel", "Object");                             // json/fittingModel
    isMemberCorrect(json["fittingModel"], "dataSets", "Object");                // json/fittingModel/dataSets
    for (const auto& dataSet : json["fittingModel"]["dataSets"].GetObject()) {  // json/fittingModel/dataSets/*
      if (verbose) cout << "info: checking dataSet "                           //
        << dataSet.name.GetString() << endl;                   //
      isMemberCorrect(dataSet.value, "exposure", "Number");                    // json/fittingModel/dataSets/*/exposure
      isMemberCorrect(dataSet.value, "components", "Object");                  // json/fittingModel/dataSets/*/components
      if (dataSet.value.HasMember("detectorResponse")) {                       // optional block:
        isMemberCorrect(dataSet.value, "detectorResponse", "Object");         // json/fittingModel/dataSets/*/detectorResponse
        for (const auto& dr : dataSet.value["detectorResponse"].GetObject()) {// json/fittingModel/dataSets/*/detectorResponse/*
          if (verbose) cout << "info: checking detectorResponse "            //
            << dr.name.GetString() << endl;                  //
          isMemberCorrect(dr.value, "responseMatrix", "Array", "String", 2); // json/fittingModel/dataSets[>/detectorResponse/<]responseMatrix
        }
      }                                                                        
      if (dataSet.value.HasMember("nadirExposurePDF")) {                       // optional block:
        isMemberCorrect(dataSet.value, "nadirExposurePDF", "Object");
        const auto& jnadir = dataSet.value["nadirExposurePDF"];
        // print jnadir elements
        for (const auto& i : jnadir.GetObject()) {
          printf("name: %s\n", i.name.GetString());
          isMemberCorrect(i.value, "pdf",                            // json/fittingModel/dataSets/*/nadirExposurePDF/pdf
              "Array", "String", 2); 
        }
      }
      for (const auto& component : dataSet.value["components"].GetObject()) {  // json/fittingModel/dataSets/*/components/*
        if (verbose) cout << "info: checking component "                      //
          << component.name.GetString() << endl;              //
        isMemberCorrect(component.value, "type", "String");                   // json/fittingModel/dataSets/*/components/*/type
        isMemberCorrect(component.value, "global", "Bool");                   // json/fittingModel/dataSets/*/components/*/global
        isMemberCorrect(component.value, "refVal", "Number");                 // json/fittingModel/dataSets/*/components/*/refVal
        isMemberCorrect(component.value, "range", "Array", "Number", 2);      // json/fittingModel/dataSets/*/components/*/range[]
        isMemberCorrect(component.value, "fitStep", "Number");                // json/fittingModel/dataSets/*/components/*/fitStep
        isMemberCorrect(component.value, "pdf", "Array", "String", 2);        // json/fittingModel/dataSets/*/components/*/pdf[]
        isMemberCorrect(component.value, "injVal", "Number");                 // json/fittingModel/dataSets/*/components/*/injVal
        isMemberCorrect(component.value, "color", "Int");                     // json/fittingModel/dataSets/*/components/*/color
        const string type = component.value["type"].GetString();              
        if ( std::strcmp(type.c_str(), "neutrino") == 0 ) {                    
          if (component.value.HasMember("oscillation")) {                       // 
            isMemberCorrect(component.value, "oscillation", "Bool");           // json/fittingModel/dataSets/*/components/*/oscillation
          }
          isMemberCorrect(component.value, "channels", "Object");     // json/fittingModel/dataSets/*/components/*/channels[]
          for (const auto& channel : component.value["channels"].GetObject()) { // json/fittingModel/dataSets/*/components/*/channels/*
            isMemberCorrect(channel.value, "responseMatrix", "String");         // json/fittingModel/dataSets/*/components/*/channels/*/responseMatrix
            isMemberCorrect(channel.value, "crossSection", "Object");           // json/fittingModel/dataSets/*/components/*/channels/*/crossSection
            const auto& jcrossSection = channel.value["crossSection"];            //
            isMemberCorrect(jcrossSection, "target", "Object");   // json/fittingModel/dataSets/*/components/*/channels/*/crossSection/target
            const rapidjson::Value& jtarget = jcrossSection["target"];
            isMemberCorrect(jtarget, "nuclides", "Array", "Number"); // json/fittingModel/dataSets/*/components/*/channels/*/crossSection/*/target/nuclides[]
            isMemberCorrect(jtarget, "atom_fractions", "Array", "Number"); // json/fittingModel/dataSets/*/components/*/channels/*/crossSection/*/target/atom_fractions[]
            isMemberCorrect(jcrossSection, "reactions", "Array", "String");       // json/fittingModel/dataSets/*/components/*/channels/*/crossSection/*/reaction
            isMemberCorrect(jcrossSection, "source", "Object");                   // json/fittingModel/dataSets/*/components/*/channels/*/crossSecrion/*/source
            isMemberCorrect(jcrossSection["source"], "neutrino", "String");       // json/fittingModel/dataSets/*/components/*/channels/*/crossSection/*/source/neutrino
          }
        } // end of neutrino block
                                                                             //
        if (component.value.HasMember("responseMatrix")) {                    // 
          isMemberCorrect(component.value, "responseMatrix", "String");      // json/fittingModel/dataSets/*/components/*/responseMatrix
        }                                                                     //
      }                                                                        //
      isMemberCorrect(dataSet.value, "projectOnAxis", "Array", "Int");         // json/fittingModel/dataSets/*/projectOnAxis[]
      isMemberCorrect(dataSet.value, "axis", "Object");                        // json/fittingModel/dataSets/*/axis
      for (const auto& axis : dataSet.value["axis"].GetObject()) {             // json/fittingModel/dataSets/*/axis/*
        if (verbose) cout << "info: checking axis "                           //
          << axis.name.GetString() << endl;              //
        isMemberCorrect(axis.value, "range", "Array", "Number", 2);           // json/fittingModel/dataSets/*/axis/*/range[]
        isMemberCorrect(axis.value, "rebin", "Int");                          // json/fittingModel/dataSets/*/axis/*/rebin
        if (axis.value.HasMember("label")) {
          isMemberCorrect(axis.value, "label", "String");                     // json/fittingModel/dataSets/*/axis/*/label
        }
        if (axis.value.HasMember("limits")) {                                   // optional block:
          isMemberCorrect(axis.value, "limits", "Array", "Number", 2);                      // json/fittingModel/dataSets/*/axis/*/xmin
        }
        if (axis.value.HasMember("nbins")) {                                    // optional block:
          isMemberCorrect(axis.value, "nbins", "Int");                         // json/fittingModel/dataSets/*/axis/*/nbin
        }
      }                                                                        // 
      if (dataSet.value.HasMember("normalizePDFInUserRange")) {                // optional block:
        isMemberCorrect(dataSet.value, "normalizePDFInUserRange", "Bool");    // json/fittingModel/dataSets/*/normalizePDFInUserRange
      }                                                                        //
    }                                                                           //
    if (json.HasMember("pulls")) {                                              // optional block:
      isMemberCorrect(json, "pulls", "Object");                                // json/pulls
      for (const auto& pull : json["pulls"].GetObject()) {                     // json/pulls/*
        if (verbose) cout << "info: checking pull "                           // 
          << pull.name.GetString() << endl;                   //
        isMemberCorrect(pull.value, "type", "String");                        // json/pulls/*/type
        isMemberCorrect(pull.value, "centroid", "Number");                    // json/pulls/*/centroid
        isMemberCorrect(pull.value, "sigma", "Number");                       // json/pulls/*/sigma
      }                                                                        //
    }
    if (json["fittingModel"].HasMember("oscillation")) {
      isMemberCorrect(json["fittingModel"], "oscillation", "Object");          // json/fittingModel/oscillation
      const auto& joscillation = json["fittingModel"]["oscillation"];          //
      isMemberCorrect(joscillation, "parameters", "Object");                   // json/fittingModel/oscillation/parameters
      for (const auto& par : joscillation["parameters"].GetObject()) {
        isMemberCorrect(par.value, "value", "Number");                         // json/fittingModel/oscillation/*/parameters/*/value
        isMemberCorrect(par.value, "range", "Array", "Number", 2);             // json/fittingModel/oscillation/*/parameters/*/range[]
        isMemberCorrect(par.value, "fitStep", "Number");                       // json/fittingModel/oscillation/*/parameters/*/fitStep
        isMemberCorrect(par.value, "fixed", "Bool");                           // json/fittingModel/oscillation/*/parameters/*/fixed
        isMemberCorrect(par.value, "refVal", "Number");                        // json/fittingModel/oscillation/*/parameters/*/refVal
      }
    }
    isMemberCorrect(json, "MinimizerSteps", "Object");                          // json/MinimizerSteps
    for (const auto& step : json["MinimizerSteps"].GetObject()) {               // json/MinimizerSteps/*
      if (verbose) cout << "info: checking step "                           //  
        << step.name.GetString() << endl;                   //
      isMemberCorrect(step.value, "method", "String");                         // json/MinimizerSteps/*
      isMemberCorrect(step.value, "resetMinuit", "Bool");                      // json/MinimizerSteps/*/resetMinuit
      isMemberCorrect(step.value, "maxCall", "Number");                        // json/MinimizerSteps/*/maxCall
      isMemberCorrect(step.value, "tollerance", "Number");                     // json/MinimizerSteps/*/tollerance
      isMemberCorrect(step.value, "verbosity", "Int");                         // json/MinimizerSteps/*/verbosity
    }                                                                           //
    if (json.HasMember("MC")) {                                                 // optional block:
      isMemberCorrect(json, "MC", "Object");                                   // json/MC
      isMemberCorrect(json["MC"], "procedure", "String");                   // json/MC/procedure
      isMemberCorrect(json["MC"], "realizations", "Number");                   // json/MC/realizations
      isMemberCorrect(json["MC"], "seed", "Int");                              // json/MC/seed
      isMemberCorrect(json["MC"], "enablePoissonFluctuations", "Bool");        // json/MC/enablePoissonFluctuations
      isMemberCorrect(json["MC"], "outputFile", "String");                     // json/MC/outputFile
    }

    return json;
  }

  inline std::string BuildMARLEYGeneratorConfig(const rapidjson::Value& jcrossSection, long int run_seed) {
    rapidjson::Document mdoc; 
    rapidjson::Document::AllocatorType& allocator = mdoc.GetAllocator(); 
    mdoc.CopyFrom(jcrossSection, allocator);
    mdoc["source"].AddMember("type", "monoenergetic", allocator);
    mdoc["source"].AddMember("energy", 10.0, allocator);

    FILE* marley_cfg_tmp; 
    char writeBuffer[65536];
    char filePathBuffer[200]; 
    std::snprintf(filePathBuffer, 200, "/tmp/solarstats_marley_config%ld.json", run_seed);
    std::string marley_cfg_path = filePathBuffer;
    marley_cfg_tmp = std::fopen(marley_cfg_path.data(), "w"); 
    rapidjson::FileWriteStream buffer(marley_cfg_tmp, writeBuffer, sizeof(writeBuffer)); 
    rapidjson::PrettyWriter<rapidjson::FileWriteStream> writer(buffer); 
    mdoc.Accept(writer);
    fclose(marley_cfg_tmp); 

    printf("Setting up marley generator with configuration:\n");
    std::ifstream tmp_file(marley_cfg_path.data());  // Apri il file

    if (!tmp_file.is_open()) {
      std::cerr << "Unable to open temporary marley configuration file" << std::endl;
      exit(EXIT_FAILURE);
    }
    std::string line;
    while (std::getline(tmp_file, line)) {
      std::cout << line << std::endl;
    }
    tmp_file.close();

    return filePathBuffer;
  }

  /* 
   * Initialize all analysis structures, i.e.: the minimizer, the PDFBuilder and
   * the statistical models composing the likelihood. 
   */
  inline MSMinimizer* InitializeAnalysis (const rapidjson::Document& json,  
      const std::string& datafileName) {

    // initialize fitter
    MSMinimizer* fitter = new MSMinimizer();

    if (json.HasMember("fittingModel") == false) {
      cerr << "error in json config file: member fittingModel not found" << endl;
      exit(1);
    }
    const auto& jfittingModel = json["fittingModel"];

    MSParameterMap oscillationParMap;
    if (jfittingModel.HasMember("oscillation")) {
      const auto& joscillation = jfittingModel["oscillation"];
      for (const auto& par : joscillation["parameters"].GetObject()) {
        auto p = new mst::MSParameter(par.name.GetString());
        p->SetFixed(par.value["fixed"].GetBool());
        p->SetRange(par.value["range"].GetArray()[0].GetDouble(),
            par.value["range"].GetArray()[1].GetDouble());
        p->SetGlobal( true );
        p->SetFitStartStep(par.value["fitStep"].GetDouble());
        p->SetFitStartValue(par.value["refVal"].GetDouble());
        p->SetOscillation();
        oscillationParMap.insert(MSParameterPair(p->GetName(), p));
      }
    }

    // loop over data sets and for each create the model and PDFBuilder
    for (const auto& dataSet : jfittingModel["dataSets"].GetObject()) {

      // Initilize analysis model and associate to them the pdfBuilder.
      // The name of the model is needed to retrieve it from the fitter
      // and to retrieve its local parameters
      MSModelTHnBMLF* mod = new MSModelTHnBMLF(dataSet.name.GetString());

      // Set the exposure
      mod->SetExposure(dataSet.value["exposure"].GetDouble());

      // Create a separate pdfBuilder for each model. The pointer of each pdfBuilder 
      // will be associated to the model.
      MSPDFBuilderTHn* pdfBuilder = new MSPDFBuilderTHn(dataSet.name.GetString());

      // Initialize hist handler //////////////////////////////////////////////
      MSTHnHandler& handler = pdfBuilder->GetHistHandler();

      // read list of projected axis from the json file 
      // (only if the block projectOnAxis is defined)
      if (dataSet.value.HasMember("projectOnAxis")) {
        std::vector<int> v;
        for (int i = 0; i < dataSet.value["projectOnAxis"].Size(); i++) 
          v.push_back((dataSet.value["projectOnAxis"].GetArray())[i].GetInt());
        handler.ProjectToAxis(v);
      }

      // set the range and binning of the axis
      if (dataSet.value.HasMember("axis")) {
        for (const auto& axis : dataSet.value["axis"].GetObject()) {
          // the stringstream is used just to convert a string into an integer
          stringstream conversion; 
          int axisID = 0;
          conversion << axis.name.GetString();
          conversion >> axisID;
          handler.SetRange(axisID, axis.value["range"][0].GetDouble(), 
              axis.value["range"][1].GetDouble());
          handler.Rebin(axisID, axis.value["rebin"].GetInt());
          if (axis.value.HasMember("label")) {
            handler.SetLabel(axisID, axis.value["label"].GetString());
          }
          if (axis.value.HasMember("limits")) {
            handler.SetLimits(axisID, 
                axis.value["limits"][0].GetDouble(), 
                axis.value["limits"][1].GetDouble());
          }
          if (axis.value.HasMember("nbins")) {
            handler.SetNbins(axisID, axis.value["nbins"].GetInt());
          }
        }
      }

      // Renormilize the histograms, in a specic range or over the full axis
      // (ecluding over- and under-flow bins)
      if (dataSet.value.HasMember("normalizePDFInUserRange"))
        handler.RespectAxisUserRange(dataSet.value["normalizePDFInUserRange"].GetBool());

      // end hist handler ////////////////////////////////////////////////////

      // Set seed of the random number generator. Note that each PSDBuilder is
      // initialized with the same seed. This will produce the same data set if
      // the parameters of the data sets are all exactly the same.
      if (json.HasMember("MC")) pdfBuilder->SetSeed(json["MC"]["seed"].GetInt());

      // Load the response matrix if present
      if (dataSet.value.HasMember("detectorResponse")) {
        for (const auto& dr : dataSet.value["detectorResponse"].GetObject()) {
          TString pathToFile (getenv("M_STATS_HIST_FIT"));
          if (pathToFile != "") pathToFile += "/";
          pathToFile += dr.value["responseMatrix"][0].GetString();
          pdfBuilder->RegisterResponseMatrix(
              handler.LoadHist(pathToFile.Data(),
                dr.value["responseMatrix"][1].GetString(),
                dr.name.GetString())); 
        }
      }

      // Load the PDF with the nadir exposure if present
      if (dataSet.value.HasMember("nadirExposurePDF")) {
        TString pathToFile (getenv("M_STATS_HIST_FIT"));
        if (pathToFile != "") pathToFile += "/";
        // print path to nadir pdf file
        const auto& jnadir = dataSet.value["nadirExposurePDF"].GetObject();

        for (const auto& i : jnadir){
          pathToFile += i.value["pdf"][0].GetString();
          THn* hnNadir = handler.LoadHist(pathToFile.Data(),
              i.value["pdf"][1].GetString(),
              "nadirExposurePDF", true);
          pdfBuilder->RegisterNadirPDF( hnNadir->Projection(0) );
          delete hnNadir;
        }
      }

      //check if handler has a nadir axis
      const auto& axes_list = handler.GetProjectID();
      bool has_nadir = false;
      THn* hnNadir = nullptr;
      for (const auto& axis_id : axes_list) {
        TString axis_name = handler.GetAxes().at(axis_id).fLabel;
        if (axis_name.Contains("nadir", TString::ECaseCompare::kIgnoreCase)) {
          has_nadir = true;
          hnNadir = pdfBuilder->BuildNadirPDF(); 
          break;
        }
      }

      // Assign the pointer to the neutrino propagator and initialize the oscillation parameters
      if (jfittingModel.HasMember("oscillation")) {
        printf("adding oscillation paramters to model %s\n", mod->GetName().data()); 
        mod->SetNeutrinoPropagator( fitter->GetNeutrinoPropagator() ); 
        for (const auto& par : oscillationParMap) {
          mod->AddParameter(par.second);
        }
      }

      // Initialize the parameters of each model, set their properties and load
      // the pdf's
      for (const auto& component : dataSet.value["components"].GetObject()) {

        auto par = new mst::MSParameter(component.name.GetString());
        par->SetGlobal(component.value["global"].GetBool());
        par->SetFixed(component.value["fixed"].GetBool());
        par->SetRange(component.value["range"].GetArray()[0].GetDouble(),
            component.value["range"].GetArray()[1].GetDouble());

        // Check if fitStep is registered and is different from zero
        if (component.value["fitStep"].GetDouble())
          par->SetFitStartStep(component.value["fitStep"].GetDouble());
        // otherwise use the range to guess to the the fit starting step
        else 
          par->SetFitStartStep( ((component.value["range"].GetArray())[1].GetDouble()
                -  (component.value["range"].GetArray())[0].GetDouble())/100.);

        par->SetFitStartValue(component.value["refVal"].GetDouble());
        par->SetInput();

        mod->AddParameter(par);

        // Load histograms for each component and possibly project it on a sub-set
        // of the axis
        TString pathToFile (getenv("M_STATS_HIST_FIT"));
        if (pathToFile != "") pathToFile += "/";
        pathToFile +=component.value["pdf"][0].GetString();

        EPDFType pdfType = get_pdf_type_code(component.value["type"].GetString());
        MSTHnPDF* pdf = nullptr;

        switch (pdfType) {
          case kComponent: 
            {
              MSTHnPDFComponent* pdf_ = new MSTHnPDFComponent(component.name.GetString());
              pdf_->SetPDFType(pdfType);
              THn* hn = nullptr;

              if (component.value.HasMember("responseMatrix")) {
                hn = handler.LoadHist( pathToFile.Data(),
                      component.value["pdf"][1].GetString(),
                      component.name.GetString(), true);
                pdf_->SetRespMatrix(pdfBuilder->GetResponseMatrix(component.value["responseMatrix"].GetString()));
              }
              else {
                hn = handler.LoadHist( pathToFile.Data(),
                      component.value["pdf"][1].GetString(),
                      component.name.GetString(),
                      true);
              }

              if ( has_nadir ) {
                printf("factorizing nadir pdf for %s\n", component.name.GetString());
                THn* hn_tmp = handler.FactorizeTHn( hn, hnNadir );
                handler.NormalizeHn( hn_tmp );
                delete hn; 
                hn = hn_tmp;

                printf("hn dimensions: %d\n", hn->GetNdimensions());
                for (int idim =0; idim < hn->GetNdimensions() ; idim++) {
                  printf("\taxis %d: %d bins - [%g, %g]\n", idim, hn->GetAxis(idim)->GetNbins(), hn->GetAxis(idim)->GetXmin(), hn->GetAxis(idim)->GetXmax());
                }
                getchar();
              }
              pdf_->SetTHn( hn );


              pdf = pdf_;
              break;
            }

          case kNeutrino: 
            {
              MSTHnPDFNeutrino* pdf_ = new MSTHnPDFNeutrino(component.name.GetString());
              pdf_->SetPDFType(pdfType);
              THn* hn = nullptr;
              if (component.value.HasMember("oscillation")) {
                pdf_->SetApplyOscillation(component.value["oscillation"].GetBool());
              }

              hn = handler.LoadHist( pathToFile.Data(),
                      component.value["pdf"][1].GetString(),
                      component.name.GetString(), true);
              if ( has_nadir ) {
                printf("factorizing nadir pdf for %s\n", component.name.GetString());
                THn* hn_tmp = handler.FactorizeTHn( hn, hnNadir );
                handler.NormalizeHn( hn_tmp );
                delete hn; 
                hn = hn_tmp;

                printf("hn dimensions: %d\n", hn->GetNdimensions());
                for (int idim =0; idim < hn->GetNdimensions() ; idim++) {
                  printf("\taxis %d: %d bins - [%g, %g]\n", idim, hn->GetAxis(idim)->GetNbins(), hn->GetAxis(idim)->GetXmin(), hn->GetAxis(idim)->GetXmax());
                }
                getchar();
              }
              pdf_->SetTHn( hn ); 

              for (const auto& jchannel : component.value["channels"].GetObject()) {
                const string channelName = jchannel.name.GetString();
                const string respMatrixName = jchannel.value["responseMatrix"].GetString();
                const string neutrino_str = jchannel.value["crossSection"]["source"]["neutrino"].GetString();
                const int neutrino_pdg = std::stoi(neutrino_str);
                MSTHnPDFNeutrino::NuIntChannel_t& channel =  
                  pdf_->AddChannel(channelName, 
                      neutrino_pdg,
                      respMatrixName); 

                const rapidjson::Value& jcrossSection = jchannel.value["crossSection"];
                const std::string gen_config = BuildMARLEYGeneratorConfig(jcrossSection, json["MC"]["seed"].GetInt());
                pdfBuilder->SetupMarleyGenerator(channelName, gen_config);

                pdfBuilder->EvaluateTotalCrossSection( channelName, pdf_ ); 
                pdf = pdf_; 
              }

              break;
            }

          case kUndefined:
            {
              std::cerr << "ERROR setting component " << component.name.GetString() << ": pdf type undefined\n";
              exit(EXIT_FAILURE);
            }

          default:
            {
              std::cerr << "ERROR setting component " << component.name.GetString() << ": pdf not set\n";
              exit(EXIT_FAILURE);
            }
        }

        pdfBuilder->RegisterPDF( pdf );
      }

      // Move the pointer of the pdfBuilder to the model
      mod->SetPDFBuilder(pdfBuilder);

      // Set data if loaded from file
      if (!datafileName.empty()) {
        std::cout << "info: loading input data from " << datafileName << std::endl;
        mod->SetDataSet( handler.BuildHist(datafileName, mod->GetName(), mod->GetName(), false) );
      }

      // Move pointer of the model to the fitter
      fitter->AddModel(mod);

      if (hnNadir) delete hnNadir;
    }

    // Add models implementing pulls if requested in the config file
    // FIXME: check if pull are present
    if (json["fittingModel"].HasMember("pulls")) {
      for (const auto& pull : json["fittingModel"]["pulls"].GetObject()) {
        // initialize and add gaussian pulls
        if (strcmp(pull.value["type"].GetString(), "gauss") == 0) {
          MSModelPullGaus* mod = new MSModelPullGaus(pull.name.GetString());
          mod->SetPullPar(pull.name.GetString());
          mod->SetGaussPar(pull.value["centroid"].GetDouble(),
              pull.value["sigma"].GetDouble());
          fitter->AddModel(mod);
          // initialize and add exponential pulls
        } else if (strcmp(pull.value["type"].GetString(), "exp") == 0) {
          MSModelPullExp* mod = new MSModelPullExp(pull.name.GetString());
          mod->SetPullPar(pull.name.GetString());
          mod->SetExpPar(pull.value["limit"].GetDouble(),
              pull.value["quantile"].GetDouble());
          fitter->AddModel(mod);
        } else {
          std::cerr << "error: pull type not implemented.\n";
          exit(1);
        }
      }
    }

    // Sync the parameters. This call is needed to finilize the initializatoin of
    // the minimizer
    fitter->SyncFitParameters();
    return fitter;
  }

  /*
   * Create data sets and automatically associate it to the models
   */
  inline bool SetDataSetFromMC (const rapidjson::Document& json, MSMinimizer* fitter) {
    // initialze the oscillation parameters (if any)
    fitter->UpdateOscillationParameters(); 
    if ( json["fittingModel"].HasMember("oscillation") ) {
      auto& oscillation_pars = fitter->GetPropagatorInputs(); 
      const auto& joscillation = json["fittingModel"]["oscillation"];
      for (const auto& par : joscillation["parameters"].GetObject()) {
        const TString parName = par.name.GetString();
        const double trueVal = par.value["value"].GetDouble();
        if (parName.Contains("Theta12", TString::ECaseCompare::kIgnoreCase)) {
          oscillation_pars.x12 = trueVal;
        } else if (parName.Contains("Theta13", TString::ECaseCompare::kIgnoreCase)) {
          oscillation_pars.x13 = trueVal;
        } else if (parName.Contains("Theta23", TString::ECaseCompare::kIgnoreCase)) {
          oscillation_pars.x23 = trueVal;
        } else if (parName.Contains("dm21", TString::ECaseCompare::kIgnoreCase)) {
          oscillation_pars.dm21 = trueVal;
        } else if (parName.Contains("dm32", TString::ECaseCompare::kIgnoreCase)) {
          oscillation_pars.dm32 = trueVal;
        } else if (parName.Contains("DeltaCP", TString::ECaseCompare::kIgnoreCase)) {
          oscillation_pars.dcp = trueVal;
        }
      }
      fitter->UpdateOscillationParameters();
    }

    // loop over the models and create a new data set for each of them
    for (const auto& model: *fitter->GetModels()) {

      // filter only models that  have a data sets (no pulls)
      const auto mod = dynamic_cast<MSModelTHnBMLF*>(model);
      if(mod == nullptr) continue;

      // get the specific pdfBuilder and reset it
      const auto pdfBuilder = mod->GetPDFBuilder();
      pdfBuilder->ResetPDF();

      // add hists to pdfBuilder with the desired rate and copute the number 
      // of counts to extract to create the data set
      double totalCounts = 0;
      for (const auto& parName: *mod->GetLocalParameters()) {
        const auto par = mod->GetParameter(parName);
        if (par->IsInput() == false) continue;
        const double trueVal = json["fittingModel"]["dataSets"][mod->GetName().c_str()]
          ["components"][parName.c_str()]["injVal"].GetDouble();

        printf("calling AddHistToPDF with par=%s, trueVal=%f and passing propagator %p\n", 
        parName.c_str(), trueVal, fitter->GetNeutrinoPropagator());
        double count_rate = pdfBuilder->AddHistToPDF(parName.c_str(), trueVal, fitter->GetNeutrinoPropagator());
        printf("count_rate = %f\n", count_rate);
        totalCounts += count_rate * mod->GetExposure();
        getchar();
      }
      // register new data set. The previous one is delete inside the model
      // Define whether to add poission fluctuation to the number of events
      printf("Generating new dataset with %f counts\n", totalCounts);
      const std::string sampling_procedure = json["MC"]["procedure"].GetString();
      mod->SetDataSet(pdfBuilder->GetMCRealizaton(
            totalCounts, 
            json["MC"]["enablePoissonFluctuations"].GetBool(), 
            sampling_procedure));
      getchar();
    }
    return true;
  }

    /*
     * Minimization of the likelihood
     */
    inline bool Minimize (const rapidjson::Document& json, MSMinimizer* fitter) {

      // Take Minuit calls from config file in the proper order
      for (const auto& step : json["MinimizerSteps"].GetObject()) {
        fitter->SetMinuitVerbosity(step.value["verbosity"].GetInt());
        fitter->Minimize(step.value["method"].GetString(),
            step.value["resetMinuit"].GetBool(),
            step.value["maxCall"].GetDouble(),
            step.value["tollerance"].GetDouble());
      }

      if (fitter->GetMinuitStatus()) {
        std::cerr << "MSMinimizer: minuit return status=" << fitter->GetMinuitStatus()
          << ", indicating problems in the convergence of the fit\n"; 
      }
      return true;
    }

    /*
     * General routine for plotting the results
     */
    inline TCanvas* GetCanvasFit (const rapidjson::Document& json, const MSMinimizer* fitter) {

      // Retrieve the canvas and define number of canvases
      int nMaxDim = 0;
      int nMaxMod = 0;
      delete gROOT->GetListOfCanvases()->FindObject("cMLF");
      // Look for the max number of dimensions 
      for (int modelNum = 0; modelNum < fitter->GetNModels(); modelNum++) {
        const auto mod = dynamic_cast<MSModelTHnBMLF*>(fitter->GetModels()->at(modelNum));
        if (mod == nullptr) continue;

        nMaxMod++;
        if (mod->GetDataSet()->GetNdimensions() > nMaxDim) 
          nMaxDim = mod->GetDataSet()->GetNdimensions();
      }
      TCanvas* cc = new TCanvas("cMLF","Max likelihood fit",400*nMaxDim, 400*nMaxMod); 
      cc->Divide(nMaxDim, nMaxMod);

      // Set general style
      gStyle->SetOptStat(0);
      gStyle->SetTitle(0);

      // index storing the number of models added to the canvas
      int m = 0;
      // loop over the models
      for (const auto& model : *fitter->GetModels())  {
        const auto mod = dynamic_cast<MSModelTHnBMLF*>(model);
        if (mod == nullptr) continue; 

        // retrieve structures to be plot
        const auto pdfBuilder   = mod->GetPDFBuilder();
        const auto dataHist     = mod->GetDataSet();
        const auto localParList = mod->GetLocalParameters();

        // loop over the dimensions
        for (int d = 0; d < mod->GetDataSet()->GetNdimensions(); d++) {
          if (dataHist->GetNdimensions()<d) break;
          // set the pad
          cc->cd(1 + m*nMaxDim + d);

          // Plot first the original data to set the plots ranges
          auto data_pr = dataHist->Projection(d);
          data_pr->SetLineWidth(0);
          data_pr->SetMinimum(0.01);
          data_pr->DrawCopy();
          gPad->SetLogy();

          // Draw best fit components
          for (auto& parName : (*localParList)) {
            pdfBuilder->ResetPDF();
            mst::MSParameter*  par = mod->GetParameter(parName.c_str());
            if (par->IsInput() == false) continue;
            pdfBuilder->AddHistToPDF(parName, mod->GetExposure()*par->GetFitBestValue(), fitter->GetNeutrinoPropagator());
            auto pdf_tmp = pdfBuilder->GetPDF("");

            auto pdf_tmp_pr = pdf_tmp->Projection(d);
            pdf_tmp_pr->SetName(parName.c_str());
            pdf_tmp_pr->SetTitle(parName.c_str());
            int color = json["fittingModel"]["dataSets"][mod->GetName().c_str()]
              ["components"][parName.c_str()]["color"].GetInt();
            pdf_tmp_pr->SetLineColor(color);
            pdf_tmp_pr->SetLineWidth(2);
            pdf_tmp_pr->DrawCopy("hist same");

            delete pdf_tmp_pr;
            delete pdf_tmp;
          }

          // Build and Draw best fit
          pdfBuilder->ResetPDF();
          for (auto&i: (*localParList)) {
            mst::MSParameter*  par = mod->GetParameter(i);
            if (par->IsInput() == false) continue;
            pdfBuilder->AddHistToPDF(i, mod->GetExposure()*par->GetFitBestValue(), fitter->GetNeutrinoPropagator());
          }
          auto tot = pdfBuilder->GetPDF("tot");
          auto tot_pr = tot->Projection(d);
          tot_pr->SetLineColor(kGray);
          tot_pr->SetLineWidth(2);
          tot_pr->DrawCopy("hist same");
          delete tot_pr;
          delete tot;

          // Plot data
          data_pr->SetMarkerColor(kBlack);
          data_pr->SetLineColor(kBlack);
          data_pr->SetMarkerStyle(20);
          data_pr->SetMarkerSize(0.7);
          data_pr->DrawCopy("PE hist same");
          delete data_pr;
        }
        // increase model index
        m++;
      }
      cc->Update();
      return cc;
    }
    
    /*
     * Build profile likelihood scan for a specific parameter
     */
    inline TGraph* Profile (const rapidjson::Document& json, MSMinimizer* fitter, 
        const string& parName, const double NLL, const int nPts) {

      printf("Building profile for parameter %s\n", parName.c_str());
      // retrieve parameter of interest (poi) from the fitter
      mst::MSParameter* poi  = fitter->GetParameter(parName.c_str());
      if (!poi) {
        std::cerr << "Profile >> error: parameter " << parName << " not found\n";
        return nullptr;
      }

      // save best fit values to restore the status of the parameters after the
      // scanning
      vector<double> fitBestValue    (fitter->GetParameterMap()->size(), -1);
      vector<double> fitBestValueErr (fitter->GetParameterMap()->size(), -1);
      {
        int parIndex = 0;
        for ( auto it : *fitter->GetParameterMap()) {
          fitBestValue.at(parIndex)    = it.second->GetFitBestValue();
          fitBestValueErr.at(parIndex) = it.second->GetFitBestValueErr();
          parIndex++;
        }
      }

      // Initialize output graph
      TGraph* gpll = new TGraph();
      // store absolute minimum of the likelihood
      double absMinNLL = std::numeric_limits<double>::max();

      // define auxiliary lambda function for profiling, which has visibility over
      // all variables defined up to now
      auto Scan = [&] (double tVal) {
        if (tVal < poi->GetRangeMin() || tVal > poi->GetRangeMax()) return false;
        poi->FixTo(tVal);

        for (const auto& step : json["MinimizerSteps"].GetObject()) {
          fitter->SetMinuitVerbosity(step.value["verbosity"].GetInt());
          fitter->Minimize(step.value["method"].GetString(),
              step.value["resetMinuit"].GetBool(),
              step.value["maxCall"].GetDouble(),
              step.value["tollerance"].GetDouble());
        }


        // extract temporary best fit value
        const double tmpMinNLL = fitter->GetMinNLL();
        gpll->SetPoint(gpll->GetN(), tVal, tmpMinNLL);
        // update absolute minimum if needed
        if (tmpMinNLL < absMinNLL) absMinNLL = tmpMinNLL;

        // check the status of minuit
        if (fitter->GetMinuitStatus()) {
          std::cerr << "Profile >> error: minuit returned failed status ["
            << fitter->GetMinuitStatus() << "]"
            << " while fitting with " << parName 
            << " fixed to " << tVal << std::endl;
        } 
        // chek if the exit conditions are met
        if (tmpMinNLL - absMinNLL <= NLL)  return true;
        else return false;
      };

      // perform actual scan
      {
        // index used to define the points to scans
        // best fit value and error
        const double poiFitBestValue    = poi->GetFitBestValue();
        const double poiFitBestValueErr = poi->GetFitBestValueErr();

        const double step = 2*poiFitBestValueErr / double(nPts);
        // scan to the right of the min
        int counter = 0;
        fitter->SyncFitParameters(true);
        while (Scan(poiFitBestValue + counter*step) && counter<10*nPts) counter++;
        // scan to the left of the min starting from -1 to not add again the best fit
        // value in the TGraph
        counter = -1;
        fitter->SyncFitParameters(true);
        while (Scan(poiFitBestValue + counter*step) && counter<10*nPts) counter--;

        // normilize profile to the absolute minimum found during while profiling
        for (int i = 0; i < gpll->GetN(); i++ ) {
          gpll->SetPoint(i, gpll->GetX()[i], (gpll->GetY()[i]-absMinNLL));
        }
      }

      // reset parameter original status
      // restore best fit values of the parameters
      {
        int parIndex = 0;
        for ( auto it : *fitter->GetParameterMap()) {
          it.second->SetFitBestValue(fitBestValue.at(parIndex));
          it.second->SetFitBestValueErr(fitBestValueErr.at(parIndex));
          parIndex++;
        }
      }
      poi->Release();

      // Set titles (this must be done after filling the TGraph. Probably it's a
      // bug of ROOT
      gpll->SetName(Form("nll_%s", parName.c_str())); 
      gpll->GetXaxis()->SetTitle(Form("%s rate [cts/100T/d]", parName.c_str())); 
      gpll->GetYaxis()->SetTitle("-LogLikelihood)"); 
      gpll->Sort();
      return gpll;
    }


    /*
     * Build profile likelihood scan for a pair of specific parameters
     */
    inline TH2D* Profile2D(const rapidjson::Document& json, MSMinimizer* fitter, 
        const string& parName1, const string& parName2, const double NLL, 
        const int nPts1, const int nPts2) 
    {

      // retrieve parameters of interest (poi) from the fitter
      mst::MSParameter* poi1  = fitter->GetParameter(parName1.c_str());
      if (!poi1) {
        std::cerr << "Profile >> error: parameter " << parName1 << " not found\n";
        return nullptr;
      }
      mst::MSParameter* poi2  = fitter->GetParameter(parName2.c_str());
      if (!poi2) {
        std::cerr << "Profile >> error: parameter " << parName2 << " not found\n";
        return nullptr;
      }
      poi1->SetRange( poi1->GetFitBestValue() - NLL*poi1->GetFitBestValueErr(),
          poi1->GetFitBestValue() + NLL*poi1->GetFitBestValueErr());
      poi2->SetRange( poi2->GetFitBestValue() - NLL*poi2->GetFitBestValueErr(),
          poi2->GetFitBestValue() + NLL*poi2->GetFitBestValueErr());

      printf("Building profile likelihood for %s and %s\n", parName1.c_str(), parName2.c_str());
      printf("parameter %s range: [%g, %g]\n", parName1.c_str(), poi1->GetRangeMin(), poi1->GetRangeMax());
      printf("parameter %s range: [%g, %g]\n", parName2.c_str(), poi2->GetRangeMin(), poi2->GetRangeMax());

      // save best fit values to restore the status of the parameters after the
      // scanning
      vector<double> fitBestValue    (fitter->GetParameterMap()->size(), -1);
      vector<double> fitBestValueErr (fitter->GetParameterMap()->size(), -1);
      {
        int parIndex = 0;
        for ( auto it : *fitter->GetParameterMap()) {
          fitBestValue.at(parIndex)    = it.second->GetFitBestValue();
          fitBestValueErr.at(parIndex) = it.second->GetFitBestValueErr();
          parIndex++;
        }
      }

      // Initialize output map
      TH2D* hpll = new TH2D("hpll", "Profiled likelihood", 
          nPts1, poi1->GetRangeMin(), poi1->GetRangeMax(), 
          nPts2, poi2->GetRangeMin(), poi2->GetRangeMax());
      // store absolute minimum of the likelihood
      double absMinNLL = std::numeric_limits<double>::max();

      // define auxiliary lambda function for profiling, which has visibility over
      // all variables defined up to now
      auto Scan = [&] (double t1Val, double t2Val) {
        if (t1Val < poi1->GetRangeMin() || t1Val > poi1->GetRangeMax()) return false;
        poi1->FixTo(t1Val);
        if (t2Val < poi2->GetRangeMin() || t2Val > poi2->GetRangeMax()) return false;
        poi2->FixTo(t2Val);

        // get bin index in the nll map 
        const int i1 = hpll->GetXaxis()->FindBin(t1Val);
        const int i2 = hpll->GetYaxis()->FindBin(t2Val);

        for (const auto& step : json["MinimizerSteps"].GetObject()) {
          fitter->SetMinuitVerbosity(step.value["verbosity"].GetInt());
          fitter->Minimize(step.value["method"].GetString(),
              step.value["resetMinuit"].GetBool(),
              step.value["maxCall"].GetDouble(),
              step.value["tollerance"].GetDouble());
        }

        // extract temporary best fit value
        const double tmpMinNLL = fitter->GetMinNLL();
        hpll->SetBinContent(i1, i2, tmpMinNLL); 
        // update absolute minimum if needed
        if (tmpMinNLL < absMinNLL) absMinNLL = tmpMinNLL;

        // check the status of minuit
        if (fitter->GetMinuitStatus()) {
          std::cerr << "Profile2D >> error: minuit returned failed status ["
            << fitter->GetMinuitStatus() << "]"
            << " while fitting with " << parName1.data()
            << " fixed to " << t1Val
            << " and " << parName2.data()
            << " fixed to " << t2Val << std::endl;
        } 
        // chek if the exit conditions are met
        if (tmpMinNLL - absMinNLL <= NLL)  return true;
        else return false;
      };

      // perform actual scan
      {
        // index used to define the points to scans
        // best fit value and error
        fitter->SyncFitParameters(true);
        
        for (int j1 = 1; j1 <= hpll->GetNbinsX(); j1++) {
          for (int j2 = 1; j2 <= hpll->GetNbinsY(); j2++) {
            const double t1Val = hpll->GetXaxis()->GetBinCenter(j1);
            const double t2Val = hpll->GetYaxis()->GetBinCenter(j2);
            Scan(t1Val, t2Val);
          }
        }

        // normilize profile to the absolute minimum found during while profiling
        const double min_pll = hpll->GetMinimum(); 
        for (int i = 0; i < hpll->GetNbinsX(); i++ ) {
          for (int j = 0; j < hpll->GetNbinsY(); j++ ) {
            hpll->SetBinContent(i, j, hpll->GetBinContent(i, j) - min_pll);
          }
        }
      }

      // reset parameter original status
      // restore best fit values of the parameters
      {
        int parIndex = 0;
        for ( auto it : *fitter->GetParameterMap()) {
          it.second->SetFitBestValue(fitBestValue.at(parIndex));
          it.second->SetFitBestValueErr(fitBestValueErr.at(parIndex));
          parIndex++;
        }
      }
      poi1->Release();
      poi2->Release();

      // Set titles (this must be done after filling the TGraph. Probably it's a
      // bug of ROOT
      hpll->SetName(Form("nll_%s_%s", parName1.c_str(), parName2.c_str())); 
      hpll->GetXaxis()->SetTitle(parName1.data()); 
      hpll->GetYaxis()->SetTitle(parName2.data()); 
      hpll->GetZaxis()->SetTitle("-#Delta LogLikelihood)");
      return hpll;
    }

    /*
     * Build Profile for each parameter of the fit
     */

    inline TCanvas* GetCanvasProfiles (const rapidjson::Document& json, MSMinimizer* fitter, 
        const double NLL, const int nPts) {

      // retrieve the canvas or initialize it
      delete gROOT->GetListOfCanvases()->FindObject("cPLL");
      TCanvas* cc = nullptr;

      if (json["MC"].HasMember("profile2D")) {
        const int nWindows = json["MC"]["profile2D"].Size();
        cc = new TCanvas("cPL2D", "profile log likelihoods 2D",
            400*ceil(nWindows/2.0), 400*ceil(nWindows/2.0));
        cc->DivideSquare( nWindows );

        int iwindow = 1;

        for (const auto& window : json["MC"]["profile2D"].GetArray()) {
          TH2D* tmp = Profile2D(json, fitter, 
              window["par1"].GetString(), 
              window["par2"].GetString(), 3.0, 
              window["nPoints1"].GetInt(), window["nPoints2"].GetInt());
          cc->cd( iwindow );
          tmp->Draw("colz");
          cc->Update();
          iwindow++;
        }
      }
      else if (json["MC"].HasMember("profile1D")) {
        const int nWindows = json["MC"]["profile1D"].Size();
        cc = new TCanvas("cPL1D", "profile log likelihoods 1D",
            400*ceil(nWindows/2.0), 400*ceil(nWindows/2.0));
        cc->DivideSquare( nWindows );

        int iwindow = 1;

        for (const auto& window : json["MC"]["profile1D"].GetArray()) {
          TGraph* tmp = Profile(json, fitter, 
              window["par"].GetString(), 
              3.0, 
              window["nPoints"].GetInt());
          cc->cd( iwindow );
          tmp->Draw("apl");
          cc->Update();
          iwindow++;
        }
      }
      return cc;
    }
  } // namespace mst

#endif // MST_MSHistFit_CXX
