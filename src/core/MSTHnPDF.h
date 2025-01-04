/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : MSPDF.h
 * @created     : Wednesday Dec 18, 2024 14:56:50 CET
 */

#ifndef MSTHnPDF_H

#define MSTHnPDF_H

#include <map>
#include <vector>

#include  "MSObject.h"
#include  "THn.h"

namespace mst {
  enum EPDFType {
    kUndefined = -1,
    kComponent = 0,
    kNeutrino = 1, 
    kResponseMatrix = 2
  };

  inline EPDFType get_pdf_type_code(const std::string& pdfType) {
    if (pdfType == "component") return kComponent;
    if (pdfType == "neutrino") return kNeutrino;
    if (pdfType == "responseMatrix") return kResponseMatrix;
    return kUndefined;
  }

  class MSTHnPDF : public MSObject {
    public:
      inline MSTHnPDF(const std::string& name = ""): MSObject(name) {}
      inline MSTHnPDF(const std::string& name, THn* thn): MSObject(name), fTHnPDF(thn) {}
      inline ~MSTHnPDF() { if (fTHnPDF) delete fTHnPDF; }
      
      inline THn* GetTHn() {return fTHnPDF;}
      inline void SetTHn(THn* thn) {fTHnPDF = thn;}
      
      inline void SetPDFType(const std::string& pdfType) {fPDFType = get_pdf_type_code(pdfType);}
      inline void SetPDFType(EPDFType pdfType) {fPDFType = pdfType;}
      
      inline EPDFType GetPDFType() const {return fPDFType;}

    protected:
      EPDFType fPDFType = {kUndefined};
      THn* fTHnPDF = {nullptr};
  };

  class MSTHnPDFComponent : public MSTHnPDF {
    public: 
      inline MSTHnPDFComponent(const std::string& name = ""): MSTHnPDF(name) {}
      inline MSTHnPDFComponent(const std::string& name, THn* thn): MSTHnPDF(name, thn) {}
      inline ~MSTHnPDFComponent() { if (fTHnPDF) delete fTHnPDF; }

      inline THn* GetRespMatrix() {return fRespMatrix;}
      inline void SetRespMatrix(THn* respMatrix) {fRespMatrix = respMatrix;}

    protected:
      THn* fRespMatrix = {nullptr};
  };

  class MSTHnPDFNeutrino : public MSTHnPDF {
    public:
      struct NuIntChannel_t {
        std::string fName = {};
        int fPDG = {0};
        double fNormalization = {0}; 
        std::string fResponeMatrix = {};
        std::vector<double> fCrossSection = {};
        int fColor = {0};

        NuIntChannel_t() {}
        
        NuIntChannel_t(const std::string& name, const int& pdg, 
            const std::string& respMatrix, int color = 0)
          : fName(name), fPDG(pdg), fResponeMatrix(respMatrix), fColor(color) {}

        inline bool operator==(const NuIntChannel_t& other) const {
          return 
            fName == other.fName && 
            fPDG == other.fPDG &&
            fResponeMatrix == other.fResponeMatrix;
        }
      };

      inline MSTHnPDFNeutrino(const std::string& name = ""): MSTHnPDF(name) {}
      inline MSTHnPDFNeutrino(const std::string& name, THn* thn): MSTHnPDF(name, thn) {}
      inline ~MSTHnPDFNeutrino() { if (fTHnPDF) delete fTHnPDF; }
      
      inline THn* GetTHn() {return fTHnPDF;}
      inline void SetTHn(THn* thn) {fTHnPDF = thn;}

      inline NuIntChannel_t& AddChannel(const NuIntChannel_t& channel) {
        fChannels.push_back(channel); return fChannels.back();}
      inline NuIntChannel_t& AddChannel(const std::string& name, 
          const int& pdg, const std::string& respMatrix, int color = 0) {
        fChannels.push_back(NuIntChannel_t(name, pdg, respMatrix, color)); 
        return fChannels.back();
      }
      inline NuIntChannel_t& GetChannel(const std::string& name) {
        for (auto& channel : fChannels) {
          if (channel.fName == name) return channel;
        }
        throw std::runtime_error("Channel not found");
      }
      
      inline std::vector<NuIntChannel_t>& GetChannels() {return fChannels;}
      inline const std::vector<NuIntChannel_t>& GetChannels() const {return fChannels;}
      
      inline bool ApplyOscillation() const {return fApplyOscillation;}
      inline void SetApplyOscillation(bool apply) {fApplyOscillation = apply;}

    protected:
      std::vector<NuIntChannel_t> fChannels;
      bool fApplyOscillation = {false};
  };
}
#endif /* end of include guard MSPDF_H */

