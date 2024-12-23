/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : MSPDF.h
 * @created     : Wednesday Dec 18, 2024 14:56:50 CET
 */

#ifndef MSTHnPDF_H

#define MSTHnPDF_H

#include  "MSObject.h"
#include  "THn.h"

namespace mst {
  class MSTHnPDF : public MSObject {
    public:
      inline MSTHnPDF(const std::string& name = ""): MSObject(name) {}
      inline MSTHnPDF(const std::string& name, THn* thn): MSObject(name), fTHnPDF(thn) {}
      inline ~MSTHnPDF() { if (fTHnPDF) delete fTHnPDF; }
      inline THn* GetTHn() {return fTHnPDF;}
      inline void SetTHn(THn* thn) {fTHnPDF = thn;}
      inline void SetRespMatrix(THn* thn) {fRespMatrix = thn;}
      inline THn* GetRespMatrix() {return fRespMatrix;}
      inline bool ApplyOscillation() const {return fApplyOscillation;}
      inline void SetApplyOscillation(bool apply) {fApplyOscillation = apply;}

    private:
      THn* fTHnPDF = {nullptr};
      THn* fRespMatrix = {nullptr};
      bool fApplyOscillation = {false};
  };
}
#endif /* end of include guard MSPDF_H */

