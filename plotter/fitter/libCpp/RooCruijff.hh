#ifndef ROO_CRUIJFF
#define ROO_CRUIJFF

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooTrace.h"

class RooCruijff : public RooAbsPdf {
public:
  RooCruijff();
  RooCruijff(const char *name, const char *title, 
             RooAbsReal& _x,
             RooAbsReal& _m0, 
             RooAbsReal& _sigmaL, RooAbsReal& _sigmaR,
             RooAbsReal& _alphaL, RooAbsReal& _alphaR);
  
  RooCruijff(const RooCruijff& other, const char* name = 0);

  virtual TObject* clone(const char* newname) const { 
    return new RooCruijff(*this,newname); }

  inline virtual ~RooCruijff() { }

protected:
  RooRealProxy x;
  RooRealProxy m0;
  RooRealProxy sigmaL;
  RooRealProxy sigmaR;
  RooRealProxy alphaL;
  RooRealProxy alphaR;

  Double_t evaluate() const;

private:
  ClassDef(RooCruijff,0)
};

#endif
