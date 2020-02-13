#include "RooDataHist.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "TH1.h"
#include "TSystem.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include <boost/python.hpp>

#include <vector>
#include <string>
#ifdef __CINT__
#pragma link C++ class std::vector<std::string>+;
#endif


using namespace RooFit;
using namespace std;

namespace py = boost::python;

class histFitter {
public:
  histFitter( TH1 *hData, std::string histname  );
  ~histFitter(void) {if( _work != 0 ) delete _work; }
  void setBkgPDF(TH1 *hBkg);
  void setWorkspace(std::vector<std::string>);
  void setOutputFile(TFile *fOut ) {_fOut = fOut;}
  void fits(std::string title = "");
  float efficiency(double xthres,std::string pdf="sigPdf");
  void useMinos(bool minos = true) {_useMinos = minos;}

  void setFitRange(double xMin,double xMax) { _xFitMin = xMin; _xFitMax = xMax; }
private:
  RooWorkspace *_work;
  std::string _histname_base;
  TFile *_fOut;
  double _nTot;
  bool _useMinos;
  double _xFitMin,_xFitMax;
};


histFitter::histFitter(TH1 *hData, std::string histname  ) : _useMinos(false) {
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  _histname_base = histname;
  
  _nTot = hData->Integral();

  _work = new RooWorkspace("w") ;
  _work->factory("x[5,30]");
  
  RooDataHist rooData("hData","hData",*_work->var("x"),hData);
  _work->import(rooData) ;
  _xFitMin = 5;
  _xFitMax = 30;
  
}

void histFitter::setBkgPDF(TH1 *hBkg) {
  RooDataHist rooData("hBkg","hBkg",*_work->var("x"),hBkg);
  _work->import(rooData) ;
}

void histFitter::setWorkspace(std::vector<std::string> workspace) {
  for( unsigned icom = 0 ; icom < workspace.size(); ++icom ) {
    _work->factory(workspace[icom].c_str());
  }

  _work->factory("HistPdf::bkgPdf(x,hBkg)");
  _work->factory(TString::Format("nSig[%f,0.5,%f]",_nTot*0.1,_nTot*1.5));
  _work->factory(TString::Format("nBkg[%f,0.5,%f]",_nTot*0.9,_nTot*1.5));
  _work->factory("SUM::pdfTot(nSig*sigPdf,nBkg*bkgPdf)");
  _work->Print();			         
}

void histFitter::fits(string title) {

  cout << " title : " << title << endl;

  
  RooAbsPdf *pdfTot = _work->pdf("pdfTot");

  // FC: seems to be better to change the actual range than using a fitRange in the fit itself (???)
  /// FC: I don't know why but the integral is done over the full range in the fit not on the reduced range
  _work->var("x")->setRange(_xFitMin,_xFitMax);
  _work->var("x")->setRange("fitRange",_xFitMin,_xFitMax);
  RooFitResult* res = pdfTot->fitTo(*_work->data("hData"),Minos(_useMinos),SumW2Error(kTRUE),Save(),Range("fitRange"));

  RooPlot *pTot = _work->var("x")->frame(5,30);
  pTot->SetTitle("selected clusters");
  
  _work->data("hData") ->plotOn( pTot );
  _work->pdf("pdfTot") ->plotOn( pTot, LineColor(kRed) );
  _work->pdf("pdfTot") ->plotOn( pTot, Components("bkgPdf"),LineColor(kBlue),LineStyle(kDashed));
  _work->pdf("pdfTot") ->plotOn( pTot, Components("sigPdf"),LineColor(kGreen+2),LineStyle(kDashed));
  _work->data("hData") ->plotOn( pTot );
  
  TCanvas c("c","c",1100,450);
  pTot->Draw();

  _fOut->cd();
  c.Write(TString::Format("%s_Canv",_histname_base.c_str()),TObject::kOverwrite);
  res->Write(TString::Format("%s_res",_histname_base.c_str()),TObject::kOverwrite);
  
}

float histFitter::efficiency(double xthres,std::string pdf) {

  RooRealVar *x = _work->var("x");
  RooAbsPdf *sigPdf = _work->pdf(pdf.c_str());

  x->setRange("my_range", xthres, _xFitMax);
  RooAbsReal* integral = sigPdf->createIntegral(*x, RooFit::NormSet(RooArgSet(*x)), RooFit::Range("my_range"));

  float effi = integral->getVal();
  return effi;
}
