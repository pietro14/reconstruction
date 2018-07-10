#include "TStyle.h"

void rootlogon() {

  TStyle *vecbosStyle = new TStyle("vecbosStyle","Style for P-TDR");

  ///////// pretty palette ///////////

  const Int_t NRGBs = 5;
  const Int_t NCont = 255;

  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  vecbosStyle->SetNumberContours(NCont);

  // For the canvas:
  vecbosStyle->SetCanvasBorderMode(0);
  vecbosStyle->SetCanvasColor(kWhite);
  vecbosStyle->SetCanvasDefH(600); //Height of canvas
  vecbosStyle->SetCanvasDefW(600); //Width of canvas
  vecbosStyle->SetCanvasDefX(0);   //POsition on screen
  vecbosStyle->SetCanvasDefY(0);

  // For the Pad:
  vecbosStyle->SetPadBorderMode(0);
  // vecbosStyle->SetPadBorderSize(Width_t size = 1);
  vecbosStyle->SetPadColor(kWhite);
  vecbosStyle->SetPadGridX(false);
  vecbosStyle->SetPadGridY(false);
  vecbosStyle->SetGridColor(0);
  vecbosStyle->SetGridStyle(3);
  vecbosStyle->SetGridWidth(1);

  // For the frame:
  vecbosStyle->SetFrameBorderMode(0);
  vecbosStyle->SetFrameBorderSize(1);
  vecbosStyle->SetFrameFillColor(0);
  vecbosStyle->SetFrameFillStyle(0);
  vecbosStyle->SetFrameLineColor(1);
  vecbosStyle->SetFrameLineStyle(1);
  vecbosStyle->SetFrameLineWidth(1);

  // set the paper & margin sizes
  vecbosStyle->SetPaperSize(20,26);
  vecbosStyle->SetPadTopMargin(0.05);
  vecbosStyle->SetPadRightMargin(0.05);
  vecbosStyle->SetPadBottomMargin(0.16);
  vecbosStyle->SetPadLeftMargin(0.12);

  // use large Times-Roman fonts
  vecbosStyle->SetTitleFont(132,"xyz");  // set the all 3 axes title font
  vecbosStyle->SetTitleFont(132," ");    // set the pad title font
  vecbosStyle->SetTitleSize(0.06,"xyz"); // set the 3 axes title size
  vecbosStyle->SetTitleSize(0.06," ");   // set the pad title size
  vecbosStyle->SetLabelFont(132,"xyz");
  vecbosStyle->SetLabelSize(0.05,"xyz");
  vecbosStyle->SetTextFont(132);
  vecbosStyle->SetTextSize(0.08);
  vecbosStyle->SetStatFont(132);

  // use bold lines and markers
  vecbosStyle->SetMarkerStyle(8);
  vecbosStyle->SetHistLineWidth(2);
  vecbosStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  //..Get rid of X error bars
  vecbosStyle->SetErrorX(0.001);

  // do not display any of the standard histogram decorations
  vecbosStyle->SetOptTitle(0);
  vecbosStyle->SetOptStat(1111111);
  vecbosStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  vecbosStyle->SetPadTickX(1);
  vecbosStyle->SetPadTickY(1);

  vecbosStyle->cd();

}
