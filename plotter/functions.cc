#include <iostream>
#include <TROOT.h>
#include <TMath.h>

float saturationFactorNLO(float density) {
  // this gives the energy in eV/ph
  float ret = 0;
  if (density<=0) // protection for the formula below
    ret = 0.85; // seems to provide some continuity
  else {
    float x = density/1.5;
    ret = (0.11 + 1.22*x)/(130283. * (1-TMath::Exp(-1*(TMath::Power(x,0.56757)/117038.))))/1.2;
  }
  // std::cout << "saturation factor = " << ret << std::endl;
  return ret;
}

float calibratedEnergy(float integral, float density) {
  // giving the energy in keV
  return saturationFactorNLO(density) * integral / 1000.; 
}
