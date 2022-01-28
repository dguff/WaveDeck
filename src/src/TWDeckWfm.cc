/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : TWDeckWfm
 * @created     : mercoledì gen 26, 2022 09:05:05 CET
 */

#include "TWDeckWfm.h"

ClassImp(TWDeckWfm)

TWDeckWfm::TWDeckWfm() : fSize(0), fWfm(0), fWfm_re(0), fWfm_im(0) {}

TWDeckWfm::TWDeckWfm(int N) : fSize(0), fWfm(0), fWfm_re(0), fWfm_im(0) 
{
  SetSize(N);
}

TWDeckWfm::TWDeckWfm(int N, double* wfm) : fSize(0), fWfm(0), fWfm_re(0), fWfm_im(0) 
{
  LoadWave(N, wfm);
}

TWDeckWfm::TWDeckWfm(int N, double* re, double* im) 
  : fSize(0), fWfm(0), fWfm_re(0), fWfm_im(0) 
{
  LoadSpectrum(N, re, im);
}

TWDeckWfm::TWDeckWfm(const TWDeckWfm& wfm) : fSize(0), fWfm(0), fWfm_re(0), fWfm_im(0) 
{
  fName = wfm.fName;
  fTitle = wfm.fTitle;
  fSize = wfm.fSize;
  SetSize(fSize);

  for (int i=0; i<fSize; i++) {
    fWfm[i] = wfm.fWfm[i];
    fWfm_re[i] = wfm.fWfm_re[i];
    fWfm_im[i] = wfm.fWfm_im[i];
  }
}

TWDeckWfm::~TWDeckWfm() {}

void TWDeckWfm::ClearWave() {
  for (auto &v : fWfm   ) v = 0.;
  for (auto &v : fWfm_re) v = 0.;
  for (auto &v : fWfm_im) v = 0.;
}

void TWDeckWfm::LoadWave(double* w) {
  for (int i=0; i<fSize; i++) {
    fWfm[i] = w[i];
  }
}

void TWDeckWfm::LoadWave(int n, double* w) {
  SetSize(n);
  LoadWave(w);
}

void TWDeckWfm::LoadSpectrum(double* re, double* im) {
  for (int i=0; i<fSize; i++) {
    fWfm_re[i] = re[i];
    fWfm_im[i] = im[i];
  }
}

void TWDeckWfm::LoadSpectrum(int N, double* re, double* im) {
  SetSize(N);
  LoadSpectrum(re, im); 
}

void TWDeckWfm::SetSize(int n) {
  fSize = n;

  fWfm   .resize(fSize, 0.);
  fWfm_re.resize(fSize, 0.);
  fWfm_im.resize(fSize, 0.);
}

std::vector<TComplex> TWDeckWfm::GetPointsComplex() {
  std::vector<TComplex> vWfmC(fSize);
  for (int i=0; i<fSize; i++) {
    vWfmC[i] = TComplex(fWfm_re[i], fWfm_im[i]);
  }

  return vWfmC;
}

double TWDeckWfm::GetSpectralDensity(int i) {
  TComplex c(fWfm_re.at(i), fWfm_im.at(i));
  return c.Rho2();
}

std::vector<double> TWDeckWfm::GetSpectralDensityPoints() {
  std::vector<double> rho2(fSize*0.5+1, 0.);
  for (int i=0; i<fSize*0.5+1; i++)
    rho2.at(i) = GetSpectralDensity(i);
  return rho2;
}


