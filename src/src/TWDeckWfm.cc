/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : TWDeckWfm
 * @created     : mercoledì gen 26, 2022 09:05:05 CET
 */

#include "TWDeckWfm.h"

ClassImp(TWDeckWfm)

TWDeckWfm::TWDeckWfm() : fSize(0) {}

TWDeckWfm::TWDeckWfm(int N) {
  SetSize(N);
}

TWDeckWfm::TWDeckWfm(int N, double* wfm) {
  LoadWave(N, wfm);
}

TWDeckWfm::TWDeckWfm(const TWDeckWfm& wfm) {
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

TWDeckWfm::~TWDeckWfm() {
  printf("Calling ~TWDeckWfm\n");
}

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

void TWDeckWfm::SetSize(int n) {
  fSize = n;

  fWfm   .resize(fSize, 0.);
  fWfm_re.resize(fSize, 0.);
  fWfm_im.resize(fSize, 0.);
}

std::vector<TComplex> TWDeckWfm::GetWfmC() {
  std::vector<TComplex> vWfmC(fSize);
  for (int i=0; i<fSize; i++) {
    vWfmC[i] = TComplex(fWfm_re[i], fWfm_im[i]);
  }

  return vWfmC;
}
