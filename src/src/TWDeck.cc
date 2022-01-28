/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : TWDeck
 * @created     : mercoledì gen 26, 2022 12:36:38 CET
 */

#include "TWDeck.h"
#include "TMath.h"
#include <utility>

ClassImp(TWDeck)

TWDeck::TWDeck() : 
  fFFT_R2C(nullptr), fFFT_C2R(nullptr), 
  fSize(1000), fFFTSize(2000)
{
  BuildFFT();
}

TWDeck::TWDeck(int n) : 
  fFFT_R2C(nullptr), fFFT_C2R(nullptr)
{
  SetSize(n);
}

TWDeck::~TWDeck()
{
  if (fFFT_R2C) {delete fFFT_R2C;}
  if (fFFT_C2R) {delete fFFT_C2R;}

  for (auto &f : fFilters) {delete f.second;}
  fFilters.clear();
}

void TWDeck::RegisterFilter(const char* filter_name, TWDeckWfmFilter* filter, wdeck::EWfmDomain kDomain) {
  TWDeckWfmFilter* pdec_filter = new TWDeckWfmFilter(*filter);
  ResizeFilter(pdec_filter);
  int size_tmp = fFFTSize;
  fFFTSize = pdec_filter->GetSize() + fSize;
  if (kDomain == wdeck::kReal) {
    FFTR2C(pdec_filter);
  } else if (kDomain == wdeck::kComplex) {
    FFTC2R(pdec_filter);
  }
  fFilters.insert(std::make_pair(filter_name, pdec_filter));

  fFFTSize = size_tmp;
}


void TWDeck::FFTR2C(TWDeckWfm* wfm) {
  fFFT_R2C->SetPoints(&wfm->GetWfm().at(0));
  fFFT_R2C->Transform();
  double xre[10000] = {0};
  double xim[10000] = {0};
  fFFT_R2C->GetPointsComplex(xre, xim);
  for (int i=0; i<fFFTSize; i++) {
    wfm->GetWfmRe()[i] = xre[i];
    wfm->GetWfmIm()[i] = xim[i];
  }

  return;
}

void TWDeck::FFTR2C(TWDeckWfm* wfm, int size) {
  int size_fft_tmp = fFFTSize;
  BuildFFT(size);
  fFFT_R2C->SetPoints(&wfm->GetWfm().at(0));
  fFFT_R2C->Transform();
  double xre[10000] = {0};
  double xim[10000] = {0};
  fFFT_R2C->GetPointsComplex(xre, xim);
  for (int i=0; i<fFFTSize; i++) {
    wfm->GetWfmRe()[i] = xre[i];
    wfm->GetWfmIm()[i] = xim[i];
  }
  // restore original FFT size
  BuildFFT(fFFTSize);
  return;
}


void TWDeck::FFTC2R(TWDeckWfm* wfm) {
  fFFT_C2R->SetPointsComplex(&wfm->GetWfmRe().at(0), &wfm->GetWfmIm().at(0));
  fFFT_C2R->Transform();
  double* v = fFFT_C2R->GetPointsReal();
  double scale_ = 1./fFFTSize;
  for (int j=0; j<wfm->GetSize(); j++)
    wfm->GetWfm().at(j) = v[j] * scale_;

  return;
}

void TWDeck::FFTC2R(TWDeckWfm* wfm, int size) {
  int size_tmp = fFFTSize;
  BuildFFT(size);
  fFFT_C2R->SetPointsComplex(&wfm->GetWfmRe().at(0), &wfm->GetWfmIm().at(0));
  fFFT_C2R->Transform();
  double* v = fFFT_C2R->GetPointsReal();
  double scale_ = 1./fFFTSize;
  for (int j=0; j<wfm->GetSize(); j++)
    wfm->GetWfm().at(j) = v[j] * scale_;
  // restore original fFFTSize
  BuildFFT(size_tmp);
  return;
}

void TWDeck::Add2Model(TWDeckWfm* wfm, TWDeckWfmModel* model) {
  if (fFFTSize != wfm->GetSize()) BuildFFT(wfm->GetSize());
  FFTR2C(wfm);
  
  model->AddWaveform(&wfm->GetWfm()[0]);
  model->AddSpectrum(&wfm->GetWfmRe()[0], &wfm->GetWfmIm()[0]);
}

void TWDeck::ResizeFilters() {
  for (auto &f : fFilters) {
    if (!f.first.Contains("Wiener")) {
      ResizeFilter(f.second);
    }
  }   
}

void TWDeck::ResizeFilter(TWDeckWfmFilter* filter) {
  int size_filter = filter->GetSize();
  if (size_filter != size_filter + fSize)
  {
    filter->GetWfm  ().resize(fSize+size_filter, 0.);
    filter->GetWfmRe().resize(fSize+size_filter, 0.);
    filter->GetWfmIm().resize(fSize+size_filter, 0.);
  }
  
  return;
}

void TWDeck::ApplyFilter(TWDeckWfm* wfm, TString filter_name) {
  if ( !(fFilters.count(filter_name) ) ) {
    printf("TWDeck::ApplyFilter ERROR: no filter named '%s' is registered. quit.\n", 
        filter_name.Data());
    return;
  }

  TWDeckWfmFilter* filter = fFilters.find(filter_name)->second;

  TWDeckWfm wfm_tmp(*wfm);
  if (fSize !=wfm_tmp.GetSize()) {
    fSize = wfm_tmp.GetSize();
    ResizeFilter(filter);
  }

  int size_tmp = filter->GetSize() + wfm_tmp.GetSize();
  fFFTSize = size_tmp;

  wfm_tmp.SetSize(size_tmp);
  
  FFTR2C(&wfm_tmp);

  double fltr_shift = filter->GetShift();
  bool apply_shift = !(fltr_shift == 0);
  for (int i = 0; i < size_tmp; i++) {
    TComplex F = TComplex(filter->GetWfmRe()[i], filter->GetWfmIm()[i]);
    TComplex W = TComplex(wfm_tmp.GetWfmRe()[i], wfm_tmp.GetWfmIm()[i]);
  
    TComplex C = W*F;
    if (apply_shift) {
      // apply shift to account for filter possible shift
      TComplex ph = TComplex(0., +TMath::TwoPi()*i*filter->GetShift()/fFFTSize);
      C = C * TComplex::Exp(ph);
    }
    wfm_tmp.GetWfmRe()[i] = C.Re();
    wfm_tmp.GetWfmIm()[i] = C.Im();
  }

  FFTC2R(&wfm_tmp);
  
  for (int i=0; i<wfm->GetSize(); i++) {
    wfm->GetWfm()[i] = wfm_tmp.GetWfm().at(i) / size_tmp;
  }

  return;
}

void TWDeck::ApplyFilter(TWDeckWfm* wfm, TWDeckWfmFilter* filter) {

  TWDeckWfm wfm_tmp(*wfm);
  if (fSize !=wfm_tmp.GetSize()) {
    fSize = wfm_tmp.GetSize();
    ResizeFilter(filter);
  }

  int size_tmp = filter->GetSize() + wfm_tmp.GetSize();
  BuildFFT(size_tmp);

  wfm_tmp.SetSize(size_tmp);
  
  FFTR2C(&wfm_tmp);
  
  double fltr_shift = filter->GetShift();
  bool apply_shift = !(fltr_shift == 0);

  for (int i = 0; i < size_tmp; i++) {
    TComplex F = TComplex(filter->GetWfmRe()[i], filter->GetWfmIm()[i]);
    TComplex W = TComplex(wfm_tmp.GetWfmRe()[i], wfm_tmp.GetWfmIm()[i]);
  
    TComplex C = W*F;

    if (apply_shift) {
      // apply shift to account for filter possible shift
      TComplex ph = TComplex(0., +TMath::TwoPi()*i*filter->GetShift()/fFFTSize);
      C = C * TComplex::Exp(ph);
    }

    wfm_tmp.GetWfmRe()[i] = C.Re();
    wfm_tmp.GetWfmIm()[i] = C.Im();
  }

  FFTC2R(&wfm_tmp);
  
  for (int i=0; i<wfm->GetSize(); i++) {
    wfm->GetWfm()[i] = wfm_tmp.GetWfm().at(i);
  }

  return;
}

void TWDeck::SetSize(int n) {
  fSize    = n;
  fFFTSize = n;

  BuildFFT(fFFTSize);
}

void TWDeck::BuildFFT(int size) {
  fFFTSize = size;
  if (fFFT_R2C) {delete fFFT_R2C; fFFT_R2C = nullptr;}
  if (fFFT_C2R) {delete fFFT_C2R; fFFT_C2R = nullptr;}

  fFFT_R2C = TVirtualFFT::FFT(1, &fFFTSize, "M R2C K");
  fFFT_C2R = TVirtualFFT::FFT(1, &fFFTSize, "M C2R K");
}
