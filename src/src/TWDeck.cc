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
  fOrigin(nullptr), fResponse(nullptr), 
  fNoiseDensity(nullptr), fH2NoiseDensity(nullptr),
  fFFT_R2C(nullptr), fFFT_C2R(nullptr), 
  fSize(100), fFFTSize(200)
{
  BuildFFT();
}

TWDeck::TWDeck(int n) : fOrigin(nullptr), fResponse(nullptr), 
  fNoiseDensity(nullptr), fH2NoiseDensity(nullptr), 
  fFFT_R2C(nullptr), fFFT_C2R(nullptr)
{
  fSize = n;
  fFFTSize = 2*n;
  BuildFFT();
}

TWDeck::~TWDeck()
{
  if (fOrigin) {delete fOrigin;}
  if (fResponse) {delete fResponse;}
  if (fNoiseDensity  ) {delete fNoiseDensity;}
  if (fH2NoiseDensity) {delete fH2NoiseDensity;}
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

void TWDeck::FFTC2R(TWDeckWfm* wfm) {
  fFFT_C2R->SetPointsComplex(&wfm->GetWfmRe().at(0), &wfm->GetWfmIm().at(0));
  fFFT_C2R->Transform();
  double* v = fFFT_C2R->GetPointsReal();
  for (int j=0; j<wfm->GetSize(); j++)
    wfm->GetWfm().at(j) = v[j];

  return;
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

  for (int i = 0; i < size_tmp; i++) {
    TComplex F = TComplex(filter->GetWfmRe()[i], filter->GetWfmIm()[i]);
    TComplex W = TComplex(wfm_tmp.GetWfmRe()[i], wfm_tmp.GetWfmIm()[i]);
  
    TComplex C = W*F;
    // apply shift to account for filter possible shift
    TComplex ph = TComplex(0., +TMath::TwoPi()*i*filter->GetShift()/fFFTSize);
    C = C * TComplex::Exp(ph);
    wfm_tmp.GetWfmRe()[i] = C.Re();
    wfm_tmp.GetWfmIm()[i] = C.Im();
  }

  FFTC2R(&wfm_tmp);
  
  for (int i=0; i<wfm->GetSize(); i++) {
    wfm->GetWfm()[i] = wfm_tmp.GetWfm().at(i) / size_tmp;
  }

  return;
}


void TWDeck::BuildFFT() {
  if (!fFFT_R2C) 
    fFFT_R2C = TVirtualFFT::FFT(1, &fFFTSize, "M R2C K");

  if (!fFFT_C2R) 
    fFFT_C2R = TVirtualFFT::FFT(1, &fFFTSize, "M C2R K");
}
