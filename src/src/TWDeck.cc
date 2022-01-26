/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : TWDeck
 * @created     : mercoledì gen 26, 2022 12:36:38 CET
 */

#include "TWDeck.h"
#include <utility>

ClassImp(TWDeck)

TWDeck::TWDeck() : 
  fOrigin(nullptr), fResponse(nullptr), 
  fNoiseDensity(nullptr), fH2NoiseDensity(nullptr),
  fFFT_R2C(nullptr), fFFT_C2R(nullptr)
{}

TWDeck::TWDeck(int n) : fOrigin(nullptr), fResponse(nullptr), 
  fNoiseDensity(nullptr), fH2NoiseDensity(nullptr), 
  fFFT_R2C(nullptr), fFFT_C2R(nullptr)
{
  fSize = n;
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

void TWDeck::RegisterFilter(const char* filter_name, TWDeckWfmFilter* filter) {
  TWDeckWfmFilter* pdec_filter = new TWDeckWfmFilter(*filter);
  ResizeFilter(pdec_filter);
  fFilters.insert(std::make_pair(filter_name, pdec_filter));
}

void TWDeck::ResizeFilters() {
  for (auto &f : fFilters) {
    if (!f.first.Contains("Wiener")) {
      ResizeFilter(f.second);
    }
  }   
}

void TWDeck::ResizeFilter(TWDeckWfmFilter* filter) {
  size_t size_filter = filter->GetSize();
  if (filter->GetWfm().size() != size_filter + fSize)
    filter->GetWfm().resize(fSize+size_filter, 0.);
  
  return;
}

void TWDeck::ApplyFilter(TWDeckWfm* wfm, TString filter_name) {
  if ( !(fFilters.count(filter_name) ) ) {
    printf("TWDeck::ApplyFilter ERROR: no filter named '%s' is registered. quit.\n", 
        filter_name.Data());
    return;
  }

  TWDeckWfmFilter* filter = fFilters.find(filter_name)->second;

  printf("TWDeck::ApplyFilter: copy input waveform...\n");
  TWDeckWfm wfm_tmp(*wfm);
  if (fSize !=wfm_tmp.GetSize()) {
    fSize = wfm_tmp.GetSize();
    ResizeFilter(filter);
  }

  int size_tmp = filter->GetSize() + wfm_tmp.GetSize();
  fFFTSize = size_tmp;

  printf("TWDeck::ApplyFilter: resize waveform to %i...\n", size_tmp);
  wfm_tmp.SetSize(size_tmp);
  
  printf("TWDeck::ApplyFilter: instance FFTs...\n");
  BuildFFT();

  printf("TWDeck::ApplyFilter: Transform wave...\n");
  fFFT_R2C->SetPoints(&wfm_tmp.GetWfm().at(0));
  fFFT_R2C->Transform();
  fFFT_R2C->GetPointsComplex(&wfm_tmp.GetWfmRe().at(0), &wfm_tmp.GetWfmIm().at(0));
  
  printf("TWDeck::ApplyFilter: Transform filter...\n");
  fFFT_R2C->SetPoints(&filter->GetWfm().at(0));
  fFFT_R2C->Transform();
  fFFT_R2C->GetPointsComplex(&filter->GetWfmRe().at(0), &filter->GetWfmIm().at(0));

  printf("TWDeck::ApplyFilter: Perform convolution...\n");
  for (int i = 0; i < size_tmp; i++) {
    TComplex F = TComplex(filter->GetWfmRe()[i], filter->GetWfmIm()[i]);
    TComplex W = TComplex(wfm_tmp.GetWfmRe()[i], wfm_tmp.GetWfmIm()[i]);

    TComplex C = W*F*TMath::Exp(TComplex(0, -filter->GetBandwidth() / fFFTSize));
    wfm_tmp.GetWfmRe()[i] = C.Re();
    wfm_tmp.GetWfmIm()[i] = C.Im();
  }

  fFFT_C2R->SetPointsComplex(&wfm_tmp.GetWfmRe().at(0), &wfm_tmp.GetWfmIm().at(0));
  fFFT_C2R->Transform();
  double* xw_tmp;
  xw_tmp = fFFT_C2R->GetPointsReal();

  for (int i=0; i<wfm->GetSize(); i++) {
    wfm->GetWfm()[i] = xw_tmp[i] / size_tmp;
  }

  printf("DONE. Ready to return.\n\n");

  return;
}


void TWDeck::BuildFFT() {
  if (!fFFT_R2C) 
    fFFT_R2C = TVirtualFFT::FFT(1, &fFFTSize, "M R2C K");

  if (!fFFT_C2R) 
    fFFT_C2R = TVirtualFFT::FFT(1, &fFFTSize, "M C2R K");
}
