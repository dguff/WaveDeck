/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : TWDeck
 * @created     : mercoledì gen 26, 2022 12:09:47 CET
 */

#ifndef TWDECK_H
             
#define TWDECK_H

#include <vector>
#include <map>

#include "TVirtualFFT.h"
#include "TGraphErrors.h"
#include "TH2D.h"
#include "TWDeckWfmFilter.h"

class TWDeck : public TObject {
  public: 
    TWDeck();
    TWDeck(int n);
    ~TWDeck();

    void ApplyFilter(TWDeckWfm* wfm, TString filter_name);
    void RegisterFilter(const char* filter_name, TWDeckWfmFilter* filter, wdeck::EWfmDomain kDomain = wdeck::kReal);
    void ResizeFilters();


  private:
    TWDeckWfm* fOrigin;
    TWDeckWfm* fResponse;
    std::map<TString, TWDeckWfmFilter*> fFilters;
    TGraph*   fNoiseDensity;
    TH2D*     fH2NoiseDensity;
    TVirtualFFT* fFFT_R2C;
    TVirtualFFT* fFFT_C2R;
    int       fSize;
    int       fFFTSize;

    void      BuildFFT();
    void      ResizeFilter(TWDeckWfmFilter* filter);
    void      FFTR2C(TWDeckWfm* wfm);
    void      FFTC2R(TWDeckWfm* wfm);
  public:
    ClassDef(TWDeck, 1)
};


#endif /* end of include guard TWDECK_H */

