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
#include "TWDeckWfmModel.h"

class TWDeck : public TObject {
  public: 
    TWDeck();
    TWDeck(int n);
    ~TWDeck();

    void Add2Model  (TWDeckWfm* wfm, TWDeckWfmModel* model);
    void ApplyFilter(TWDeckWfm* wfm, TString filter_name);
    void ApplyFilter(TWDeckWfm* wfm, TWDeckWfmFilter* filter);

    void BuildFFT(int size = 1000);
    void FFTR2C(TWDeckWfm* wfm);
    void FFTR2C(TWDeckWfm* wfm, int size);
    void FFTC2R(TWDeckWfm* wfm);
    void FFTC2R(TWDeckWfm* wfm, int size);

    void RegisterFilter(const char* filter_name, TWDeckWfmFilter* filter, wdeck::EWfmDomain kDomain = wdeck::kReal);
    void ResizeFilters();
    void SetSize(int n);


  private:
    std::map<TString, TWDeckWfmFilter*> fFilters;
    TVirtualFFT* fFFT_R2C;
    TVirtualFFT* fFFT_C2R;
    int       fSize;
    int       fFFTSize;

    void      ResizeFilter(TWDeckWfmFilter* filter);
  public:
    ClassDef(TWDeck, 1)
};


#endif /* end of include guard TWDECK_H */

