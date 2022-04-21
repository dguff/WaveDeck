/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : TWDeck.h
 * @created     : mercoledì gen 26, 2022 12:09:47 CET
 */

#ifndef TWDECK_H
             
#define TWDECK_H

#include <vector>
#include <map>

#include "TVirtualFFT.h"
#include "TGraphErrors.h"
#include "TH2D.h"
#include "TWDeckWfm.h"
#include "TWDeckWfmFilter.h"
#include "TWDeckWfmModel.h"

/**
 * @brief Waveform convolution/deconvolution interface
 *
 * The TWDeck class is the cornerstone of the WaveDeck package, 
 * providing a (hopefully convenient) interface to the FFT utilities 
 * used by the ROOT package. 
 * This class operates on waveforms described by TWDeckWfm and derived classes
 * to perform Fourier transform to/from the real and complex domain. 
 * 
 * Furthermore, one can use this class to add samples to a given model
 * filling both the waveform and spectrum models.
 * 
 * A single instance of TWDeck can register several filters that can be applied 
 * to a waveform.
 */
class TWDeck : public TObject {
  public: 
    //! Empty constructor
    TWDeck();
    //! Create an instance of TWDeck with size `n`
    TWDeck(int n);
    //! Destructor
    ~TWDeck();

    //! Add the waveform `wfm` to the `model`
    void Add2Model  (TWDeckWfm* wfm, TWDeckWfmModel* model);
    //! Apply the filter registered under the given `filter_name` to the waveform `wfm`
    void ApplyFilter(TWDeckWfm* wfm, TString filter_name, bool padding);
    //! Apply the filter `filter` to the waveform `wfm`
    void ApplyFilter(TWDeckWfm* wfm, TWDeckWfmFilter* filter, bool padding);

    //! Build the FFT plans to be used for waveform Fourier Transform
    void BuildFFT(int size = 1000);
    //! Perform a FFT from Real to Complex on `wfm`
    void FFTR2C(TWDeckWfm* wfm);
    //! Perform a FFT from Real to Complex on `wfm`
    void FFTR2C(TWDeckWfm* wfm, int size);
    //! Perform a FFT from Complex to Real on `wfm`
    void FFTC2R(TWDeckWfm* wfm);
    //! Perform a FFT from Complex to Real on `wfm`
    void FFTC2R(TWDeckWfm* wfm, int size);

    //! Produce sythetic waveform from the model with the same power spectral density
    TWDeckWfm* Produce(TWDeckWfmModel* model);

    //! Register a filter in the filter list
    void RegisterFilter(const char* filter_name, TWDeckWfmFilter* filter, bool padding);
    //! Get filter from the register
    TWDeckWfmFilter* GetFilter(const char* filter_name);
    //! Get wavedeck size
    int GetSize() {return fSize;}
    //! Get wavedeck FFT size
    int GetFFTSize() {return fFFTSize;}
    //! Set wavedeck size
    void SetSize(int n);
    //! Resize the filter
    void ResizeFilter(TWDeckWfmFilter* filter, bool padding, wdeck::EWfmDomain kDomain);


  protected:
    std::map<TString, TWDeckWfmFilter*> fFilters; //!< Filters register
    TVirtualFFT* fFFT_R2C;      //!< FFT_R2C engine
    TVirtualFFT* fFFT_C2R;      //!< FFT_C2R engine
    int       fSize;            //!< wavedeck size
    int       fFFTSize;         //!< wavedeck FFT size

  public:
    ClassDef(TWDeck, 1)
};


#endif /* end of include guard TWDECK_H */

