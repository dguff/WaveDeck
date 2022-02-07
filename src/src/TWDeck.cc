/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : TWDeck
 * @created     : mercoledì gen 26, 2022 12:36:38 CET
 */

#include "TWDeck.h"
#include "TMath.h"
#include "TWDeckWfm.h"
#include "TRandom3.h"
#include <utility>

ClassImp(TWDeck)

/**
 * @details Initialize the WaveDeck size to 1000 (FFT size 2000) 
 * and build the FFT interfaces
 */
TWDeck::TWDeck() : 
  fFFT_R2C(nullptr), fFFT_C2R(nullptr), 
  fSize(1), fFFTSize(1)
{
  BuildFFT();
}

/**
 * @details Initialize the WaveDeck to `n` and built the 
 * FFT interfaces
 *
 * @param n WaveDeck size
 */
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

/**
 * @details Register the `filter` under the key `filter_name` in the list 
 * #fFilters. The filters can be recalled from TWDeck by their `filter_name`.
 * When storing a filter, the used must specify if this in the 
 * time (#wdeck::kReal) or complex (#wdeck::kComplex) domain, so that the missing 
 * representation can be computed. 
 * If the `padding` option is set true, the filter's values containers 
 * (i.e. TWDeckWfmFilter#fWfm, TWDeckWfmFilter#fWfm_re, TWDeckWfmFilter#fWfm_im)
 * are "padded" with zeroes up to size `filter_size + wavedeck_size`. This is 
 * important if the user wants to avoid circular problems in the waveform 
 * convolution/deconvolution
 *
 * @param filter_name Filter name (key)
 * @param filter Filter
 * @param padding Apply waveform padding
 * @param kDomain
 */
void TWDeck::RegisterFilter(const char* filter_name, TWDeckWfmFilter* filter, bool padding, wdeck::EWfmDomain kDomain ) {
  TWDeckWfmFilter* pdec_filter = new TWDeckWfmFilter(*filter);
  ResizeFilter(pdec_filter, padding, kDomain);
  int size_tmp = fFFTSize;
  int size = fSize;
  if (padding) size += filter->GetSize();
  BuildFFT(size);
  if (kDomain == wdeck::kReal) {
    FFTR2C(pdec_filter);
  } else if (kDomain == wdeck::kComplex) {
    FFTC2R(pdec_filter);
  }
  fFilters.insert(std::make_pair(filter_name, pdec_filter));

  BuildFFT(size_tmp);
}


/**
 * @details Perform a FFT from Real to complex of 
 * the TWDeckWfm `wfm` assuming the size given by #fSize.
 *
 * @param wfm Waveform to be transformed
 */
void TWDeck::FFTR2C(TWDeckWfm* wfm) {
  fFFT_R2C->SetPoints(&wfm->GetWfm().at(0));
  fFFT_R2C->Transform();
  double xre[100000] = {0};
  double xim[100000] = {0};
  fFFT_R2C->GetPointsComplex(xre, xim);
  for (int i=0; i<fFFTSize; i++) {
    wfm->GetWfmRe()[i] = xre[i];
    wfm->GetWfmIm()[i] = xim[i];
  }

  return;
}

/**
 * @details Perform a FFT from Real to complex of 
 * the TWDeckWfm `wfm` with the specified `size`.
 * If `size` does not match the one currently used by #fFFT_R2C, 
 * the FFT engine is delete and recreated with the correct size. 
 * The original #fFFTSize is then restored after the transformation
 * has been performed. 
 *
 * @param wfm Waveform to be transformed
 * @param size Transform size
 */
void TWDeck::FFTR2C(TWDeckWfm* wfm, int size) {
  int size_fft_tmp = fFFTSize;
  BuildFFT(size);
  fFFT_R2C->SetPoints(&wfm->GetWfm().at(0));
  fFFT_R2C->Transform();
  double xre[100000] = {0};
  double xim[100000] = {0};
  fFFT_R2C->GetPointsComplex(xre, xim);
  for (int i=0; i<fFFTSize; i++) {
    wfm->GetWfmRe()[i] = xre[i];
    wfm->GetWfmIm()[i] = xim[i];
  }
  // restore original FFT size
  BuildFFT(size_fft_tmp);
  return;
}

/**
 * @details Perform a FFT from Complex to Real of 
 * the TWDeckWfm `wfm` assuming the size given by #fSize.
 * The output is scaled for the inverse of the size of the 
 * transform to account for the implementation of the FFT in fftw3.
 *
 * @param wfm Waveform to be transformed
 */
void TWDeck::FFTC2R(TWDeckWfm* wfm) {
  fFFT_C2R->SetPointsComplex(&wfm->GetWfmRe().at(0), &wfm->GetWfmIm().at(0));
  fFFT_C2R->Transform();
  double* v = fFFT_C2R->GetPointsReal();
  double scale_ = 1./fFFTSize;
  for (int j=0; j<wfm->GetSize(); j++)
    wfm->GetWfm().at(j) = v[j] * scale_;

  return;
}

/**
 * @details Perform a FFT from Complex to Real of 
 * the TWDeckWfm `wfm` with the specified `size`.
 * If `size` does not match the one currently used by #fFFT_C2R, 
 * the FFT engine is delete and recreated with the correct size. 
 * The original #fFFTSize is then restored after the transformation
 * has been performed. 
 * The output is scaled for the inverse of the size of the 
 * transform to account for the implementation of the FFT in fftw3.
 *
 * @param wfm Waveform to be transformed
 * @param size Transform size
 */
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

/**
 * @details Add the sample `wfm` to the `model` recording both the 
 * waveform and the spectrum.
 */
void TWDeck::Add2Model(TWDeckWfm* wfm, TWDeckWfmModel* model) {
  if (fFFTSize != wfm->GetSize()) BuildFFT(wfm->GetSize());
  FFTR2C(wfm);
  
  model->AddWaveform(&wfm->GetWfm  ()[0]);
  model->AddSpectrum(&wfm->GetWfmRe()[0], &wfm->GetWfmIm()[0]);
}

/**
 * @details Resize the `filter` containers' size. 
 * If the option `padding` is true the size will be made equal to 
 * `filter_size + waveform_size`, otherwise the WaveDeck size #fSize
 * is imposed
 *
 * @param filter Filter to be resized
 * @param padding Apply padding to filter containers
 */
void TWDeck::ResizeFilter(TWDeckWfmFilter* filter, bool padding, wdeck::EWfmDomain kDomain) {
  int size_filter = filter->GetSize();
  int size_ = fSize;
  if (padding) {
    size_ += size_filter;
  }

  if (size_filter != size_)
  {
    filter->GetWfm  ().resize(size_, 0.);
    filter->GetWfmRe().resize(size_, 0.);
    filter->GetWfmIm().resize(size_, 0.);

    BuildFFT(size_);

    // re-compute the fourier transform
    if (kDomain == wdeck::kReal) {
      FFTR2C(filter);
    } else if (kDomain == wdeck::kComplex) {
      FFTC2R(filter);
    }
  }
  
  return;
}

/**
 * @details Apply the filter registered under `filter_name` in fFilters to 
 * the waveform `wfm`. If the flag `padding` is active, both the filters 
 * and the waveform are resize to have size `filter_size + waveform_size` to 
 * avoid issues due to the boundary conditions of the FFT algorithms.
 *
 * @param wfm
 * @param filter_name
 * @param padding
 */
void TWDeck::ApplyFilter(TWDeckWfm* wfm, TString filter_name, bool padding) {
  if ( !(fFilters.count(filter_name) ) ) {
    printf("TWDeck::ApplyFilter ERROR: no filter named '%s' is registered. quit.\n", 
        filter_name.Data());
    return;
  }

  TWDeckWfmFilter* filter = fFilters.find(filter_name)->second;

  TWDeckWfm wfm_tmp(*wfm);
  if (fSize !=wfm_tmp.GetSize()) {
    fSize = wfm_tmp.GetSize();
    ResizeFilter(filter, padding, filter->GetOriginDomain());
  }

  int size_tmp = filter->GetSize();
  if (padding) size_tmp += wfm_tmp.GetSize();
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

/**
 * @details Apply the `filter` to the waveform `wfm`. 
 * If the option `padding` is set, then both the waveform and the filter
 * are padded with zeroes up to a size given by `filter_size + waveform_size`.
 * Otherwise, the WaveDeck waveform is set as the waveform size and the 
 * filter size is adjusted without extra padding. 
 *
 * @param wfm Waveform to be filtered
 * @param filter Filter to be applied
 * @param padding Apply zero padding
 */
void TWDeck::ApplyFilter(TWDeckWfm* wfm, TWDeckWfmFilter* filter, bool padding) {

  TWDeckWfm wfm_tmp(*wfm);
  if (fSize !=wfm_tmp.GetSize()) {
    fSize = wfm_tmp.GetSize();
  }
  ResizeFilter(filter, padding, filter->GetOriginDomain());

  int size_tmp = filter->GetSize();
  if (padding) size_tmp += wfm_tmp.GetSize();
  BuildFFT(size_tmp);

  wfm_tmp.SetSize(size_tmp);
  
  printf("FFT size is %i\n", fFFTSize);
  FFTR2C(&wfm_tmp);
  
  double fltr_shift = filter->GetShift();
  bool apply_shift = !(fltr_shift == 0);

  int nloop = 0.5*size_tmp+1; 
  if (padding) nloop = size_tmp;


  for (int i = 0; i < nloop; i++) {
    TComplex F = TComplex(filter->GetWfmRe()[i], filter->GetWfmIm()[i]);
    TComplex W = TComplex(wfm_tmp.GetWfmRe()[i], wfm_tmp.GetWfmIm()[i]);
  
    TComplex C = W*F;
    //printf("[%i] %g +i(%g) = [[F]%g +i(%g)] * [[W]%g +i(%g)]\n", 
       //i, C.Re(), C.Im(), F.Re(), F.Im(), W.Re(), W.Im() );

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

void TWDeck::BuildFFT(int size) 
{
  fFFTSize = size;

  if (fFFT_R2C) {delete fFFT_R2C; fFFT_R2C = nullptr;}
  if (fFFT_C2R) {delete fFFT_C2R; fFFT_C2R = nullptr;}

  fFFT_R2C = TVirtualFFT::FFT(1, &fFFTSize, "M R2C K");
  fFFT_C2R = TVirtualFFT::FFT(1, &fFFTSize, "M C2R K");
}

/**
 * @details This method produces a waveform with same power spectral 
 * density of the model. To do so, unitary white noise is filtered with 
 * with a filter having the same spectral density of the `model` and a constant
 * phase. 
 *
 * @return 
 */
TWDeckWfm* TWDeck::Produce(TWDeckWfmModel* model) 
{
  printf("TWDeck::Produce\n");
  const int size = model->GetSize();
  double ytime[size]; double fre[size]; double fim[size];

  // create white noise waveform
  for (int i=0; i<size; i++) ytime[i] = gRandom->Gaus(0, 1.);
  TWDeckWfm* ywfm = new TWDeckWfm(size, ytime);
  SetSize(size);

  // create filter from model
  // note that a factor 1/size must be considered to account for 
  // fftw scaling convention
  auto vpsd = model->GetSpectralDensityPoints();
  for (int j=0; j<size; j++) {
    fre[j] = TMath::Sqrt(vpsd.at(j)/size); 
    fim[j] = 0.;
  }
  TWDeckWfmFilter filter(size, fre, fim);

  //apply the filter
  ApplyFilter(ywfm, &filter, false);

  return ywfm;
}

