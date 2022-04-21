/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : TWDeckWfmModel.cc
 * @created     : venerdì gen 28, 2022 10:22:43 CET
 */

#include "TWDeckWfmModel.h"

ClassImp(TWDeckWfmModel)

TWDeckWfmModel::TWDeckWfmModel() : 
  fNSampleWave(0), fNSampleSpectrum(0), 
  fSpectralDensityHist(nullptr), fWaveDensityHist(nullptr) {}

  /**
   * @details Create a new waveform model of size `N`
   *
   * @param N waveform model size
   */
TWDeckWfmModel::TWDeckWfmModel(int N) : 
  fNSampleWave(0), fNSampleSpectrum(0), 
  fSpectralDensityHist(nullptr), fWaveDensityHist(nullptr)  
{
  SetSize(N);
}

/**
 * @details Create a new waveform model of size `N` and 
 * initialize it with the values in `data`
 *
 * @param N waveform model size
 * @param data first waveform values entry
 */
TWDeckWfmModel::TWDeckWfmModel(int N, double* data) :
  fNSampleWave(0),  fNSampleSpectrum(0), 
  fSpectralDensityHist(nullptr), fWaveDensityHist(nullptr)  
{
  SetSize(N);
  AddWaveform(data);
}

/**
 * @details Copy constructor (patrially reimplemented from 
 * TWDeckWfm::TWDeckWfm(const TWDeckWfm&))
 */
TWDeckWfmModel::TWDeckWfmModel(const TWDeckWfmModel& model) : TWDeckWfm(model)
{
  fSpectralDensity.resize(fSize);
  for (int i=0; i<fSize; i++) {
    fSpectralDensity[i] = model.fSpectralDensity[i];
  }
  fNSampleWave         = model.fNSampleWave;
  fNSampleSpectrum     = model.fNSampleSpectrum;
  fSpectralDensityHist = (TH2D*)model.fSpectralDensityHist->Clone();
  fWaveDensityHist     = (TH2D*)model.fWaveDensityHist->Clone();
}

TWDeckWfmModel::~TWDeckWfmModel()
{
  if (fSpectralDensityHist) {delete fSpectralDensityHist;}
  if (fWaveDensityHist    ) {delete fWaveDensityHist;}
}

/**
 * @details Resize the waveform model to have size `n`
 *
 * @param n waveform model size
 */
void TWDeckWfmModel::SetSize(int n)  {
  TWDeckWfm::SetSize(n);
  fSpectralDensity.resize(n, 0.);
}

/**
 * @details Create the waveform density histogram #fWaveDensityHist.
 * In the x-axis the histogram has a number of bins equal to the waveform 
 * model size, while the y-axis is divided into 100 bins. 
 * The y-axis range is computed based on the values provided in `data`
 *
 * @param data waveform values
 */
void TWDeckWfmModel::BuildWaveDensity(double* data) {
  double ymax = *std::max_element(data, data+fSize);
  double ymin = *std::min_element(data, data+fSize);

  double yax_min = 0.; double yax_max = 0.;
  (ymax < 0.) ? yax_max = 0.7*yax_max : yax_max = 1.3*ymax; 
  (ymin < 0.) ? yax_min = 1.3*ymin : yax_min = 0.7*ymin;
  std::vector<double> xbins = linspace(0., (double)fSize, fSize+1);
  std::vector<double> ybins = linspace(yax_min, yax_max, 101);
  fWaveDensityHist = new TH2D(
      Form("%s_wdensity", fName.Data()), 
      Form("%s wfm density:Time [ticks]:Amplitude", fTitle.Data()), 
      fSize, &xbins[0], 100, &ybins[0]);
  return;
}

/**
 * @details Create the spectral density histogram #fSpectralDensityHist. 
 * Because of the symmetry of the Fourier transform, 
 * the histogram x-axis is divided into a number of bins equal to #fSize*0.5+1, 
 * while the y-axis is divided into 50 bins equally spaced in log scale. 
 * The y-axis range is determined based on the value of the Fourier 
 * cofficients provided by the `xre` and `xim` arrays.
 *
 * @param xre Real parts of the Fourier coefficients
 * @param xim Imaginary parts of the Fourier coefficients
 */
void TWDeckWfmModel::BuildSpectralDensity(double* xre, double* xim) {

  std::vector<double> vmag(fSize, 0.);
  for (int i=0; i<fSize; i++) {
    TComplex c_(xre[i], xim[i]);
    vmag.at(i) = c_.Rho2();
  }

  double ymax = *(std::max_element(vmag.begin(), vmag.end()));
  double ymin = *(std::min_element(vmag.begin(), vmag.end()));

  if (ymin == 0) ymin = 1e-3;
  std::vector<double> xbins = linspace(0., (double)fSize*0.5, fSize*0.5+2);
  std::vector<double> ybins_ex = linspace(0.8*log10(ymin), 1.2*log10(ymax), 51);
  std::vector<double> ybins(ybins_ex);
  for (auto &v : ybins) v = TMath::Power(10., v);

  fSpectralDensityHist = new TH2D(
      Form("%s_sdensity", fName.Data()), 
      Form("%s spectral density", fTitle.Data()), 
      fSize*0.5+1, &xbins[0], 50, &ybins[0]);
  return;
}


/**
 * @details Add the waveform contained in `data` to the model: 
 * the function updates the average waveform #fWfm and fills the 
 * waveform density histogram.
 * In case this is the first waveform sample, the waveform density histogram 
 * is created by calling the TWDeckWfmModel::BuildWaveDensity function. 
 *
 * @param data waveform sample
 */
void TWDeckWfmModel::AddWaveform(double* data) {
  if (fNSampleWave==0) BuildWaveDensity(data);

  fNSampleWave++;
  int ii=0;
  for (auto &v : fWfm) {
    v = ((fNSampleWave-1)*v + data[ii]) / fNSampleWave;
    fWaveDensityHist->Fill(ii, data[ii]);
    ++ii;
  }
  return;
}

/**
 * @details Add the Fourier coefficients contained in the `re` and `im`
 * arrays to the model: the average value of the Fourier coefficients, as well 
 * as the average spectral density is updated and the spectral density histogram 
 * is filled. 
 * If this is the first spectrum sample, the spectral density histogram is created 
 * by means of the TWDeckWfmModel::BuildSpectralDensity function.
 *
 * @param re
 * @param im
 */
void TWDeckWfmModel::AddSpectrum(double* re, double* im) 
{
  if (fNSampleSpectrum==0) BuildSpectralDensity(re, im);

  fNSampleSpectrum++;
  for (int j=0; j<fSize; j++) {
    TComplex c_(re[j], im[j]);
    fWfm_re[j] = ((fNSampleSpectrum-1)*fWfm_re[j] + re[j]) / fNSampleSpectrum;
    fWfm_im[j] = ((fNSampleSpectrum-1)*fWfm_im[j] + im[j]) / fNSampleSpectrum;

    fSpectralDensity.at(j) = ((fNSampleSpectrum-1)*fSpectralDensity.at(j) + c_.Rho2()) / fNSampleSpectrum;
    fSpectralDensityHist->Fill(j, c_.Rho2());
  }
  return;
}

