/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : TWDeckWfmModel
 * @created     : venerdì gen 28, 2022 10:22:43 CET
 */

#include "TWDeckWfmModel.h"

ClassImp(TWDeckWfmModel)

TWDeckWfmModel::TWDeckWfmModel() : 
  fNSampleWave(0), fNSampleSpectrum(0), 
  fSpectralDensityHist(nullptr), fWaveDensityHist(nullptr) {}

TWDeckWfmModel::TWDeckWfmModel(int N) : 
  fNSampleWave(0), fNSampleSpectrum(0), 
  fSpectralDensityHist(nullptr), fWaveDensityHist(nullptr)  
{
  SetSize(N);
}

TWDeckWfmModel::TWDeckWfmModel(int N, double* data) :
  fNSampleWave(0),  fNSampleSpectrum(0), 
  fSpectralDensityHist(nullptr), fWaveDensityHist(nullptr)  
{
  SetSize(N);
  AddWaveform(data);
}

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

void TWDeckWfmModel::SetSize(int n)  {
  TWDeckWfm::SetSize(n);
  fSpectralDensity.resize(n, 0.);
}

void TWDeckWfmModel::BuildWaveDensity(double* data) {
  double ymax = *std::max_element(data, data+fSize);
  double ymin = *std::min_element(data, data+fSize);

  std::vector<double> xbins = linspace(0., (double)fSize, fSize+1);
  std::vector<double> ybins = linspace(1.3*ymin, 1.3*ymax, 101);
  fWaveDensityHist = new TH2D(
      Form("%s_wdensity", fName.Data()), 
      Form("%s wfm density:Time [ticks]:Amplitude", fTitle.Data()), 
      fSize, &xbins[0], 100, &ybins[0]);
  return;
}

void TWDeckWfmModel::BuildSpectralDensity(double* xre, double* xim) {

  std::vector<double> vmag(fSize, 0.);
  for (int i=0; i<fSize; i++) {
    TComplex c_(xre[i], xim[i]);
    vmag.at(i) = c_.Rho2();
  }

  double ymax = *(std::max_element(vmag.begin(), vmag.end()));
  double ymin = *(std::min_element(vmag.begin(), vmag.end()));

  if (ymin == 0) ymin = 1e-3;
  std::vector<double> xbins = linspace(0., (double)fSize*0.5, fSize*0.5+1);
  std::vector<double> ybins_ex = linspace(0.8*log10(ymin), 1.2*log10(ymax), 51);
  std::vector<double> ybins(ybins_ex);
  for (auto &v : ybins) v = TMath::Power(10., v);

  fSpectralDensityHist = new TH2D(
      Form("%s_sdensity", fName.Data()), 
      Form("%s spectral density", fTitle.Data()), 
      fSize*0.5, &xbins[0], 50, &ybins[0]);
  return;
}


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
