/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : example_wdeck_wiener
 * @created     : martedì gen 25, 2022 11:47:54 CET
 */

#include <stdio.h>
#include <iostream>

#include "TApplication.h"
#include "TMath.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TRandom3.h"

#include "TWDeckUtils.h"
#include "TWDeckWfm.h"
#include "TWDeckWfmFilter.h"
#include "TWDeckWfmModel.h"
#include "TWDeck.h"

/**
 * @brief Impulse response 
 *
 * @param x
 *
 * @return 
 */
double spe_response(double x) {
  double y = 0.;
  if (x>0.) y = 20*TMath::Exp(-x/2.0);
  return y;
}

int example_wdeck_wiener(int n_p = 2) {

  gStyle->SetPalette(kSunset);

  const int size = 1024;
  const double t0 = 0.0;
  const double t1 = 204.8;
  double dt = (t1-t0) / size;
  const double noise_rms = 1.0;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -     
  // B U I L D   W A V E F O R M

  double xv  [size] = {0}; // waveform array
  std::vector<double> xt = linspace(t0, t1, size);
  const int nph_true = gRandom->Poisson(n_p);   // true number of p.e.
  std::vector<double> tph_true(nph_true, 0.); // true time of p.e. 
  for (auto &t : tph_true) t = gRandom->Rndm()*100 + 10.;
  
  for (int it=0; it<size; it++) {
    for (auto &tp : tph_true) {
      double tt = xt[it] - tp;
      xv[it] += spe_response(tt);
    }
    xv[it] += gRandom->Gaus(0, noise_rms);
  }

  // plot waveform
  gStyle->SetOptTitle(0);
  TCanvas* cWaveform = new TCanvas("cWaveform", "waveform", 0, 0, 1400, 600);
  cWaveform->SetGrid(1, 1); cWaveform->SetTicks(1, 1);
  cWaveform->SetRightMargin(0.05); cWaveform->SetTopMargin(0.05);
  TGraph* gW = new TGraph(size, &xt.at(0), xv);
  gW->SetNameTitle("gWfm", "Original wfm;Time [#mus];Amplitude [a.u.]");
  gW->SetMarkerStyle(20);
  gW->SetMarkerColorAlpha(kGray+2, 0.5);
  gW->Draw("awp");

  //TWDeckWfm* wfm = new TWDeckWfm(size, xv);
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -     
  // B U I L D   W I E N E R   F I L T E R
  //

  TWDeck wdeck(size);

  //create a model for the noise
  TWDeckWfmModel* wm_noise = new TWDeckWfmModel(size);
  wm_noise->SetNameTitle("noise", "Gaussian noise");
  int n_noise_wfms = 1e3;
  double noise_control[size] = {0};
  for (int jn=0; jn<n_noise_wfms; jn++) {
    for (int it=0; it<size; it++) noise_control[it] = gRandom->Gaus(0, noise_rms);
    TWDeckWfm noise_tmp(size, noise_control);
    wdeck.Add2Model(&noise_tmp, wm_noise);
  }

  TCanvas* cNoiseDensity = new TCanvas("cNoiseDensity", "Noise Density plots", 1000, 600);
  cNoiseDensity->Divide(2, 1);
  cNoiseDensity->cd(1);
  wm_noise->GetWavefmDensityHist()  ->Draw("colz");
  cNoiseDensity->cd(2);
  gPad->SetLogy(1);
  wm_noise->GetSpectralDensityHist()->DrawClone("colz");

  // create template for the single p.e. response
  double spe_xv[size] = {0};
  double spe_offset = 5.;
  for (int i=0; i<100; i++) spe_xv[i] = spe_response(i-spe_offset);
  TWDeckWfm* wfm_spe = new TWDeckWfm(size, spe_xv);
  TGraph gspe(size, &xt[0], spe_xv); 
  double spe_integral = g_integral(&gspe, 0, 100);

  // create template for the delta (here is Gaussian with σ = 0.5*δt)
  double delta_xv[size] = {0};
  double delta_offset = 5.;
  for (int i=0; i<100; i++) 
    delta_xv[i] = spe_integral*TMath::Gaus(i, delta_offset, 0.5*dt, true);
  TWDeckWfm* wfm_delta = new TWDeckWfm(size, delta_xv);

  // assemble the Wiener filter
  wdeck.FFTR2C(wfm_spe);
  wdeck.FFTR2C(wfm_delta);

  TWDeckWfmFilter* wiener = new TWDeckWfmFilter(size);
  for (int i=0; i<0.5*size+1; i++) {
    TComplex cw = 
      TComplex::Conjugate(wfm_spe->GetPointComplex(i))*wfm_delta->GetSpectralDensity(i);
    cw /= 
      (wfm_spe->GetSpectralDensity(i)*wfm_delta->GetSpectralDensity(i) + 
       wm_noise->GetSpectralDensity(i));
    wiener->GetWfmRe()[i] = cw.Re();
    wiener->GetWfmIm()[i] = cw.Im();
  }
  
  return 0;
}

int main(int argc, char* argv[])
{
  int n_pulses = 3;
  if (argc > 1)
    n_pulses = std::atoi(argv[1]);

  TApplication* tapp = new TApplication("example_wdeck_wiener", &argc, argv);
  
  example_wdeck_wiener(n_pulses);
  
  tapp->Run();
    
  return 0;
}

