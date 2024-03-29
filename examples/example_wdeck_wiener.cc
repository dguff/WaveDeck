/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : example_wdeck_wiener.cc
 * @created     : martedì gen 25, 2022 11:47:54 CET
 *
 * \page wiener Example of Wiener deconvolution
 *
 * This script shows how to use **WaveDeck** to assemble a Wiener filter
 * and to use it to perform a waveform deconvolution. 
 * The Wiener filter \f$G(\omega)\f$ is defined as 
 * \f[
 * G(\omega) = \frac{H^*(\omega)S(\omega)}{|H(\omega)|^2S(\omega) + N(\omega)}
 * \f]
 * where \f$H(\omega)\f$ is the Fourier transform of the impulse response function
 * while \f$S(\omega)\f$ and \f$N(\omega)\f$ denote the spectral density of the 
 * original signal and the _mean_ spectral density of the noise respectively. 
 *
 * The original signal shape is defined as a very thin Gaussian with standard 
 * deviation given by 1 tick. 
 * The original waveform, the impulse response template, the original signal 
 * shape are treated via the TWDeck interface to compute 
 * their Fourier transform.
 * In particular, TWDeck is used to build the noise model described by the 
 * TWDeckWfmModel class using a sample of simulated noise waveforms. 
 *
 * Finally, we compute the Wiener filter and we use TWDeck to perform the waveform 
 * deconvolution. 
 *
 * ![wiener deconvolution example](example_wiener.png)
 *
 * \include example_wdeck_wiener.cc
 */

#include <stdio.h>
#include <iostream>

#include "TApplication.h"
#include "TMath.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TRandom3.h"

#include "TWDeckPATH.h"
#include "TWDeckUtils.h"
#include "TWDeckWfm.h"
#include "TWDeckWfmFilter.h"
#include "TWDeckWfmModel.h"
#include "TWDeck.h"


int example_wdeck_wiener(int n_p = 2, bool do_display = true) {

  gStyle->SetPalette(kSunset);

  const int size = 10000;
  const double t0 =  0.0;
  const double t1 = 50.0;
  double dt = (t1-t0) / size;
  const double noise_rms = 0.1;
  // D E F I N E   S P E   R E S P O N S E
  TGraph* gspe_origin = new TGraph(Form("%s/examples/spe_generic_model.txt", WDECK_PROJECT_DIR));

  auto spe_response = [gspe_origin](double x) {
    double y = 0;
    if (x>0 && x < 9) y = gspe_origin->Eval(x+1.1);
    return y;
  };

  new TCanvas();
  gspe_origin->Draw("awpl");


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -     
  // B U I L D   W A V E F O R M

  double xv  [size] = {0}; // waveform array
  std::vector<double> xt = linspace(t0, t1, size);
  std::vector<double> xtick = linspace(0.5, 0.5+(size), size+1);
  const int nph_true = gRandom->Poisson(n_p);   // true number of p.e.
  std::vector<double> tph_true(nph_true, 0.); // true time of p.e. 
  for (auto &t : tph_true) t = gRandom->Rndm()*30 + 5.;
  
  for (int it=0; it<size; it++) {
    for (auto &tp : tph_true) {
      double tt = xt[it] - tp;
      xv[it] += spe_response(tt);
    }
    xv[it] += gRandom->Gaus(0, noise_rms);
  }

  TWDeckWfm* wfm_origin = new TWDeckWfm(size, xv);
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -     
  // B U I L D   W I E N E R   F I L T E R
  //

  // - - - - - - - - - - - - - - - - - - - - - - - - - - create a WaveDeck
  TWDeck wdeck(size);

  // - - - - - - - - - - - - - - - - - - - -  create a model for the noise
  TWDeckWfmModel* wm_noise = new TWDeckWfmModel(size);
  wm_noise->SetNameTitle("noise", "Gaussian noise");

  int n_noise_wfms = 1e3; // number of noise waveforms used to build the model
  double noise_control[size] = {0};
  for (int jn=0; jn<n_noise_wfms; jn++) {
    for (int it=0; it<size; it++) noise_control[it] = gRandom->Gaus(0, noise_rms);
    TWDeckWfm noise_tmp(size, noise_control);
    wdeck.Add2Model(&noise_tmp, wm_noise);
  }

  // - - - - - - - - - - - - -  create a template for the s.p.e. response
  double spe_xv[size] = {0};
  double spe_offset = 5.;
  for (int i=0; i<size; i++) spe_xv[i] = spe_response(xt[i]-spe_offset);
  TWDeckWfm* wfm_spe = new TWDeckWfm(size, spe_xv);
  //
  // - - - - - - - - - - - - - - - - - - -  create template for the delta 
  double delta_xv[size] = {0};
  delta_xv[0] = 1.0; 
  TWDeckWfm* wfm_delta = new TWDeckWfm(size, delta_xv);

  // - - - - - - - - - - Use WaveDeck to compute the FFT of the templates
  wdeck.FFTR2C(wfm_spe);
  wdeck.FFTR2C(wfm_delta);

  // - - - - - - - - - - - - - - - - - - - - - Assemble the Wiener filter
  TWDeckWfmFilter* wiener = new TWDeckWfmFilter(size);
  wiener->SetOriginDomain(wdeck::kComplex);
  for (int i=0; i<size; i++) {
    TComplex h  = TComplex::Conjugate(wfm_spe->GetPointComplex(i));
    double   H2 = wfm_spe->GetSpectralDensity(i);
    double   S2 = wfm_delta->GetSpectralDensity(i); 
    double   N2 = wm_noise->GetSpectralDensity(i);
    TComplex cw = h*S2 / (H2*S2 + N2);
    wiener->GetWfmRe()[i] = cw.Re();
    wiener->GetWfmIm()[i] = cw.Im();
  }
  // compensate for the offset in the response template
  wiener->SetShift(-spe_offset/dt);
  
  // - - - - - - - - Apply the filter to a clone of the original waveform 
  TWDeckWfm* wfm_filtered = new TWDeckWfm(*wfm_origin);
  wdeck.ApplyFilter(wfm_filtered, wiener, false);
  for (auto &v : wfm_filtered->GetWfm()) v /= dt;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -     
  // P L O T S
  //
  // - - - - - - - - - - - - - - - - Plot original and post-convolution waveforms
  gStyle->SetOptTitle(0);
  TCanvas* cWaveform = nullptr;
  TCanvas* cNoiseDensity = nullptr;
  TCanvas* cSpectralDensity = nullptr;
  if (do_display) {
    cWaveform = new TCanvas("cWaveform", "waveform", 0, 0, 1400, 600);
    cWaveform->SetGrid(1, 1); cWaveform->SetTicks(1, 1);
    cWaveform->SetRightMargin(0.05); cWaveform->SetTopMargin(0.05);
    TGraph* gWfm_origin = new TGraph(size, &xt.at(0), &wfm_origin->GetWfm()[0]);
    gWfm_origin->SetNameTitle("gWfm", "Original wfm;Time [#mus];Amplitude [a.u.]");
    gWfm_origin->SetMarkerStyle(20);
    gWfm_origin->SetMarkerColorAlpha(kGray+2, 0.5);
    gWfm_origin->Draw("awp");
    TGraph* gW_filtered = new TGraph(size, &xt[0], &wfm_filtered->GetWfm()[0]);
    gW_filtered->SetLineColor(kRed+1);
    gW_filtered->SetLineWidth(2);
    gW_filtered->Draw("l)");


    // - - - - - - - - - - - - - -  Plot noise model waveform and spectral density
    cNoiseDensity = new TCanvas("cNoiseDensity", "Noise Density plots", 1000, 600);
    cNoiseDensity->Divide(2, 1);
    cNoiseDensity->cd(1);
    wm_noise->GetWavefmDensityHist()  ->Draw("colz");
    TGraph* gwmodel_wave = new TGraph(size, &xtick[0], &wm_noise->GetWfm()[0]);
    gwmodel_wave->SetLineWidth(3); gwmodel_wave->SetLineWidth(kMagenta+7);
    gwmodel_wave->Draw("l");
    cNoiseDensity->cd(2);
    gPad->SetLogy(1);
    wm_noise->GetSpectralDensityHist()->DrawClone("colz");
    std::vector<double> wnoise_sd(0.5*size+1, 0.);
    int iw=0;
    for (auto &v : wnoise_sd) {v = wm_noise->GetSpectralDensity(iw); ++iw;}
    TGraph* gwmodel_sptd = new TGraph(0.5*size+1, &xtick[0], &wnoise_sd[0]);
    gwmodel_sptd->SetLineWidth(3); gwmodel_sptd->SetLineWidth(kMagenta+7);
    gwmodel_sptd->Draw("l");

    // - - - - - - - - - - - - - -  Plot the spectral density of the spe response,
    //                                       noise model and original pulse shape
    cSpectralDensity = new TCanvas("cSpectralDensity", "Spectral Density", 
        0, 0, 800, 600);
    cSpectralDensity->SetTicks(1, 1);
    cSpectralDensity->SetGrid(1, 1);
    cSpectralDensity->SetLogy(1);
    std::vector<double> xDdelta = wfm_delta->GetSpectralDensityPoints();
    std::vector<double> xDnoise = wm_noise ->GetSpectralDensityPoints();
    std::vector<double> xDspe   = wfm_spe  ->GetSpectralDensityPoints();
    std::vector<double> xDwiener= wiener   ->GetSpectralDensityPoints();
    TGraph* gDdelta = new TGraph(1+0.5*size, &xtick[0], &xDdelta[0]);
    TGraph* gDnoise = new TGraph(1+0.5*size, &xtick[0], &xDnoise[0]);
    TGraph* gDspe   = new TGraph(1+0.5*size, &xtick[0], &xDspe  [0]);
    TGraph* gDwiener= new TGraph(1+0.5*size, &xtick[0], &xDwiener[0]);
    gDdelta ->SetName("gDdelta" ); gDdelta ->SetLineColor(kBlue);
    gDspe   ->SetName("gDspe"   ); gDspe   ->SetLineColor(kRed+1);
    gDnoise ->SetName("gDnoise" ); gDnoise ->SetLineColor(kGray+2);
    gDwiener->SetName("gDwiener"); gDwiener->SetLineColor(kGreen+2);

    gDspe   ->Draw("awl");
    gDdelta ->Draw(  "l");
    gDnoise ->Draw(  "l");
    gDwiener->Draw(  "l");
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

