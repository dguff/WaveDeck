/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : example_wdeck_waveproducer.cc
 * @created     : venerdì feb 04, 2022 16:04:31 CET
 *
 *
 * \page producer Example of synthetic waveform production
 * 
 * Source file: `examples/example_wdeck_waveproducer.C`
 *
 * This script shows the basic usage of **WaveDeck** for the production of
 * a synthetic waveform. 
 * The first part of the code shows how to convolve the time series of the 
 * pulses with the impulse response function in order to obtain a waveform. 
 * 
 * ![example of impulse response convolution](example_producer_wave.png)
 * 
 * The second part of the script completes the production of a simulated waveform 
 * by producing a noise signal according to a given power spectral density 
 * (see for a more detailed documentation).
 *
 * ![example of synthetic noise](example_producer_noise.png)
 *
 * \include example_wdeck_waveproducer.cc
 *  
 */

#include <iostream>
#include <cstdlib>
#include <vector>
#include "TStyle.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TApplication.h"

#include "TF1.h"
#include "TWDeck.h"
#include "TWDeckWfm.h"
#include "TWDeckWfmModel.h"
#include "TWDeckUtils.h"


std::vector<double> build_spe_response(int size, double t0, double t1);

int example_wdeck_waveproducer()
{
  const bool use_batch_mode = (std::getenv("CI")) ? true : false;
  if (use_batch_mode) {
    gROOT->SetBatch(kTRUE);
  }

  //- - - - - - - - - - - - - - - - - - - - - - - -environment setup
  const int size = 2500;
  const double t0 = 0.;
  const double t1 = 1e4;
  double Dt = t1-t0;       //!< Time interval
  double dt = Dt / size;   //!< Δt between two samples
  printf("waveform size: %i; t0 = %.2f, t1 = %.2f -> dt = %g\n", 
      size, t0, t1, dt);
  double ddt = 1;        //!< desired time resolution limit (1 ns)
  int    ssize = Dt / ddt;
  std::vector<double> xt  = linspace(t0, t1, size+1); 
  std::vector<double> xtt = linspace(t0, t1, ssize+1); 
  std::vector<double> xf = linspace(0., size/Dt, size+1);
  printf("waveform super-size: %i; t0 = %.2f, t1 = %.2f -> dt = %g\n", 
      ssize, t0, t1, ddt);
  //- - - - - - - - - - - - - - - - - - - - define original waveform 
  // original signal waveform
  TH1D* hOrigin = new TH1D("hs", "Original waveform", ssize, &xtt[0]);
  // create a δ-pulse in i = 600, 800 and 1500;
  hOrigin->Fill( 600);
  hOrigin->Fill( 800);
  hOrigin->Fill(1500);
   
  // define spe reponse
  std::vector<double> xh = build_spe_response(ssize, t0, t1);

  //- - - - - - - - - - - - - - - - - - - - - - - instance wavedeck
  TWDeck wdeck(ssize);
  double xs[ssize]; 
  for (int i=0; i<ssize; i++) xs[i] = hOrigin->GetBinContent(i);
  TWDeckWfm*       wy = new TWDeckWfm(ssize, xs);
  printf("building filter...\n");
  TWDeckWfmFilter* wh = new TWDeckWfmFilter(ssize, &xh[0]);
  wh->SetShift(500/ddt);
  printf("computing response FFT...\n");
  printf("filter size: %i - wdeck size: %i\n", 
      wh->GetSize(), wdeck.GetFFTSize());
  wdeck.FFTR2C(wh);

  printf("applying response...\n");
  wdeck.ApplyFilter(wy, wh, true);

  //- - - - - - - - - - -Create an arbitrary noise spectral density
  //- - - - - - - - - - - - - - - - - - (e.g. a Lorentian function)
  printf("Noise model...\n");
  TWDeckWfmModel* noise_model = new TWDeckWfmModel(size); 
  noise_model->SetOriginDomain(wdeck::kComplex);
  double nre_[size] = {0}; double nim_[size] = {0.};
  for (int j=0; j<size; j++) {
    nre_[j] = sqrt(TMath::CauchyDist(j, 0.25*size, 0.05*size));
  }
  noise_model->AddSpectrum(nre_, nim_);

  TWDeckWfm* wnoise = wdeck.Produce(noise_model);
  wdeck.SetSize(size);
  wdeck.FFTR2C(wnoise);
  
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - -plots 
  TGraph* gs = new TGraph(ssize, &xtt[0], &xs[0]);
  TGraph* gh = new TGraph(ssize, &xtt[0], &xh[0]);
  TGraph* gy = new TGraph(ssize, &xtt[0], &wy->GetWfm()[0]);

  gStyle->SetOptStat(0);
  TCanvas* cOrigins = new TCanvas("cOrigins", "original signals", 0, 0, 800, 600);
  cOrigins->SetTicks(1, 1);
  gh->SetNameTitle("gh", "Response function;Time [ns];Amplitude [a.u.]");
  //gs->SetNameTitle("gs", "Original signal;Time [ns];Amplitude [a.u.]");
  gy->SetNameTitle("gy", "Convolved signal;Time [ns];Amplitude [a.u.]");
  hOrigin->SetTitle("Original signal;Time [ns]; Amplitude [a.u.]");
  hOrigin->SetLineWidth(2); hOrigin->SetLineColor(kAzure-4);
  hOrigin->Draw("hist");
  hOrigin->GetYaxis()->SetRangeUser(-3., +3);
  gh->SetLineStyle(7);
  gh->SetLineWidth(2); gh->SetLineColor(kGray+1) ; gh->Draw("l");
  //gs->SetLineWidth(2); gs->SetLineColor(kAzure+4); gs->Draw("l");
  gy->SetLineWidth(2); gy->SetLineColor(kBlue+1) ; gy->Draw("l");
  TLegend* leg_sig = new TLegend(0.6, 0.88, 0.88, 0.6);
  leg_sig->AddEntry(hOrigin, hOrigin->GetTitle(), "l");
  leg_sig->AddEntry(gh, gh->GetTitle(), "l");
  leg_sig->AddEntry(gy, gy->GetTitle(), "l");
  leg_sig->SetTextFont(43); leg_sig->SetTextSize(22);
  leg_sig->Draw();

  std::vector<double> xn = wnoise->GetWfm();
  std::vector<double> xN = wnoise->GetSpectralDensityPoints();
  std::vector<double> xM = noise_model->GetSpectralDensityPoints();
  TGraph* gn    = new TGraph(size, &xt[0], &xn[0]);
  TGraph* gnpds = new TGraph(0.5*size+1, &xf[0], &xN[0]);
  TGraph* gmpds = new TGraph(0.5*size+1, &xf[0], &xM[0]);
  TCanvas* cNoise = new TCanvas("cnoise", "noise plots", 0, 0, 800, 600);
  cNoise->Divide(1, 2);
  cNoise->cd(1); 
  gn->SetNameTitle("gnoise", "Synthetic noise;Time [ns];Amplitude [a.u.]");
  gn->SetLineColorAlpha(kGray+2, 0.5);
  gn->Draw("awl");

  cNoise->cd(2);
  gnpds->SetNameTitle("gnPDS", "Synthetic noise;Frequency [a.u.];Power Spectral Density");
  gmpds->SetNameTitle("gmPDS", "Noise Model;Frequency [a.u.];Power Spectral Density");
  gnpds->SetLineWidth(2); gnpds->SetLineColorAlpha(kGray+2, 0.5); gnpds->Draw("awl");
  gmpds->SetLineWidth(3); gmpds->SetLineColorAlpha(kRed +2, 1.0); gmpds->Draw("l");
  TLegend* leg_noise = new TLegend(0.6, 0.88, 0.88, 0.6);
  leg_noise->AddEntry(gnpds, gnpds->GetTitle(), "l");
  leg_noise->AddEntry(gmpds, gmpds->GetTitle(), "l");
  leg_noise->Draw();

  if (use_batch_mode) {
    cOrigins->SaveAs("example_2_producer_wave.png");
    cNoise->SaveAs("example_2_producer_noise.png");
  }

  return 0;
}

std::vector<double> build_spe_response(int size, double t0, double t1) {
  std::vector<double> xh(size, 0.);
  std::vector<double> xt = linspace(t0, t1, size+1);
  double Dt = t1 -t0;

  // define nominal spe response
  auto l_spe = [](double* x, double* p) {
    double B = p[0];
    double C = p[1];
    double A = -B -C;
    double t0 = p[2];
    double t_rise = p[3];
    double t_fall = p[4];
    double t_reco = p[5];

    double y = 0.;
    double t = x[0] - t0;
    if (t > 0) 
      y = A*exp(-t/t_rise) + B*exp(-t/t_fall) + C*exp(-t/t_reco);
    return y;
  };

  TF1 fspe("fspe", l_spe, t0, t1, 6);
  TF1 fgaus("fgaus", "TMath::Gaus(x, 200, [0], true)", t0, t1);

  TF1 fspeg("fspeg", "CONV(fspe, fgaus)", t0 - 0.1*Dt, t1 + 0.1*Dt);
  fspeg.SetNpx(2*size);
  fspeg.SetParameters(8, -4, 500., 20, 300, 800, 50);

  int it = 0;
  for (auto &t : xt) {
    xh[it] = fspeg.Eval(t);
    ++it;
  }

  return xh;
}

int main(int argc, char** argv) {
  const bool use_batch_mode = 
    (std::getenv("CI") || std::getenv("GITHUB_ACTIONS")) ? true : false;

  TApplication* tapp = new TApplication("example_wdeck_waveproducer", &argc, argv);
  example_wdeck_waveproducer();
  if (!use_batch_mode) tapp->Run();

  return 0;
}
