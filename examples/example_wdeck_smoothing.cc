/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : example_wdeck_smoothing.cc
 * @created     : martedì gen 25, 2022 11:47:54 CET
 *
 * \page Examples
 * # Examples
 *
 * The following scripts showcase a few of the WaveDeck functionalities and 
 * offers commented examples that hopefully can help integrating WaveDeck
 * into the user code
 *
 * - \subpage convolution
 * - \subpage wiener
 * - \subpage producer
 *
 *
 * \page convolution Example of waveform convolution/smoothing
 *
 * Source file: `examples/example_wdeck_smoothing.cc`
 *
 * This script shows how to use **WaveDeck** to perform waveform convolution 
 * using simple filters defined in the time domain. 
 *
 * In the code we produce a sysntetic waveform assuming white noise and 
 * two pulses with an exponential decay. Using the facilitated constructor
 * of TWDeckWfmFilter::TWDeckWfmFilter(int,wdeck::EFltrShape,double) 
 * we implement a few pre-defined filters to be applied to the waveform. 
 * In addition we also build a simple custom filter. 
 *
 * The **WaveDeck** interface is then used to perform the waveform 
 * consolution with the above filters. 
 *
 * Note that the filters defined in the constructor of 
 * TWDeckWfmFilter::TWDeckWfmFilter(int,wdeck::EFltrShape,double) 
 * are defined as starting at their first point, which will naturally 
 * induce a time _shift_ in the convolution output. To compensate this 
 * behaviour, we introduced a TWDeckWfmFilter::fShift member that 
 * produces an additional _phase shift_ in the convolution which 
 * correct this effect.
 *
 * ![convolution example output](example_convolution.png)
 *
 *
 * \include example_wdeck_smoothing.cc
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
#include "TWDeckWfmFilter.h"
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
  if (x>0.) y = TMath::Exp(-x/2.0);
  return y;
}

int example_wdeck_smoothing(int n_p = 2) {

  const int size = 1024;
  const double t0 = 0.0;
  const double t1 = 204.8;

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
      xv[it] += 17*spe_response(tt);
    }
    xv[it] += gRandom->Gaus(0, 1.0);
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

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -     
  // B U I L D   F I L T E R S
  //
  // In case one needs multiple filters, it is more convenient to have 
  // them all of the same size.
  double bandwidth =  6; // in ticks

  const int n_fltr = 5;
  TWDeckWfmFilter* filter[n_fltr] = {0};
  wdeck::EFltrShape fltr_shape[4] = {
    wdeck::kSqr,   // square box (moving average)
    wdeck::kTrngl, // triangle
    wdeck::kGauss, // Gaussian
    wdeck::kDiff}; // differential
  const char* fltr_name[n_fltr] = {"square", "triangle", "gauss", "diff", "custom"};
  Color_t     fltr_cols[n_fltr] = {kRed+1, kGreen+2, kAzure-4, kOrange+7, kMagenta+1};

  // create pre-defined filters (kSqr, kTrng, kGauss, kDiff)
  for (int i=0; i<n_fltr-1; i++)
    filter[i] = new TWDeckWfmFilter(size, fltr_shape[i], bandwidth);

  // create a custom filter
  double xf_cstm[size] = {0.};
  double scale = 2*TMath::InvPi() / TMath::Sq(bandwidth);
  for (int i=0; i<=2*bandwidth; i++) 
    xf_cstm[i] = sqrt(TMath::Sq(bandwidth) - TMath::Sq(i-bandwidth))*scale ;
  filter[4] = new TWDeckWfmFilter(size, xf_cstm);
  filter[4]->SetBandwidth(bandwidth);
  filter[4]->SetShift(bandwidth);


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -     
  // I N S T A N C E   W D E C K
  TWDeck wdeck(size);
  for (int i=0; i<n_fltr; i++) {
    wdeck.RegisterFilter(fltr_name[i], filter[i], true);
  }
  
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -     
  // A P P L Y   F I L T E R S
  // create copies of the original waveform and apply the different filters
  TWDeckWfm*       wfm_filtered[n_fltr] = {0};
  TGraph*          gw_filtered [n_fltr] = {0};
  cWaveform->cd();
  for (int i=0; i<n_fltr; i++) {
    wfm_filtered[i] = new TWDeckWfm(size, xv);
    wdeck.ApplyFilter(wfm_filtered[i], fltr_name[i], true);

    gw_filtered[i] = new TGraph(size, &xt[0], &wfm_filtered[i]->GetWfm()[0]);
    gw_filtered[i]->SetLineWidth(3);
    gw_filtered[i]->SetLineColor(fltr_cols[i]);
    gw_filtered[i]->Draw("l");
    gw_filtered[i]->SetNameTitle(Form("gW%s", fltr_name[i]), 
        Form("%s filter:Time [#mus];Amplitude [a.u.]", fltr_name[i]));
  }
  
  
  TPad* pFilter = new TPad("cFilter", "filters", 0.6, 0.4, 0.92, 0.92);
  pFilter->Draw();
  pFilter->cd();

  pFilter->SetGrid(1, 1);
  pFilter->SetTicks(1, 1);
  
  TGraph* gfltr[5] = {0};
  std::vector<double> yfltr_tmp[5];
  double offset    = 30; // in ticks
  for (int i=0; i<5; i++) {
    int idiff = filter[i]->GetShift() - offset;
    yfltr_tmp[i] = std::vector<double>(filter[i]->GetWfm());
    if (idiff < 0) {
      for (int j=0; j<TMath::Abs(idiff); j++)
        yfltr_tmp[i].insert(yfltr_tmp[i].begin(), 0.);
    } else if (idiff > 0) {
      yfltr_tmp[i].erase(yfltr_tmp[i].begin(), yfltr_tmp[i].begin()+idiff );
    }

    gfltr[i] = new TGraph(60, &xt[0], &yfltr_tmp[i].at(0));
    gfltr[i]->SetLineWidth(3);
    gfltr[i]->SetLineColor(fltr_cols[i]);
    gfltr[i]->GetYaxis()->SetRangeUser(-0.3, +0.3);

    if (!i) gfltr[i]->Draw("awl");
    else    gfltr[i]->Draw("l");
  }


  return 0;
}

int main(int argc, char* argv[])
{
  int n_pulses = 3;
  if (argc > 1)
    n_pulses = std::atoi(argv[1]);

  TApplication* tapp = new TApplication("example_wdeck_smoothing", &argc, argv);
  
  example_wdeck_smoothing(n_pulses);
  
  tapp->Run();
    
  return 0;
}

