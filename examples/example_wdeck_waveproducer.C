/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : example_wdeck_waveproducer
 * @created     : venerdì feb 04, 2022 16:04:31 CET
 */

#include <iostream>
#include <vector>
#include "TStyle.h"
#include "TCanvas.h"
#include "TRandom3.h"

#include "TF1.h"
#include "TWDeck.h"
#include "TWDeckUtils.h"


std::vector<double> build_spe_response(int size, double t0, double t1);

int example_wdeck_waveproducer()
{
  //- - - - - - - - - - - - - - - - - - - - - - environment setup
  const int size = 2500;
  const double t0 = 0.;
  const double t1 = 1e4;
  double Dt = t1-t0;       //!< Time interval
  double dt = Dt / size;   //!< Δt between two samples
  printf("waveform size: %i; t0 = %.2f, t1 = %.2f -> dt = %g\n", 
      size, t0, t1, dt);
  double ddt = 1;        //!< desired time resolution limit (1 ns)
  int    ssize = Dt / ddt;
  std::vector<double> xt = linspace(0., 1e4, ssize+1); 
  printf("waveform super-size: %i; t0 = %.2f, t1 = %.2f -> dt = %g\n", 
      ssize, t0, t1, ddt);
  // original signal waveform
  TH1D* hOrigin = new TH1D("hs", "Original waveform", ssize, &xt[0]);
  // create a δ-pulse in i = 600, 800 and 1500;
  hOrigin->Fill(3000, 1./ddt);
   
  // define spe reponse
  std::vector<double> xh = build_spe_response(ssize, t0, t1);

  //- - - - - - - - - - - - - - - - - - - - - - - instance wavedeck
  TWDeck wdeck(ssize);
  double xs[ssize]; 
  for (int i=0; i<ssize; i++) xs[i] = hOrigin->GetBinContent(i);
  TWDeckWfm*       wy = new TWDeckWfm(ssize, xs);
  TWDeckWfmFilter* wh = new TWDeckWfmFilter(ssize, &xh[0]);
  wh->SetShift(500/ddt);
  wdeck.FFTR2C(wh);

  wdeck.ApplyFilter(wy, wh, true);

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - -plots 
  TGraph* gs = new TGraph(ssize, &xt[0], &xs[0]);
  TGraph* gh = new TGraph(ssize, &xt[0], &xh[0]);
  TGraph* gy = new TGraph(ssize, &xt[0], &wy->GetWfm()[0]);

  TCanvas* cOrigins = new TCanvas("cOrigins", "original signals", 0, 0, 800, 600);
  cOrigins->SetTicks(1, 1);
  gh->SetNameTitle("gh", "Response function;Time [ns];Amplitude [a.u.]");
  gs->SetNameTitle("gs", "Original signal;Time [ns];Amplitude [a.u.]");
  gy->SetNameTitle("gy", "Convolved signal;Time [ns];Amplitude [a.u.]");
  hOrigin->Draw("hist");
  hOrigin->GetYaxis()->SetRangeUser(-1., +3);
  gh->SetLineWidth(2); gh->SetLineColor(kGray+1); gh->Draw("l");
  gs->SetLineWidth(2); gs->SetLineColor(kAzure+4); gs->Draw("l");
  gy->SetLineWidth(2); gy->SetLineColor(kBlue+1); gy->Draw("l");

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
