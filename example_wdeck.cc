/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : example_wdeck
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
#include "TWDeckWfmFilter.h"
#include "TWDeck.h"

double spe_response(double x) {
  double y = 0.;
  if (x>0.) y = TMath::Exp(-x/2.0);
  return y;
}

int example_wdeck(int n_p = 2) {

  const int size = 1024;
  const double dt = 0.2;
  const double t0 = 0.0;
  const double t1 = 204.8;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -     
  // B U I L D   W A V E F O R M

  double xv  [size] = {0}; // waveform array
  std::vector<double> xt = linspace(t0, t1, size);

  const int nph_true = gRandom->Poisson(3);   // true number of p.e.
  std::vector<double> tph_true(nph_true, 0.); // true time of p.e. 
  for (auto &t : tph_true) t = gRandom->Rndm()*100 + 10.;
  
  for (int it=0; it<size; it++) {
    for (auto &tp : tph_true) {
      double tt = xt[it] - tp;
      xv[it] += 17*spe_response(tt);
    }
    xv[it] += gRandom->Gaus(0, 1.0);
  }

  TGraph* gW = new TGraph(size, &xt.at(0), xv);
  gW->Draw("awpl");

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -     
  // B U I L D   F I L T E R
  // simple squared filter (i.e. moving average)
  double bandwidth = 5; // in ticks
  double xf_sqr[size]  = {0};
  for (int i=0; i<bandwidth*2; i++) xf_sqr[i] = 0.5/bandwidth;
  TWDeckWfmFilter* square_f = new TWDeckWfmFilter(size, xf_sqr);
  square_f->SetBandwidth(bandwidth);
  // simple Gaussian filter 

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -     
  // I N S T A N C E   P D E C
  TWDeck wdeck(size);
  wdeck.RegisterFilter("square", square_f);
  
  // apply square filter
  TWDeckWfm* wfm = new TWDeckWfm(size, xv);
  //printf("wfm size is %i\n", wfm->GetSize());
  wdeck.ApplyFilter(wfm, "square");

  TGraph* gs = new TGraph(size, &xt.at(0), &wfm->GetWfm().at(0));
  gs->SetLineColor(kRed);

  gs->Draw("l");
  

  

  return 0;
}

int main()
{
    
    return 0;
}

