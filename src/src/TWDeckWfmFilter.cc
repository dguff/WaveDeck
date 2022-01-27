/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : TWDeckWfmFilter
 * @created     : mercoledì gen 26, 2022 12:28:34 CET
 */

#include "TWDeckWfmFilter.h"

ClassImp(TWDeckWfmFilter)

TWDeckWfmFilter::TWDeckWfmFilter() : TWDeckWfm()
{ }

TWDeckWfmFilter::TWDeckWfmFilter(int N, double* x) : TWDeckWfm(N, x) {}

TWDeckWfmFilter::TWDeckWfmFilter(int N, wdeck::EFltrShape kShape, double bw) {
  SetSize(N);
  if (kShape==wdeck::kSqr) {
    if (N <= bw) SetSize(2*bw);
    double h = 1./bw;
    for (int i=0; i<bw; i++) fWfm.at(i) = h;
    SetShift(0.5*bw);
    SetBandwidth(bw);
  } else if (kShape == wdeck::kTrngl) {
    double h = 1./bw;
    double slope = h/bw;
    if (N <= bw) SetSize(3*bw);
    for (int i=0; i<2*bw; i++) 
      (i <= bw) ? fWfm.at(i) = i*slope : fWfm.at(i) = h-slope*(i-bw);
    SetBandwidth(bw);
    SetShift(bw);
  } else if (kShape == wdeck::kGauss) {
    double offset = 4*bw;
    if (N <= 8*bw) SetSize(8*bw);
    for (int i=0; i<8*bw; i++)
      fWfm.at(i) = TMath::Gaus(i, offset, bw, true);
    SetBandwidth(bw);
    SetShift(offset);
  } else if (kShape == wdeck::kDiff) {
    if (N <= 2*bw) SetSize(3*bw);
    double h = 1./bw;
    for (int i=0; i<2*bw; i++) 
      (i<bw) ? fWfm.at(i) = h : fWfm.at(i) = -h;
    SetShift(bw);
    SetBandwidth(bw);
  } else {
    printf("TWDeckWfmFilter::TWDeckWfmFilter(int N = %i, EFltrShape = %i) ERROR:", 
        N, kShape);
    printf("Unknown filter shape.\n");
  }
}

TWDeckWfmFilter::TWDeckWfmFilter(const TWDeckWfmFilter& filter) : TWDeckWfm(filter) 
{
  fBandwidth = filter.fBandwidth;
  fShift     = filter.fShift;
}

TWDeckWfmFilter::~TWDeckWfmFilter() {}
