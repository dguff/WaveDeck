/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : TWDeckWfmFilter.cc
 * @created     : mercoledì gen 26, 2022 12:28:34 CET
 */

#include "TWDeckWfmFilter.h"

ClassImp(TWDeckWfmFilter)

/**
 * @details Inherited from TWDeckWfm::TWDeckWfm()
 */
TWDeckWfmFilter::TWDeckWfmFilter() : TWDeckWfm()
{ }

/**
 * @details Inherited from TWDeckWfm::TWDeckWfm(int)
 *
 * @param N filter waveform size
 */
TWDeckWfmFilter::TWDeckWfmFilter(int N) : TWDeckWfm(N) {}

/**
 * @details Inherited from TWDeckWfm::TWDeckWfm(int, double*)
 *
 * @param N filter waveform size
 * @param x filter waveform values
 */
TWDeckWfmFilter::TWDeckWfmFilter(int N, double* x) : TWDeckWfm(N, x) {}

/**
 * @details Inherited from TWDeckWfm::TWDeckWfm(int, double*, double*)
 *
 * @param N filter size
 * @param re Real parts of filter Fourier components
 * @param im Imaginary parts of filter Fourier components
 */
TWDeckWfmFilter::TWDeckWfmFilter(int N, double* re, double* im) 
  : TWDeckWfm(N, re, im) {}

/**
 * @details Facilitated constructor implementing some 
 * pre-defined commonly used filters. The filters 
 * currently implemented in the WaveDeck package are
 * called via the dedicated enumerator #wdeck::EFltrShape.
 * In this particular case, the bandwidth (expressed in ticks)
 * indicates the different filters' width using the following conventions:
 * - #wdeck::kSqr: bandwidth = **half** of the square gate length
 * - #wdeck::kTrngl: bandwidth = **half** of the triangle base
 * - #wdeck::kGauss: bandiwdth = **sigma** of the Gaussian
 * - #wdeck::kDiff: bdnwidth = **half** of the flat smoothing gate 
 *
 * Note that all the above mentioned filters are normalized. 
 * @param N Filter waveform size
 * @param kShape Filter shape 
 * @param bw Filter bandwidth
 */
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

/**
 * @details Copy constructor partially reimplemented from 
 * TWDeckWfm::TWDeckWfm(const TWDeckWfm&)
 */
TWDeckWfmFilter::TWDeckWfmFilter(const TWDeckWfmFilter& filter) : TWDeckWfm(filter) 
{
  fBandwidth = filter.fBandwidth;
  fShift     = filter.fShift;
}

TWDeckWfmFilter::~TWDeckWfmFilter() {}
