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

TWDeckWfmFilter::TWDeckWfmFilter(const TWDeckWfmFilter& filter) : TWDeckWfm(filter) 
{
  fBandWindth = filter.fBandWindth;
}

TWDeckWfmFilter::~TWDeckWfmFilter() {
  printf("Calling ~TWDeckWfmFilter\n");
}
