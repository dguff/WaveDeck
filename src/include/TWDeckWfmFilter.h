/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : TWDeckWfmFilter
 * @created     : mercoledì gen 26, 2022 12:25:23 CET
 */

#ifndef TPDECWFMFILTER_H

#define TPDECWFMFILTER_H

#include "TNamed.h"
#include "TWDeckWfm.h"

namespace wdeck {
  enum EFltrShape  {kSqr=0, kTrngl = 1, kGauss = 2, kDiff = 3};
}

class TWDeckWfmFilter : public TWDeckWfm {
  public : 
    TWDeckWfmFilter();
    TWDeckWfmFilter(int N, double* x);
    TWDeckWfmFilter(int N, wdeck::EFltrShape kShape, double bw);
    TWDeckWfmFilter(const TWDeckWfmFilter& filter);
    ~TWDeckWfmFilter();

    void SetBandwidth(double x) {fBandwidth = x;}
    double GetBandwidth() {return fBandwidth;}
    void SetShift(double s) {fShift = s;}
    double GetShift() {return fShift;}

  protected:
    double fBandwidth;
    double fShift;

  public:
    ClassDef(TWDeckWfmFilter, 1)
};


#endif /* end of include guard TPDECWFMFILTER_H */

