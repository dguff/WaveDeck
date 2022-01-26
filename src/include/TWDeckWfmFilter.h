/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : TWDeckWfmFilter
 * @created     : mercoledì gen 26, 2022 12:25:23 CET
 */

#ifndef TPDECWFMFILTER_H

#define TPDECWFMFILTER_H

#include "TNamed.h"
#include "TWDeckWfm.h"

class TWDeckWfmFilter : public TWDeckWfm {
  public : 
    TWDeckWfmFilter();
    TWDeckWfmFilter(int N, double* x);
    TWDeckWfmFilter(const TWDeckWfmFilter& filter);
    ~TWDeckWfmFilter();

    void SetBandwidth(double x) {fBandWindth = x;}
    double GetBandwidth() {return fBandWindth;}
    void SetShift(double s) {fShift = s;}
    double GetShift() {return fShift;}

  protected:
    double fBandWindth;
    double fShift;

  public:
    ClassDef(TWDeckWfmFilter, 1)
};


#endif /* end of include guard TPDECWFMFILTER_H */

