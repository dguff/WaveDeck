/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : TWDeckWfm
 * @created     : martedì gen 25, 2022 11:55:43 CET
 */

#ifndef TWDECKWFM_H

#define TWDECKWFM_H

#include <iostream>
#include <vector>
#include "TNamed.h"
#include "TGraph.h"
#include "TComplex.h"

namespace wdeck{
  enum EWfmDomain {kReal = 0, kComplex = 1};
}

class TWDeckWfm : public TNamed {
  public:
    TWDeckWfm();
    TWDeckWfm(int N);
    TWDeckWfm(int N, double* wfm);
    TWDeckWfm(int N, double* re, double* im);
    TWDeckWfm(const TWDeckWfm& wfm);
    ~TWDeckWfm();
    
    void                   ClearWave();

    void                   LoadSpectrum(int n, double* re, double* im);
    void                   LoadSpectrum(double* re, double* im);
    void                   LoadWave(int n, double* w);
    void                   LoadWave(double* w);

    int                    GetSize  () {return fSize;}
    std::vector<double>&   GetWfm   () {return fWfm;}
    std::vector<double>&   GetWfmRe () {return fWfm_re;}
    std::vector<double>&   GetWfmIm () {return fWfm_im;}
    inline TComplex        
      GetPointComplex(int i) {return TComplex(fWfm_re.at(i), fWfm_im.at(i));}
    inline double         
      GetPointRead(int i) {return fWfm.at(i);}
    double GetSpectralDensity(int i);

    std::vector<TComplex>  GetPointsComplex();
    std::vector<double>&   GetPointsReal() {return fWfm;}
    void                   SetSize(int n);

  protected:
    int                 fSize;
    std::vector<double> fWfm;
    std::vector<double> fWfm_re;
    std::vector<double> fWfm_im;

  public:
    ClassDef(TWDeckWfm, 1);
};


#endif /* end of include guard TWDECKWFM_H */

