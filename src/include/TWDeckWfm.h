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

class TWDeckWfm : public TNamed {
  public:
    TWDeckWfm();
    TWDeckWfm(int N);
    TWDeckWfm(int N, double* wfm);
    TWDeckWfm(const TWDeckWfm& wfm);
    ~TWDeckWfm();
    
    void                   ClearWave();

    void                   LoadWave(int n, double* w);
    void                   LoadWave(double* w);

    int                    GetSize  () {return fSize;}
    std::vector<double>&   GetWfm   () {return fWfm;}
    std::vector<double>&   GetWfmRe () {return fWfm_re;}
    std::vector<double>&   GetWfmIm () {return fWfm_im;}
    std::vector<TComplex>  GetWfmC  ();
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

