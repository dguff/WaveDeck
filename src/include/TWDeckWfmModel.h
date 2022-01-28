/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : TWDeckWfmModel
 * @created     : venerdì gen 28, 2022 10:17:52 CET
 */

#ifndef TWDECKWFMMODEL_H

#define TWDECKWFMMODEL_H


#include "TWDeckUtils.h"
#include "TH2D.h"
#include "TNamed.h"
#include "TWDeckWfm.h"

class TWDeckWfmModel : public TWDeckWfm {
  public:
    TWDeckWfmModel();
    TWDeckWfmModel(int size);
    TWDeckWfmModel(int size, double *data);
    TWDeckWfmModel(const TWDeckWfmModel& model);
    ~TWDeckWfmModel();

    void AddWaveform(double* data);
    void AddSpectrum(double* re, double* im);

    TH2D* GetSpectralDensityHist() {return fSpectralDensity;}
    TH2D* GetWavefmDensityHist()   {return fWaveDensity;}

    int   GetNSampleWave() {return fNSampleWave;}
    int   GetNSampleSpectrum() {return fNSampleSpectrum;}

  protected:
    int   fNSampleWave;
    int   fNSampleSpectrum;
    TH2D* fSpectralDensity;
    TH2D* fWaveDensity;

    void  BuildWaveDensity(double* data);
    void  BuildSpectralDensity(double* xre, double* xim);

  public:
    ClassDef(TWDeckWfmModel, 1);
};



#endif /* end of include guard TWDECKWFMMODEL_H */

