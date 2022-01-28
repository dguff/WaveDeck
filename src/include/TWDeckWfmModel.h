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
#include <vector>

class TWDeckWfmModel : public TWDeckWfm {
  public:
    TWDeckWfmModel();
    TWDeckWfmModel(int size);
    TWDeckWfmModel(int size, double *data);
    TWDeckWfmModel(const TWDeckWfmModel& model);
    ~TWDeckWfmModel();

    void AddWaveform(double* data);
    void AddSpectrum(double* re, double* im);

    TH2D* GetSpectralDensityHist() {return fSpectralDensityHist;}
    TH2D* GetWavefmDensityHist()   {return fWaveDensityHist;}

    int   GetNSampleWave() {return fNSampleWave;}
    int   GetNSampleSpectrum() {return fNSampleSpectrum;}
    double GetSpectralDensity(int i) {return fSpectralDensity.at(i);}
    std::vector<double> GetSpectralDensityPoints() {return fSpectralDensity;}

    void  SetSize(int n);
  protected:
    int   fNSampleWave;
    int   fNSampleSpectrum;
    TH2D* fSpectralDensityHist;
    TH2D* fWaveDensityHist;
    std::vector<double> fSpectralDensity;

    void  BuildWaveDensity(double* data);
    void  BuildSpectralDensity(double* xre, double* xim);

  public:
    ClassDef(TWDeckWfmModel, 1);
};



#endif /* end of include guard TWDECKWFMMODEL_H */

