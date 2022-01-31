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

/**
 * @brief Derived waveform class for statistical waveform model building
 *
 * The TWDeckWfmModel class, based on TWDeckWfm, is intended for creating 
 * an easily and accessible container for building waveform models based 
 * on a sample of data. 
 * Single waveform are added to the model via the TWDeckWfmModel::AddWaveform 
 * (or TWDeckWfmModel::AddSpectrum, depending if one wants to add a waveform or
 * its Fourier transform), and the average value of the waveform if automatically 
 * updated. 
 *
 * While the model does not store each individual waveform used for building it, 
 * it contains two density plots (one for the waveform itself, one for the 
 * spectral density), so that the user can later check out the consistency 
 * of the sample. 
 */
class TWDeckWfmModel : public TWDeckWfm {
  public:
    //! Empty constructor 
    TWDeckWfmModel();
    //! Create a waveform model with size `N`
    TWDeckWfmModel(int N);
    //! Create a waveform model with size `N` and initialize it with `data`
    TWDeckWfmModel(int N, double *data);
    //! Copy constructor
    TWDeckWfmModel(const TWDeckWfmModel& model);
    //! Destructor
    ~TWDeckWfmModel();

    //! Add a new waveform to the model
    void AddWaveform(double* data);
    //! Add a new sample to the spectrum model
    void AddSpectrum(double* re, double* im);

    //! Get the model spectral density histogram
    inline TH2D* GetSpectralDensityHist() {return fSpectralDensityHist;}
    //! Get the model waveform density histogram
    inline TH2D* GetWavefmDensityHist()   {return fWaveDensityHist;}

    //! Get the model waveform sample size
    inline int   GetNSampleWave() {return fNSampleWave;}
    //! Get the model spectrum sample size
    inline int   GetNSampleSpectrum() {return fNSampleSpectrum;}
    //! Return the i-th Fourier component spectral density
    inline double GetSpectralDensity(int i) {return fSpectralDensity.at(i);}
    //! Return the Fourier components spectral density
    inline std::vector<double> GetSpectralDensityPoints() {return fSpectralDensity;}

    //! Set waveform model size
    void  SetSize(int n);
  protected:
    int   fNSampleWave;                         //!< Waveform model sample size
    int   fNSampleSpectrum;                     //!< Spectrum model sample size
    TH2D* fSpectralDensityHist;                 //!< Spectrum density histogram
    TH2D* fWaveDensityHist;                     //!< Waveform density histogram
    std::vector<double> fSpectralDensity;       //!< Values of Fourier components spectral density

    //! Create the waveform density histogram
    void  BuildWaveDensity(double* data);
    //! Create the spectrum density histogram
    void  BuildSpectralDensity(double* xre, double* xim);

  public:
    ClassDef(TWDeckWfmModel, 1);
};



#endif /* end of include guard TWDECKWFMMODEL_H */

