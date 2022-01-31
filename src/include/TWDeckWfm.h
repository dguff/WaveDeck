/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : TWDeckWfm
 * @created     : martedì gen 25, 2022 11:55:43 CET
 *
 */

#ifndef TWDECKWFM_H

#define TWDECKWFM_H

#include <iostream>
#include <vector>
#include "TNamed.h"
#include "TGraph.h"
#include "TComplex.h"

//! WaveDeck enumerators' namespace
namespace wdeck{
  /**
   * @brief Indicates real or complex domain
   */
  enum EWfmDomain {kReal = 0, kComplex = 1};
}

/**
 * @brief Base waveform object for WaveDeck
 * Base class representing a waveform to be processed by WaveDeck, 
 * containing its representation both in the real and in the 
 * complex domain   
 *
 */
class TWDeckWfm : public TNamed {
  public:
    //! Empty constructor
    TWDeckWfm();
    //!Create and empty waveform with size `N`
    TWDeckWfm(int N);
    //! Create a real N-size waveform from array wfm
    TWDeckWfm(int N, double* wfm);
    //! Create a waveform from its Fourier components
    TWDeckWfm(int N, double* re, double* im);
    //! Copy constructor
    TWDeckWfm(const TWDeckWfm& wfm);
    //! Destructor
    ~TWDeckWfm();
    
    //! Reset waveform values to 0
    void ClearWave();
    //! Set waveform Fourier components
    void LoadSpectrum(int n, double* re, double* im);
    //! Set waveform Fourier components
    void LoadSpectrum(double* re, double* im);
    //! Set waveform values
    void LoadWave(int n, double* w);
    //! Set waveform values
    void LoadWave(double* w);

    //! Get waveform size
    int                    GetSize  () {return fSize;}
    /**
     * @brief Get waveform values
     *
     * @return reference to the vector of waveform values
     */
    std::vector<double>&   GetWfm   () {return fWfm;}
    /**
     * @brief Get real part of waveform Fourier coefficient
     *
     * @return reference to the vector of waveform FFT Re
     */
    std::vector<double>&   GetWfmRe () {return fWfm_re;}
    /**
     * @brief Get imaginary part of waveform Fourier coefficient
     *
     * @return reference to the vector of waveform FFT Im
     */
    std::vector<double>&   GetWfmIm () {return fWfm_im;}
    
    //! Get the `i`-th Fourier coefficient
    inline TComplex        
      GetPointComplex(int i) {return TComplex(fWfm_re.at(i), fWfm_im.at(i));}
    //! Get the i-th waveform point
    inline double         
      GetPointReal(int i) {return fWfm.at(i);}
    //! Get spectral density    
    double GetSpectralDensity(int i);
    //! Get vector of spectral density values    
    std::vector<double>    GetSpectralDensityPoints();
    //! Get vector of Fourier coefficients    
    std::vector<TComplex>  GetPointsComplex();
    //! Get reference to the vector of waveform values 
    std::vector<double>&   GetPointsReal() {return fWfm;}
    //! Set waveform size
    void                   SetSize(int n);

  protected:
    int                 fSize;   //!< Waveform size
    std::vector<double> fWfm;    //!< Waveform values 
    std::vector<double> fWfm_re; //!< Waveform real part containers
    std::vector<double> fWfm_im; //!< Waveform imaginary part containers

  public:
    ClassDef(TWDeckWfm, 1);
};


#endif /* end of include guard TWDECKWFM_H */

