/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : TWDeckWfmFilter.h
 * @created     : mercoledì gen 26, 2022 12:25:23 CET
 */

#ifndef TPDECWFMFILTER_H

#define TPDECWFMFILTER_H

#include "TNamed.h"
#include "TWDeckWfm.h"

namespace wdeck {
  /**
   * @brief Predefined filter shapes
   */
  enum EFltrShape  {
    kSqr=0,     //!< Squared filter
    kTrngl = 1, //!< Triangular filter
    kGauss = 2, //!< Gaussian filter
    kDiff = 3   //!< Smoothed differential filter
  };
}

/**
 * @brief Derived waveform class used to implement WaveDeck filters
 *
 * The #TWDeckWfmFilter class, based on the #TWDeckWfm waveform class,
 * is used to treat filters in the context of the WaveDeck interface. 
 *
 * In addition to the #TWDeckWfm attributes and memebers, this class
 * present two additional attributes, namely the filter's bandwidth 
 * and a possible shift.
 *
 * The class comes with a facilitated constructor implementing 
 * few pre-defined simple filters that could be used for waveform 
 * smoothing through the #TWDeck interface.
 */
class TWDeckWfmFilter : public TWDeckWfm {
  public : 
    //! Empty constructor
    TWDeckWfmFilter();
    //! Create an empty filter of size N
    TWDeckWfmFilter(int N);
    //! Create a filter in the real domain of size `N` and values `x`
    TWDeckWfmFilter(int N, double* x);
    //! Create a filter in the complex domain of size `N`
    TWDeckWfmFilter(int N, double* re, double* im);
    //! Facilitated constructor for predefined filters
    TWDeckWfmFilter(int N, wdeck::EFltrShape kShape, double bw);
    //! Copy constructor
    TWDeckWfmFilter(const TWDeckWfmFilter& filter);
    //! Filter destructor
    ~TWDeckWfmFilter();

    //! Set filter bandwidth
    void SetBandwidth(double x) {fBandwidth = x;}
    //! Get filter bandwidth
    double GetBandwidth() {return fBandwidth;}
    //! Set filter shift
    void SetShift(double s) {fShift = s;}
    //! Get filter shift
    double GetShift() {return fShift;}

  protected:
    double fBandwidth;  //!< Filter bandwidth
    double fShift;      //!< Filter shift

  public:
    ClassDef(TWDeckWfmFilter, 1)
};


#endif /* end of include guard TPDECWFMFILTER_H */

