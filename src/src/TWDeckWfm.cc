/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : TWDeckWfm
 * @created     : mercoledì gen 26, 2022 09:05:05 CET
 */

#include "TWDeckWfm.h"

ClassImp(TWDeckWfm)

TWDeckWfm::TWDeckWfm() : fSize(0), fWfm(0), fWfm_re(0), fWfm_im(0) {}

/**
 * @details Create an empty waveform of size \a N
 *
 * @param N waveform size
 */
TWDeckWfm::TWDeckWfm(int N) : fSize(0), fWfm(0), fWfm_re(0), fWfm_im(0) 
{
  SetSize(N);
}

/**
 * @details Create a real waveform of size \a N from the values 
 * stored in the \a wfm array
 *
 * @param N waveform size
 * @param wfm waveform array values
 */
TWDeckWfm::TWDeckWfm(int N, double* wfm) : fSize(0), fWfm(0), fWfm_re(0), fWfm_im(0) 
{
  LoadWave(N, wfm);
}


/**
 * @details Create a waveform passing its complex
 * representation as \a N pairs from the #re and #im
 * arrays as their real and imaginary part.
 *
 * @param N waveform size
 * @param re real parts of the Fourier components
 * @param im imaginary parts of the Fourier components
 */
TWDeckWfm::TWDeckWfm(int N, double* re, double* im) 
  : fSize(0), fWfm(0), fWfm_re(0), fWfm_im(0) 
{
  LoadSpectrum(N, re, im);
}

TWDeckWfm::TWDeckWfm(const TWDeckWfm& wfm) : fSize(0), fWfm(0), fWfm_re(0), fWfm_im(0) 
{
  fName = wfm.fName;
  fTitle = wfm.fTitle;
  fSize = wfm.fSize;
  SetSize(fSize);

  for (int i=0; i<fSize; i++) {
    fWfm[i] = wfm.fWfm[i];
    fWfm_re[i] = wfm.fWfm_re[i];
    fWfm_im[i] = wfm.fWfm_im[i];
  }
}

TWDeckWfm::~TWDeckWfm() {}

/**
 * @details Reset all waveform values, both real and complex, to 0.
 */
void TWDeckWfm::ClearWave() {
  for (auto &v : fWfm   ) v = 0.;
  for (auto &v : fWfm_re) v = 0.;
  for (auto &v : fWfm_im) v = 0.;
}

/**
 * Set the real values of the waveform. 
 * Note that since the size of the array \a w is not specified, 
 * the TCWDeckWfm size #fSize is assumed
 *
 * @param w waveform values array
 */
void TWDeckWfm::LoadWave(double* w) {
  for (int i=0; i<fSize; i++) {
    fWfm[i] = w[i];
  }
}

/**
 * Set the real values of waveform of size n from the array w
 *
 * @param n waveform size
 * @param w waveform values array
 */
void TWDeckWfm::LoadWave(int n, double* w) {
  SetSize(n);
  LoadWave(w);
}

/**
 * Set the waveform Fourier components. 
 * Note that since the size of the real and imaginary parts arrays is not
 * specified, the TWDeckWfm size #fSize is assumed
 *
 * @param re Real part array
 * @param im Imaginary part array
 */
void TWDeckWfm::LoadSpectrum(double* re, double* im) {
  for (int i=0; i<fSize; i++) {
    fWfm_re[i] = re[i];
    fWfm_im[i] = im[i];
  }
}

/**
 * @param n Array size
 * @param re Real part array
 * @param im Imaginary part array
 */
void TWDeckWfm::LoadSpectrum(int N, double* re, double* im) {
  SetSize(N);
  LoadSpectrum(re, im); 
}

/**
 * @details resize the waveform containers (both real and complex)
 * to have size n
 *
 * @param n new waveform size
 */
void TWDeckWfm::SetSize(int n) {
  fSize = n;

  fWfm   .resize(fSize, 0.);
  fWfm_re.resize(fSize, 0.);
  fWfm_im.resize(fSize, 0.);
}

std::vector<TComplex> TWDeckWfm::GetPointsComplex() {
  std::vector<TComplex> vWfmC(fSize);
  for (int i=0; i<fSize; i++) {
    vWfmC[i] = TComplex(fWfm_re[i], fWfm_im[i]);
  }

  return vWfmC;
}

/**
 * @details Returns the spectral density (i.e. the squared modulus) of the 
 * i-th Fourier components
 */
double TWDeckWfm::GetSpectralDensity(int i) {
  TComplex c(fWfm_re.at(i), fWfm_im.at(i));
  return c.Rho2();
}

/**
 * @details Returns a vector containing the values of the spectral 
 * density (i.e. the squared modulus) of the waveform Fourier components
 */
std::vector<double> TWDeckWfm::GetSpectralDensityPoints() {
  std::vector<double> rho2(fSize*0.5+1, 0.);
  for (int i=0; i<fSize*0.5+1; i++)
    rho2.at(i) = GetSpectralDensity(i);
  return rho2;
}


