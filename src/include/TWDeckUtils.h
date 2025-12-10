/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : TWDeckUtils.h
 * @created     : mercoledì gen 26, 2022 11:37:35 CET
 *
 * This file includes a few public utility functions to be 
 * used in the context of WaveDeck (or outside)
 */

#ifndef TWDECKUTILS_H

#define TWDECKUTILS_H

#include <vector>
#include "TF1.h"
#include "TGraph.h"

/**
 * @brief Returns a std::vector of N equally spaced points between a and b
 *
 * @tparam T number type 
 * @param a min value
 * @param b max value
 * @param N number of points
 *
 * @return 
 */
template <typename T>
inline std::vector<T> linspace(T a, T b, size_t N) {
  T h = (b - a) / static_cast<T>(N-1);
  std::vector<T> xs(N);
  typename std::vector<T>::iterator x;
  T val;
  for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
    *x = val;
  return xs;
}


/**
 * @brief Compute the integral of a TGraph object
 *
 * @param g Target `TGraph`
 * @param x0 Integration lower boundary
 * @param x1 Integration upper boundary
 *
 * @return 
 */
inline double g_integral(TGraph* g, double x0, double x1) {
  auto fc_gintegral = [g](double *x, double* p) {
    return p[0]*g->Eval(x[0]);
  };
  TF1 f("f", fc_gintegral, x0, x1, 1);
  f.SetNpx(g->GetN());
  f.SetParameter(0, 1.);
  return f.Integral(x0, x1, 1e-1);
}

/**
 * @brief Scale all the Y-points of a `TGraph`
 *
 * @param g targer `TGraph`
 * @param c scale factor
 */
inline void g_scale_Y(TGraph* g, double c) {
  for (int i=0; i<g->GetN(); i++) {
    g->GetY()[i] *= c;
  }
  return;
}

/**
 * @brief Scale all the X-points of a `TGraph`
 *
 * @param g target `TGraph`
 * @param c scale factor
 */
inline void g_scale_X(TGraph* g, double c) {
  for (int i=0; i<g->GetN(); i++) {
    g->GetX()[i] *= c;
  }
  return;
}

/**
 * @brief Find the value of x corresponding to a given value of y in the range [x0, x1]
 *
 * @param g target `TGraph`
 * @param y target y value
 * @param x0 lower bound
 * @param x1 upper bound
 *
 * @return 
 */
inline double g_find_x(TGraph* g, double y, double x0, double x1) {
  auto fc_g= [g](double *x, double* p) {
    return p[0]*g->Eval(x[0]);
  };
  TF1 f("f", fc_g, x0, x1, 1);
  f.SetNpx(g->GetN());
  f.SetParameter(0, 1.);

  return f.GetX(y, x0, x1, 1e-2);
}
#endif /* end of include guard TWDECKUTILS_H */

