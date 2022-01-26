/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : TWDeckUtils
 * @created     : mercoledì gen 26, 2022 11:37:35 CET
 */

#ifndef TWDECKUTILS_H

#define TWDECKUTILS_H

#include <vector>

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



#endif /* end of include guard TWDECKUTILS_H */

